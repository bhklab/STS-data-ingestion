from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote

import requests
from sqlalchemy import inspect, select, text
from sqlalchemy.dialects.mysql import insert as mysql_insert
from sqlalchemy.orm import Session

from .seeding_coordinator_engine import alchemy_engine
from ..models.tables import Base, PreClinicalCellLine, PreClinicalCellLineInfo


DEFAULT_ANNOTATIONDB_BASE_URL = "https://annotationdb.bhklab.ca"
DEFAULT_ANNOTATIONDB_ENDPOINT = "/cell_line/many"
DEFAULT_ONCOTREE_BASE_URL = "https://oncotree.mskcc.org"
DISEASE_LIST_DELIMITER = "|"
DEFAULT_BATCH_SIZE = 250
DEFAULT_AUDIT_DIR = Path("extraction/data/proc/preclinical/cell_line_info")
DEFAULT_RETURNED_CSV = DEFAULT_AUDIT_DIR / "cell_line_accessions_returned_annotationdb.csv"
DEFAULT_NOT_FOUND_CSV = DEFAULT_AUDIT_DIR / "cell_line_accessions_not_found_annotationdb.csv"
DEFAULT_ONCOTREE_UNRESOLVED_CSV = DEFAULT_AUDIT_DIR / "cell_line_info_oncotree_unresolved.csv"

# Top-level fields returned by AnnotationDB for /cell_line/many, excluding the
# nested `diseases` list which is flattened separately.
CELL_LINE_INFO_FIELDS: tuple[str, ...] = (
    "accession",
    "cell_line_name",
    "category",
    "date",
    "age_at_sampling",
    "sex_of_cell",
    "hierarchy",
    "cell_type",
    "derived_from_site",
    "donor_information",
    "doubling_time",
    "genome_ancestry",
    "hla_typing",
    "microsatellite_instability",
    "omics",
    "part_of",
    "population",
    "sequence_variation",
    "anecdotal",
    "biotechnology",
    "discontinued",
    "group_col",
    "misspelling",
    "registration",
    "virology",
    "caution",
    "characteristics",
    "karyotypic_information",
    "problematic_cell_line",
    "transformant",
    "miscellaneous",
    "from_col",
    "genetic_integration",
    "knockout_cell",
    "selected_for_resistance_to",
)

CREATE_ONLY_TABLES = [PreClinicalCellLineInfo.__table__]
ONCOTREE_LEVELS = '"1,2,3,4,5,6,7"'


def clean_value(value: Any) -> Any | None:
    if value is None:
        return None

    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None

    if isinstance(value, str):
        value = value.strip()
        if value == "" or value.upper() in {"NA", "N/A", "NS", "NAN", "NONE", "NULL"}:
            return None
        return value

    return value


def clean_str(value: Any) -> str | None:
    value = clean_value(value)
    if value is None:
        return None
    return str(value)


def normalize_accession(value: Any) -> str | None:
    value = clean_str(value)
    if value is None:
        return None
    return value.strip()


def ordered_unique(values: Iterable[str | None]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        value = clean_str(value)
        if value is None:
            continue
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def chunked(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def ensure_accession_index(engine: Any, *, verbose: bool = True) -> None:
    """Ensure MySQL has an index on pre_clinical_cell_line.accession.

    SQLAlchemy create_all() does not add a new index to an existing table. The
    FK from pre_clinical_cell_line_info.accession to pre_clinical_cell_line.accession
    requires the referenced column to be indexed in MySQL.
    """
    with engine.begin() as conn:
        rows = conn.execute(
            text(
                """
                SELECT COUNT(*) AS n
                FROM information_schema.statistics
                WHERE table_schema = DATABASE()
                  AND table_name = 'pre_clinical_cell_line'
                  AND index_name = 'ix_pc_cell_line_accession'
                """
            )
        ).mappings().one()

        if int(rows["n"]) == 0:
            if verbose:
                print("Creating index ix_pc_cell_line_accession on pre_clinical_cell_line(accession)")
            conn.execute(text("ALTER TABLE pre_clinical_cell_line ADD INDEX ix_pc_cell_line_accession (accession)"))


def create_required_tables(engine: Any, *, verbose: bool = True) -> None:
    ensure_accession_index(engine, verbose=verbose)

    existing_tables = set(inspect(engine).get_table_names())
    if "pre_clinical_cell_line" not in existing_tables:
        raise RuntimeError(
            "pre_clinical_cell_line does not exist. Seed the dataset cell-line tables before seeding cell-line info."
        )

    Base.metadata.create_all(bind=engine, tables=CREATE_ONLY_TABLES)


def get_unique_cell_line_accessions(session: Session) -> list[str]:
    rows = session.execute(
        select(PreClinicalCellLine.accession)
        .where(PreClinicalCellLine.accession.is_not(None))
        .distinct()
        .order_by(PreClinicalCellLine.accession)
    ).scalars()

    accessions = [normalize_accession(row) for row in rows]
    return [accession for accession in accessions if accession]


def request_json_with_retries(
    *,
    url: str,
    params: dict[str, Any] | list[tuple[str, Any]] | None,
    timeout_seconds: int,
    max_retries: int,
    retry_sleep_seconds: float,
) -> Any:
    last_error: Exception | None = None

    for attempt in range(max_retries + 1):
        try:
            response = requests.get(
                url,
                params=params,
                headers={"accept": "application/json"},
                timeout=timeout_seconds,
            )
            if response.status_code == 200:
                return response.json()

            snippet = response.text[:1000]
            raise RuntimeError(
                f"HTTP {response.status_code} for {response.url}. Response snippet: {snippet}"
            )
        except Exception as exc:  # requests and JSON failures
            last_error = exc
            if attempt >= max_retries:
                break
            time.sleep(retry_sleep_seconds)

    raise RuntimeError(f"Request failed after {max_retries + 1} attempts for {url}: {last_error}")


def request_annotationdb_cell_line_batch(
    *,
    accessions: list[str],
    annotationdb_base_url: str,
    annotationdb_endpoint: str,
    timeout_seconds: int,
    max_retries: int,
    retry_sleep_seconds: float,
) -> list[dict[str, Any]]:
    if not accessions:
        return []

    url = annotationdb_base_url.rstrip("/") + "/" + annotationdb_endpoint.strip("/")
    data = request_json_with_retries(
        url=url,
        params={"cell_lines": ",".join(accessions), "format": "json"},
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        retry_sleep_seconds=retry_sleep_seconds,
    )

    if data is None:
        return []
    if isinstance(data, list):
        return [row for row in data if isinstance(row, dict)]
    if isinstance(data, dict):
        # Defensive support in case AnnotationDB wraps results in a key.
        for key in ("data", "results", "cell_lines"):
            value = data.get(key)
            if isinstance(value, list):
                return [row for row in value if isinstance(row, dict)]
    raise RuntimeError(f"Unexpected AnnotationDB response shape for batch starting {accessions[0]}: {type(data)}")


class OncoTreeResolver:
    def __init__(
        self,
        *,
        base_url: str = DEFAULT_ONCOTREE_BASE_URL,
        timeout_seconds: int = 30,
        max_retries: int = 2,
        retry_sleep_seconds: float = 1.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout_seconds = timeout_seconds
        self.max_retries = max_retries
        self.retry_sleep_seconds = retry_sleep_seconds
        self._name_cache: dict[str, list[dict[str, Any]]] = {}
        self._code_cache: dict[str, list[dict[str, Any]]] = {}

    def _request_oncotree_json_or_none(
        self,
        *,
        url: str,
        params: dict[str, Any] | None,
    ) -> Any | None:
        """Request OncoTree JSON.

        OncoTree uses HTTP 404 for a valid "no tumor types found" result.
        For this seeder, 404 is not an error and should not be retried: it
        means "this description/code did not resolve; try the next disease
        description if one exists."
        """
        last_error: Exception | None = None

        for attempt in range(self.max_retries + 1):
            try:
                response = requests.get(
                    url,
                    params=params,
                    headers={"accept": "application/json"},
                    timeout=self.timeout_seconds,
                )

                if response.status_code == 200:
                    return response.json()

                if response.status_code == 404:
                    return None

                snippet = response.text[:1000]
                raise RuntimeError(
                    f"HTTP {response.status_code} for {response.url}. Response snippet: {snippet}"
                )
            except Exception as exc:
                last_error = exc
                if attempt >= self.max_retries:
                    break
                time.sleep(self.retry_sleep_seconds)

        raise RuntimeError(
            f"OncoTree request failed after {self.max_retries + 1} attempts for {url}: {last_error}"
        )

    def search_by_name(self, disease_description: str) -> list[dict[str, Any]]:
        key = disease_description.strip().lower()
        if key in self._name_cache:
            return self._name_cache[key]

        url = f"{self.base_url}/api/tumorTypes/search/name/{quote(disease_description, safe='')}"
        data = self._request_oncotree_json_or_none(
            url=url,
            params={"levels": ONCOTREE_LEVELS},
        )

        # No initial disease-name hit: do not recurse by code; move to the next
        # disease description for the same cell line.
        if data is None:
            self._name_cache[key] = []
            return []

        rows = data if isinstance(data, list) else []
        rows = [row for row in rows if isinstance(row, dict)]
        self._name_cache[key] = rows
        return rows

    def search_by_code(self, code: str) -> list[dict[str, Any]]:
        key = code.strip().upper()
        if key in self._code_cache:
            return self._code_cache[key]

        url = f"{self.base_url}/api/tumorTypes/search/code/{quote(code, safe='')}"
        data = self._request_oncotree_json_or_none(
            url=url,
            params={"levels": ONCOTREE_LEVELS},
        )

        # Parent-code miss after a valid name hit means this particular
        # disease-description traversal failed. The caller can still try the
        # next disease description.
        if data is None:
            self._code_cache[key] = []
            return []

        rows = data if isinstance(data, list) else []
        rows = [row for row in rows if isinstance(row, dict)]
        self._code_cache[key] = rows
        return rows

    @staticmethod
    def pick_name_hit(rows: list[dict[str, Any]], disease_description: str) -> dict[str, Any] | None:
        if not rows:
            return None

        target = disease_description.strip().lower()
        exact = [row for row in rows if str(row.get("name", "")).strip().lower() == target]
        if exact:
            return exact[0]

        # Prefer the deepest/most specific hit if exact matching is unavailable.
        def level_value(row: dict[str, Any]) -> int:
            try:
                return int(row.get("level") or -1)
            except (TypeError, ValueError):
                return -1

        return sorted(rows, key=level_value, reverse=True)[0]

    @staticmethod
    def pick_code_hit(rows: list[dict[str, Any]], code: str) -> dict[str, Any] | None:
        if not rows:
            return None

        target = code.strip().upper()
        exact = [row for row in rows if str(row.get("code", "")).strip().upper() == target]
        return exact[0] if exact else rows[0]

    def resolve_level_2_from_description(self, disease_description: str) -> dict[str, str | None]:
        """Resolve one disease description to OncoTree level 2.

        The first request is by disease name. If a hit is returned, parent codes
        are queried recursively with the code endpoint until `level == 2`.
        """
        disease_description = disease_description.strip()
        if not disease_description:
            return {
                "first_level": None,
                "second_level": None,
                "matched_description": disease_description,
                "matched_code": None,
                "level_2_code": None,
                "reason": "empty disease description",
            }

        current = self.pick_name_hit(self.search_by_name(disease_description), disease_description)
        if current is None:
            return {
                "first_level": None,
                "second_level": None,
                "matched_description": disease_description,
                "matched_code": None,
                "level_2_code": None,
                "reason": "no OncoTree name hit",
            }

        visited_codes: set[str] = set()
        while current is not None:
            code = clean_str(current.get("code"))
            parent = clean_str(current.get("parent"))
            name = clean_str(current.get("name"))

            try:
                level = int(current.get("level"))
            except (TypeError, ValueError):
                return {
                    "first_level": None,
                    "second_level": None,
                    "matched_description": disease_description,
                    "matched_code": code,
                    "level_2_code": None,
                    "reason": f"invalid OncoTree level for code {code}",
                }

            if level == 2:
                return {
                    "first_level": parent.replace("_", " ") if parent else None,
                    "second_level": name,
                    "matched_description": disease_description,
                    "matched_code": code,
                    "level_2_code": code,
                    "reason": None,
                }

            if level < 2:
                return {
                    "first_level": parent.replace("_", " ") if parent else None,
                    "second_level": name,
                    "matched_description": disease_description,
                    "matched_code": code,
                    "level_2_code": None,
                    "reason": f"reached level {level} before level 2",
                }

            if not parent:
                return {
                    "first_level": None,
                    "second_level": None,
                    "matched_description": disease_description,
                    "matched_code": code,
                    "level_2_code": None,
                    "reason": f"missing parent for code {code}",
                }

            parent_key = parent.upper()
            if parent_key in visited_codes:
                return {
                    "first_level": None,
                    "second_level": None,
                    "matched_description": disease_description,
                    "matched_code": code,
                    "level_2_code": None,
                    "reason": f"cycle detected at parent {parent}",
                }
            visited_codes.add(parent_key)

            current = self.pick_code_hit(self.search_by_code(parent), parent)

        return {
            "first_level": None,
            "second_level": None,
            "matched_description": disease_description,
            "matched_code": None,
            "level_2_code": None,
            "reason": "OncoTree code traversal failed",
        }

    def resolve_first_success(self, disease_descriptions: list[str]) -> dict[str, str | None]:
        if not disease_descriptions:
            return {
                "first_level": None,
                "second_level": None,
                "matched_description": None,
                "matched_code": None,
                "level_2_code": None,
                "reason": "no disease descriptions",
            }

        last_result: dict[str, str | None] | None = None
        for description in disease_descriptions:
            result = self.resolve_level_2_from_description(description)
            if result.get("second_level"):
                return result
            last_result = result

        return last_result or {
            "first_level": None,
            "second_level": None,
            "matched_description": None,
            "matched_code": None,
            "level_2_code": None,
            "reason": "no disease descriptions resolved",
        }


def flatten_diseases(record: dict[str, Any]) -> tuple[str | None, str | None, list[str]]:
    diseases = record.get("diseases")
    if not isinstance(diseases, list):
        return None, None, []

    disease_ids: list[str] = []
    disease_descriptions: list[str] = []
    for disease in diseases:
        if not isinstance(disease, dict):
            continue
        disease_id = clean_str(disease.get("id"))
        disease_description = clean_str(disease.get("description"))
        if disease_id:
            disease_ids.append(disease_id)
        if disease_description:
            disease_descriptions.append(disease_description)

    return (
        DISEASE_LIST_DELIMITER.join(disease_ids) if disease_ids else None,
        DISEASE_LIST_DELIMITER.join(disease_descriptions) if disease_descriptions else None,
        ordered_unique(disease_descriptions),
    )


def build_cell_line_info_row(
    record: dict[str, Any],
    *,
    oncotree_resolver: OncoTreeResolver | None,
    skip_oncotree: bool,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    accession = normalize_accession(record.get("accession"))
    if accession is None:
        return None, None

    disease_ids, disease_descriptions, disease_description_list = flatten_diseases(record)

    oncotree_result: dict[str, str | None] = {
        "first_level": None,
        "second_level": None,
        "matched_description": None,
        "matched_code": None,
        "level_2_code": None,
        "reason": "OncoTree lookup skipped" if skip_oncotree else None,
    }

    if not skip_oncotree:
        if oncotree_resolver is None:
            raise RuntimeError("oncotree_resolver cannot be None when skip_oncotree is False")
        oncotree_result = oncotree_resolver.resolve_first_success(disease_description_list)

    row = {field: clean_str(record.get(field)) for field in CELL_LINE_INFO_FIELDS}
    row["accession"] = accession
    row["disease_ids"] = disease_ids
    row["disease_descriptions"] = disease_descriptions
    row["first_level"] = clean_str(oncotree_result.get("first_level"))
    row["second_level"] = clean_str(oncotree_result.get("second_level"))

    unresolved = None
    if not row["second_level"]:
        unresolved = {
            "accession": accession,
            "cell_line_name": row.get("cell_line_name"),
            "disease_descriptions": disease_descriptions,
            "matched_description": oncotree_result.get("matched_description"),
            "matched_code": oncotree_result.get("matched_code"),
            "reason": oncotree_result.get("reason"),
        }

    return row, unresolved


def upsert_cell_line_info(session: Session, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return

    table = PreClinicalCellLineInfo.__table__
    stmt = mysql_insert(table).values(rows)
    update_cols = {
        column.name: stmt.inserted[column.name]
        for column in table.columns
        if column.name != "accession"
    }
    stmt = stmt.on_duplicate_key_update(**update_cols)
    session.execute(stmt)


def seed_cell_line_info_from_annotationdb(
    *,
    batch_size: int,
    annotationdb_base_url: str,
    annotationdb_endpoint: str,
    oncotree_base_url: str,
    audit_dir: Path,
    timeout_seconds: int,
    max_retries: int,
    retry_sleep_seconds: float,
    continue_on_api_error: bool,
    skip_oncotree: bool,
    limit: int | None,
) -> None:
    engine = alchemy_engine()
    create_required_tables(engine)

    returned_rows: list[dict[str, Any]] = []
    not_found_rows: list[dict[str, Any]] = []
    unresolved_oncotree_rows: list[dict[str, Any]] = []

    with Session(engine) as session:
        accessions = get_unique_cell_line_accessions(session)
        if limit is not None:
            accessions = accessions[:limit]

        print(f"Unique non-null cell-line accessions found: {len(accessions)}")
        if not accessions:
            return

        oncotree_resolver = None if skip_oncotree else OncoTreeResolver(
            base_url=oncotree_base_url,
            timeout_seconds=timeout_seconds,
            max_retries=max_retries,
            retry_sleep_seconds=retry_sleep_seconds,
        )

        for batch_number, batch in enumerate(chunked(accessions, batch_size), start=1):
            print(f"Querying AnnotationDB cell-line batch {batch_number}: {len(batch)} accessions")

            try:
                records = request_annotationdb_cell_line_batch(
                    accessions=batch,
                    annotationdb_base_url=annotationdb_base_url,
                    annotationdb_endpoint=annotationdb_endpoint,
                    timeout_seconds=timeout_seconds,
                    max_retries=max_retries,
                    retry_sleep_seconds=retry_sleep_seconds,
                )
            except Exception as exc:
                if not continue_on_api_error:
                    raise
                print(f"WARNING: AnnotationDB batch failed for first accession {batch[0]}: {exc}")
                not_found_rows.extend({"accession": accession, "reason": f"api_error: {exc}"} for accession in batch)
                continue

            returned_accessions = {
                normalize_accession(record.get("accession"))
                for record in records
                if normalize_accession(record.get("accession"))
            }

            missing = sorted(set(batch) - returned_accessions)
            not_found_rows.extend({"accession": accession, "reason": "not returned by AnnotationDB"} for accession in missing)

            rows_to_insert: list[dict[str, Any]] = []
            for record in records:
                try:
                    row, unresolved = build_cell_line_info_row(
                        record,
                        oncotree_resolver=oncotree_resolver,
                        skip_oncotree=skip_oncotree,
                    )
                except Exception as exc:
                    accession = normalize_accession(record.get("accession"))
                    if not continue_on_api_error:
                        raise
                    print(f"WARNING: Failed to prepare accession {accession}: {exc}")
                    unresolved_oncotree_rows.append(
                        {
                            "accession": accession,
                            "cell_line_name": clean_str(record.get("cell_line_name")),
                            "disease_descriptions": None,
                            "matched_description": None,
                            "matched_code": None,
                            "reason": f"prepare_error: {exc}",
                        }
                    )
                    continue

                if row is None:
                    continue
                rows_to_insert.append(row)
                returned_rows.append(
                    {
                        "accession": row.get("accession"),
                        "cell_line_name": row.get("cell_line_name"),
                        "first_level": row.get("first_level"),
                        "second_level": row.get("second_level"),
                    }
                )
                if unresolved is not None:
                    unresolved_oncotree_rows.append(unresolved)

            upsert_cell_line_info(session, rows_to_insert)
            session.commit()
            print(f"Seeded cell-line info rows from batch {batch_number}: {len(rows_to_insert)}")

    audit_dir.mkdir(parents=True, exist_ok=True)
    returned_csv = audit_dir / DEFAULT_RETURNED_CSV.name
    not_found_csv = audit_dir / DEFAULT_NOT_FOUND_CSV.name
    unresolved_csv = audit_dir / DEFAULT_ONCOTREE_UNRESOLVED_CSV.name

    write_csv(
        returned_csv,
        returned_rows,
        ["accession", "cell_line_name", "first_level", "second_level"],
    )
    write_csv(not_found_csv, not_found_rows, ["accession", "reason"])
    write_csv(
        unresolved_csv,
        unresolved_oncotree_rows,
        ["accession", "cell_line_name", "disease_descriptions", "matched_description", "matched_code", "reason"],
    )

    print(f"Returned AnnotationDB audit CSV: {returned_csv}")
    print(f"Not-found AnnotationDB audit CSV: {not_found_csv}")
    print(f"Unresolved OncoTree audit CSV: {unresolved_csv}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Seed pre_clinical_cell_line_info from AnnotationDB cell_line/many, "
            "then resolve disease descriptions to OncoTree level 2."
        )
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--annotationdb-base-url", default=DEFAULT_ANNOTATIONDB_BASE_URL)
    parser.add_argument("--annotationdb-endpoint", default=DEFAULT_ANNOTATIONDB_ENDPOINT)
    parser.add_argument("--oncotree-base-url", default=DEFAULT_ONCOTREE_BASE_URL)
    parser.add_argument("--audit-dir", type=Path, default=DEFAULT_AUDIT_DIR)
    parser.add_argument("--timeout-seconds", type=int, default=30)
    parser.add_argument("--max-retries", type=int, default=2)
    parser.add_argument("--retry-sleep-seconds", type=float, default=1.0)
    parser.add_argument("--continue-on-api-error", action="store_true")
    parser.add_argument("--skip-oncotree", action="store_true")
    parser.add_argument("--limit", type=int, default=None)

    args = parser.parse_args()

    if args.batch_size <= 0:
        raise ValueError("--batch-size must be positive")

    seed_cell_line_info_from_annotationdb(
        batch_size=args.batch_size,
        annotationdb_base_url=args.annotationdb_base_url,
        annotationdb_endpoint=args.annotationdb_endpoint,
        oncotree_base_url=args.oncotree_base_url,
        audit_dir=args.audit_dir,
        timeout_seconds=args.timeout_seconds,
        max_retries=args.max_retries,
        retry_sleep_seconds=args.retry_sleep_seconds,
        continue_on_api_error=args.continue_on_api_error,
        skip_oncotree=args.skip_oncotree,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()