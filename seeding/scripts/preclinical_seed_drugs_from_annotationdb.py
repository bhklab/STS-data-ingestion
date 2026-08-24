from __future__ import annotations

import argparse
import csv
import math
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import requests
from sqlalchemy import select
from sqlalchemy.dialects.mysql import insert as mysql_insert
from sqlalchemy.orm import Session

from .seeding_coordinator_engine import alchemy_engine
from ..models.tables import Base, PreClinicalDrug, PreClinicalTreatmentResponse


DEFAULT_ANNOTATIONDB_BASE_URL = "https://annotationdb.bhklab.ca"
DEFAULT_ANNOTATIONDB_ENDPOINTS = (
    # AnnotationDB expects repeated `compound=` query parameters.
    # Example:
    #   /compound/many?compound=2244&compound=53355&format=json&mechanism=true
    "/compound/many",
    # Older/fallback candidates kept for compatibility if the deployed route changes.
    "/compounds/compound/many",
    "/compounds/many",
)
DEFAULT_BATCH_SIZE = 250
DEFAULT_AUDIT_DIR = Path("extraction/data/proc/preclinical/drugs")
DEFAULT_NOT_FOUND_AUDIT_CSV = DEFAULT_AUDIT_DIR / "drug_cids_not_found_annotationdb.csv"
DEFAULT_RETURNED_CSV = DEFAULT_AUDIT_DIR / "drug_cids_returned_annotationdb.csv"


DRUG_TOP_LEVEL_FIELDS: tuple[str, ...] = (
    "cid",
    "title",
    "mapped_name",
    "molecule_chembl_id",
    "molecule_chembl_id_from_synonyms",
    "molecular_formula",
    "molecular_weight",
    "smiles",
    "connectivity_smiles",
    "inchi",
    "inchikey",
    "iupac_name",
    "xlogp",
    "exact_mass",
    "monoisotopic_mass",
    "tpsa",
    "complexity",
    "charge",
    "h_bond_donor_count",
    "h_bond_acceptor_count",
    "rotatable_bond_count",
    "heavy_atom_count",
    "isotope_atom_count",
    "atom_stereo_count",
    "defined_atom_stereo_count",
    "undefined_atom_stereo_count",
    "bond_stereo_count",
    "defined_bond_stereo_count",
    "undefined_bond_stereo_count",
    "covalent_unit_count",
    "volume_3d",
    "x_steric_quadrupole_3d",
    "y_steric_quadrupole_3d",
    "z_steric_quadrupole_3d",
    "feature_count_3d",
    "feature_acceptor_count_3d",
    "feature_donor_count_3d",
    "feature_anion_count_3d",
    "feature_cation_count_3d",
    "feature_ring_count_3d",
    "feature_hydrophobe_count_3d",
    "conformer_model_rmsd_3d",
    "effective_rotor_count_3d",
    "conformer_count_3d",
    "fingerprint_2d",
    "patent_count",
    "patent_family_count",
    "literature_count",
    "annotation_types",
    "annotation_type_count",
    "chembl_max_phase",
    "drug_like",
    "fda_approval",
    "date_added",
    "atc_code",
)

INT_FIELDS = {
    "charge",
    "h_bond_donor_count",
    "h_bond_acceptor_count",
    "rotatable_bond_count",
    "heavy_atom_count",
    "isotope_atom_count",
    "atom_stereo_count",
    "defined_atom_stereo_count",
    "undefined_atom_stereo_count",
    "bond_stereo_count",
    "defined_bond_stereo_count",
    "undefined_bond_stereo_count",
    "covalent_unit_count",
    "feature_count_3d",
    "feature_acceptor_count_3d",
    "feature_donor_count_3d",
    "feature_anion_count_3d",
    "feature_cation_count_3d",
    "feature_ring_count_3d",
    "feature_hydrophobe_count_3d",
    "conformer_count_3d",
    "patent_count",
    "patent_family_count",
    "literature_count",
    "annotation_type_count",
    "chembl_max_phase",
}

FLOAT_FIELDS = {
    "xlogp",
    "tpsa",
    "complexity",
    "volume_3d",
    "x_steric_quadrupole_3d",
    "y_steric_quadrupole_3d",
    "z_steric_quadrupole_3d",
    "conformer_model_rmsd_3d",
    "effective_rotor_count_3d",
}

BOOL_FIELDS = {
    "molecule_chembl_id_from_synonyms",
    "drug_like",
    "fda_approval",
}


CREATE_ONLY_TABLES = [PreClinicalDrug.__table__]


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


def clean_cid(value: Any) -> str | None:
    value = clean_value(value)
    if value is None:
        return None

    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]

    return text or None


def clean_int(value: Any) -> int | None:
    value = clean_value(value)
    if value is None:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def clean_float(value: Any) -> float | None:
    value = clean_value(value)
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(out):
        return None
    return out


def clean_bool(value: Any) -> bool | None:
    value = clean_value(value)
    if value is None:
        return None
    if isinstance(value, bool):
        return value

    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return None


def clean_datetime(value: Any) -> datetime | None:
    value = clean_value(value)
    if value is None:
        return None

    text = str(value).strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"

    try:
        return datetime.fromisoformat(text)
    except ValueError:
        return None


def chunked(values: list[str], batch_size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), batch_size):
        yield values[start : start + batch_size]


def create_required_tables(engine: Any) -> None:
    Base.metadata.create_all(bind=engine, tables=CREATE_ONLY_TABLES)


def query_unique_treatment_response_cids(session: Session) -> list[str]:
    rows = session.execute(
        select(PreClinicalTreatmentResponse.cid)
        .where(PreClinicalTreatmentResponse.cid.is_not(None))
        .distinct()
    ).scalars()

    cids = sorted({cid for cid in (clean_cid(row) for row in rows) if cid is not None})
    return cids


def build_annotationdb_many_params(cids: list[str]) -> list[tuple[str, str]]:
    """Build the AnnotationDB /compound/many query-string shape.

    Expected API form:
        /compound/many?compound=2244&compound=53355&format=json&mechanism=true

    Keep mechanism=true because the drugs table flattens the first mechanism object.
    """
    params: list[tuple[str, str]] = [("compound", str(cid)) for cid in cids]
    params.extend(
        [
            ("format", "json"),
            ("bioassay", "false"),
            ("mechanism", "true"),
            ("toxicity", "false"),
            ("golden_bioassay", "false"),
        ]
    )
    return params


def request_annotationdb_batch(
    *,
    cids: list[str],
    base_url: str,
    endpoints: tuple[str, ...],
    timeout: float,
    max_retries: int,
    retry_sleep_seconds: float,
) -> list[dict[str, Any]]:
    last_error: Exception | None = None
    base_url = base_url.rstrip("/")

    query_param_variants: tuple[Any, ...] = (
        # Correct AnnotationDB route shape from the docs/example.
        build_annotationdb_many_params(cids),

        # Fallback shapes retained only in case another deployed AnnotationDB route
        # accepts cids/cid. These still include flags so the response contains mechanisms.
        [("cids", cid) for cid in cids]
        + [
            ("format", "json"),
            ("bioassay", "false"),
            ("mechanism", "true"),
            ("toxicity", "false"),
            ("golden_bioassay", "false"),
        ],
        [("cid", cid) for cid in cids]
        + [
            ("format", "json"),
            ("bioassay", "false"),
            ("mechanism", "true"),
            ("toxicity", "false"),
            ("golden_bioassay", "false"),
        ],
        {
            "cids": ",".join(cids),
            "format": "json",
            "bioassay": "false",
            "mechanism": "true",
            "toxicity": "false",
            "golden_bioassay": "false",
        },
        {
            "cid": ",".join(cids),
            "format": "json",
            "bioassay": "false",
            "mechanism": "true",
            "toxicity": "false",
            "golden_bioassay": "false",
        },
    )

    for endpoint in endpoints:
        url = f"{base_url}/{endpoint.lstrip('/')}"

        for params in query_param_variants:
            for attempt in range(1, max_retries + 1):
                try:
                    response = requests.get(url, params=params, timeout=timeout)

                    if response.status_code == 404:
                        last_error = RuntimeError(
                            f"404 from {url}; tried params beginning with {list(params)[:3] if isinstance(params, list) else params}"
                        )
                        # Try the next endpoint/param shape.
                        break

                    if response.status_code >= 400:
                        snippet = response.text[:500].replace("\n", " ")
                        raise RuntimeError(
                            f"HTTP {response.status_code} from {url}: {snippet}"
                        )

                    payload = response.json()

                    if isinstance(payload, list):
                        return payload

                    if isinstance(payload, dict):
                        # FastAPI routes sometimes wrap list responses under a key.
                        for key in ("data", "results", "compounds", "items"):
                            maybe_records = payload.get(key)
                            if isinstance(maybe_records, list):
                                return maybe_records

                    raise ValueError(
                        f"Unexpected AnnotationDB response shape from {url}: {type(payload).__name__}"
                    )

                except Exception as exc:  # noqa: BLE001 - include HTTP/JSON errors in retry loop.
                    last_error = exc
                    if attempt < max_retries:
                        time.sleep(retry_sleep_seconds)
                    else:
                        break

    if last_error is not None:
        raise RuntimeError(f"AnnotationDB lookup failed for batch starting {cids[0]}: {last_error}")

    raise RuntimeError(f"AnnotationDB lookup failed for batch starting {cids[0]}.")


def normalize_annotationdb_record(record: dict[str, Any]) -> dict[str, Any] | None:
    cid = clean_cid(record.get("cid"))
    if cid is None:
        return None

    row: dict[str, Any] = {"cid": cid}

    for field in DRUG_TOP_LEVEL_FIELDS:
        if field == "cid":
            continue

        value = record.get(field)
        if field in INT_FIELDS:
            row[field] = clean_int(value)
        elif field in FLOAT_FIELDS:
            row[field] = clean_float(value)
        elif field in BOOL_FIELDS:
            row[field] = clean_bool(value)
        elif field == "date_added":
            row[field] = clean_datetime(value)
        else:
            row[field] = clean_str(value)

    mechanisms = record.get("mechanisms")
    first_mechanism = mechanisms[0] if isinstance(mechanisms, list) and mechanisms else {}
    if not isinstance(first_mechanism, dict):
        first_mechanism = {}

    row["mechanism_molecule_chembl_id"] = clean_str(first_mechanism.get("molecule_chembl_id"))
    row["mechanism_parent_molecule_chembl_id"] = clean_str(
        first_mechanism.get("parent_molecule_chembl_id")
    )
    row["mechanism_action_type"] = clean_str(first_mechanism.get("action_type"))
    row["mechanism_of_action"] = clean_str(first_mechanism.get("mechanism_of_action"))

    return row


def write_audit_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def upsert_drugs(session: Session, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return

    stmt = mysql_insert(PreClinicalDrug.__table__).values(rows)
    update_cols = {
        col.name: stmt.inserted[col.name]
        for col in PreClinicalDrug.__table__.columns
        if col.name != "cid"
    }
    stmt = stmt.on_duplicate_key_update(**update_cols)
    session.execute(stmt)


def seed_drugs_from_annotationdb(
    *,
    batch_size: int,
    annotationdb_base_url: str,
    annotationdb_endpoints: tuple[str, ...],
    timeout: float,
    max_retries: int,
    retry_sleep_seconds: float,
    audit_csv: Path,
    returned_csv: Path,
    continue_on_api_error: bool,
    limit: int | None,
) -> None:
    engine = alchemy_engine()
    create_required_tables(engine)

    all_missing_audit_rows: list[dict[str, Any]] = []
    all_returned_audit_rows: list[dict[str, Any]] = []
    total_inserted_or_updated = 0

    with Session(engine) as session:
        cids = query_unique_treatment_response_cids(session)
        if limit is not None:
            cids = cids[:limit]

        print(f"Unique treatment-response CIDs to query: {len(cids)}")

        for batch_number, batch in enumerate(chunked(cids, batch_size), start=1):
            print(f"Querying AnnotationDB batch {batch_number}: {len(batch)} CIDs")

            try:
                records = request_annotationdb_batch(
                    cids=batch,
                    base_url=annotationdb_base_url,
                    endpoints=annotationdb_endpoints,
                    timeout=timeout,
                    max_retries=max_retries,
                    retry_sleep_seconds=retry_sleep_seconds,
                )
            except Exception as exc:  # noqa: BLE001 - optionally audit and continue.
                if not continue_on_api_error:
                    raise

                reason = f"annotationdb_request_failed: {exc}"
                all_missing_audit_rows.extend(
                    {"cid": cid, "reason": reason, "batch_number": batch_number}
                    for cid in batch
                )
                continue

            normalized_rows = [normalize_annotationdb_record(record) for record in records]
            normalized_rows = [row for row in normalized_rows if row is not None]

            returned_cids = {row["cid"] for row in normalized_rows}
            requested_cids = set(batch)
            missing_cids = sorted(requested_cids - returned_cids)

            all_returned_audit_rows.extend(
                {"cid": row["cid"], "title": row.get("title"), "batch_number": batch_number}
                for row in normalized_rows
            )

            all_missing_audit_rows.extend(
                {"cid": cid, "reason": "not_returned_by_annotationdb", "batch_number": batch_number}
                for cid in missing_cids
            )

            upsert_drugs(session, normalized_rows)
            session.commit()
            total_inserted_or_updated += len(normalized_rows)

            print(
                "  Returned:",
                len(normalized_rows),
                "Missing:",
                len(missing_cids),
            )

    write_audit_csv(
        audit_csv,
        all_missing_audit_rows,
        fieldnames=["cid", "reason", "batch_number"],
    )
    write_audit_csv(
        returned_csv,
        all_returned_audit_rows,
        fieldnames=["cid", "title", "batch_number"],
    )

    print(f"Seeded/updated drugs rows: {total_inserted_or_updated}")
    print(f"Wrote missing-CID audit CSV: {audit_csv}")
    print(f"Wrote returned-CID audit CSV: {returned_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Query unique CIDs from pre_clinical_treatment_response, fetch compound "
            "annotations from AnnotationDB in batches, and seed the drugs table."
        )
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--annotationdb-base-url", default=DEFAULT_ANNOTATIONDB_BASE_URL)
    parser.add_argument(
        "--annotationdb-endpoint",
        action="append",
        default=None,
        help=(
            "AnnotationDB endpoint path. Can be passed multiple times. Defaults to "
            "/compound/many with older fallback candidates."
        ),
    )
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--retry-sleep-seconds", type=float, default=2.0)
    parser.add_argument("--audit-csv", type=Path, default=DEFAULT_NOT_FOUND_AUDIT_CSV)
    parser.add_argument("--returned-csv", type=Path, default=DEFAULT_RETURNED_CSV)
    parser.add_argument(
        "--continue-on-api-error",
        action="store_true",
        help="Write failed batches to the audit CSV instead of raising immediately.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional debug limit on number of unique CIDs queried.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.batch_size <= 0:
        raise ValueError("--batch-size must be positive")
    if args.batch_size > 250:
        raise ValueError("--batch-size must be <= 250 for AnnotationDB")

    endpoints = tuple(args.annotationdb_endpoint or DEFAULT_ANNOTATIONDB_ENDPOINTS)

    seed_drugs_from_annotationdb(
        batch_size=args.batch_size,
        annotationdb_base_url=args.annotationdb_base_url,
        annotationdb_endpoints=endpoints,
        timeout=args.timeout,
        max_retries=args.max_retries,
        retry_sleep_seconds=args.retry_sleep_seconds,
        audit_csv=args.audit_csv,
        returned_csv=args.returned_csv,
        continue_on_api_error=args.continue_on_api_error,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()