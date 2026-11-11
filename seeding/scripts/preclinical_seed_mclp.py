from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from sqlalchemy import func, select
from sqlalchemy.dialects.mysql import insert as mysql_insert
from sqlalchemy.orm import Session

from .seeding_coordinator_engine import alchemy_engine
from ..models.tables import (
    Base,
    PreClinicalCellLine,
    PreClinicalDataset,
    PreClinicalGene,
    PreClinicalRppa,
    PreClinicalSample,
)


DEFAULT_DATA_DIR = Path("extraction/data/proc/preclinical/MCLP")
DEFAULT_DATASET_NAME = "MCLP"
SAMPLE_ID_PREFIX = "MCLP_"
DEFAULT_CHUNK_SIZE = 100_000


# ---------------------------------------------------------------------------
# Required CSV columns
# ---------------------------------------------------------------------------

REQUIRED_CELL_LINE_COLUMNS = {
    "cell_line_name",
}

REQUIRED_SAMPLE_COLUMNS = {
    "cell_line_name",
}

REQUIRED_RPPA_COLUMNS = {
    "feature_id",
    "sample_id",
    "gene_symbol_primary",
    "gene_id",
    "value",
}


# ---------------------------------------------------------------------------
# Cleaning helpers
# ---------------------------------------------------------------------------

def clean_value(value: Any) -> Any | None:
    """
    Convert missing / NA-like values to None and strip whitespace
    from string values.
    """
    if pd.isna(value):
        return None

    if isinstance(value, str):
        value = value.strip()

        if (
            value == ""
            or value.upper()
            in {
                "NA",
                "N/A",
                "NS",
                "NAN",
                "NONE",
                "NULL",
            }
        ):
            return None

        return value

    return value


def clean_str(value: Any) -> str | None:
    """Return a cleaned string or None."""
    value = clean_value(value)

    if value is None:
        return None

    return str(value)


def clean_int(value: Any) -> int | None:
    """
    Convert a value to integer when possible.

    Example:
        46 -> 46
        "46" -> 46
    """
    value = clean_value(value)

    if value is None:
        return None

    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def clean_float(value: Any) -> float | None:
    """
    Convert a value to a finite float.

    NaN and Inf are converted to None so they are never sent to MySQL.
    """
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


def clean_gene_id(value: Any) -> str | None:
    """
    Clean MCLP gene IDs.

    Steps:
      1. Remove whitespace.
      2. If multiple IDs are separated by ///, keep the first.
      3. Remove the Ensembl version suffix after the first period.

    Examples:

        ENSG00000160710.18
            -> ENSG00000160710

        ENSG00000142208.18///ENSG00000105221.18
            -> ENSG00000142208
    """
    gene_id = clean_str(value)

    if gene_id is None:
        return None

    # Keep first ID when multiple IDs exist.
    gene_id = gene_id.split("///", 1)[0].strip()

    # Remove Ensembl version suffix.
    gene_id = gene_id.split(".", 1)[0].strip()

    return gene_id or None


def clean_gene_symbol(value: Any) -> str | None:
    """
    Clean gene symbols.

    If multiple gene symbols are separated by ///, keep only the first.

    Example:

        AKT1///AKT2///AKT3
            -> AKT1
    """
    symbol = clean_str(value)

    if symbol is None:
        return None

    symbol = symbol.split("///", 1)[0].strip()

    return symbol or None


def normalize_sample_id(value: Any) -> str | None:
    """
    Ensure all MCLP samples have an MCLP_ prefix.

    Examples:

        SW872
            -> MCLP_SW872

        MCLP_SW872
            -> MCLP_SW872
    """
    sample_id = clean_str(value)

    if sample_id is None:
        return None

    if sample_id.startswith(SAMPLE_ID_PREFIX):
        return sample_id

    return f"{SAMPLE_ID_PREFIX}{sample_id}"


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def read_csv(path: Path) -> pd.DataFrame:
    """Read one of the smaller MCLP CSV files."""
    if not path.exists():
        raise FileNotFoundError(
            f"Missing required MCLP CSV: {path}"
        )

    return pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
        na_values=[],
    )


def iter_csv_chunks(
    path: Path,
    *,
    chunksize: int,
) -> Iterable[pd.DataFrame]:
    """
    Read a large CSV in chunks.

    Used for RPPA so the complete file does not need to be held
    in memory at once.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Missing required MCLP CSV: {path}"
        )

    yield from pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
        na_values=[],
        chunksize=chunksize,
    )


def require_columns(
    df: pd.DataFrame,
    required: set[str],
    path: Path,
) -> None:
    """Verify that an extracted CSV has the expected columns."""
    missing = sorted(
        required - set(df.columns)
    )

    if missing:
        raise ValueError(
            f"{path} is missing required columns: "
            f"{', '.join(missing)}. "
            f"Found columns: {', '.join(df.columns)}"
        )


def get_sample_id_column(
    sample_df: pd.DataFrame,
    sample_path: Path,
) -> str:
    """
    Allow the sample ID column to use one of the standard names
    used by the extraction scripts.
    """
    for candidate in (
        "id",
        "sample_id",
        "sampleid",
    ):
        if candidate in sample_df.columns:
            return candidate

    raise ValueError(
        f"{sample_path} must contain one of: "
        "id, sample_id, sampleid. "
        f"Found columns: {', '.join(sample_df.columns)}"
    )


def chunked(
    items: list[Any],
    size: int,
) -> Iterable[list[Any]]:
    """Break a Python list into smaller insertion batches."""
    for start in range(
        0,
        len(items),
        size,
    ):
        yield items[start:start + size]


# ---------------------------------------------------------------------------
# Database setup
# ---------------------------------------------------------------------------

def create_required_tables(engine) -> None:
    """
    Create only the tables needed for MCLP if they do not already exist.
    """
    Base.metadata.create_all(
        bind=engine,
        tables=[
            PreClinicalDataset.__table__,
            PreClinicalCellLine.__table__,
            PreClinicalSample.__table__,
            PreClinicalGene.__table__,
            PreClinicalRppa.__table__,
        ],
    )


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

def get_or_create_dataset(
    session: Session,
    *,
    dataset_name: str,
) -> PreClinicalDataset:
    """
    Find the existing MCLP dataset row.

    If it does not exist, create a minimal dataset entry containing
    the dataset name.
    """
    dataset = session.scalar(
        select(
            PreClinicalDataset
        ).where(
            func.lower(
                PreClinicalDataset.name
            )
            == dataset_name.lower()
        )
    )

    if dataset is not None:
        return dataset

    stmt = (
        mysql_insert(
            PreClinicalDataset.__table__
        )
        .values(
            name=dataset_name
        )
        .prefix_with("IGNORE")
    )

    session.execute(stmt)
    session.flush()

    dataset = session.scalar(
        select(
            PreClinicalDataset
        ).where(
            func.lower(
                PreClinicalDataset.name
            )
            == dataset_name.lower()
        )
    )

    if dataset is None:
        raise RuntimeError(
            f"Could not create or find dataset "
            f"{dataset_name!r}."
        )

    return dataset


# ---------------------------------------------------------------------------
# Cell lines
# ---------------------------------------------------------------------------

def seed_cell_lines(
    session: Session,
    *,
    dataset_id: int,
    data_dir: Path,
) -> set[str]:
    """
    Insert MCLP cell lines.

    Existing rows are skipped by INSERT IGNORE.

    The database's logical uniqueness constraint is:

        (cell_line_name, dataset_id)
    """
    path = (
        data_dir
        / "pre_clinical_cell_line.csv"
    )

    df = read_csv(path)

    require_columns(
        df,
        REQUIRED_CELL_LINE_COLUMNS,
        path,
    )

    df = df.copy()

    df["cell_line_name"] = (
        df["cell_line_name"]
        .map(clean_str)
    )

    df = df[
        df["cell_line_name"].notna()
    ]

    df = df.drop_duplicates(
        subset=["cell_line_name"],
        keep="first",
    )

    rows = [
        {
            "dataset_id": dataset_id,
            "cell_line_name": clean_str(
                row.get("cell_line_name")
            ),
            "tissueid": clean_str(
                row.get("tissueid")
            ),
            "mod_tissueid": clean_str(
                row.get("mod_tissueid")
            ),
            "accession": clean_str(
                row.get("accession")
            ),
            "category": clean_str(
                row.get("category")
            ),
            "sex": clean_str(
                row.get("sex")
            ),
            "age": clean_int(
                row.get("age")
            ),
        }
        for row
        in df.to_dict(
            orient="records"
        )
    ]

    inserted = 0

    if rows:
        for row_chunk in chunked(
            rows,
            5_000,
        ):
            stmt = (
                mysql_insert(
                    PreClinicalCellLine.__table__
                )
                .prefix_with("IGNORE")
            )

            result = session.execute(
                stmt,
                row_chunk,
            )

            if (
                result.rowcount
                and result.rowcount > 0
            ):
                inserted += result.rowcount

    session.flush()

    valid_cell_lines = set(
        session.scalars(
            select(
                PreClinicalCellLine.cell_line_name
            ).where(
                PreClinicalCellLine.dataset_id
                == dataset_id
            )
        )
    )

    print(
        f"MCLP cell lines: "
        f"{len(rows)} candidates, "
        f"{inserted} inserted, "
        f"{len(rows) - inserted} "
        "already present/skipped."
    )

    return valid_cell_lines


# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

def seed_samples(
    session: Session,
    *,
    dataset_id: int,
    data_dir: Path,
    valid_cell_lines: set[str],
) -> set[str]:
    """
    Insert MCLP samples.

    Sample IDs are prefixed with MCLP_ because
    pre_clinical_sample.id is a global primary key.

    Existing primary keys are skipped.

    If an existing sample ID belongs to a different dataset,
    it is excluded so RPPA cannot accidentally attach to that
    dataset's sample.
    """
    path = (
        data_dir
        / "pre_clinical_sample.csv"
    )

    df = read_csv(path)

    require_columns(
        df,
        REQUIRED_SAMPLE_COLUMNS,
        path,
    )

    sample_id_col = (
        get_sample_id_column(
            df,
            path,
        )
    )

    df = df.copy()

    df["__id"] = (
        df[sample_id_col]
        .map(normalize_sample_id)
    )

    df["cell_line_name"] = (
        df["cell_line_name"]
        .map(clean_str)
    )

    df = df[
        df["__id"].notna()
        & df["cell_line_name"].notna()
    ]

    df = df.drop_duplicates(
        subset=["__id"],
        keep="first",
    )

    # Check that every sample references an MCLP cell line.
    missing_cell_lines = sorted(
        set(
            df["cell_line_name"]
        )
        - valid_cell_lines
    )

    if missing_cell_lines:
        preview = ", ".join(
            missing_cell_lines[:20]
        )

        raise ValueError(
            "Some MCLP samples reference cell lines "
            "that are not present in "
            "pre_clinical_cell_line for "
            f"dataset_id={dataset_id}. "
            f"Examples: {preview}"
        )

    candidate_ids = (
        df["__id"]
        .tolist()
    )

    # Determine whether any sample ID already exists
    # under another dataset.
    existing_dataset_by_sample: dict[
        str,
        int,
    ] = {}

    for id_chunk in chunked(
        candidate_ids,
        5_000,
    ):
        existing_dataset_by_sample.update(
            dict(
                session.execute(
                    select(
                        PreClinicalSample.id,
                        PreClinicalSample.dataset_id,
                    ).where(
                        PreClinicalSample.id.in_(
                            id_chunk
                        )
                    )
                ).all()
            )
        )

    cross_dataset_conflicts = {
        sample_id
        for (
            sample_id,
            existing_dataset_id,
        )
        in existing_dataset_by_sample.items()
        if existing_dataset_id
        != dataset_id
    }

    if cross_dataset_conflicts:
        preview = ", ".join(
            sorted(
                cross_dataset_conflicts
            )[:20]
        )

        print(
            "Warning: skipping "
            f"{len(cross_dataset_conflicts)} "
            "sample IDs that already belong "
            "to another dataset. "
            f"Examples: {preview}"
        )

        df = df[
            ~df["__id"].isin(
                cross_dataset_conflicts
            )
        ]

    rows = [
        {
            "id": row["__id"],
            "dataset_id": dataset_id,
            "cell_line_name": (
                row["cell_line_name"]
            ),
        }
        for row
        in df.to_dict(
            orient="records"
        )
    ]

    inserted = 0

    if rows:
        for row_chunk in chunked(
            rows,
            5_000,
        ):
            stmt = (
                mysql_insert(
                    PreClinicalSample.__table__
                )
                .prefix_with("IGNORE")
            )

            result = session.execute(
                stmt,
                row_chunk,
            )

            if (
                result.rowcount
                and result.rowcount > 0
            ):
                inserted += result.rowcount

    session.flush()

    valid_sample_ids = set(
        session.scalars(
            select(
                PreClinicalSample.id
            ).where(
                PreClinicalSample.dataset_id
                == dataset_id
            )
        )
    )

    print(
        f"MCLP samples: "
        f"{len(rows)} candidates after "
        "conflict filtering, "
        f"{inserted} inserted, "
        f"{len(rows) - inserted} "
        "already present/skipped."
    )

    return valid_sample_ids


# ---------------------------------------------------------------------------
# Genes
# ---------------------------------------------------------------------------

def seed_genes_from_rppa(
    session: Session,
    *,
    data_dir: Path,
    chunksize: int,
) -> set[str]:
    """
    Derive genes directly from pre_clinical_rppa.csv.

    Mapping:

        gene_id
            -> pre_clinical_gene.id

        gene_symbol_primary
            -> pre_clinical_gene.name

    Existing gene IDs are skipped with INSERT IGNORE.

    Example RPPA source:

        gene_symbol_primary:
            AKT1

        gene_id:
            ENSG00000142208

    becomes:

        pre_clinical_gene:
            id   = ENSG00000142208
            name = AKT1
    """
    path = (
        data_dir
        / "pre_clinical_rppa.csv"
    )

    total_candidates = 0
    total_inserted = 0

    # Prevent repeated insert attempts when the same gene occurs
    # in many RPPA chunks.
    seen_gene_ids: set[str] = set()

    for chunk_index, df in enumerate(
        iter_csv_chunks(
            path,
            chunksize=chunksize,
        ),
        start=1,
    ):
        require_columns(
            df,
            REQUIRED_RPPA_COLUMNS,
            path,
        )

        df = df.copy()

        df["gene_id"] = (
            df["gene_id"]
            .map(clean_gene_id)
        )

        df["gene_symbol_primary"] = (
            df["gene_symbol_primary"]
            .map(clean_gene_symbol)
        )

        # A gene_id is required for pre_clinical_gene.
        df = df[
            df["gene_id"].notna()
        ]

        # One row per gene within this chunk.
        df = df.drop_duplicates(
            subset=["gene_id"],
            keep="first",
        )

        # Do not repeatedly try the same gene in later chunks.
        df = df[
            ~df["gene_id"].isin(
                seen_gene_ids
            )
        ]

        rows = [
            {
                "id": row["gene_id"],
                "name": row[
                    "gene_symbol_primary"
                ],
            }
            for row
            in df.to_dict(
                orient="records"
            )
        ]

        inserted = 0

        if rows:
            for row_chunk in chunked(
                rows,
                5_000,
            ):
                stmt = (
                    mysql_insert(
                        PreClinicalGene.__table__
                    )
                    .prefix_with("IGNORE")
                )

                result = session.execute(
                    stmt,
                    row_chunk,
                )

                if (
                    result.rowcount
                    and result.rowcount > 0
                ):
                    inserted += result.rowcount

            seen_gene_ids.update(
                row["id"]
                for row in rows
            )

        total_candidates += len(rows)
        total_inserted += inserted

        print(
            f"MCLP genes chunk "
            f"{chunk_index}: "
            f"{len(rows)} candidate genes, "
            f"{inserted} inserted, "
            f"{len(rows) - inserted} "
            "already present/skipped."
        )

    session.flush()

    print(
        "Finished MCLP gene seeding: "
        f"{total_candidates} unique "
        "candidate genes, "
        f"{total_inserted} inserted, "
        f"{total_candidates - total_inserted} "
        "already present/skipped."
    )

    return seen_gene_ids


# ---------------------------------------------------------------------------
# RPPA
# ---------------------------------------------------------------------------

def seed_rppa(
    session: Session,
    *,
    data_dir: Path,
    valid_sample_ids: set[str],
    chunksize: int,
) -> None:
    """
    Seed MCLP RPPA data.

    RPPA is loaded in chunks because the extracted assay can be large.

    Database uniqueness is:

        (sample_id, feature_id)

    Therefore rerunning this script safely skips RPPA rows already
    present in the database.
    """
    path = (
        data_dir
        / "pre_clinical_rppa.csv"
    )

    total_candidates = 0
    total_inserted = 0
    total_skipped_invalid = 0
    total_skipped_missing_sample = 0

    for chunk_index, df in enumerate(
        iter_csv_chunks(
            path,
            chunksize=chunksize,
        ),
        start=1,
    ):
        require_columns(
            df,
            REQUIRED_RPPA_COLUMNS,
            path,
        )

        df = df.copy()

        # Clean all RPPA fields.
        df["sample_id"] = (
            df["sample_id"]
            .map(normalize_sample_id)
        )

        df["feature_id"] = (
            df["feature_id"]
            .map(clean_str)
        )

        df["gene_symbol_primary"] = (
            df["gene_symbol_primary"]
            .map(clean_gene_symbol)
        )

        df["gene_id"] = (
            df["gene_id"]
            .map(clean_gene_id)
        )

        df["__value"] = (
            df["value"]
            .map(clean_float)
        )

        # sample_id, feature_id and numeric value are required.
        # gene_id / gene_symbol_primary are allowed to be NULL by
        # the RPPA database model.
        before_required_filter = len(df)

        df = df[
            df["sample_id"].notna()
            & df["feature_id"].notna()
            & df["__value"].notna()
        ]

        total_skipped_invalid += (
            before_required_filter
            - len(df)
        )

        # Only upload RPPA rows whose sample belongs to MCLP.
        before_sample_filter = len(df)

        df = df[
            df["sample_id"].isin(
                valid_sample_ids
            )
        ]

        total_skipped_missing_sample += (
            before_sample_filter
            - len(df)
        )

        # Match the database unique constraint.
        df = df.drop_duplicates(
            subset=[
                "sample_id",
                "feature_id",
            ],
            keep="first",
        )

        rows = [
            {
                "sample_id": (
                    row["sample_id"]
                ),
                "feature_id": (
                    row["feature_id"]
                ),
                "gene_symbol_primary": (
                    row[
                        "gene_symbol_primary"
                    ]
                ),
                "gene_id": (
                    row["gene_id"]
                ),
                "value": float(
                    row["__value"]
                ),
            }
            for row
            in df.to_dict(
                orient="records"
            )
        ]

        inserted = 0

        if rows:
            stmt = (
                mysql_insert(
                    PreClinicalRppa.__table__
                )
                .prefix_with("IGNORE")
            )

            result = session.execute(
                stmt,
                rows,
            )

            if (
                result.rowcount
                and result.rowcount > 0
            ):
                inserted = result.rowcount

        # Commit each large RPPA chunk separately.
        session.commit()

        total_candidates += len(rows)
        total_inserted += inserted

        print(
            f"MCLP RPPA chunk "
            f"{chunk_index}: "
            f"{len(rows)} candidates, "
            f"{inserted} inserted, "
            f"{len(rows) - inserted} "
            "already present/skipped."
        )

    print(
        "Finished MCLP RPPA seeding: "
        f"{total_candidates} candidate rows, "
        f"{total_inserted} inserted, "
        f"{total_candidates - total_inserted} "
        "existing duplicates skipped. "
        f"Invalid/missing required values skipped: "
        f"{total_skipped_invalid}. "
        "Rows whose sample was not valid "
        f"for MCLP skipped: "
        f"{total_skipped_missing_sample}."
    )


# ---------------------------------------------------------------------------
# Main seeding workflow
# ---------------------------------------------------------------------------

def seed_mclp(
    *,
    data_dir: Path,
    dataset_name: str,
    chunksize: int,
) -> None:
    """
    Seed the MCLP preclinical dataset.

    Order matters because of foreign-key relationships:

        dataset
          ↓
        cell lines
          ↓
        samples
          ↓
        genes
          ↓
        RPPA
    """
    engine = alchemy_engine()

    create_required_tables(
        engine
    )

    with Session(engine) as session:

        # ---------------------------------------------------------------
        # Dataset
        # ---------------------------------------------------------------

        dataset = get_or_create_dataset(
            session,
            dataset_name=dataset_name,
        )

        session.commit()

        print(
            f"Using dataset_id="
            f"{dataset.id} "
            f"for dataset_name="
            f"{dataset.name}"
        )

        # ---------------------------------------------------------------
        # Cell lines
        # ---------------------------------------------------------------

        valid_cell_lines = seed_cell_lines(
            session,
            dataset_id=dataset.id,
            data_dir=data_dir,
        )

        session.commit()

        # ---------------------------------------------------------------
        # Samples
        # ---------------------------------------------------------------

        valid_sample_ids = seed_samples(
            session,
            dataset_id=dataset.id,
            data_dir=data_dir,
            valid_cell_lines=valid_cell_lines,
        )

        session.commit()

        # ---------------------------------------------------------------
        # Genes derived from RPPA
        # ---------------------------------------------------------------

        seed_genes_from_rppa(
            session,
            data_dir=data_dir,
            chunksize=chunksize,
        )

        session.commit()

        # ---------------------------------------------------------------
        # RPPA
        # ---------------------------------------------------------------

        seed_rppa(
            session,
            data_dir=data_dir,
            valid_sample_ids=valid_sample_ids,
            chunksize=chunksize,
        )

    print(
        "Finished MCLP preclinical seeding: "
        "cell lines, samples, genes, and RPPA."
    )


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Seed MCLP cell-line, sample, gene, "
            "and RPPA preclinical data."
        )
    )

    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=(
            "Directory containing "
            "pre_clinical_cell_line.csv, "
            "pre_clinical_sample.csv, and "
            "pre_clinical_rppa.csv. "
            f"Default: {DEFAULT_DATA_DIR}"
        ),
    )

    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help=(
            "Dataset name in datasets table. "
            f"Default: {DEFAULT_DATASET_NAME}"
        ),
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help=(
            "RPPA CSV rows processed per chunk. "
            f"Default: {DEFAULT_CHUNK_SIZE}"
        ),
    )

    return parser


def main() -> None:
    args = (
        build_parser()
        .parse_args()
    )

    seed_mclp(
        data_dir=args.data_dir,
        dataset_name=args.dataset_name,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()