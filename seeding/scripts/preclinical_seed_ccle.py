from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd
from sqlalchemy import delete, select
from sqlalchemy.dialects.mysql import insert as mysql_insert
from sqlalchemy.orm import Session

from .seeding_coordinator_engine import alchemy_engine
from ..models.tables import (
    Base,
    PreClinicalCellLine,
    PreClinicalCopyNumberVariation,
    PreClinicalDataset,
    PreClinicalGene,
    PreClinicalMicroarray,
    PreClinicalRnaSeq,
    PreClinicalSample,
    PreClinicalTreatmentResponse,
)


DEFAULT_DATA_DIR = Path("extraction/data/proc/preclinical/CCLE")
DEFAULT_DATASET_NAME = "CCLE"
SAMPLE_ID_PREFIX = "CCLE_"
DEFAULT_CHUNK_SIZE = 100_000


REQUIRED_CELL_LINE_COLUMNS = {
    "cell_line_name",
    "tissueid",
    "mod_tissueid",
    "accession",
    "category",
    "sex",
    "age",
}

REQUIRED_SAMPLE_COLUMNS = {"cell_line_name"}

REQUIRED_TREATMENT_COLUMNS = {
    "cell_line_name",
    "treatment_id",
    "ic50_recomputed",
    "acc_recomputed",
    "mechanism_of_action",
}

REQUIRED_GENE_COLUMNS = {"id", "name"}
REQUIRED_EXPRESSION_COLUMNS = {"sample_id", "gene_id", "expression_value"}
REQUIRED_CNV_COLUMNS = {"sample_id", "gene_id", "value"}


def ccle_linear_cnv_to_log2(value: float | None) -> float | None:
    """Convert CCLE linear copy-number ratio to gene-level log2 copy-number value."""
    if value is None:
        return None

    # CCLE CNV values are linear and centered around ~1.0. The DB convention for
    # this project is log2 copy-number value, where ~0 is copy-neutral. Do not
    # floor positive values; invalid non-positive values are stored as NULL.
    if value <= 0:
        return None

    out = math.log2(value)
    if not math.isfinite(out):
        return None
    return out


MOLECULAR_LOAD_PLAN = (
    {
        "label": "RNA-seq",
        "filename": "pre_clinical_rna_seq.csv",
        "model": PreClinicalRnaSeq,
        "value_column": "expression_value",
        "required_columns": REQUIRED_EXPRESSION_COLUMNS,
        "value_transform": None,
    },
    {
        "label": "microarray",
        "filename": "pre_clinical_microarray.csv",
        "model": PreClinicalMicroarray,
        "value_column": "expression_value",
        "required_columns": REQUIRED_EXPRESSION_COLUMNS,
        "value_transform": None,
    },
    {
        "label": "copy-number variation",
        "filename": "pre_clinical_copy_number_variation.csv",
        "model": PreClinicalCopyNumberVariation,
        "value_column": "value",
        "required_columns": REQUIRED_CNV_COLUMNS,
        "value_transform": ccle_linear_cnv_to_log2,
    },
)


# Intentionally excluded for now:
#   pre_clinical_mutation.csv -> PreClinicalMutation


def clean_value(value: Any) -> Any | None:
    if pd.isna(value):
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

    # MySQL/PyMySQL cannot insert NaN, Inf, or -Inf into FLOAT/DOUBLE columns.
    if not math.isfinite(out):
        return None

    return out


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required CCLE CSV: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])


def iter_csv_chunks(path: Path, *, chunksize: int) -> Iterable[pd.DataFrame]:
    if not path.exists():
        raise FileNotFoundError(f"Missing required CCLE CSV: {path}")

    yield from pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
        na_values=[],
        chunksize=chunksize,
    )


def require_columns(df: pd.DataFrame, required: set[str], path: Path) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(missing)}. "
            f"Found columns: {', '.join(df.columns)}"
        )


def get_sample_id_column(sample_df: pd.DataFrame, sample_path: Path) -> str:
    # Current extraction scripts may write id; the database model uses sampleid.
    for candidate in ("sampleid", "id", "sample_id"):
        if candidate in sample_df.columns:
            return candidate

    raise ValueError(
        f"{sample_path} must contain one of: sampleid, id, sample_id. "
        f"Found columns: {', '.join(sample_df.columns)}"
    )


def ensure_prefixed_sample_ids(sample_ids: list[Any], *, auto_prefix: bool) -> list[str | None]:
    cleaned: list[str | None] = []
    unprefixed: list[str] = []

    for raw_sample_id in sample_ids:
        sample_id = clean_str(raw_sample_id)
        if sample_id is None:
            cleaned.append(None)
            continue

        if sample_id.startswith(SAMPLE_ID_PREFIX):
            cleaned.append(sample_id)
        elif auto_prefix:
            cleaned.append(f"{SAMPLE_ID_PREFIX}{sample_id}")
        else:
            unprefixed.append(sample_id)
            cleaned.append(sample_id)

    if unprefixed:
        preview = ", ".join(unprefixed[:10])
        raise ValueError(
            "CCLE sample IDs must already be prefixed before database insert. "
            f"Found unprefixed sample IDs, for example: {preview}. "
            "Fix the extractor output or rerun with --auto-prefix-samples."
        )

    return cleaned


def validate_final_tables_model() -> None:
    if not hasattr(PreClinicalSample, "sampleid"):
        raise RuntimeError(
            "tables.py must define PreClinicalSample.sampleid as the primary key."
        )
    if not hasattr(PreClinicalSample, "dataset_id"):
        raise RuntimeError("tables.py must define PreClinicalSample.dataset_id.")
    if not hasattr(PreClinicalSample, "cell_line_name"):
        raise RuntimeError("tables.py must define PreClinicalSample.cell_line_name.")
    if not hasattr(PreClinicalTreatmentResponse, "dataset_id"):
        raise RuntimeError("tables.py must define PreClinicalTreatmentResponse.dataset_id.")
    if not hasattr(PreClinicalTreatmentResponse, "cell_line_name"):
        raise RuntimeError("tables.py must define PreClinicalTreatmentResponse.cell_line_name.")
    for col in ("tissueid", "mod_tissueid"):
        if not hasattr(PreClinicalCellLine, col):
            raise RuntimeError(f"tables.py must define PreClinicalCellLine.{col}.")


def create_required_tables(engine) -> None:
    Base.metadata.create_all(
        bind=engine,
        tables=[
            PreClinicalDataset.__table__,
            PreClinicalCellLine.__table__,
            PreClinicalSample.__table__,
            PreClinicalTreatmentResponse.__table__,
            PreClinicalGene.__table__,
            PreClinicalRnaSeq.__table__,
            PreClinicalMicroarray.__table__,
            PreClinicalCopyNumberVariation.__table__,
        ],
    )


def chunked(items: list[Any], size: int) -> Iterable[list[Any]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def delete_existing_dataset(session: Session, dataset_name: str) -> None:
    dataset = session.scalar(
        select(PreClinicalDataset).where(PreClinicalDataset.name == dataset_name)
    )
    if dataset is None:
        return

    sample_ids = list(
        session.scalars(
            select(PreClinicalSample.sampleid).where(
                PreClinicalSample.dataset_id == dataset.id
            )
        )
    )

    cell_line_ids = list(
        session.scalars(
            select(PreClinicalCellLine.id).where(
                PreClinicalCellLine.dataset_id == dataset.id
            )
        )
    )

    if sample_ids:
        for sample_id_chunk in chunked(sample_ids, 10_000):
            session.execute(
                delete(PreClinicalRnaSeq).where(
                    PreClinicalRnaSeq.sample_id.in_(sample_id_chunk)
                )
            )
            session.execute(
                delete(PreClinicalMicroarray).where(
                    PreClinicalMicroarray.sample_id.in_(sample_id_chunk)
                )
            )
            session.execute(
                delete(PreClinicalCopyNumberVariation).where(
                    PreClinicalCopyNumberVariation.sample_id.in_(sample_id_chunk)
                )
            )

    session.execute(
        delete(PreClinicalTreatmentResponse).where(
            PreClinicalTreatmentResponse.dataset_id == dataset.id
        )
    )
    session.execute(
        delete(PreClinicalSample).where(PreClinicalSample.dataset_id == dataset.id)
    )

    if cell_line_ids:
        session.execute(
            delete(PreClinicalCellLine).where(PreClinicalCellLine.id.in_(cell_line_ids))
        )

    session.execute(delete(PreClinicalDataset).where(PreClinicalDataset.id == dataset.id))
    session.flush()


def get_or_create_dataset(session: Session, dataset_name: str) -> PreClinicalDataset:
    dataset = session.scalar(
        select(PreClinicalDataset).where(PreClinicalDataset.name == dataset_name)
    )
    if dataset is not None:
        return dataset

    dataset = PreClinicalDataset(name=dataset_name)
    session.add(dataset)
    session.flush()
    return dataset


def seed_cell_lines(session: Session, *, dataset_id: int, data_dir: Path) -> dict[str, int]:
    cell_line_path = data_dir / "pre_clinical_cell_line.csv"
    cell_line_df = read_csv(cell_line_path)
    require_columns(cell_line_df, REQUIRED_CELL_LINE_COLUMNS, cell_line_path)

    cell_line_df = cell_line_df.copy()
    cell_line_df["cell_line_name"] = cell_line_df["cell_line_name"].map(clean_str)
    cell_line_df = cell_line_df[cell_line_df["cell_line_name"].notna()]
    cell_line_df = cell_line_df.drop_duplicates(subset=["cell_line_name"], keep="first")

    if cell_line_df.empty:
        raise ValueError(f"No usable cell lines found in {cell_line_path}")

    rows = [
        PreClinicalCellLine(
            dataset_id=dataset_id,
            cell_line_name=clean_str(row.get("cell_line_name")),
            tissueid=clean_str(row.get("tissueid")),
            mod_tissueid=clean_str(row.get("mod_tissueid")),
            accession=clean_str(row.get("accession")),
            category=clean_str(row.get("category")),
            sex=clean_str(row.get("sex")),
            age=clean_int(row.get("age")),
        )
        for row in cell_line_df.to_dict(orient="records")
    ]

    session.add_all(rows)
    session.flush()

    lookup = dict(
        session.execute(
            select(PreClinicalCellLine.cell_line_name, PreClinicalCellLine.id).where(
                PreClinicalCellLine.dataset_id == dataset_id
            )
        ).all()
    )
    print(f"Seeded CCLE cell lines: {len(lookup)}")
    return lookup


def seed_samples(
    session: Session,
    *,
    dataset_id: int,
    data_dir: Path,
    cell_line_name_to_id: dict[str, int],
    auto_prefix_samples: bool,
) -> set[str]:
    sample_path = data_dir / "pre_clinical_sample.csv"
    sample_df = read_csv(sample_path)
    require_columns(sample_df, REQUIRED_SAMPLE_COLUMNS, sample_path)

    sample_id_col = get_sample_id_column(sample_df, sample_path)
    sample_df = sample_df.copy()
    sample_df["__sampleid"] = ensure_prefixed_sample_ids(
        sample_df[sample_id_col].tolist(),
        auto_prefix=auto_prefix_samples,
    )
    sample_df["cell_line_name"] = sample_df["cell_line_name"].map(clean_str)
    sample_df = sample_df[
        sample_df["__sampleid"].notna() & sample_df["cell_line_name"].notna()
    ]
    sample_df = sample_df.drop_duplicates(subset=["__sampleid"], keep="first")

    missing_cell_lines = sorted(
        set(sample_df["cell_line_name"]) - set(cell_line_name_to_id.keys())
    )
    if missing_cell_lines:
        preview = ", ".join(missing_cell_lines[:20])
        raise ValueError(
            "Some CCLE sample rows reference cell lines that were not loaded into "
            f"pre_clinical_cell_line for dataset_id={dataset_id}. Examples: {preview}"
        )

    rows = [
        PreClinicalSample(
            sampleid=clean_str(row.get("__sampleid")),
            dataset_id=dataset_id,
            cell_line_name=clean_str(row.get("cell_line_name")),
        )
        for row in sample_df.to_dict(orient="records")
    ]

    session.add_all(rows)
    session.flush()

    sample_ids = set(sample_df["__sampleid"])
    print(f"Seeded CCLE samples: {len(rows)}")
    return sample_ids


def seed_treatment_response(
    session: Session,
    *,
    dataset_id: int,
    data_dir: Path,
    cell_line_name_to_id: dict[str, int],
) -> None:
    treatment_path = data_dir / "pre_clinical_treatment_response.csv"
    treatment_df = read_csv(treatment_path)
    require_columns(treatment_df, REQUIRED_TREATMENT_COLUMNS, treatment_path)

    treatment_df = treatment_df.copy()
    treatment_df["cell_line_name"] = treatment_df["cell_line_name"].map(clean_str)
    treatment_df["treatment_id"] = treatment_df["treatment_id"].map(clean_str)
    treatment_df = treatment_df[
        treatment_df["cell_line_name"].notna() & treatment_df["treatment_id"].notna()
    ]

    missing_cell_lines = sorted(
        set(treatment_df["cell_line_name"]) - set(cell_line_name_to_id.keys())
    )
    if missing_cell_lines:
        preview = ", ".join(missing_cell_lines[:20])
        raise ValueError(
            "Some CCLE treatment-response rows reference cell lines that were not "
            f"loaded into pre_clinical_cell_line for dataset_id={dataset_id}. "
            f"Examples: {preview}"
        )

    treatment_df = treatment_df.drop_duplicates(
        subset=["cell_line_name", "treatment_id"],
        keep="first",
    )

    rows = [
        {
            "dataset_id": dataset_id,
            "cell_line_name": clean_str(row.get("cell_line_name")),
            "treatment_id": clean_str(row.get("treatment_id")),
            "ic50_recomputed": clean_float(row.get("ic50_recomputed")),
            "acc_recomputed": clean_float(row.get("acc_recomputed")),
            "mechanism_of_action": clean_str(row.get("mechanism_of_action")),
        }
        for row in treatment_df.to_dict(orient="records")
    ]

    if rows:
        session.bulk_insert_mappings(PreClinicalTreatmentResponse, rows)
        session.flush()

    print(f"Seeded CCLE treatment responses: {len(rows)}")


def seed_genes(session: Session, *, data_dir: Path) -> set[str]:
    gene_path = data_dir / "pre_clinical_gene.csv"
    gene_df = read_csv(gene_path)
    require_columns(gene_df, REQUIRED_GENE_COLUMNS, gene_path)

    gene_df = gene_df.copy()
    gene_df["id"] = gene_df["id"].map(clean_str)
    gene_df["name"] = gene_df["name"].map(clean_str)
    gene_df = gene_df[gene_df["id"].notna()]
    gene_df = gene_df.drop_duplicates(subset=["id"], keep="first")

    rows = [
        {
            "id": clean_str(row.get("id")),
            "name": clean_str(row.get("name")),
        }
        for row in gene_df.to_dict(orient="records")
    ]
    rows = [row for row in rows if row["id"] is not None]

    if rows:
        for row_chunk in chunked(rows, 5_000):
            stmt = mysql_insert(PreClinicalGene.__table__).values(row_chunk)
            stmt = stmt.on_duplicate_key_update(name=stmt.inserted.name)
            session.execute(stmt)
            session.flush()

    gene_ids = set(gene_df["id"])
    print(f"Seeded/updated CCLE genes: {len(gene_ids)}")
    return gene_ids


def get_sample_ids_for_dataset(session: Session, dataset_id: int) -> set[str]:
    return set(
        session.scalars(
            select(PreClinicalSample.sampleid).where(
                PreClinicalSample.dataset_id == dataset_id
            )
        )
    )


def get_gene_ids(session: Session) -> set[str]:
    return set(session.scalars(select(PreClinicalGene.id)))


def validate_molecular_samples_prefixed(sample_ids: set[str], *, label: str) -> None:
    unprefixed = sorted(s for s in sample_ids if not s.startswith(SAMPLE_ID_PREFIX))
    if unprefixed:
        preview = ", ".join(unprefixed[:10])
        raise ValueError(
            f"CCLE {label} contains unprefixed sample_id values. Examples: {preview}. "
            "Fix the extractor output before loading."
        )


def seed_molecular_file(
    session: Session,
    *,
    data_dir: Path,
    filename: str,
    label: str,
    model: type,
    value_column: str,
    required_columns: set[str],
    valid_sample_ids: set[str],
    valid_gene_ids: set[str],
    chunksize: int,
    value_transform: Callable[[float | None], float | None] | None = None,
) -> None:
    path = data_dir / filename
    total_insert_candidates = 0
    total_skipped_missing_value = 0

    for chunk_index, chunk_df in enumerate(iter_csv_chunks(path, chunksize=chunksize), start=1):
        require_columns(chunk_df, required_columns, path)

        chunk_df = chunk_df.copy()
        chunk_df["sample_id"] = chunk_df["sample_id"].map(clean_str)
        chunk_df["gene_id"] = chunk_df["gene_id"].map(clean_str)
        chunk_df["__value"] = chunk_df[value_column].map(clean_float)

        if value_transform is not None:
            chunk_df["__value"] = chunk_df["__value"].map(value_transform)

        chunk_df = chunk_df[chunk_df["sample_id"].notna() & chunk_df["gene_id"].notna()]

        chunk_sample_ids = set(chunk_df["sample_id"])
        validate_molecular_samples_prefixed(chunk_sample_ids, label=label)

        missing_sample_ids = sorted(chunk_sample_ids - valid_sample_ids)
        if missing_sample_ids:
            preview = ", ".join(missing_sample_ids[:20])
            raise ValueError(
                f"CCLE {label} has sample_id values missing from pre_clinical_sample. "
                f"Examples: {preview}"
            )

        chunk_gene_ids = set(chunk_df["gene_id"])
        missing_gene_ids = sorted(chunk_gene_ids - valid_gene_ids)
        if missing_gene_ids:
            preview = ", ".join(missing_gene_ids[:20])
            raise ValueError(
                f"CCLE {label} has gene_id values missing from pre_clinical_gene. "
                f"Examples: {preview}"
            )

        before_value_filter = len(chunk_df)
        chunk_df = chunk_df[chunk_df["__value"].notna()]
        total_skipped_missing_value += before_value_filter - len(chunk_df)

        chunk_df = chunk_df.drop_duplicates(subset=["sample_id", "gene_id"], keep="first")

        rows = [
            {
                "sample_id": row["sample_id"],
                "gene_id": row["gene_id"],
                value_column: float(row["__value"]),
            }
            for row in chunk_df.to_dict(orient="records")
        ]

        if rows:
            stmt = mysql_insert(model.__table__).prefix_with("IGNORE")
            session.execute(stmt, rows)
            session.flush()

        total_insert_candidates += len(rows)
        print(
            f"Seeded CCLE {label} chunk {chunk_index}: "
            f"{len(rows)} candidate rows"
        )

    print(
        f"Finished CCLE {label}: {total_insert_candidates} candidate rows. "
        f"Skipped non-finite/missing values: {total_skipped_missing_value}"
    )


def seed_ccle(
    *,
    data_dir: Path,
    dataset_name: str,
    replace: bool,
    auto_prefix_samples: bool,
    chunksize: int,
) -> None:
    validate_final_tables_model()
    engine = alchemy_engine()
    create_required_tables(engine)

    with Session(engine) as session:
        if replace:
            print(f"Replacing existing dataset rows for {dataset_name}")
            delete_existing_dataset(session, dataset_name)
            session.commit()

        dataset = get_or_create_dataset(session, dataset_name)
        session.flush()
        print(f"Using dataset_id={dataset.id} for dataset_name={dataset.name}")

        cell_line_name_to_id = seed_cell_lines(
            session,
            dataset_id=dataset.id,
            data_dir=data_dir,
        )
        seed_samples(
            session,
            dataset_id=dataset.id,
            data_dir=data_dir,
            cell_line_name_to_id=cell_line_name_to_id,
            auto_prefix_samples=auto_prefix_samples,
        )
        seed_treatment_response(
            session,
            dataset_id=dataset.id,
            data_dir=data_dir,
            cell_line_name_to_id=cell_line_name_to_id,
        )
        seed_genes(session, data_dir=data_dir)
        session.commit()

        valid_sample_ids = get_sample_ids_for_dataset(session, dataset.id)
        valid_gene_ids = get_gene_ids(session)

        for plan in MOLECULAR_LOAD_PLAN:
            seed_molecular_file(
                session,
                data_dir=data_dir,
                filename=plan["filename"],
                label=plan["label"],
                model=plan["model"],
                value_column=plan["value_column"],
                required_columns=plan["required_columns"],
                valid_sample_ids=valid_sample_ids,
                valid_gene_ids=valid_gene_ids,
                chunksize=chunksize,
                value_transform=plan["value_transform"],
            )
            session.commit()

    print("Finished CCLE preclinical seeding. Mutation upload was skipped.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Seed the CCLE preclinical dataset from extracted CSVs. Mutation upload is intentionally skipped."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=f"Directory containing CCLE extracted CSVs. Default: {DEFAULT_DATA_DIR}",
    )
    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help=f"Dataset name to create/use in pre_clinical_dataset. Default: {DEFAULT_DATASET_NAME}",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Delete existing CCLE dataset rows before reloading.",
    )
    parser.add_argument(
        "--auto-prefix-samples",
        action="store_true",
        help=(
            "Automatically add CCLE_ to sample IDs if the sample CSV is not already prefixed. "
            "By default, the loader fails on unprefixed sample IDs."
        ),
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help=f"Rows per molecular CSV chunk. Default: {DEFAULT_CHUNK_SIZE}",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    seed_ccle(
        data_dir=args.data_dir,
        dataset_name=args.dataset_name,
        replace=args.replace,
        auto_prefix_samples=args.auto_prefix_samples,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()