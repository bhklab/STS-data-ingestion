from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import pandas as pd
from sqlalchemy import delete, select
from sqlalchemy.orm import Session

from .seeding_coordinator_engine import alchemy_engine
from ..models.tables import (
    Base,
    PreClinicalCellLine,
    PreClinicalDataset,
    PreClinicalSample,
    PreClinicalTreatmentResponse,
)


DEFAULT_DATA_DIR = Path("extraction/data/proc/preclinical/CTRPv2")
DEFAULT_DATASET_NAME = "CTRPv2"
SAMPLE_ID_PREFIX = "CTRP_"


REQUIRED_CELL_LINE_COLUMNS = {
    "cell_line_name",
    "tissueid",
    "mod_tissueid",
    "accession",
    "category",
    "sex",
    "age",
}

REQUIRED_SAMPLE_COLUMNS = {
    "cell_line_name",
}

REQUIRED_TREATMENT_COLUMNS = {
    "cell_line_name",
    "treatment_id",
    "ic50_recomputed",
    "acc_recomputed",
    "mechanism_of_action",
}


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
    # Treat non-finite assay values as missing so they are stored as NULL.
    if not math.isfinite(out):
        return None

    return out


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required CTRPv2 CSV: {path}")

    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])


def require_columns(df: pd.DataFrame, required: set[str], path: Path) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(missing)}. "
            f"Found columns: {', '.join(df.columns)}"
        )


def get_sample_id_column(sample_df: pd.DataFrame, sample_path: Path) -> str:
    # Newer desired extraction output can use sampleid.
    # Existing extractors may still write id; this loader maps it to DB sampleid.
    for candidate in ("sampleid", "id", "sample_id"):
        if candidate in sample_df.columns:
            return candidate

    raise ValueError(
        f"{sample_path} must contain one of: sampleid, id, sample_id. "
        f"Found columns: {', '.join(sample_df.columns)}"
    )


def ensure_prefixed_sample_ids(sample_ids: list[str], *, auto_prefix: bool) -> list[str]:
    cleaned: list[str] = []
    unprefixed: list[str] = []

    for sample_id in sample_ids:
        sample_id = clean_str(sample_id)
        if sample_id is None:
            cleaned.append(None)  # type: ignore[arg-type]
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
            "CTRPv2 sample IDs must already be prefixed before database insert. "
            f"Found unprefixed sample IDs, for example: {preview}. "
            "Fix the extractor output or rerun with --auto-prefix-samples."
        )

    return cleaned


def validate_final_tables_model() -> None:
    if not hasattr(PreClinicalSample, "sampleid"):
        raise RuntimeError(
            "tables.py must define PreClinicalSample.sampleid as the primary key. "
            "Your current model still appears to use PreClinicalSample.id."
        )

    if not hasattr(PreClinicalSample, "dataset_id"):
        raise RuntimeError("tables.py must define PreClinicalSample.dataset_id.")

    if not hasattr(PreClinicalSample, "cell_line_name"):
        raise RuntimeError(
            "tables.py must define PreClinicalSample.cell_line_name. "
            "This loader expects samples to reference cell lines by "
            "(cell_line_name, dataset_id), not cell_line_id."
        )

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
        ],
    )


def delete_existing_dataset(session: Session, dataset_name: str) -> None:
    dataset = session.scalar(
        select(PreClinicalDataset).where(PreClinicalDataset.name == dataset_name)
    )

    if dataset is None:
        return

    cell_line_ids = list(
        session.scalars(
            select(PreClinicalCellLine.id).where(
                PreClinicalCellLine.dataset_id == dataset.id
            )
        )
    )

    session.execute(
        delete(PreClinicalTreatmentResponse).where(
            PreClinicalTreatmentResponse.dataset_id == dataset.id
        )
    )

    session.execute(
        delete(PreClinicalSample).where(
            PreClinicalSample.dataset_id == dataset.id
        )
    )

    if cell_line_ids:
        session.execute(
            delete(PreClinicalCellLine).where(
                PreClinicalCellLine.id.in_(cell_line_ids)
            )
        )

    session.execute(
        delete(PreClinicalDataset).where(PreClinicalDataset.id == dataset.id)
    )
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


def seed_cell_lines(
    session: Session,
    *,
    dataset_id: int,
    data_dir: Path,
) -> dict[str, int]:
    cell_line_path = data_dir / "pre_clinical_cell_line.csv"
    cell_line_df = read_csv(cell_line_path)
    require_columns(cell_line_df, REQUIRED_CELL_LINE_COLUMNS, cell_line_path)

    cell_line_df = cell_line_df.copy()
    cell_line_df["cell_line_name"] = cell_line_df["cell_line_name"].map(clean_str)
    cell_line_df = cell_line_df[cell_line_df["cell_line_name"].notna()]
    cell_line_df = cell_line_df.drop_duplicates(subset=["cell_line_name"], keep="first")

    if cell_line_df.empty:
        raise ValueError(f"No usable cell lines found in {cell_line_path}")

    rows = []
    for row in cell_line_df.to_dict(orient="records"):
        rows.append(
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
        )

    session.add_all(rows)
    session.flush()

    lookup = dict(
        session.execute(
            select(PreClinicalCellLine.cell_line_name, PreClinicalCellLine.id).where(
                PreClinicalCellLine.dataset_id == dataset_id
            )
        ).all()
    )

    print(f"Seeded CTRPv2 cell lines: {len(lookup)}")
    return lookup


def seed_samples(
    session: Session,
    *,
    dataset_id: int,
    data_dir: Path,
    cell_line_name_to_id: dict[str, int],
    auto_prefix_samples: bool,
) -> None:
    sample_path = data_dir / "pre_clinical_sample.csv"
    sample_df = read_csv(sample_path)
    require_columns(sample_df, REQUIRED_SAMPLE_COLUMNS, sample_path)

    sample_id_col = get_sample_id_column(sample_df, sample_path)

    sample_df = sample_df.copy()
    sample_df["__sampleid"] = ensure_prefixed_sample_ids(
        [clean_str(x) for x in sample_df[sample_id_col].tolist()],
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
            "Some CTRPv2 sample rows reference cell lines that were not loaded into "
            f"pre_clinical_cell_line for dataset_id={dataset_id}. Examples: {preview}"
        )

    rows = []
    for row in sample_df.to_dict(orient="records"):
        rows.append(
            PreClinicalSample(
                sampleid=clean_str(row.get("__sampleid")),
                dataset_id=dataset_id,
                cell_line_name=clean_str(row.get("cell_line_name")),
            )
        )

    session.add_all(rows)
    session.flush()

    print(f"Seeded CTRPv2 samples: {len(rows)}")


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
            "Some CTRPv2 treatment-response rows reference cell lines that were not "
            f"loaded into pre_clinical_cell_line for dataset_id={dataset_id}. "
            f"Examples: {preview}"
        )

    # The extractor should already summarize duplicates, but this protects the DB
    # unique constraint and keeps the first non-empty row ordering from the CSV.
    treatment_df = treatment_df.drop_duplicates(
        subset=["cell_line_name", "treatment_id"],
        keep="first",
    )

    rows = []
    for row in treatment_df.to_dict(orient="records"):
        rows.append(
            PreClinicalTreatmentResponse(
                dataset_id=dataset_id,
                cell_line_name=clean_str(row.get("cell_line_name")),
                treatment_id=clean_str(row.get("treatment_id")),
                ic50_recomputed=clean_float(row.get("ic50_recomputed")),
                acc_recomputed=clean_float(row.get("acc_recomputed")),
                mechanism_of_action=clean_str(row.get("mechanism_of_action")),
            )
        )

    session.add_all(rows)
    session.flush()

    print(f"Seeded CTRPv2 treatment responses: {len(rows)}")


def seed_ctrpv2(
    *,
    data_dir: Path,
    dataset_name: str,
    replace: bool,
    auto_prefix_samples: bool,
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

        session.commit()

    print("Finished CTRPv2 preclinical seeding.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Seed the CTRPv2 preclinical dataset from extracted CSVs."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=f"Directory containing CTRPv2 extracted CSVs. Default: {DEFAULT_DATA_DIR}",
    )
    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help=f"Dataset name to create/use in pre_clinical_dataset. Default: {DEFAULT_DATASET_NAME}",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Delete existing CTRPv2 dataset rows before reloading.",
    )
    parser.add_argument(
        "--auto-prefix-samples",
        action="store_true",
        help=(
            "Automatically add CTRPv2_ to sample IDs if the sample CSV is not already prefixed. "
            "By default, the loader fails on unprefixed sample IDs."
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    seed_ctrpv2(
        data_dir=args.data_dir,
        dataset_name=args.dataset_name,
        replace=args.replace,
        auto_prefix_samples=args.auto_prefix_samples,
    )


if __name__ == "__main__":
    main()