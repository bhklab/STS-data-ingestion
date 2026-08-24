from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd
from sqlalchemy import delete, func, select
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


DEFAULT_DATA_DIR = Path("extraction/data/proc/preclinical/PRISM")
DEFAULT_DATASET_NAME = "PRISM"
SAMPLE_ID_PREFIX = "PRISM_"
DEFAULT_DATASET_METADATA_CSV = Path("extraction/data/raw/preclinical/combined_datasets.csv")
DEFAULT_CHUNK_SIZE = 100_000
LOAD_RNA_SEQ = False
LOAD_MICROARRAY = False
LOAD_CNV = False
RNA_TRANSFORM = "none"
CNV_TRANSFORM = "none"


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
    "cid",
    "ic50_recomputed",
    "acc_recomputed",
    "mechanism_of_action",
}

REQUIRED_GENE_COLUMNS = {"id", "name"}
REQUIRED_MOLECULAR_COLUMNS = {"sample_id", "gene_id", "value"}


def log2_tpm_plus_pseudocount_to_tpm(value: float | None) -> float | None:
    """Convert log2(TPM + 0.001) values back to TPM."""
    if value is None:
        return None

    out = max((2 ** value) - 0.001, 0.0)
    if not math.isfinite(out):
        return None
    return out


def linear_cnv_to_log2(value: float | None) -> float | None:
    """Convert positive linear CNV/copy-ratio values to log2 values."""
    if value is None:
        return None
    if value <= 0:
        return None

    out = math.log2(value)
    if not math.isfinite(out):
        return None
    return out


def get_value_transform(transform_name: str) -> Callable[[float | None], float | None] | None:
    if transform_name == "log2_tpm_plus_pseudocount_to_tpm":
        return log2_tpm_plus_pseudocount_to_tpm
    if transform_name == "linear_cnv_to_log2":
        return linear_cnv_to_log2
    if transform_name in {"", "none", "None", "null"}:
        return None
    raise ValueError(f"Unknown transform name: {transform_name}")


def build_molecular_load_plan() -> tuple[dict[str, Any], ...]:
    plan: list[dict[str, Any]] = []

    if LOAD_RNA_SEQ:
        plan.append(
            {
                "label": "RNA-seq",
                "filename": "pre_clinical_rna_seq.csv",
                "model": PreClinicalRnaSeq,
                "required_columns": REQUIRED_MOLECULAR_COLUMNS,
                "value_transform": get_value_transform(RNA_TRANSFORM),
            }
        )

    if LOAD_MICROARRAY:
        plan.append(
            {
                "label": "microarray",
                "filename": "pre_clinical_microarray.csv",
                "model": PreClinicalMicroarray,
                "required_columns": REQUIRED_MOLECULAR_COLUMNS,
                "value_transform": None,
            }
        )

    if LOAD_CNV:
        plan.append(
            {
                "label": "copy-number variation",
                "filename": "pre_clinical_copy_number_variation.csv",
                "model": PreClinicalCopyNumberVariation,
                "required_columns": REQUIRED_MOLECULAR_COLUMNS,
                "value_transform": get_value_transform(CNV_TRANSFORM),
            }
        )

    return tuple(plan)


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


def clean_gene_id(value: Any) -> str | None:
    gene_id = clean_str(value)
    if gene_id is None:
        return None
    return gene_id.split(".", 1)[0]


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


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required PRISM CSV: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])


def iter_csv_chunks(path: Path, *, chunksize: int) -> Iterable[pd.DataFrame]:
    if not path.exists():
        raise FileNotFoundError(f"Missing required PRISM CSV: {path}")

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
    for candidate in ("id", "sampleid", "sample_id"):
        if candidate in sample_df.columns:
            return candidate

    raise ValueError(
        f"{sample_path} must contain one of: id, sampleid, sample_id. "
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
            f"PRISM sample IDs must already be prefixed before database insert. "
            f"Found unprefixed sample IDs, for example: {preview}. "
            "Fix the extractor output or rerun with --auto-prefix-samples."
        )

    return cleaned


def validate_final_tables_model() -> None:
    if PreClinicalDataset.__tablename__ != "datasets":
        raise RuntimeError("tables.py must map PreClinicalDataset to the datasets table.")
    for attr in (
        "name",
        "version",
        "software",
        "link",
        "publication",
        "PMID",
        "description",
        "key_study_findings",
        "clinical",
    ):
        if not hasattr(PreClinicalDataset, attr):
            raise RuntimeError(f"tables.py must define PreClinicalDataset.{attr}.")
    if not hasattr(PreClinicalSample, "id"):
        raise RuntimeError("tables.py must define PreClinicalSample.id as the primary key.")
    if not hasattr(PreClinicalSample, "dataset_id"):
        raise RuntimeError("tables.py must define PreClinicalSample.dataset_id.")
    if not hasattr(PreClinicalSample, "cell_line_name"):
        raise RuntimeError("tables.py must define PreClinicalSample.cell_line_name.")
    if not hasattr(PreClinicalTreatmentResponse, "cid"):
        raise RuntimeError("tables.py must define PreClinicalTreatmentResponse.cid.")
    for model in (PreClinicalRnaSeq, PreClinicalMicroarray, PreClinicalCopyNumberVariation):
        if not hasattr(model, "value"):
            raise RuntimeError(f"tables.py must define {model.__name__}.value.")


def create_required_tables(engine) -> None:
    tables = [
        PreClinicalDataset.__table__,
        PreClinicalCellLine.__table__,
        PreClinicalSample.__table__,
        PreClinicalTreatmentResponse.__table__,
    ]

    if LOAD_RNA_SEQ or LOAD_MICROARRAY or LOAD_CNV:
        tables.append(PreClinicalGene.__table__)
    if LOAD_RNA_SEQ:
        tables.append(PreClinicalRnaSeq.__table__)
    if LOAD_MICROARRAY:
        tables.append(PreClinicalMicroarray.__table__)
    if LOAD_CNV:
        tables.append(PreClinicalCopyNumberVariation.__table__)

    Base.metadata.create_all(bind=engine, tables=tables)


def chunked(items: list[Any], size: int) -> Iterable[list[Any]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def find_dataset_by_name(session: Session, dataset_name: str) -> PreClinicalDataset | None:
    return session.scalar(
        select(PreClinicalDataset).where(
            func.lower(PreClinicalDataset.name) == dataset_name.lower()
        )
    )


def load_dataset_metadata(dataset_name: str, metadata_csv: Path) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "name": dataset_name,
        "version": None,
        "software": None,
        "link": None,
        "publication": None,
        "PMID": None,
        "description": None,
        "key_study_findings": None,
        "clinical": None,
    }

    if not metadata_csv.exists():
        print(
            f"Dataset metadata CSV not found at {metadata_csv}. "
            f"Seeding datasets row with name={dataset_name!r} only."
        )
        return metadata

    df = pd.read_csv(metadata_csv, dtype=str, keep_default_na=False, na_values=[])
    if "name" not in df.columns:
        raise ValueError(f"{metadata_csv} must contain a 'name' column.")

    match = df[df["name"].str.lower() == dataset_name.lower()]
    if match.empty:
        raise ValueError(
            f"Could not find dataset {dataset_name!r} in {metadata_csv}. "
            f"Available names: {', '.join(df['name'].tolist())}"
        )

    row = match.iloc[0].to_dict()
    metadata.update(
        {
            "name": clean_str(row.get("name")) or dataset_name,
            "version": clean_str(row.get("version")),
            "software": clean_str(row.get("software")),
            "link": clean_str(row.get("link")),
            "publication": clean_str(row.get("publication")),
            "PMID": clean_str(row.get("PMID")),
            "description": clean_str(row.get("description")),
            "key_study_findings": clean_str(row.get("key study findings")),
            "clinical": clean_bool(row.get("clinical")),
        }
    )
    return metadata


def delete_existing_dataset(session: Session, dataset_name: str) -> None:
    dataset = find_dataset_by_name(session, dataset_name)
    if dataset is None:
        return

    sample_ids = list(
        session.scalars(
            select(PreClinicalSample.id).where(PreClinicalSample.dataset_id == dataset.id)
        )
    )

    if sample_ids:
        for sample_id_chunk in chunked(sample_ids, 10_000):
            if LOAD_RNA_SEQ:
                session.execute(
                    delete(PreClinicalRnaSeq).where(
                        PreClinicalRnaSeq.sample_id.in_(sample_id_chunk)
                    )
                )
            if LOAD_MICROARRAY:
                session.execute(
                    delete(PreClinicalMicroarray).where(
                        PreClinicalMicroarray.sample_id.in_(sample_id_chunk)
                    )
                )
            if LOAD_CNV:
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
    session.execute(delete(PreClinicalSample).where(PreClinicalSample.dataset_id == dataset.id))
    session.execute(delete(PreClinicalCellLine).where(PreClinicalCellLine.dataset_id == dataset.id))
    session.execute(delete(PreClinicalDataset).where(PreClinicalDataset.id == dataset.id))
    session.flush()


def get_or_create_dataset(
    session: Session,
    *,
    dataset_name: str,
    metadata_csv: Path,
) -> PreClinicalDataset:
    metadata = load_dataset_metadata(dataset_name, metadata_csv)
    dataset = find_dataset_by_name(session, metadata["name"])

    if dataset is None:
        dataset = PreClinicalDataset(**metadata)
        session.add(dataset)
        session.flush()
        return dataset

    for key, value in metadata.items():
        setattr(dataset, key, value)
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
    print(f"Seeded PRISM cell lines: {len(lookup)}")
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
    sample_df["__id"] = ensure_prefixed_sample_ids(
        sample_df[sample_id_col].tolist(),
        auto_prefix=auto_prefix_samples,
    )
    sample_df["cell_line_name"] = sample_df["cell_line_name"].map(clean_str)
    sample_df = sample_df[sample_df["__id"].notna() & sample_df["cell_line_name"].notna()]
    sample_df = sample_df.drop_duplicates(subset=["__id"], keep="first")

    missing_cell_lines = sorted(set(sample_df["cell_line_name"]) - set(cell_line_name_to_id.keys()))
    if missing_cell_lines:
        preview = ", ".join(missing_cell_lines[:20])
        raise ValueError(
            f"Some PRISM sample rows reference cell lines that were not loaded into "
            f"pre_clinical_cell_line for dataset_id={dataset_id}. Examples: {preview}"
        )

    rows = [
        PreClinicalSample(
            id=clean_str(row.get("__id")),
            dataset_id=dataset_id,
            cell_line_name=clean_str(row.get("cell_line_name")),
        )
        for row in sample_df.to_dict(orient="records")
    ]

    session.add_all(rows)
    session.flush()

    sample_ids = set(sample_df["__id"])
    print(f"Seeded PRISM samples: {len(rows)}")
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
    treatment_df["cid"] = treatment_df["cid"].map(clean_str)
    treatment_df = treatment_df[
        treatment_df["cell_line_name"].notna() & treatment_df["treatment_id"].notna()
    ]

    missing_cell_lines = sorted(set(treatment_df["cell_line_name"]) - set(cell_line_name_to_id.keys()))
    if missing_cell_lines:
        preview = ", ".join(missing_cell_lines[:20])
        raise ValueError(
            f"Some PRISM treatment-response rows reference cell lines that were not "
            f"loaded into pre_clinical_cell_line for dataset_id={dataset_id}. "
            f"Examples: {preview}"
        )

    treatment_df = treatment_df.drop_duplicates(subset=["cell_line_name", "treatment_id"], keep="first")

    rows = [
        {
            "dataset_id": dataset_id,
            "cell_line_name": clean_str(row.get("cell_line_name")),
            "treatment_id": clean_str(row.get("treatment_id")),
            "cid": clean_str(row.get("cid")),
            "ic50_recomputed": clean_float(row.get("ic50_recomputed")),
            "acc_recomputed": clean_float(row.get("acc_recomputed")),
            "mechanism_of_action": clean_str(row.get("mechanism_of_action")),
        }
        for row in treatment_df.to_dict(orient="records")
    ]

    if rows:
        session.bulk_insert_mappings(PreClinicalTreatmentResponse, rows)
        session.flush()

    print(f"Seeded PRISM treatment responses: {len(rows)}")


def seed_genes(session: Session, *, data_dir: Path) -> set[str]:
    gene_path = data_dir / "pre_clinical_gene.csv"
    gene_df = read_csv(gene_path)
    require_columns(gene_df, REQUIRED_GENE_COLUMNS, gene_path)

    gene_df = gene_df.copy()
    gene_df["id"] = gene_df["id"].map(clean_gene_id)
    gene_df["name"] = gene_df["name"].map(clean_str)
    gene_df = gene_df[gene_df["id"].notna()]
    gene_df = gene_df.drop_duplicates(subset=["id"], keep="first")

    rows = [
        {"id": clean_gene_id(row.get("id")), "name": clean_str(row.get("name"))}
        for row in gene_df.to_dict(orient="records")
    ]
    rows = [row for row in rows if row["id"] is not None]

    if rows:
        for row_chunk in chunked(rows, 5_000):
            stmt = mysql_insert(PreClinicalGene.__table__).prefix_with("IGNORE")
            session.execute(stmt, row_chunk)
            session.flush()

    gene_ids = set(gene_df["id"])
    print(
        f"Inserted new PRISM genes where absent; existing gene IDs were skipped. "
        f"Candidate gene IDs: {len(gene_ids)}"
    )
    return gene_ids


def get_sample_ids_for_dataset(session: Session, dataset_id: int) -> set[str]:
    return set(session.scalars(select(PreClinicalSample.id).where(PreClinicalSample.dataset_id == dataset_id)))


def get_gene_ids(session: Session) -> set[str]:
    return set(session.scalars(select(PreClinicalGene.id)))


def validate_molecular_samples_prefixed(sample_ids: set[str], *, label: str) -> None:
    unprefixed = sorted(s for s in sample_ids if not s.startswith(SAMPLE_ID_PREFIX))
    if unprefixed:
        preview = ", ".join(unprefixed[:10])
        raise ValueError(
            f"PRISM {label} contains unprefixed sample_id values. Examples: {preview}. "
            "Fix the extractor output before loading."
        )


def seed_molecular_file(
    session: Session,
    *,
    data_dir: Path,
    filename: str,
    label: str,
    model: type,
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
        chunk_df["gene_id"] = chunk_df["gene_id"].map(clean_gene_id)
        chunk_df["__value"] = chunk_df["value"].map(clean_float)

        if value_transform is not None:
            chunk_df["__value"] = chunk_df["__value"].map(value_transform)

        chunk_df = chunk_df[chunk_df["sample_id"].notna() & chunk_df["gene_id"].notna()]

        chunk_sample_ids = set(chunk_df["sample_id"])
        validate_molecular_samples_prefixed(chunk_sample_ids, label=label)

        missing_sample_ids = sorted(chunk_sample_ids - valid_sample_ids)
        if missing_sample_ids:
            preview = ", ".join(missing_sample_ids[:20])
            raise ValueError(
                f"PRISM {label} has sample_id values missing from pre_clinical_sample. "
                f"Examples: {preview}"
            )

        chunk_gene_ids = set(chunk_df["gene_id"])
        missing_gene_ids = sorted(chunk_gene_ids - valid_gene_ids)
        if missing_gene_ids:
            preview = ", ".join(missing_gene_ids[:20])
            raise ValueError(
                f"PRISM {label} has gene_id values missing from pre_clinical_gene. "
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
                "value": float(row["__value"]),
            }
            for row in chunk_df.to_dict(orient="records")
        ]

        if rows:
            stmt = mysql_insert(model.__table__).prefix_with("IGNORE")
            session.execute(stmt, rows)
            session.flush()

        total_insert_candidates += len(rows)
        print(f"Seeded PRISM {label} chunk {chunk_index}: {len(rows)} candidate rows")

    print(
        f"Finished PRISM {label}: {total_insert_candidates} candidate rows. "
        f"Skipped non-finite/missing values: {total_skipped_missing_value}"
    )


def seed_dataset(
    *,
    data_dir: Path,
    dataset_name: str,
    dataset_metadata_csv: Path,
    replace: bool,
    auto_prefix_samples: bool,
    chunksize: int,
) -> None:
    validate_final_tables_model()
    engine = alchemy_engine()
    create_required_tables(engine)

    with Session(engine) as session:
        if replace:
            print(f"Replacing existing dataset rows for {PRISM}")
            delete_existing_dataset(session, dataset_name)
            session.commit()

        dataset = get_or_create_dataset(
            session,
            dataset_name=dataset_name,
            metadata_csv=dataset_metadata_csv,
        )
        session.flush()
        print(f"Using dataset_id={dataset.id} for dataset_name={dataset.name}")

        cell_line_name_to_id = seed_cell_lines(session, dataset_id=dataset.id, data_dir=data_dir)
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

        molecular_plan = build_molecular_load_plan()
        if molecular_plan:
            seed_genes(session, data_dir=data_dir)
        session.commit()

        if molecular_plan:
            valid_sample_ids = get_sample_ids_for_dataset(session, dataset.id)
            valid_gene_ids = get_gene_ids(session)

            for plan in molecular_plan:
                seed_molecular_file(
                    session,
                    data_dir=data_dir,
                    filename=plan["filename"],
                    label=plan["label"],
                    model=plan["model"],
                    required_columns=plan["required_columns"],
                    valid_sample_ids=valid_sample_ids,
                    valid_gene_ids=valid_gene_ids,
                    chunksize=chunksize,
                    value_transform=plan["value_transform"],
                )
                session.commit()

    print(f"Finished PRISM preclinical seeding. Mutation upload was skipped.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=f"Seed the PRISM preclinical dataset from extracted CSVs. Mutation upload is intentionally skipped."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=f"Directory containing PRISM extracted CSVs. Default: {DEFAULT_DATA_DIR}",
    )
    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help=f"Dataset name to create/use in datasets. Default: {DEFAULT_DATASET_NAME}",
    )
    parser.add_argument(
        "--dataset-metadata-csv",
        type=Path,
        default=DEFAULT_DATASET_METADATA_CSV,
        help=(
            "CSV containing dataset metadata fields: name, version, software, link, "
            "publication, PMID, description, key study findings, clinical. "
            f"Default: {DEFAULT_DATASET_METADATA_CSV}"
        ),
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help=f"Delete existing PRISM dataset rows before reloading.",
    )
    parser.add_argument(
        "--auto-prefix-samples",
        action="store_true",
        help=(
            f"Automatically add {SAMPLE_ID_PREFIX} to sample IDs if the sample CSV is not already prefixed. "
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
    seed_dataset(
        data_dir=args.data_dir,
        dataset_name=args.dataset_name,
        dataset_metadata_csv=args.dataset_metadata_csv,
        replace=args.replace,
        auto_prefix_samples=args.auto_prefix_samples,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()