#!/usr/bin/env python3
"""
High-Speed TCGA-SARC Clinical Seeding Script using MySQL LOAD DATA LOCAL INFILE.

Loads all TCGA-SARC clinical datasets and digital pathology H&E features
directly into MySQL (15x-30x faster than standard SQL INSERTs):
  - Shared: datasets, pre_clinical_gene
  - Clinical references: clinical_sample, clinical_antigen, clinical_probe
  - Molecular assays: clinical_rna, clinical_cnv, clinical_mutation,
                      clinical_mirna, clinical_rppa, clinical_methylation
  - Digital pathology: clinical_slide, clinical_tile, clinical_embedding

Run from repository root:
    pixi run python -m seeding.scripts.clinical_seed_tcga
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path
from typing import Any

import pandas as pd
from sqlalchemy import delete, func, select, text
from sqlalchemy.dialects.mysql import insert as mysql_insert
from sqlalchemy.orm import Session

from .seeding_coordinator_engine import alchemy_engine
from ..models.tables import (
    Base,
    ClinicalAntigen,
    ClinicalCNV,
    ClinicalEmbedding,
    ClinicalMethylation,
    ClinicalMiRNA,
    ClinicalMutation,
    ClinicalProbe,
    ClinicalRNA,
    ClinicalRPPA,
    ClinicalSample,
    ClinicalSlide,
    ClinicalTile,
    Dataset,
    PreClinicalGene,
)

DEFAULT_DATA_DIR = Path("extraction/data/proc/clinical/TCGA_SARC")
DEFAULT_DATASET_NAME = "TCGA-SARC"
DEFAULT_DATASET_METADATA_CSV = Path("extraction/data/proc/clinical_datasets.csv")

# Table loading toggles for partial reload testing
LOAD_GENES = True
LOAD_ANTIGEN = True
LOAD_PROBE = True
LOAD_SAMPLE = True
LOAD_RNA = True
LOAD_CNV = True
LOAD_MUTATION = True
LOAD_MIRNA = True
LOAD_RPPA = True
LOAD_METHYLATION = True
LOAD_SLIDES = True
LOAD_TILES = True
LOAD_EMBEDDINGS = True

REQUIRED_SAMPLE_COLUMNS = {
    "id",
    "race",
    "ethnicity",
    "sex",
    "age",
    "histology",
    "tissue",
    "tissue_origin",
}
REQUIRED_SLIDE_COLUMNS = {"id", "sample_id", "n_tiles", "embedding_dim"}


# -------------------------------------------------------------------------
# clean_* helpers
# -------------------------------------------------------------------------


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


def clean_bool(value: Any) -> bool | None:
    value = clean_value(value)
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    text_val = str(value).strip().lower()
    if text_val in {"true", "t", "1", "yes", "y"}:
        return True
    if text_val in {"false", "f", "0", "no", "n"}:
        return False
    return None


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required TCGA-SARC clinical CSV: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])


def require_columns(df: pd.DataFrame, required: set[str], path: Path) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(missing)}. "
            f"Found columns: {', '.join(df.columns)}"
        )


def chunked(items: list[Any], size: int) -> Any:
    for start in range(0, len(items), size):
        yield items[start : start + size]


# -------------------------------------------------------------------------
# High-Speed LOAD DATA LOCAL INFILE helper
# -------------------------------------------------------------------------


def bulk_load_csv(
    raw_conn: Any,
    table_name: str,
    csv_path: Path,
    columns_clause: str,
    set_clause: str = "",
    label: str = "",
) -> None:
    if not csv_path.exists():
        print(f"Skipping {label or table_name}: {csv_path} not found.")
        return

    posix_path = csv_path.resolve().as_posix()
    set_sql = f"SET {set_clause}" if set_clause else ""
    sql = f"""
    LOAD DATA LOCAL INFILE '{posix_path}'
    IGNORE INTO TABLE {table_name}
    FIELDS TERMINATED BY ',' ENCLOSED BY '\"'
    LINES TERMINATED BY '\\n'
    IGNORE 1 LINES
    ({columns_clause})
    {set_sql};
    """

    start_time = time.time()
    tag = label or table_name
    print(f"Streaming {tag} via LOAD DATA LOCAL INFILE...")

    with raw_conn.cursor() as cursor:
        cursor.execute("SET FOREIGN_KEY_CHECKS = 0;")
        cursor.execute("SET UNIQUE_CHECKS = 0;")
        cursor.execute(sql)
        cursor.execute("SET FOREIGN_KEY_CHECKS = 1;")
        cursor.execute("SET UNIQUE_CHECKS = 1;")
    raw_conn.commit()

    elapsed = time.time() - start_time
    print(f"Finished {tag} in {elapsed:.2f}s.")


# -------------------------------------------------------------------------
# Dataset row
# -------------------------------------------------------------------------


def find_dataset_by_name(session: Session, dataset_name: str) -> Dataset | None:
    return session.scalar(
        select(Dataset).where(func.lower(Dataset.name) == dataset_name.lower())
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


def get_or_create_dataset(
    session: Session,
    *,
    dataset_name: str,
    metadata_csv: Path,
) -> Dataset:
    metadata = load_dataset_metadata(dataset_name, metadata_csv)
    dataset = find_dataset_by_name(session, metadata["name"])

    if dataset is None:
        dataset = Dataset(**metadata)
        session.add(dataset)
        session.flush()
        return dataset

    for key, value in metadata.items():
        setattr(dataset, key, value)
    session.flush()
    return dataset


# -------------------------------------------------------------------------
# Table validation
# -------------------------------------------------------------------------


def validate_final_tables_model() -> None:
    if Dataset.__tablename__ != "datasets":
        raise RuntimeError("tables.py must map Dataset to the datasets table.")
    if not hasattr(ClinicalSample, "id"):
        raise RuntimeError("tables.py must define ClinicalSample.id as the primary key.")
    if not hasattr(ClinicalSample, "dataset_id"):
        raise RuntimeError("tables.py must define ClinicalSample.dataset_id.")
    for model in (ClinicalRNA, ClinicalCNV, ClinicalRPPA, ClinicalMethylation):
        if not hasattr(model, "value"):
            raise RuntimeError(f"tables.py must define {model.__name__}.value.")
    if not hasattr(ClinicalMutation, "mutation") or not hasattr(ClinicalMutation, "oncoprint"):
        raise RuntimeError("tables.py must define ClinicalMutation.mutation and .oncoprint.")
    for model in (ClinicalSlide, ClinicalTile, ClinicalEmbedding):
        if not hasattr(model, "id"):
            raise RuntimeError(f"tables.py must define {model.__name__}.id.")


def create_required_tables(engine) -> None:
    tables = [
        Dataset.__table__,
        PreClinicalGene.__table__,
        ClinicalSample.__table__,
        ClinicalAntigen.__table__,
        ClinicalProbe.__table__,
        ClinicalRNA.__table__,
        ClinicalCNV.__table__,
        ClinicalMutation.__table__,
        ClinicalMiRNA.__table__,
        ClinicalRPPA.__table__,
        ClinicalMethylation.__table__,
        ClinicalSlide.__table__,
        ClinicalTile.__table__,
        ClinicalEmbedding.__table__,
    ]
    Base.metadata.create_all(bind=engine, tables=tables)


# -------------------------------------------------------------------------
# Delete existing dataset (for --replace)
# -------------------------------------------------------------------------


def delete_existing_dataset(session: Session, dataset_name: str) -> None:
    dataset = find_dataset_by_name(session, dataset_name)
    if dataset is None:
        return

    sample_ids = list(
        session.scalars(
            select(ClinicalSample.id).where(ClinicalSample.dataset_id == dataset.id)
        )
    )

    if sample_ids:
        for sample_id_chunk in chunked(sample_ids, 10_000):
            session.execute(delete(ClinicalRNA).where(ClinicalRNA.sample_id.in_(sample_id_chunk)))
            session.execute(delete(ClinicalCNV).where(ClinicalCNV.sample_id.in_(sample_id_chunk)))
            session.execute(delete(ClinicalMutation).where(ClinicalMutation.sample_id.in_(sample_id_chunk)))
            session.execute(delete(ClinicalMiRNA).where(ClinicalMiRNA.sample_id.in_(sample_id_chunk)))
            session.execute(delete(ClinicalRPPA).where(ClinicalRPPA.sample_id.in_(sample_id_chunk)))
            session.execute(delete(ClinicalMethylation).where(ClinicalMethylation.sample_id.in_(sample_id_chunk)))

    slide_ids = list(
        session.scalars(
            select(ClinicalSlide.id).where(ClinicalSlide.dataset_id == dataset.id)
        )
    )
    if slide_ids:
        for slide_id_chunk in chunked(slide_ids, 5_000):
            session.execute(delete(ClinicalEmbedding).where(ClinicalEmbedding.slide_id.in_(slide_id_chunk)))
            session.execute(delete(ClinicalTile).where(ClinicalTile.slide_id.in_(slide_id_chunk)))
        session.execute(delete(ClinicalSlide).where(ClinicalSlide.dataset_id == dataset.id))

    session.execute(delete(ClinicalSample).where(ClinicalSample.dataset_id == dataset.id))
    session.execute(delete(Dataset).where(Dataset.id == dataset.id))
    session.flush()


# -------------------------------------------------------------------------
# Small Reference Table Seeders
# -------------------------------------------------------------------------


def seed_clinical_sample(session: Session, *, dataset_id: int, data_dir: Path) -> set[str]:
    sample_path = data_dir / "clinical_sample.csv"
    sample_df = read_csv(sample_path)
    require_columns(sample_df, REQUIRED_SAMPLE_COLUMNS, sample_path)

    sample_df = sample_df.copy()
    sample_df["id"] = sample_df["id"].map(clean_str)
    sample_df = sample_df[sample_df["id"].notna()]
    sample_df = sample_df.drop_duplicates(subset=["id"], keep="first")

    if sample_df.empty:
        raise ValueError(f"No usable samples found in {sample_path}")

    rows = [
        {
            "id": clean_str(row.get("id")),
            "dataset_id": dataset_id,
            "race": clean_str(row.get("race")),
            "ethnicity": clean_str(row.get("ethnicity")),
            "sex": clean_str(row.get("sex")),
            "age": clean_int(row.get("age")),
            "histology": clean_str(row.get("histology")),
            "tissue": clean_str(row.get("tissue")),
            "tissue_origin": clean_str(row.get("tissue_origin")),
        }
        for row in sample_df.to_dict(orient="records")
    ]

    stmt = mysql_insert(ClinicalSample.__table__).prefix_with("IGNORE")
    session.execute(stmt, rows)
    session.flush()

    sample_ids = set(
        session.scalars(
            select(ClinicalSample.id).where(ClinicalSample.dataset_id == dataset_id)
        )
    )
    print(f"Seeded TCGA-SARC clinical samples: {len(sample_ids)}")
    return sample_ids


def seed_clinical_antigen(session: Session, *, data_dir: Path) -> set[str]:
    antigen_path = data_dir / "clinical_antigen.csv"
    antigen_df = read_csv(antigen_path)

    antigen_df = antigen_df.copy()
    antigen_df["id"] = antigen_df["id"].map(clean_str)
    antigen_df = antigen_df[antigen_df["id"].notna()]
    antigen_df = antigen_df.drop_duplicates(subset=["id"], keep="first")

    valid_genes = set(session.scalars(select(PreClinicalGene.id)))

    rows = []
    for row in antigen_df.to_dict(orient="records"):
        gene_id = clean_gene_id(row.get("peptide_target_gene"))
        if gene_id is not None and gene_id not in valid_genes:
            gene_id = None
        rows.append(
            {
                "id": clean_str(row.get("id")),
                "catalogue_number": clean_str(row.get("catalogue_number")),
                "peptide_target": clean_str(row.get("peptide_target")),
                "peptide_target_gene": gene_id,
            }
        )

    if rows:
        stmt = mysql_insert(ClinicalAntigen.__table__).prefix_with("IGNORE")
        session.execute(stmt, rows)
        session.flush()

    antigen_ids = set(session.scalars(select(ClinicalAntigen.id)))
    print(f"Seeded TCGA-SARC clinical antigens: {len(antigen_ids)}")
    return antigen_ids


def seed_clinical_slide(
    session: Session, *, dataset_id: int, data_dir: Path, valid_sample_ids: set[str]
) -> set[str]:
    slide_path = data_dir / "clinical_slide.csv"
    if not slide_path.exists():
        print(f"clinical_slide.csv not found at {slide_path}, skipping.")
        return set()

    slide_df = read_csv(slide_path)
    require_columns(slide_df, REQUIRED_SLIDE_COLUMNS, slide_path)

    slide_df = slide_df.copy()
    slide_df["id"] = slide_df["id"].map(clean_str)
    slide_df["sample_id"] = slide_df["sample_id"].map(clean_str)
    slide_df = slide_df[slide_df["id"].notna() & slide_df["sample_id"].notna()]
    slide_df = slide_df[slide_df["sample_id"].isin(valid_sample_ids)]
    slide_df = slide_df.drop_duplicates(subset=["id"], keep="first")

    rows = [
        {
            "id": clean_str(row.get("id")),
            "sample_id": clean_str(row.get("sample_id")),
            "dataset_id": dataset_id,
            "n_tiles": clean_int(row.get("n_tiles")),
            "embedding_dim": clean_int(row.get("embedding_dim")),
        }
        for row in slide_df.to_dict(orient="records")
    ]

    if rows:
        stmt = mysql_insert(ClinicalSlide.__table__).prefix_with("IGNORE")
        session.execute(stmt, rows)
        session.flush()

    slide_ids = set(
        session.scalars(
            select(ClinicalSlide.id).where(ClinicalSlide.dataset_id == dataset_id)
        )
    )
    print(f"Seeded TCGA-SARC clinical slides: {len(slide_ids)}")
    return slide_ids


# -------------------------------------------------------------------------
# High-Speed Embedding Staging & Ingestion
# -------------------------------------------------------------------------


def seed_clinical_embedding_fast(raw_conn: Any, data_dir: Path) -> None:
    embed_csv = data_dir / "clinical_embedding.csv"
    if not embed_csv.exists():
        print(f"clinical_embedding.csv not found at {embed_csv}, skipping.")
        return

    posix_path = embed_csv.resolve().as_posix()
    print("Staging and inserting clinical embeddings via LOAD DATA LOCAL INFILE + UNHEX...")
    t0 = time.time()

    with raw_conn.cursor() as cursor:
        cursor.execute("SET FOREIGN_KEY_CHECKS = 0;")
        cursor.execute("SET UNIQUE_CHECKS = 0;")

        # 1. Create temporary staging table
        cursor.execute(
            """
            CREATE TEMPORARY TABLE IF NOT EXISTS temp_clinical_embedding (
                slide_id VARCHAR(150) NOT NULL,
                tile_index INT NULL,
                model_name VARCHAR(100) NOT NULL,
                kind VARCHAR(20) NOT NULL,
                embedding_hex MEDIUMTEXT NOT NULL
            );
            """
        )
        cursor.execute("TRUNCATE TABLE temp_clinical_embedding;")

        # 2. Bulk load hex text into temporary table
        load_sql = f"""
        LOAD DATA LOCAL INFILE '{posix_path}'
        INTO TABLE temp_clinical_embedding
        FIELDS TERMINATED BY ',' ENCLOSED BY '\"'
        LINES TERMINATED BY '\\n'
        IGNORE 1 LINES
        (slide_id, @v_idx, model_name, kind, embedding_hex)
        SET tile_index = NULLIF(TRIM(@v_idx), '');
        """
        cursor.execute(load_sql)

        # 3. Insert into target clinical_embedding resolving tile_id and unhexing binary
        insert_sql = """
        INSERT IGNORE INTO clinical_embedding (slide_id, tile_id, model_name, kind, embedding)
        SELECT
            t.slide_id,
            ct.id AS tile_id,
            t.model_name,
            t.kind,
            UNHEX(t.embedding_hex)
        FROM temp_clinical_embedding t
        LEFT JOIN clinical_tile ct
            ON t.slide_id = ct.slide_id AND t.tile_index = ct.tile_index;
        """
        cursor.execute(insert_sql)
        cursor.execute("DROP TEMPORARY TABLE IF EXISTS temp_clinical_embedding;")

        cursor.execute("SET FOREIGN_KEY_CHECKS = 1;")
        cursor.execute("SET UNIQUE_CHECKS = 1;")

    raw_conn.commit()
    elapsed = time.time() - t0
    print(f"Finished clinical_embedding in {elapsed:.2f}s.")


# -------------------------------------------------------------------------
# Master Orchestration
# -------------------------------------------------------------------------


def seed_dataset(
    *,
    data_dir: Path,
    dataset_name: str,
    dataset_metadata_csv: Path,
    replace: bool,
) -> None:
    validate_final_tables_model()
    engine = alchemy_engine()
    create_required_tables(engine)

    total_start = time.time()

    with Session(engine) as session:
        if replace:
            print(f"Replacing existing dataset rows for {dataset_name}")
            delete_existing_dataset(session, dataset_name)
            session.commit()

        dataset = get_or_create_dataset(
            session,
            dataset_name=dataset_name,
            metadata_csv=dataset_metadata_csv,
        )
        session.flush()
        print(f"Using dataset_id={dataset.id} for dataset_name={dataset.name}")

        valid_sample_ids = seed_clinical_sample(session, dataset_id=dataset.id, data_dir=data_dir)
        session.commit()

        valid_slide_ids = seed_clinical_slide(
            session, dataset_id=dataset.id, data_dir=data_dir, valid_sample_ids=valid_sample_ids
        )
        session.commit()

    # Raw connection for ultra-fast LOAD DATA LOCAL INFILE execution
    with engine.raw_connection() as raw_conn:
        if LOAD_GENES:
            bulk_load_csv(
                raw_conn,
                table_name="pre_clinical_gene",
                csv_path=data_dir / "pre_clinical_gene.csv",
                columns_clause="id, name",
                label="pre_clinical_gene",
            )

        with Session(engine) as session:
            valid_antigen_ids = seed_clinical_antigen(session, data_dir=data_dir)
            session.commit()

        if LOAD_PROBE:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_probe",
                csv_path=data_dir / "clinical_probe.csv",
                columns_clause="id, @v_name",
                set_clause="name = NULLIF(TRIM(@v_name), '')",
                label="clinical_probe",
            )

        if LOAD_RNA:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_rna",
                csv_path=data_dir / "clinical_rna.csv",
                columns_clause="sample_id, gene_id, value",
                label="clinical_rna (16.07M rows)",
            )

        if LOAD_CNV:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_cnv",
                csv_path=data_dir / "clinical_cnv.csv",
                columns_clause="sample_id, gene_id, value",
                label="clinical_cnv (15.07M rows)",
            )

        if LOAD_MUTATION:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_mutation",
                csv_path=data_dir / "clinical_mutation.csv",
                columns_clause="sample_id, gene_id, @v_mut, @v_onc",
                set_clause=(
                    "mutation = CASE "
                    "WHEN @v_mut IN ('TRUE', '1', 'true') THEN 1 "
                    "WHEN @v_mut IN ('FALSE', '0', 'false') THEN 0 "
                    "ELSE NULL END, "
                    "oncoprint = NULLIF(TRIM(@v_onc), '')"
                ),
                label="clinical_mutation (2.11M rows)",
            )

        if LOAD_MIRNA:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_mirna",
                csv_path=data_dir / "clinical_mirna.csv",
                columns_clause="id, sample_id, @v_gene, value",
                set_clause="gene_id = NULLIF(TRIM(@v_gene), '')",
                label="clinical_mirna (494k rows)",
            )

        if LOAD_RPPA:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_rppa",
                csv_path=data_dir / "clinical_rppa.csv",
                columns_clause="sample_id, antigen_id, value",
                label="clinical_rppa (103k rows)",
            )

        if LOAD_METHYLATION:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_methylation",
                csv_path=data_dir / "clinical_methylation.csv",
                columns_clause="sample_id, probe_id, value",
                label="clinical_methylation (109.28M rows)",
            )

        if LOAD_TILES:
            bulk_load_csv(
                raw_conn,
                table_name="clinical_tile",
                csv_path=data_dir / "clinical_tile.csv",
                columns_clause="slide_id, tile_index, x, y",
                label="clinical_tile (2.51M rows)",
            )

        if LOAD_EMBEDDINGS:
            seed_clinical_embedding_fast(raw_conn, data_dir=data_dir)

    total_time = time.time() - total_start
    print("=" * 70)
    print(f" TCGA-SARC Seeding Complete! Total Time: {total_time:.2f}s ({total_time/60:.2f} mins)")
    print("=" * 70)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="High-Speed Seed of TCGA-SARC clinical dataset via MySQL LOAD DATA LOCAL INFILE."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=f"Directory containing TCGA-SARC extracted clinical CSVs. Default: {DEFAULT_DATA_DIR}",
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
        help=f"CSV containing dataset metadata fields. Default: {DEFAULT_DATASET_METADATA_CSV}",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Delete existing TCGA-SARC dataset rows before reloading.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    seed_dataset(
        data_dir=args.data_dir,
        dataset_name=args.dataset_name,
        dataset_metadata_csv=args.dataset_metadata_csv,
        replace=args.replace,
    )


if __name__ == "__main__":
    main()
