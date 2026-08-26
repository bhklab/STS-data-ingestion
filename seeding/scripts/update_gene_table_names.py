from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd
from sqlalchemy import select, text
from sqlalchemy.orm import Session

try:
    from .seeding_coordinator_engine import alchemy_engine
    from ..models.tables import PreClinicalGene
except ImportError:
    # Allows running as:
    #   pixi run python seeding/scripts/update_gene_table_names.py
    from seeding.scripts.seeding_coordinator_engine import alchemy_engine
    from seeding.models.tables import PreClinicalGene


BIOMART_GENE_FILE_PATH = Path("extraction/data/raw/preclinical/biomart_genes.csv")
DEFAULT_AUDIT_DIR = Path("extraction/data/raw/preclinical")


def clean_ensembl_gene_id(value: Any) -> str | None:
    """
    Normalize Ensembl IDs by removing version/postfix suffixes.

    Examples:
      ENSG00000002586.20_PAR_Y -> ENSG00000002586
      ENSG00000123456.7        -> ENSG00000123456
      ENSG00000123456          -> ENSG00000123456
    """
    if value is None or pd.isna(value):
        return None

    value = str(value).strip()
    if not value:
        return None

    return value.split(".", 1)[0]


def clean_gene_name(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None

    value = str(value).strip()
    if not value or value.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return None

    return value


def qualified_table_name(schema: str | None, table: str) -> str:
    if schema:
        return f"`{schema}`.`{table}`"
    return f"`{table}`"


def fetch_current_genes() -> pd.DataFrame:
    engine = alchemy_engine()

    with Session(engine) as session:
        rows = session.execute(
            select(
                PreClinicalGene.id.label("gene_id"),
                PreClinicalGene.name.label("current_gene_name"),
            )
        ).all()

    current_df = pd.DataFrame(rows, columns=["gene_id", "current_gene_name"])

    if current_df.empty:
        return current_df

    current_df["gene_id"] = current_df["gene_id"].map(clean_ensembl_gene_id)
    current_df["current_gene_name"] = current_df["current_gene_name"].map(clean_gene_name)
    current_df = current_df.dropna(subset=["gene_id"])
    current_df = current_df.drop_duplicates(subset=["gene_id"], keep="first")

    return current_df


def read_biomart_genes(path: Path) -> pd.DataFrame:
    biomart_df = pd.read_csv(path, dtype=str)

    required_cols = {"Gene stable ID", "Gene name"}
    missing_cols = required_cols - set(biomart_df.columns)

    if missing_cols:
        raise ValueError(
            f"Missing required Biomart columns: {sorted(missing_cols)}. "
            f"Found columns: {list(biomart_df.columns)}"
        )

    biomart_df = biomart_df.rename(
        columns={
            "Gene stable ID": "gene_id",
            "Gene name": "new_gene_name",
        }
    )

    biomart_df["gene_id"] = biomart_df["gene_id"].map(clean_ensembl_gene_id)
    biomart_df["new_gene_name"] = biomart_df["new_gene_name"].map(clean_gene_name)

    biomart_df = biomart_df.dropna(subset=["gene_id", "new_gene_name"])

    # If Biomart has duplicate gene IDs, keep the first non-null name.
    biomart_df = biomart_df.drop_duplicates(subset=["gene_id"], keep="first")

    return biomart_df[["gene_id", "new_gene_name"]]


def build_update_tables(
    current_df: pd.DataFrame,
    biomart_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    merged = current_df.merge(
        biomart_df,
        on="gene_id",
        how="inner",
    )

    if merged.empty:
        empty_updates = pd.DataFrame(
            columns=["gene_id", "old_gene_name", "new_gene_name"]
        )
        empty_correct = pd.DataFrame(
            columns=["gene_id", "gene_name"]
        )
        empty_missing = pd.DataFrame(
            columns=["gene_id", "current_gene_name"]
        )
        return empty_updates, empty_correct, empty_missing

    current_names = merged["current_gene_name"].fillna("")
    new_names = merged["new_gene_name"].fillna("")

    changed_mask = current_names != new_names

    updates_df = merged.loc[
        changed_mask,
        ["gene_id", "current_gene_name", "new_gene_name"],
    ].rename(
        columns={
            "current_gene_name": "old_gene_name",
        }
    )

    correct_df = merged.loc[
        ~changed_mask,
        ["gene_id", "current_gene_name"],
    ].rename(
        columns={
            "current_gene_name": "gene_name",
        }
    )

    biomart_gene_ids = set(biomart_df["gene_id"])
    missing_from_biomart_df = current_df.loc[
        ~current_df["gene_id"].isin(biomart_gene_ids),
        ["gene_id", "current_gene_name"],
    ]

    return updates_df, correct_df, missing_from_biomart_df


def bulk_update_gene_names(
    updates_df: pd.DataFrame,
    *,
    schema: str | None = None,
    chunksize: int = 5000,
    dry_run: bool = False,
) -> int:
    if updates_df.empty:
        return 0

    if dry_run:
        return len(updates_df)

    engine = alchemy_engine()
    gene_table = qualified_table_name(schema, "pre_clinical_gene")

    update_rows = updates_df[["gene_id", "new_gene_name"]].rename(
        columns={
            "gene_id": "id",
            "new_gene_name": "name",
        }
    )

    payload = update_rows.to_dict(orient="records")

    with engine.begin() as conn:
        conn.execute(text("DROP TEMPORARY TABLE IF EXISTS tmp_gene_name_updates"))

        conn.execute(
            text(
                """
                CREATE TEMPORARY TABLE tmp_gene_name_updates (
                    id VARCHAR(255) PRIMARY KEY,
                    name VARCHAR(255) NOT NULL
                ) ENGINE=InnoDB
                """
            )
        )

        insert_stmt = text(
            """
            INSERT INTO tmp_gene_name_updates (id, name)
            VALUES (:id, :name)
            """
        )

        for start in range(0, len(payload), chunksize):
            chunk = payload[start : start + chunksize]
            conn.execute(insert_stmt, chunk)

        result = conn.execute(
            text(
                f"""
                UPDATE {gene_table} AS g
                INNER JOIN tmp_gene_name_updates AS u
                    ON g.id = u.id
                SET g.name = u.name
                """
            )
        )

        conn.execute(text("DROP TEMPORARY TABLE IF EXISTS tmp_gene_name_updates"))

    return int(result.rowcount or 0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bulk-update pre_clinical_gene.name from Biomart gene annotations."
    )

    parser.add_argument(
        "--biomart-csv",
        type=Path,
        default=BIOMART_GENE_FILE_PATH,
        help="Path to Biomart CSV with columns 'Gene stable ID' and 'Gene name'.",
    )

    parser.add_argument(
        "--audit-dir",
        type=Path,
        default=DEFAULT_AUDIT_DIR,
        help="Directory where audit CSVs will be written.",
    )

    parser.add_argument(
        "--schema",
        type=str,
        default=None,
        help=(
            "Optional MySQL schema/database name. "
            "Example: --schema sts_portal_beta. "
            "If omitted, uses the selected DB from your SQLAlchemy connection."
        ),
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=5000,
        help="Number of rows per bulk insert into the temporary update table.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Compute updates and write audit files, but do not update the database.",
    )

    args = parser.parse_args()

    args.audit_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading current gene table from database...")
    current_df = fetch_current_genes()
    print(f"Current DB genes: {len(current_df):,}")

    print(f"Reading Biomart genes from: {args.biomart_csv}")
    biomart_df = read_biomart_genes(args.biomart_csv)
    print(f"Biomart genes: {len(biomart_df):,}")

    updates_df, correct_df, missing_from_biomart_df = build_update_tables(
        current_df=current_df,
        biomart_df=biomart_df,
    )

    updated_genes_path = args.audit_dir / "updated_genes.csv"
    already_correct_path = args.audit_dir / "already_correct_naming.csv"
    missing_biomart_path = args.audit_dir / "genes_missing_from_biomart.csv"

    updates_df.to_csv(updated_genes_path, index=False)
    correct_df.to_csv(already_correct_path, index=False)
    missing_from_biomart_df.to_csv(missing_biomart_path, index=False)

    print(f"Genes needing update: {len(updates_df):,}")
    print(f"Already correct genes: {len(correct_df):,}")
    print(f"DB genes missing from Biomart CSV: {len(missing_from_biomart_df):,}")

    updated_count = bulk_update_gene_names(
        updates_df,
        schema=args.schema,
        chunksize=args.chunksize,
        dry_run=args.dry_run,
    )

    if args.dry_run:
        print(f"Dry run complete. Would update {updated_count:,} gene names.")
    else:
        print(f"Updated {updated_count:,} gene names in the database.")

    print(f"Wrote audit files:")
    print(f"  {updated_genes_path}")
    print(f"  {already_correct_path}")
    print(f"  {missing_biomart_path}")


if __name__ == "__main__":
    main()