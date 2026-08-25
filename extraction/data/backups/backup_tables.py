#!/usr/bin/env python3
"""
Database Table Backup Utility for sts_portal_beta
=================================================

Exports all (or specified) tables from the MySQL database into CSV files.
Supports chunked reading and streaming to handle large tables (millions of rows)
without running into memory constraints.

Usage:
    # Run backup with default settings (saves to a timestamped folder inside extraction/data/backups):
    python backup_tables.py

    # Or run with pixi:
    pixi run python extraction/data/backups/backup_tables.py

    # Options:
    python backup_tables.py --help
    python backup_tables.py --no-timestamp            # Save CSVs directly in extraction/data/backups
    python backup_tables.py --compress                # Save as compressed .csv.gz files
    python backup_tables.py --tables datasets drugs   # Backup specific tables only
    python backup_tables.py --chunksize 100000        # Custom batch size
"""

import argparse
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from urllib.parse import quote_plus

import pandas as pd
from dotenv import find_dotenv, load_dotenv
from sqlalchemy import create_engine, inspect, text


def find_and_load_env(custom_env_path: str | None = None) -> None:
    """Locate and load .env file from custom path, workspace root, or parent directories."""
    if custom_env_path and os.path.exists(custom_env_path):
        load_dotenv(custom_env_path, override=True)
        return

    # Check project root by walking up from current script directory
    script_dir = Path(__file__).resolve().parent
    for parent in [script_dir, *script_dir.parents]:
        env_candidate = parent / ".env"
        if env_candidate.exists():
            load_dotenv(env_candidate, override=True)
            return

    # Fallback to dotenv find_dotenv
    dotenv_file = find_dotenv(usecwd=True)
    if dotenv_file:
        load_dotenv(dotenv_file, override=True)


def get_db_engine(db_name: str | None = None):
    """Create SQLAlchemy engine using environment variables."""
    user = os.getenv("DATABASE_USER")
    password = os.getenv("DATABASE_PASS", "")
    host = os.getenv("DATABASE_IP")
    port = os.getenv("PORT", "3306")
    db = db_name or os.getenv("SELECTED_DB", "sts_portal_beta")

    missing = []
    if not user:
        missing.append("DATABASE_USER")
    if not host:
        missing.append("DATABASE_IP")
    if not db:
        missing.append("SELECTED_DB")

    if missing:
        raise ValueError(
            f"Missing required database environment variables: {', '.join(missing)}. "
            f"Please check your .env file."
        )

    password_cleaned = quote_plus(password)
    connection_url = f"mysql+pymysql://{user}:{password_cleaned}@{host}:{port}/{db}"
    return create_engine(connection_url, echo=False)


def format_bytes(size_in_bytes: int) -> str:
    """Format bytes to human-readable string."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_in_bytes < 1024.0:
            return f"{size_in_bytes:.2f} {unit}"
        size_in_bytes /= 1024.0
    return f"{size_in_bytes:.2f} PB"


def format_duration(seconds: float) -> str:
    """Format duration in seconds to human-readable string."""
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes = int(seconds // 60)
    rem_seconds = seconds % 60
    return f"{minutes}m {rem_seconds:.1f}s"


def backup_table(
    engine,
    table_name: str,
    output_dir: Path,
    chunksize: int = 50000,
    compress: bool = False,
) -> dict:
    """Export a single table to CSV in chunks."""
    ext = ".csv.gz" if compress else ".csv"
    output_file = output_dir / f"{table_name}{ext}"
    compression_arg = "gzip" if compress else None

    start_time = time.time()

    # Get total row count
    with engine.connect() as conn:
        count_res = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`")).scalar()
        total_rows = int(count_res) if count_res is not None else 0

    print(f"\n[+] Backing up '{table_name}' ({total_rows:,} rows)...")

    if total_rows == 0:
        # Write empty DataFrame with headers if table has 0 rows
        with engine.connect() as conn:
            df_empty = pd.read_sql_query(
                text(f"SELECT * FROM `{table_name}` LIMIT 0"),
                con=conn,
            )
            df_empty.to_csv(
                output_file,
                index=False,
                compression=compression_arg,
            )
        elapsed = time.time() - start_time
        file_size = output_file.stat().st_size if output_file.exists() else 0
        print(f"    Saved empty table -> {output_file.name} ({format_duration(elapsed)})")
        return {
            "table": table_name,
            "rows": 0,
            "file": output_file.name,
            "file_path": str(output_file.resolve()),
            "size_bytes": file_size,
            "size_formatted": format_bytes(file_size),
            "duration_seconds": round(elapsed, 2),
            "status": "success",
        }

    rows_written = 0
    chunk_count = 0

    # Ensure existing file is removed before appending
    if output_file.exists():
        output_file.unlink()

    with engine.connect().execution_options(stream_results=True) as conn:
        for chunk in pd.read_sql_query(
            text(f"SELECT * FROM `{table_name}`"),
            con=conn,
            chunksize=chunksize,
        ):
            is_first_chunk = (chunk_count == 0)
            chunk.to_csv(
                output_file,
                mode="a",
                header=is_first_chunk,
                index=False,
                compression=compression_arg,
            )
            rows_written += len(chunk)
            chunk_count += 1

            pct = (rows_written / total_rows * 100) if total_rows > 0 else 100.0
            print(
                f"    -> Exported {rows_written:,}/{total_rows:,} rows ({pct:.1f}%)",
                end="\r",
                flush=True,
            )

    elapsed = time.time() - start_time
    file_size = output_file.stat().st_size if output_file.exists() else 0
    print(
        f"    Saved {rows_written:,} rows -> {output_file.name} "
        f"({format_bytes(file_size)}, {format_duration(elapsed)})"
    )

    return {
        "table": table_name,
        "rows": rows_written,
        "file": output_file.name,
        "file_path": str(output_file.resolve()),
        "size_bytes": file_size,
        "size_formatted": format_bytes(file_size),
        "duration_seconds": round(elapsed, 2),
        "status": "success",
    }


def run_backup(
    output_base_dir: Path | None = None,
    use_timestamp_dir: bool = True,
    tables_to_backup: list[str] | None = None,
    chunksize: int = 50000,
    compress: bool = False,
    db_name: str | None = None,
    env_file: str | None = None,
) -> None:
    """Execute the full database backup."""
    total_start_time = time.time()
    find_and_load_env(env_file)

    script_default_dir = Path(__file__).resolve().parent
    base_dir = output_base_dir or script_default_dir

    db = db_name or os.getenv("SELECTED_DB", "sts_portal_beta")
    timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")

    if use_timestamp_dir:
        backup_dir = base_dir / f"{db}_backup_{timestamp_str}"
    else:
        backup_dir = base_dir

    backup_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(" STS Portal Database Backup")
    print("=" * 70)
    print(f"Database:      {db}")
    print(f"Target Dir:    {backup_dir.resolve()}")
    print(f"Chunk Size:    {chunksize:,} rows")
    print(f"Compression:   {'gzip (.csv.gz)' if compress else 'None (.csv)'}")
    print(f"Start Time:    {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    engine = get_db_engine(db_name=db)
    inspector = inspect(engine)
    available_tables = sorted(inspector.get_table_names())

    if not available_tables:
        print(f"No tables found in database '{db}'.")
        return

    if tables_to_backup:
        target_tables = [t for t in tables_to_backup if t in available_tables]
        missing = [t for t in tables_to_backup if t not in available_tables]
        if missing:
            print(f"[!] Warning: Specified tables not found: {', '.join(missing)}")
    else:
        target_tables = available_tables

    print(f"Tables to backup ({len(target_tables)}): {', '.join(target_tables)}")

    results = []
    failed = []

    for tbl in target_tables:
        try:
            res = backup_table(
                engine=engine,
                table_name=tbl,
                output_dir=backup_dir,
                chunksize=chunksize,
                compress=compress,
            )
            results.append(res)
        except Exception as e:
            print(f"\n[!] ERROR backing up table '{tbl}': {e}")
            failed.append({"table": tbl, "error": str(e), "status": "failed"})

    total_elapsed = time.time() - total_start_time
    total_rows = sum(r["rows"] for r in results)
    total_bytes = sum(r["size_bytes"] for r in results)

    # Save manifest
    manifest = {
        "database": db,
        "created_at": datetime.now().isoformat(),
        "timestamp_local": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "total_tables": len(target_tables),
        "successful_tables": len(results),
        "failed_tables": len(failed),
        "total_rows": total_rows,
        "total_size_bytes": total_bytes,
        "total_size_formatted": format_bytes(total_bytes),
        "total_duration_seconds": round(total_elapsed, 2),
        "tables": results,
        "failures": failed,
    }

    manifest_file = backup_dir / "backup_manifest.json"
    with open(manifest_file, "w") as f:
        json.dump(manifest, f, indent=2)

    print("\n" + "=" * 70)
    print(" Backup Summary")
    print("=" * 70)
    print(f"Status:             {'COMPLETED' if not failed else 'COMPLETED WITH ERRORS'}")
    print(f"Tables Backed Up:   {len(results)}/{len(target_tables)}")
    print(f"Total Rows Saved:   {total_rows:,}")
    print(f"Total Backup Size:  {format_bytes(total_bytes)}")
    print(f"Elapsed Time:       {format_duration(total_elapsed)}")
    print(f"Manifest File:      {manifest_file.resolve()}")
    print(f"Backup Location:    {backup_dir.resolve()}")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Backup all tables from sts_portal_beta database to CSV files."
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=str,
        default=None,
        help="Target directory to store backups (default: extraction/data/backups/)",
    )
    parser.add_argument(
        "--no-timestamp",
        action="store_true",
        help="Do not create a timestamped subdirectory; save CSVs directly in output directory",
    )
    parser.add_argument(
        "--tables",
        "-t",
        nargs="+",
        default=None,
        help="Specify one or more table names to backup instead of all tables",
    )
    parser.add_argument(
        "--chunksize",
        "-c",
        type=int,
        default=50000,
        help="Batch size (number of rows) for reading/writing tables (default: 50000)",
    )
    parser.add_argument(
        "--compress",
        "-z",
        action="store_true",
        help="Compress CSV output using gzip (.csv.gz)",
    )
    parser.add_argument(
        "--db",
        type=str,
        default=None,
        help="Database name (defaults to SELECTED_DB from .env or 'sts_portal_beta')",
    )
    parser.add_argument(
        "--env-file",
        type=str,
        default=None,
        help="Custom path to .env file",
    )

    args = parser.parse_args()

    output_base_dir = Path(args.output_dir).resolve() if args.output_dir else None

    run_backup(
        output_base_dir=output_base_dir,
        use_timestamp_dir=not args.no_timestamp,
        tables_to_backup=args.tables,
        chunksize=args.chunksize,
        compress=args.compress,
        db_name=args.db,
        env_file=args.env_file,
    )


if __name__ == "__main__":
    main()
