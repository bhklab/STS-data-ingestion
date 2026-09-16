#!/usr/bin/env python3
"""
Extract slide-related CSVs from TCGA-SARC WSI feature files (.pt) for loading
into the clinical schema defined in seeding/models/tables.py:
  1. clinical_slide.csv      -> ClinicalSlide (slide metadata)
  2. clinical_tile.csv       -> ClinicalTile (tile coordinates: x, y)
  3. clinical_embedding.csv  -> ClinicalEmbedding (slide-level mean and tile embeddings)

Run from repository root:
    pixi run python extraction/scripts/clinical/extract_tcga_slides.py

Only input and output file locations can be configured via CLI.
"""

from __future__ import annotations

import argparse
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch

DEFAULT_PT_DIR = Path("/Users/mattbocc/uhn/TCGA_SARC/features/TCGA_SARC_pt")
DEFAULT_SAMPLE_CSV = Path("extraction/data/proc/clinical/TCGA_SARC/clinical_sample.csv")
DEFAULT_OUT_DIR = Path("extraction/data/proc/clinical/TCGA_SARC")
MODEL_NAME = "prov-gigapath"

PT_SUFFIX = ".csv.pt"
REQUIRED_SAMPLE_COLUMNS = {"id"}


# -------------------------------------------------------------------------
# CLI args (only input and output file locations)
# -------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract clinical slide, tile, and embedding CSVs from TCGA-SARC .pt feature files."
    )
    parser.add_argument("--pt-dir", type=Path, default=DEFAULT_PT_DIR, help="Path to .pt feature files directory")
    parser.add_argument("--sample-csv", type=Path, default=DEFAULT_SAMPLE_CSV, help="Path to clinical_sample.csv")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for CSVs")
    return parser.parse_args()


# -------------------------------------------------------------------------
# Sample resolution
# -------------------------------------------------------------------------


def sample_type_prefix(sample_id: str) -> str:
    """Patient + sample-type, dropping the vial letter, e.g.
    'TCGA-3B-A9HI-01A' -> 'TCGA-3B-A9HI-01' (15 chars)."""
    return sample_id[:-1]


def slide_barcode_from_filename(pt_path: Path) -> str:
    name = pt_path.name
    if not name.endswith(PT_SUFFIX):
        raise ValueError(f"Expected filename ending in {PT_SUFFIX!r}, got {name!r}")
    return name[:-len(PT_SUFFIX)]


def build_sample_lookup(sample_csv: Path) -> dict[str, list[str]]:
    """Maps sample-type prefix -> sorted list of matching clinical_sample.id values."""
    sample_df = pd.read_csv(sample_csv, dtype=str)
    missing = REQUIRED_SAMPLE_COLUMNS - set(sample_df.columns)
    if missing:
        raise ValueError(f"{sample_csv} is missing expected columns: {missing}")

    lookup: dict[str, list[str]] = defaultdict(list)
    for sample_id in sample_df["id"].dropna():
        lookup[sample_type_prefix(sample_id)].append(sample_id)
    for prefix in lookup:
        lookup[prefix].sort()
    return lookup


def resolve_sample_id(
    slide_barcode: str, lookup: dict[str, list[str]]
) -> tuple[str | None, bool]:
    """Returns (resolved clinical_sample.id or None, was_ambiguous)."""
    # slide_barcode e.g. "TCGA-3B-A9HI-01Z-00-DX1" -> prefix "TCGA-3B-A9HI-01" (15 chars)
    prefix = slide_barcode[:15]
    candidates = lookup.get(prefix, [])
    if not candidates:
        return None, False
    # Deterministic tie-break: alphabetically-first vial letter.
    return candidates[0], len(candidates) > 1


# -------------------------------------------------------------------------
# Extraction
# -------------------------------------------------------------------------


def extract_all_slide_data(
    pt_dir: Path,
    sample_lookup: dict[str, list[str]],
    out_dir: Path,
) -> None:
    pt_paths = sorted(pt_dir.glob(f"*{PT_SUFFIX}"))
    print(f"Found {len(pt_paths)} .pt files in {pt_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    slide_csv_path = out_dir / "clinical_slide.csv"
    tile_csv_path = out_dir / "clinical_tile.csv"
    embed_csv_path = out_dir / "clinical_embedding.csv"

    if tile_csv_path.exists():
        tile_csv_path.unlink()
    if embed_csv_path.exists():
        embed_csv_path.unlink()

    slide_rows: list[dict] = []
    unresolved: list[str] = []
    ambiguous: list[tuple[str, str]] = []

    total_tiles_count = 0
    total_embeds_count = 0

    first_tile_write = True
    first_embed_write = True

    start_time = time.time()

    for idx, pt_path in enumerate(pt_paths, start=1):
        slide_id = slide_barcode_from_filename(pt_path)
        slide_barcode = slide_id.split(".")[0]

        sample_id, was_ambiguous = resolve_sample_id(slide_barcode, sample_lookup)
        if sample_id is None:
            unresolved.append(slide_id)
            continue
        if was_ambiguous:
            ambiguous.append((slide_id, sample_id))

        tensors = torch.load(pt_path, map_location="cpu", weights_only=False)
        tile_embeds_t = tensors["tile_embeds"]
        coords_t = tensors.get("coords")

        n_tiles, embedding_dim = tile_embeds_t.shape

        # 1. Slide Metadata (ClinicalSlide)
        slide_rows.append(
            {
                "id": slide_id,
                "sample_id": sample_id,
                "n_tiles": int(n_tiles),
                "embedding_dim": int(embedding_dim),
            }
        )

        # 2. Tile Coordinates (ClinicalTile)
        if coords_t is not None:
            coords_np = coords_t.numpy().astype(np.float32)
            n_coords = len(coords_np)
            tile_df = pd.DataFrame(
                {
                    "slide_id": [slide_id] * n_coords,
                    "tile_index": np.arange(n_coords, dtype=np.int32),
                    "x": coords_np[:, 0],
                    "y": coords_np[:, 1],
                }
            )
            tile_df.to_csv(
                tile_csv_path,
                mode="a",
                header=first_tile_write,
                index=False,
            )
            first_tile_write = False
            total_tiles_count += n_coords

        # 3. Embeddings (ClinicalEmbedding: slide_mean + tile embeddings)
        embed_rows = []
        tile_embeds_np = tile_embeds_t.numpy().astype(np.float32)

        # Slide Mean Embedding (kind="slide_mean", tile_index=None)
        mean_vec = tile_embeds_np.mean(axis=0)
        embed_rows.append(
            {
                "slide_id": slide_id,
                "tile_index": None,
                "model_name": MODEL_NAME,
                "kind": "slide_mean",
                "embedding": mean_vec.tobytes().hex(),
            }
        )

        # Tile-level Embeddings (kind="tile", tile_index=0..n_tiles-1)
        for t_idx in range(n_tiles):
            embed_rows.append(
                {
                    "slide_id": slide_id,
                    "tile_index": t_idx,
                    "model_name": MODEL_NAME,
                    "kind": "tile",
                    "embedding": tile_embeds_np[t_idx].tobytes().hex(),
                }
            )

        embed_df = pd.DataFrame(embed_rows)
        embed_df.to_csv(
            embed_csv_path,
            mode="a",
            header=first_embed_write,
            index=False,
        )
        first_embed_write = False
        total_embeds_count += len(embed_rows)

        if idx % 25 == 0 or idx == len(pt_paths):
            elapsed = time.time() - start_time
            print(
                f"  Processed {idx}/{len(pt_paths)} slides ({idx/len(pt_paths)*100:.1f}%) "
                f"[{elapsed:.1f}s]...",
                end="\r",
                flush=True,
            )

    print()

    # Save clinical_slide.csv
    slide_df = pd.DataFrame(slide_rows, columns=["id", "sample_id", "n_tiles", "embedding_dim"])
    slide_df = slide_df.sort_values("id").reset_index(drop=True)
    slide_df.to_csv(slide_csv_path, index=False)

    if unresolved:
        print(
            f"WARNING: {len(unresolved)} slide(s) could not be resolved to a "
            "clinical_sample.id and were skipped:"
        )
        for s_id in unresolved:
            print(f"  - {s_id}")

    if ambiguous:
        print(
            f"WARNING: {len(ambiguous)} slide(s) matched more than one "
            "clinical_sample.id with the same sample type; resolved to the "
            "alphabetically-first vial letter -- review these by hand:"
        )
        for s_id, sam_id in ambiguous:
            print(f"  - {s_id} -> {sam_id}")

    total_time = time.time() - start_time
    print("=" * 70)
    print(" Slide Extraction Complete")
    print("=" * 70)
    print(f"Wrote {slide_csv_path}: {len(slide_df)} rows")
    print(f"Wrote {tile_csv_path}: {total_tiles_count:,} rows")
    print(f"Wrote {embed_csv_path}: {total_embeds_count:,} rows")
    print(f"Total Elapsed Time: {total_time:.2f}s")
    print("=" * 70)


def main() -> None:
    args = parse_args()

    if not args.pt_dir.is_dir():
        raise SystemExit(f"--pt-dir does not exist or is not a directory: {args.pt_dir}")
    if not args.sample_csv.is_file():
        raise SystemExit(
            f"--sample-csv not found: {args.sample_csv}. "
            "Run extract_tcga.R first to produce clinical_sample.csv."
        )

    sample_lookup = build_sample_lookup(args.sample_csv)

    extract_all_slide_data(
        pt_dir=args.pt_dir,
        sample_lookup=sample_lookup,
        out_dir=args.out_dir,
    )


if __name__ == "__main__":
    main()
