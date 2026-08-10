#!/usr/bin/env python3
"""Convert PLINK PED/MAP to GenoCore CSV format.

Usage:
    python plink2genocore.py --ped input.ped --map input.map -o output.csv
"""

import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import sys
import os


def plink_to_genocore(ped_path, map_path, output_csv, batch_size=10000):
    for f in [ped_path, map_path]:
        if not os.path.exists(f):
            print(f"Error: file not found {f}")
            sys.exit(1)
        if os.path.getsize(f) == 0:
            print(f"Error: file {f} is empty")
            sys.exit(1)

    map_df = pd.read_csv(
        map_path, sep=r"\s+", header=None,
        names=["chr", "snp_id", "cm", "pos"],
        usecols=["snp_id"],
        dtype={"snp_id": str}
    )
    snp_ids = map_df["snp_id"].tolist()
    total_snps = len(snp_ids)
    print(f"SNPs loaded: {total_snps}")
    if total_snps == 0:
        print("Error: no SNP IDs in MAP file")
        sys.exit(1)

    print("Loading PED file (may take a while for large files)...")
    ped_df = pd.read_csv(
        ped_path, sep=r"\s+", header=None,
        low_memory=False, dtype=str
    )
    ped_cols = ped_df.shape[1]
    expected_gt_cols = total_snps
    actual_gt_cols = ped_cols - 6
    if actual_gt_cols != expected_gt_cols:
        print(f"Warning: genotype column mismatch (map={expected_gt_cols}, ped={actual_gt_cols})")
        if actual_gt_cols > expected_gt_cols:
            ped_df = ped_df.iloc[:, :6 + expected_gt_cols]
        else:
            for i in range(actual_gt_cols, expected_gt_cols):
                ped_df[f"extra_{i}"] = "."

    sample_names = ped_df.iloc[:, 1].tolist()
    print(f"Samples loaded: {len(sample_names)}")
    if len(sample_names) == 0:
        print("Error: no samples in PED file")
        sys.exit(1)

    gt_columns = ped_df.iloc[:, 6:]
    gt_transposed = gt_columns.T
    gt_transposed.columns = sample_names
    gt_transposed.index = snp_ids[:len(gt_transposed)]
    print(f"Genotype matrix: {gt_transposed.shape} (SNPs x samples)")

    gt_transposed = gt_transposed.replace(".", "NA").astype(str)
    gt_transposed = gt_transposed.dropna(how="all")
    print(f"After filtering empty rows: {len(gt_transposed)} SNPs")

    header = [""] + sample_names
    with open(output_csv, "w", encoding="utf-8", newline="", buffering=1024*1024) as f:
        f.write(",".join(header) + "\n")
        current_snps = len(gt_transposed)
        for i in tqdm(range(0, current_snps, batch_size), desc="Writing GenoCore CSV"):
            end_idx = min(i + batch_size, current_snps)
            batch = gt_transposed.iloc[i:end_idx]
            batch_lines = []
            for snp_id, row in batch.iterrows():
                row_vals = [snp_id] + row.fillna("NA").tolist()
                batch_lines.append(",".join(row_vals))
            if batch_lines:
                f.write("\n".join(batch_lines) + "\n")

    print(f"Done: {os.path.abspath(output_csv)}")
    print("First 5 rows preview:")
    try:
        preview = pd.read_csv(output_csv, nrows=5, header=None)
        print(preview)
    except Exception:
        print("Preview failed but file was written successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert PLINK PED/MAP to GenoCore CSV format."
    )
    parser.add_argument("--ped", required=True, help="PLINK PED file path")
    parser.add_argument("--map", required=True, help="PLINK MAP file path")
    parser.add_argument("--output", "-o", required=True, help="Output CSV for GenoCore")
    parser.add_argument("--batch-size", type=int, default=10000,
                        help="Batch size for writing (default: 10000)")
    args = parser.parse_args()
    plink_to_genocore(args.ped, args.map, args.output, args.batch_size)
