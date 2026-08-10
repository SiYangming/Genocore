#!/usr/bin/env python3
"""Convert VCF to GenoCore CSV format (single-threaded, streaming).

Usage:
    python vcf2genocore_csv.py --vcf input.vcf.gz -o output.csv
"""

import pysam
import pandas as pd
import sys
import argparse
from tqdm import tqdm


GT_MAP = {
    "0/0": 0,
    "0/1": 1, "1/0": 1,
    "1/1": 2,
    "./.": "NA", ".": "NA"
}


def vcf_to_genocore_csv(vcf_path, output_csv):
    vcf = pysam.VariantFile(vcf_path, "r")
    samples = list(vcf.header.samples)
    print(f"Samples detected: {len(samples)}")

    csv_header = [""] + samples

    with open(output_csv, "w", encoding="utf-8", newline="") as f:
        pd.DataFrame([csv_header]).to_csv(f, index=False, header=False)

        for record in tqdm(vcf.fetch(), desc="Converting SNPs"):
            if not record.id:
                continue

            snp_id = record.id
            sample_genotypes = []
            for sample in samples:
                gt = record.samples[sample].get("GT", "./.")
                if isinstance(gt, tuple):
                    gt = "/".join(map(str, gt))
                geno = GT_MAP.get(gt, "NA")
                sample_genotypes.append(geno)

            row_data = [snp_id] + sample_genotypes
            pd.DataFrame([row_data]).to_csv(
                f, index=False, header=False, mode="a", na_rep="NA"
            )

    vcf.close()
    print(f"Done: {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert VCF to GenoCore CSV format (single-threaded)."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--output", "-o", required=True, help="Output CSV for GenoCore")
    args = parser.parse_args()

    try:
        pysam.VariantFile(args.vcf)
    except FileNotFoundError:
        print(f"Error: file not found {args.vcf}")
        sys.exit(1)

    vcf_to_genocore_csv(args.vcf, args.output)
