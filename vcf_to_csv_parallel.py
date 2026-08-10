#!/usr/bin/env python3
"""Convert VCF to GenoCore CSV format (multi-threaded, parallel by chunks).

Usage:
    python vcf_to_csv_parallel.py --vcf input.vcf.gz -o output.csv -w 64
"""

import pysam
import sys
import csv
import argparse
import multiprocessing as mp
import os
from tqdm import tqdm
from itertools import groupby


GT_MAP = {
    "0/0": 0,
    "0/1": 1, "1/0": 1,
    "1/1": 2,
    "./.": "NA", ".": "NA"
}


def process_chunk(task):
    regions, samples, gt_map, temp_file, vcf_file = task
    vcf = pysam.VariantFile(vcf_file, "r")
    with open(temp_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        for chrom, start, end in regions:
            for record in vcf.fetch(chrom, start, end):
                if not record.id:
                    continue
                snp_id = record.id
                genotypes = []
                for sample in samples:
                    gt = record.samples[sample].get("GT")
                    if gt is None or None in gt:
                        geno = "NA"
                    else:
                        gt_str = "/".join(map(str, gt))
                        geno = gt_map.get(gt_str, "NA")
                    genotypes.append(geno)
                writer.writerow([snp_id] + genotypes)
    vcf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert VCF to GenoCore CSV format (parallel)."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--output", "-o", required=True, help="Output CSV for GenoCore")
    parser.add_argument("--workers", "-w", type=int, default=None,
                        help="Number of parallel workers (default: CPU count)")
    args = parser.parse_args()

    # Validate VCF
    try:
        with pysam.VariantFile(args.vcf, "r") as vcf:
            samples = list(vcf.header.samples)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(f"Samples detected: {len(samples)}")
    header = [""] + samples

    # Scan variant positions
    print("Step 1: scanning VCF for variant positions...")
    variant_positions = []
    with pysam.VariantFile(args.vcf, "r") as vcf:
        for record in tqdm(vcf.fetch(), desc="Scanning"):
            if record.id:
                variant_positions.append((record.chrom, record.pos))

    total = len(variant_positions)
    print(f"Found {total:,} valid SNP sites")
    if total == 0:
        print("No valid sites, exiting")
        sys.exit(0)

    # Chunking
    num_workers = args.workers or mp.cpu_count()
    chunk_size = max(1, total // num_workers)
    chunks = []
    for i in range(0, total, chunk_size):
        chunk = variant_positions[i:i + chunk_size]
        regions = []
        for chrom, group in groupby(chunk, key=lambda x: x[0]):
            g = list(group)
            start = g[0][1]
            end = g[-1][1] + 1
            regions.append((chrom, start, end))
        chunks.append(regions)

    temp_files = [f"temp_chunk_{i:04d}.csv" for i in range(len(chunks))]

    print(f"Step 2: converting with {num_workers} workers...")
    tasks = [(regions, samples, GT_MAP, temp_file, args.vcf)
             for regions, temp_file in zip(chunks, temp_files)]

    with mp.Pool(num_workers) as pool:
        list(tqdm(pool.imap_unordered(process_chunk, tasks),
                  total=len(tasks), desc="Parallel conversion"))

    # Merge
    print("Step 3: merging temporary files...")
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for tf in temp_files:
            with open(tf, "r", encoding="utf-8") as src:
                f.writelines(src)
            os.remove(tf)

    print(f"Done: {args.output}")
    print(f"Total sites: {total:,}  Workers: {num_workers}")
