#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import tempfile


def read_genome_order(genome_file):
    order = {}
    with open(genome_file) as handle:
        for index, line in enumerate(handle):
            contig = line.rstrip().split("\t")[0]
            order[contig] = index
    return order


def combine_beds(targets, genome_order, combined_bed):
    combined = []
    for feature_type, bed_file, _output in targets:
        with open(bed_file) as handle:
            for line in handle:
                if not line.strip():
                    continue
                chrom, start, end, feature = line.rstrip().split("\t")[:4]
                combined.append((genome_order[chrom], int(start), int(end), chrom, start, end, feature_type, feature))

    combined.sort(key=lambda row: row[:3])
    with open(combined_bed, "w") as handle:
        for _order, _start, _end, chrom, start, end, feature_type, feature in combined:
            handle.write(f"{chrom}\t{start}\t{end}\t{feature_type}|{feature}\n")


def split_coverage(lines, targets):
    handles = {feature_type: open(output, "w") for feature_type, _bed, output in targets}
    try:
        for line in lines:
            fields = line.rstrip().split("\t")
            feature_type, feature = fields[3].split("|", 1)
            handles[feature_type].write("\t".join(fields[:3] + [feature, fields[-1]]) + "\n")
    finally:
        for handle in handles.values():
            handle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--target", nargs=3, action="append", metavar=("TYPE", "BED", "OUTPUT"), required=True)
    args = parser.parse_args()

    genome_order = read_genome_order(args.genome)
    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        combined_bed = tmp.name

    try:
        combine_beds(args.target, genome_order, combined_bed)
        process = subprocess.Popen(
            ["bedtools", "coverage", "-a", combined_bed, "-b", args.bam, "-g", args.genome, "-mean", "-sorted"],
            text=True,
            stdout=subprocess.PIPE,
        )
        split_coverage(process.stdout, args.target)
        returncode = process.wait()
        if returncode != 0:
            sys.exit(returncode)
    finally:
        os.unlink(combined_bed)
