#!/usr/bin/env python
# -*- coding: latin-1 -*-
import argparse
import itertools
import math
import os
import resource
import tempfile
from os.path import basename


def get_open_file_budget():
    try:
        soft_limit, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    except (OSError, ValueError):
        return 100
    if soft_limit == resource.RLIM_INFINITY:
        return 200
    return max(1, min(200, soft_limit - 50))


def choose_batch_size(file_count, requested_batch_size=None):
    budget = get_open_file_budget()
    if requested_batch_size is not None:
        return max(1, min(requested_batch_size, budget))

    batch_size = min(200, budget)
    if file_count:
        min_batch_for_merge = math.ceil(file_count / max(1, budget - 1))
        batch_size = max(batch_size, min_batch_for_merge)
    return max(1, min(batch_size, budget))


def parse_cov_line(path, line_number, line):
    fields = line.rstrip().split("\t")
    if len(fields) < 5:
        raise ValueError(f"{path}:{line_number} has fewer than five columns")
    return fields[3], str(float(fields[4]))


def write_feature_template(first_cov_file):
    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        tmp.write("cov\n")
        with open(first_cov_file) as handle:
            for line_number, line in enumerate(handle, start=1):
                feature, _cov = parse_cov_line(first_cov_file, line_number, line)
                tmp.write(feature + "\n")
        return tmp.name


def write_batch(batch_files, batch_samples, feature_file):
    handles = [open(path) for path in batch_files]
    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        tmp.write("\t".join(batch_samples) + "\n")
        try:
            with open(feature_file) as feature_handle:
                next(feature_handle)
                rows = itertools.zip_longest(feature_handle, *handles)
                for line_number, row in enumerate(rows, start=1):
                    expected_feature = row[0]
                    if expected_feature is None:
                        extra_files = [path for path, value in zip(batch_files, row[1:]) if value is not None]
                        raise ValueError(
                            "Coverage file has more rows than the first coverage file: " + ", ".join(extra_files)
                        )

                    expected_feature = expected_feature.rstrip("\n")
                    cov_values = []
                    for path, line in zip(batch_files, row[1:]):
                        if line is None:
                            raise ValueError(f"{path} has fewer rows than the first coverage file")
                        feature, cov = parse_cov_line(path, line_number, line)
                        if feature != expected_feature:
                            raise ValueError(
                                f"{path}:{line_number} has feature {feature!r}, expected {expected_feature!r}"
                            )
                        cov_values.append(cov)
                    tmp.write("\t".join(cov_values) + "\n")
        finally:
            for handle in handles:
                handle.close()
        return tmp.name


def write_output(feature_file, batch_files, output):
    handles = [open(feature_file)] + [open(path) for path in batch_files]
    try:
        with open(output, "w") as out:
            for row in itertools.zip_longest(*handles):
                if any(value is None for value in row):
                    raise ValueError("Temporary coverage batch files have different row counts")
                out.write("\t".join(value.rstrip("\n") for value in row) + "\n")
    finally:
        for handle in handles:
            handle.close()


def collate_coverage(cov_files, get_sample_name, output, batch_size=None):
    sample_to_cov = {}
    for path in cov_files:
        sample = get_sample_name(path)
        if sample in sample_to_cov:
            raise ValueError(f"Duplicate sample name derived from coverage files: {sample}")
        sample_to_cov[sample] = path

    sorted_samples = sorted(sample_to_cov.keys())
    sorted_cov_files = [sample_to_cov[sample] for sample in sorted_samples]
    batch_size = choose_batch_size(len(sorted_cov_files), batch_size)
    temp_files = []
    try:
        feature_file = write_feature_template(sorted_cov_files[0])
        temp_files.append(feature_file)
        batch_outputs = []
        for index in range(0, len(sorted_cov_files), batch_size):
            batch_files = sorted_cov_files[index:index + batch_size]
            batch_samples = sorted_samples[index:index + batch_size]
            batch_output = write_batch(batch_files, batch_samples, feature_file)
            batch_outputs.append(batch_output)
            temp_files.append(batch_output)
        write_output(feature_file, batch_outputs, output)
    finally:
        for path in temp_files:
            if os.path.exists(path):
                os.unlink(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", required=True, help="output file name")
    parser.add_argument("-s", required=False,
                        help="suffix to remove to get the name of the sample", default="")
    parser.add_argument("-l", nargs="+", required=True,
                        help="list of coverage files")
    parser.add_argument("--batch-size", type=int, help="maximum number of coverage files to open at once")
    args = parser.parse_args()

    if args.s:
        get_sample_name = lambda x: basename(x).replace(args.s, "")
    else:
        get_sample_name = lambda x: ".".join(basename(x).split('.')[:-2])

    collate_coverage(args.l, get_sample_name, args.o, args.batch_size)
