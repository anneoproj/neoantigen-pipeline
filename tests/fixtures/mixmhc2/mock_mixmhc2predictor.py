#!/usr/bin/env python3

import argparse
import csv


def parse_pairs(input_path, output_path):
    pairs = []
    with open(input_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            hla, peptide = row[:2]
            pairs.append((hla.strip(), peptide.strip()))

    with open(output_path, "w") as handle:
        handle.write("HLA,Peptide,Score,PercentileRank\n")
        for index, (hla, peptide) in enumerate(pairs, start=1):
            score = 0.1 * index
            rank = 100 - index
            handle.write(f"{hla},{peptide},{score:.3f},{rank:.3f}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    parse_pairs(args.input, args.output)


if __name__ == "__main__":
    main()
