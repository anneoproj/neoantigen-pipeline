import argparse
import os
from pathlib import Path

import pandas as pd

from scoring.hla import normalize_hla, normalize_peptide


def _split_alleles(value):
    if pd.isna(value):
        return []
    items = []
    for chunk in str(value).replace(",", ";").replace("/", ";").split(";"):
        chunk = chunk.strip()
        if chunk:
            items.append(normalize_hla(chunk))
    return [item for item in dict.fromkeys(items) if item]


def _write_batches(merged_df, output_dir, chunk_size=2):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    peptides = merged_df["Peptide"].map(normalize_peptide).dropna().unique().tolist()
    if not peptides:
        return

    hlas = merged_df["HLA"].dropna().map(normalize_hla).dropna().unique().tolist()
    if not hlas:
        return

    chunk_size = int(chunk_size)
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")

    chunk_count = (len(peptides) + chunk_size - 1) // chunk_size
    width = max(5, len(str(chunk_count)))
    for hla in hlas:
        for index, start in enumerate(range(0, len(peptides), chunk_size), start=1):
            batch = peptides[start : start + chunk_size]
            batch_name = f"HLA-{hla}__chunk_{index:0{width}d}.txt"
            batch_path = output_dir / batch_name
            with open(batch_path, "w") as handle:
                for peptide in batch:
                    handle.write(f"{peptide}\n")


def prepare_neopep(input_file, output_dir=".", chunk_size=1000, chunksize=3):
    df = pd.read_csv(input_file, sep=None, engine="python", dtype=str)
    if df.empty:
        raise ValueError(f"No data in {input_file}")

    peptide_col = "Peptide" if "Peptide" in df.columns else "mutant_seq"
    if peptide_col not in df.columns:
        raise ValueError("Input must contain Peptide or mutant_seq")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    expanded_rows = []
    for _, row in df.iterrows():
        hla_values = []
        if "mutant_best_alleles" in df.columns:
            hla_values.extend(_split_alleles(row.get("mutant_best_alleles")))
        if "mutant_other_significant_alleles_netMHCpan" in df.columns:
            hla_values.extend(
                _split_alleles(row.get("mutant_other_significant_alleles_netMHCpan"))
            )
        hla_values = [hla for hla in dict.fromkeys(hla_values) if hla]
        if not hla_values:
            continue

        for hla in hla_values:
            data = row.to_dict()
            data["HLA"] = hla
            data["Peptide"] = normalize_peptide(row[peptide_col])
            expanded_rows.append(data)

    if not expanded_rows:
        raise ValueError("No valid HLA expansion entries in neopep input")

    expanded = pd.DataFrame(expanded_rows).drop_duplicates()
    expanded = expanded.dropna(how="all", axis=1)

    expanded_path = output_dir / "Neopep_data_org_expanded.csv"
    expanded.to_csv(expanded_path, index=False)

    mhcflurry_input = output_dir / "Neopep_data_org_mhcflurry_input.csv"
    mhcflurry_input_df = pd.DataFrame({
        "HLA": expanded["HLA"],
        "Peptide": expanded["Peptide"],
    })
    mhcflurry_input_df.to_csv(mhcflurry_input, index=False)

    batches_dir = output_dir / "Neopep_data_org_netmhcpan_batches"
    _write_batches(expanded, batches_dir, chunk_size=chunk_size)
    _write_batches(expanded, batches_dir, chunk_size=chunksize)

    return expanded_path, mhcflurry_input, batches_dir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_dir", default=".")
    parser.add_argument("--chunk_size", type=int, default=2)
    parser.add_argument("--chunksize", type=int, default=3)
    args = parser.parse_args()
    prepare_neopep(
        input_file=args.input_file,
        output_dir=args.output_dir,
        chunk_size=args.chunk_size,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()
