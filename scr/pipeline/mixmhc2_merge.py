import argparse

import pandas as pd

from scoring.hla import normalize_hla, normalize_peptide
from scoring.mixmhc2 import MIXMHC2_COLUMNS, parse_mixmhc2_output


def _peptide_column(df):
    if "Peptide" in df.columns:
        return "Peptide"
    if "Tumor_Peptide" in df.columns:
        return "Tumor_Peptide"
    raise ValueError("Input CSV must contain Peptide or Tumor_Peptide")


def add_mixmhc2_scores(df, prediction_file):
    predictions = parse_mixmhc2_output(prediction_file)
    if predictions.empty:
        df = df.copy()
        for column in MIXMHC2_COLUMNS:
            df[column] = pd.NA
        return df

    df = df.copy()
    pep_col = _peptide_column(df)
    if "HLA" not in df.columns:
        raise ValueError("Input CSV must contain HLA column")

    df["_merge_hla"] = df["HLA"].map(normalize_hla)
    df["_merge_peptide"] = df[pep_col].map(normalize_peptide)
    predictions["_merge_hla"] = predictions["HLA"].map(normalize_hla)
    predictions["_merge_peptide"] = predictions["Peptide"].map(normalize_peptide)
    predictions = predictions.drop_duplicates(["_merge_hla", "_merge_peptide"])

    merged = pd.merge(
        df,
        predictions[["_merge_hla", "_merge_peptide", *MIXMHC2_COLUMNS]],
        on=["_merge_hla", "_merge_peptide"],
        how="left",
    )

    for column in MIXMHC2_COLUMNS:
        merged[column] = pd.to_numeric(merged[column], errors="coerce")
    return merged.drop(columns=["_merge_hla", "_merge_peptide"])


def process_file(input_file, prediction_file, output_file):
    df = pd.read_csv(input_file)
    merged = add_mixmhc2_scores(df, prediction_file)
    merged.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--prediction_file", required=True)
    parser.add_argument("--output_file", required=True)

    args = parser.parse_args()
    process_file(
        input_file=args.input_file,
        prediction_file=args.prediction_file,
        output_file=args.output_file,
    )


if __name__ == "__main__":
    main()
