import argparse
import pandas as pd

from scoring.hla import normalize_hla, normalize_peptide
from scoring.mhcflurry import MHCFLURRY_COLUMNS, parse_mhcflurry_output


def _first_column(df, names):
    normalized = {name.lower(): name for name in df.columns}
    for candidate in names:
        if candidate in normalized:
            return normalized[candidate]
        candidate = candidate.lower()
        for name in df.columns:
            if name.lower() == candidate:
                return name
    return None


def _add_prediction_scores(df, prediction_file, source_column, required_column):
    predictions = pd.read_csv(prediction_file) if prediction_file else pd.DataFrame()
    if predictions.empty:
        out = df.copy()
        out[required_column] = pd.NA
        return out

    hla_col = _first_column(predictions, ["HLA", "hla", "allele"])
    peptide_col = _first_column(
        predictions,
        ["Peptide", "peptide", "Tumor_Peptide", "mutant_seq"],
    )
    value_col = _first_column(predictions, [source_column, source_column.lower()])
    if value_col is None:
        out = df.copy()
        out[required_column] = pd.NA
        return out

    predictions = predictions[[hla_col, peptide_col, value_col]].copy()
    predictions = predictions.rename(columns={hla_col: "HLA", peptide_col: "Peptide", value_col: required_column})
    predictions["HLA"] = predictions["HLA"].map(normalize_hla)
    predictions["Peptide"] = predictions["Peptide"].map(normalize_peptide)
    predictions = predictions.drop_duplicates(["HLA", "Peptide"])

    out = df.copy()
    out["HLA_norm"] = out["HLA"].map(normalize_hla)
    out["Peptide_norm"] = out["Peptide"].map(normalize_peptide)
    predictions["HLA_norm"] = predictions["HLA"]
    predictions["Peptide_norm"] = predictions["Peptide"]
    merged = out.merge(
        predictions[["HLA_norm", "Peptide_norm", required_column]],
        on=["HLA_norm", "Peptide_norm"],
        how="left",
    )
    merged = merged.drop(columns=["HLA_norm", "Peptide_norm"])
    merged[required_column] = pd.to_numeric(merged[required_column], errors="coerce")
    return merged


def add_mhcnuggets_scores(df, prediction_file):
    return _add_prediction_scores(df, prediction_file, "MHCnuggets_IC50", "MHCnuggets_IC50")


def add_pssmhcpan_scores(df, prediction_file):
    return _add_prediction_scores(df, prediction_file, "PSSMHCpan_IC50", "PSSMHCpan_IC50")


def add_mhcflurry_scores(df, mhcflurry_file):
    predictions = parse_mhcflurry_output(mhcflurry_file)
    if predictions.empty:
        for column in MHCFLURRY_COLUMNS:
            df[column] = pd.NA
        return df

    out = df.copy()
    out["HLA_norm"] = out["HLA"].map(normalize_hla)
    out["Peptide_norm"] = out["Peptide"].map(normalize_peptide)
    predictions = predictions.copy()
    predictions["HLA_norm"] = predictions["HLA"].map(normalize_hla)
    predictions["Peptide_norm"] = predictions["Peptide"].map(normalize_peptide)
    predictions = predictions.drop_duplicates(["HLA_norm", "Peptide_norm"])

    merged = out.merge(
        predictions[["HLA_norm", "Peptide_norm", *MHCFLURRY_COLUMNS]],
        on=["HLA_norm", "Peptide_norm"],
        how="left",
    )
    return merged.drop(columns=["HLA_norm", "Peptide_norm"])


def process_existing_file(input_file, mhcflurry_file, mhcnuggets_file, pssmhcpan_file, output_file):
    df = pd.read_csv(input_file)
    df = add_mhcflurry_scores(df, mhcflurry_file)
    df = add_mhcnuggets_scores(df, mhcnuggets_file)
    df = add_pssmhcpan_scores(df, pssmhcpan_file)
    df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--mhcflurry_file", required=True)
    parser.add_argument("--mhcnuggets_file", required=True)
    parser.add_argument("--pssmhcpan_file", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()
    process_existing_file(
        input_file=args.input_file,
        mhcflurry_file=args.mhcflurry_file,
        mhcnuggets_file=args.mhcnuggets_file,
        pssmhcpan_file=args.pssmhcpan_file,
        output_file=args.output_file,
    )


if __name__ == "__main__":
    main()
