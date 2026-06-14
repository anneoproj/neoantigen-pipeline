import argparse
import os
import re
from pathlib import Path

import pandas as pd

from scoring.expression_score import add_expression_score
from scoring.hla import normalize_hla, normalize_peptide
from scoring.mhcflurry import MHCFLURRY_COLUMNS, parse_mhcflurry_output
from scoring.mixmhcpred import MIXMHCPRED_COLUMNS, parse_mixmhcpred_output
from scoring.netchop import netchop_score_for_peptide, parse_netchop
from scoring.netmhcpan import parse_netmhcpan_output
from scoring.peptide_score import add_peptide_score


def peptide_column(df):
    if "Tumor_Peptide" in df.columns:
        return "Tumor_Peptide"
    if "Peptide" in df.columns:
        return "Peptide"
    raise ValueError("Input CSV must contain Peptide or Tumor_Peptide")


def transcript_id(value):
    if pd.isna(value):
        return None
    return str(value).split(".")[0]


def netchop_index(netchop_dir):
    index = {}
    for root, _, files in os.walk(netchop_dir, followlinks=True):
        if "netchop.txt" not in files:
            continue
        path = Path(root) / "netchop.txt"
        match = re.search(r"(ENST\d+)", str(path.parent))
        if match:
            index[match.group(1)] = parse_netchop(path)
    return index


def expression_lookup(expression_file):
    df = pd.read_csv(expression_file)
    if df.empty:
        return pd.DataFrame(columns=["Gene", "gene_abundance"])

    first_row = df.iloc[0]
    rows = []
    for gene, value in first_row.items():
        if str(gene).startswith("Unnamed:"):
            continue
        rows.append((gene, value))

    expr = pd.DataFrame(rows, columns=["Gene", "gene_abundance"])
    expr["gene_abundance"] = pd.to_numeric(expr["gene_abundance"], errors="coerce").fillna(0.0)
    return expr


def add_netchop_scores(df, netchop_dir):
    pep_col = peptide_column(df)
    scores_by_transcript = netchop_index(netchop_dir)
    scores = []

    for _, row in df.iterrows():
        enst = transcript_id(row.get("Transcript_ID"))
        netchop_scores = scores_by_transcript.get(enst, {})
        peptide = str(row.get(pep_col, "")).replace("-", "")
        length = row.get("Peptide_Length", len(peptide))
        start = row.get("Tumor_start_position")
        scores.append(netchop_score_for_peptide(start, length, netchop_scores))

    df = df.copy()
    if "NetChop_score" in df.columns:
        df = df.drop(columns=["NetChop_score"])
    df["NetChop_Tumor"] = scores
    return df


def parse_netmhcpan_files(netmhcpan_files):
    frames = []
    for netmhcpan_file in netmhcpan_files:
        with open(netmhcpan_file) as handle:
            df = parse_netmhcpan_output(handle.read())
        if not df.empty:
            frames.append(df)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)


def add_ic50_scores(df, netmhcpan_files):
    pep_col = peptide_column(df)
    df_netmhc = parse_netmhcpan_files(netmhcpan_files)

    if df_netmhc.empty:
        df = df.copy()
        df["IC50"] = pd.NA
        return df

    df = df.copy()
    df_netmhc = df_netmhc.copy()

    df["_merge_hla"] = df["HLA"].map(normalize_hla)
    df["_merge_peptide"] = df[pep_col].map(normalize_peptide)
    df_netmhc["_merge_hla"] = df_netmhc["MHC"].map(normalize_hla)
    df_netmhc["_merge_peptide"] = df_netmhc["Peptide"].map(normalize_peptide)
    df_netmhc = df_netmhc.drop_duplicates(["_merge_hla", "_merge_peptide"])

    merged = pd.merge(
        df,
        df_netmhc[["_merge_hla", "_merge_peptide", "Affinity(nM)"]],
        on=["_merge_hla", "_merge_peptide"],
        how="left",
    )

    merged = merged.rename(columns={"Affinity(nM)": "IC50"})
    merged["IC50"] = pd.to_numeric(merged["IC50"], errors="coerce")
    return merged.drop(columns=["_merge_hla", "_merge_peptide"])


def add_prediction_scores(df, predictions, prediction_columns):
    pep_col = peptide_column(df)
    df = df.copy()

    for column in prediction_columns:
        if column in df.columns:
            df = df.drop(columns=[column])

    if predictions.empty:
        for column in prediction_columns:
            df[column] = pd.NA
        return df

    predictions = predictions.copy()
    df["_merge_hla"] = df["HLA"].map(normalize_hla)
    df["_merge_peptide"] = df[pep_col].map(normalize_peptide)
    predictions["_merge_hla"] = predictions["HLA"].map(normalize_hla)
    predictions["_merge_peptide"] = predictions["Peptide"].map(normalize_peptide)
    predictions = predictions.drop_duplicates(["_merge_hla", "_merge_peptide"])

    available_columns = [
        column for column in prediction_columns if column in predictions.columns
    ]
    merged = pd.merge(
        df,
        predictions[["_merge_hla", "_merge_peptide", *available_columns]],
        on=["_merge_hla", "_merge_peptide"],
        how="left",
    )

    for column in prediction_columns:
        if column not in merged.columns:
            merged[column] = pd.NA
        merged[column] = pd.to_numeric(merged[column], errors="coerce")

    return merged.drop(columns=["_merge_hla", "_merge_peptide"])


def add_mhcflurry_scores(df, mhcflurry_file):
    return add_prediction_scores(
        df,
        parse_mhcflurry_output(mhcflurry_file),
        MHCFLURRY_COLUMNS,
    )


def add_mixmhcpred_scores(df, mixmhcpred_file):
    return add_prediction_scores(
        df,
        parse_mixmhcpred_output(mixmhcpred_file),
        MIXMHCPRED_COLUMNS,
    )


def add_mhcnuggets_scores(df, prediction_file):
    predictions = pd.read_csv(prediction_file)
    if predictions.empty:
        return add_prediction_scores(
            df,
            predictions,
            ["MHCnuggets_IC50"],
        )

    if "HLA" not in predictions.columns:
        for col in predictions.columns:
            if col.lower() == "hla":
                predictions = predictions.rename(columns={col: "HLA"})
                break
    if "Peptide" not in predictions.columns:
        for col in predictions.columns:
            if col.lower() == "peptide":
                predictions = predictions.rename(columns={col: "Peptide"})
                break

    if "MHCnuggets_IC50" not in predictions.columns:
        for col in predictions.columns:
            if col.lower() in {"mhcnuggets_ic50", "mhcnuggets"}:
                predictions = predictions.rename(columns={col: "MHCnuggets_IC50"})
                break

    return add_prediction_scores(
        df,
        predictions,
        ["MHCnuggets_IC50"],
    )


def add_pssmhcpan_scores(df, prediction_file):
    predictions = pd.read_csv(prediction_file)
    if predictions.empty:
        return add_prediction_scores(
            df,
            predictions,
            ["PSSMHCpan_IC50"],
        )

    if "HLA" not in predictions.columns:
        for col in predictions.columns:
            if col.lower() == "hla":
                predictions = predictions.rename(columns={col: "HLA"})
                break
    if "Peptide" not in predictions.columns:
        for col in predictions.columns:
            if col.lower() == "peptide":
                predictions = predictions.rename(columns={col: "Peptide"})
                break

    if "PSSMHCpan_IC50" not in predictions.columns:
        for col in predictions.columns:
            if col.lower() in {"pssmhcpan_ic50", "pssmhcpan"}:
                predictions = predictions.rename(columns={col: "PSSMHCpan_IC50"})
                break

    return add_prediction_scores(
        df,
        predictions,
        ["PSSMHCpan_IC50"],
    )


def add_expression_scores(df, expression_file):
    expr = expression_lookup(expression_file)
    df = df.copy()

    if "gene_abundance" in df.columns:
        df = df.drop(columns=["gene_abundance"])
    if "ExpressionScore" in df.columns:
        df = df.drop(columns=["ExpressionScore"])

    merged = df.merge(expr, on="Gene", how="left")
    merged["gene_abundance"] = merged["gene_abundance"].fillna(0.0)
    return add_expression_score(merged)


def move_columns_after(df, columns, anchor):
    existing = [column for column in columns if column in df.columns]
    if not existing or anchor not in df.columns:
        return df

    remaining = [column for column in df.columns if column not in existing]
    anchor_index = remaining.index(anchor) + 1
    ordered = remaining[:anchor_index] + existing + remaining[anchor_index:]
    return df[ordered]


def process_file(
    input_file,
    netchop_dir,
    netmhcpan_files,
    mhcflurry_file,
    mixmhcpred_file,
    expression_file,
    output_file,
):
    df = pd.read_csv(input_file)
    df = add_netchop_scores(df, netchop_dir)
    df = add_ic50_scores(df, netmhcpan_files)
    df = add_mhcflurry_scores(df, mhcflurry_file)
    df = add_mixmhcpred_scores(df, mixmhcpred_file)
    df = add_expression_scores(df, expression_file)
    df = add_peptide_score(df)
    df = move_columns_after(df, [*MHCFLURRY_COLUMNS, *MIXMHCPRED_COLUMNS], "IC50")
    df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--netchop_dir", required=True)
    parser.add_argument("--netmhcpan_files", nargs="+", required=True)
    parser.add_argument("--mhcflurry_file", required=True)
    parser.add_argument("--mixmhcpred_file", required=True)
    parser.add_argument("--expression_file", required=True)
    parser.add_argument("--output_file", required=True)

    args = parser.parse_args()
    process_file(
        input_file=args.input_file,
        netchop_dir=args.netchop_dir,
        netmhcpan_files=args.netmhcpan_files,
        mhcflurry_file=args.mhcflurry_file,
        mixmhcpred_file=args.mixmhcpred_file,
        expression_file=args.expression_file,
        output_file=args.output_file,
    )


if __name__ == "__main__":
    main()
