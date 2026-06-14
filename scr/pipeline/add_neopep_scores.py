import argparse
import pandas as pd

from scoring.hla import normalize_hla, normalize_peptide
from scoring.mhcflurry import parse_mhcflurry_output, MHCFLURRY_COLUMNS
from scoring.netmhcpan import parse_netmhcpan_output


def _first_column(df, names):
    normalized = {name.lower(): name for name in df.columns}
    for candidate in names:
        lc = candidate.lower()
        if lc in normalized:
            return normalized[lc]
        for name in df.columns:
            if name.lower() == lc:
                return name
    return None


def _read_netmhcpan_files(netmhcpan_files):
    frames = []
    for path in netmhcpan_files:
        with open(path) as handle:
            frame = parse_netmhcpan_output(handle.read())
        if not frame.empty:
            frame = frame.rename(columns={"MHC": "HLA"})[["HLA", "Peptide", "Affinity(nM)"]]
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["HLA", "Peptide", "Affinity(nM)"])
    return pd.concat(frames, ignore_index=True)


def _merge_scores(df, predictions, source_col, target_col):
    if predictions.empty:
        df[target_col] = pd.NA
        return df

    out = df.copy()
    out["HLA_norm"] = out["HLA"].map(normalize_hla)
    out["Peptide_norm"] = out["Peptide"].map(normalize_peptide)
    predictions = predictions.copy()
    if "HLA" not in predictions.columns or source_col not in predictions.columns:
        out[target_col] = pd.NA
        return out.drop(columns=["HLA_norm", "Peptide_norm"])

    predictions["HLA_norm"] = predictions["HLA"].map(normalize_hla)
    predictions["Peptide_norm"] = predictions["Peptide"].map(normalize_peptide)
    predictions = predictions.rename(columns={source_col: target_col})
    predictions = predictions.drop_duplicates(["HLA_norm", "Peptide_norm"])
    merged = out.merge(
        predictions[["HLA_norm", "Peptide_norm", target_col]],
        on=["HLA_norm", "Peptide_norm"],
        how="left",
    )
    merged[target_col] = pd.to_numeric(merged[target_col], errors="coerce")
    return merged.drop(columns=["HLA_norm", "Peptide_norm"])


def add_mhcnuggets_scores(df, prediction_file):
    predictions = pd.read_csv(prediction_file)
    return _merge_scores(
        df,
        predictions.rename(columns=_column_renames(predictions, [
            ("MHCnuggets_IC50", "MHCnuggets_IC50"),
            ("ic50", "MHCnuggets_IC50"),
        ]),
        source_col="MHCnuggets_IC50",
        target_col="MHCnuggets_IC50",
    )


def add_pssmhcpan_scores(df, prediction_file):
    predictions = pd.read_csv(prediction_file)
    return _merge_scores(
        df,
        predictions.rename(columns=_column_renames(predictions, [
            ("PSSMHCpan_IC50", "PSSMHCpan_IC50"),
            ("ic50", "PSSMHCpan_IC50"),
        ]),
        source_col="PSSMHCpan_IC50",
        target_col="PSSMHCpan_IC50",
    )


def _column_renames(df, pairs):
    normalized = {column.lower(): column for column in df.columns}
    mapping = {}
    for source, target in pairs:
        source_lower = source.lower()
        if source_lower in normalized:
            mapping[normalized[source_lower]] = target
    return mapping


def _add_mhcflurry(df, mhcflurry_file):
    out = df.copy()
    predictions = parse_mhcflurry_output(mhcflurry_file)
    if predictions.empty:
        for column in MHCFLURRY_COLUMNS:
            out[column] = pd.NA
        return out

    out["HLA_norm"] = out["HLA"].map(normalize_hla)
    out["Peptide_norm"] = out["Peptide"].map(normalize_peptide)
    predictions["HLA_norm"] = predictions["HLA"].map(normalize_hla)
    predictions["Peptide_norm"] = predictions["Peptide"].map(normalize_peptide)
    merged = out.merge(
        predictions[["HLA_norm", "Peptide_norm", *MHCFLURRY_COLUMNS]],
        on=["HLA_norm", "Peptide_norm"],
        how="left",
    )
    return merged.drop(columns=["HLA_norm", "Peptide_norm"])


def add_ic50_scores(df, netmhcpan_files):
    predictions = _read_netmhcpan_files(netmhcpan_files)
    if predictions.empty:
        df["IC50"] = pd.NA
        return df
    predictions = predictions.rename(columns={"Affinity(nM)": "IC50"})
    return _merge_scores(df, predictions, "IC50", "IC50")


def process_neopep_file(input_file, netmhcpan_files, mhcflurry_file, mhcnuggets_file, pssmhcpan_file, output_file):
    df = pd.read_csv(input_file)
    if "Peptide" not in df.columns and "mutant_seq" in df.columns:
        df["Peptide"] = df["mutant_seq"]

    if "HLA" not in df.columns:
        raise ValueError("Input CSV must contain HLA column")

    df = add_ic50_scores(df, netmhcpan_files)
    df = _add_mhcflurry(df, mhcflurry_file)
    df = add_mhcnuggets_scores(df, mhcnuggets_file)
    df = add_pssmhcpan_scores(df, pssmhcpan_file)
    df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--netmhcpan_files", nargs="+", required=True)
    parser.add_argument("--mhcflurry_file", required=True)
    parser.add_argument("--mhcnuggets_file", required=True)
    parser.add_argument("--pssmhcpan_file", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()
    process_neopep_file(
        input_file=args.input_file,
        netmhcpan_files=args.netmhcpan_files,
        mhcflurry_file=args.mhcflurry_file,
        mhcnuggets_file=args.mhcnuggets_file,
        pssmhcpan_file=args.pssmhcpan_file,
        output_file=args.output_file,
    )


if __name__ == "__main__":
    main()
