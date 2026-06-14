import argparse
import re
from pathlib import Path

import pandas as pd

from scoring.hla import hla_to_mixmhc2, normalize_hla, normalize_peptide


def _normalize_column(name):
    return re.sub(r"\W+", "", str(name).strip().lower())


def _load_table(input_file):
    input_path = Path(input_file)
    suffix = input_path.suffix.lower()

    if suffix == ".txt":
        with open(input_file) as handle:
            lines = [line.strip() for line in handle if line.strip()]

        if not lines:
            return pd.DataFrame(columns=["Peptide"])

        sample_text = "\n".join(lines[:5])
        if "," in sample_text or "\t" in sample_text or ";" in sample_text:
            return pd.read_csv(input_file, sep=None, engine="python", dtype=str)
        return pd.DataFrame({"Peptide": lines}, dtype=str)

    return pd.read_csv(input_file, sep=None, engine="python", dtype=str)


def _read_hla_file(hla_file):
    hla = []
    with open(hla_file) as handle:
        for line in handle:
            for value in re.split(r"[,;\s]+", line.strip()):
                normalized = normalize_hla(value)
                if normalized:
                    hla.append(normalized)
    return list(dict.fromkeys(hla))


def _pick_column(columns, candidates):
    normalized_candidates = {
        re.sub(r"\W+", "", str(candidate).strip().lower())
        for candidate in candidates
    }
    for candidate in candidates:
        normalized_candidate = re.sub(r"\W+", "", str(candidate).strip().lower())
        if normalized_candidate in columns:
            return columns[normalized_candidate]
    for normalized, original in columns.items():
        if normalized in normalized_candidates:
            return original
    return None


def _split_hlas(value):
    if pd.isna(value):
        return []

    text = str(value).strip()
    if not text:
        return []

    return list(dict.fromkeys([
        normalize_hla(part)
        for part in re.split(r"[;,/\s]+", text)
        if normalize_hla(part)
    ]))


def _ensure_input_source_column(row, prepared_columns, input_source):
    row = row.copy()
    if "InputSource" not in prepared_columns:
        row["InputSource"] = input_source
        return row
    if pd.isna(row.get("InputSource")) or str(row.get("InputSource", "")).strip() == "":
        row["InputSource"] = input_source
        return row

    return row


def prepare_mixmhc2_inputs(
    input_file,
    output_file,
    pairs_file,
    hla_dir,
    default_sample=None,
    input_source="custom_peptides",
):
    if not hla_dir:
        raise ValueError("hla_dir is required")

    source = _load_table(input_file)
    if source.empty:
        raise ValueError(f"No data found in {input_file}")

    columns = { _normalize_column(name): name for name in source.columns }
    peptide_column = _pick_column(
        columns,
        {"peptide", "tumorpeptide", "mutant_seq", "sequence"},
    )
    if not peptide_column:
        raise ValueError(f"Input must contain a Peptide column: {input_file}")

    sample_column = _pick_column(
        columns,
        {"sample", "sampleid", "sample_id", "patientid", "patient"},
    )
    hla_column = _pick_column(columns, {"hla", "hla_allele", "h2", "mhc"})
    source_columns = set(source.columns)

    prepared_sample = default_sample.strip() if default_sample else ""
    prepared = []
    sample_cache = {}

    for _, row in source.iterrows():
        peptide = normalize_peptide(row.get(peptide_column))
        if not peptide:
            continue

        row_sample = row.get(sample_column) if sample_column else prepared_sample
        sample = str(row_sample).strip() if not pd.isna(row_sample) else prepared_sample
        if not sample:
            sample = prepared_sample
        if sample and sample == "nan":
            sample = ""

        if sample and sample not in sample_cache:
            if hla_dir:
                hla_file = Path(hla_dir) / f"{sample}_hla.txt"
                if hla_file.exists():
                    sample_cache[sample] = _read_hla_file(hla_file)

        hlas_from_row = _split_hlas(row.get(hla_column)) if hla_column else []

        if not hlas_from_row:
            if not sample:
                raise ValueError(
                    f"Sample is required for row without HLA: {input_file}"
                )

            hlas_from_row = sample_cache.get(sample)
            if not hlas_from_row:
                hla_file = Path(hla_dir) / f"{sample}_hla.txt"
                if not hla_file.exists():
                    raise FileNotFoundError(f"Missing HLA file: {hla_file}")
                hlas_from_row = _read_hla_file(hla_file)
                sample_cache[sample] = hlas_from_row

        for hla in hlas_from_row:
            normalized = normalize_hla(hla)
            if not normalized:
                continue

            data = row.to_dict()
            data["Peptide"] = peptide
            data["HLA"] = normalized
            data["Sample"] = sample if sample else data.get("Sample")
            data = _ensure_input_source_column(data, source_columns, input_source)
            prepared.append(data)

    if not prepared:
        raise ValueError("No valid peptide records after preparation")

    prepared_df = pd.DataFrame(prepared)
    if "Peptide" in prepared_df.columns:
        prepared_df["Peptide"] = prepared_df["Peptide"].astype(str)
    prepared_df = prepared_df.drop_duplicates()

    pairs_df = prepared_df[["HLA", "Peptide"]].drop_duplicates()
    pairs_df["HLA"] = pairs_df["HLA"].map(hla_to_mixmhc2)
    pairs_df = pairs_df[pairs_df["HLA"].astype(bool)]

    if pairs_df.empty:
        raise ValueError("No valid HLA-Peptide pairs prepared")

    prepared_df.to_csv(output_file, index=False)
    pairs_df.to_csv(pairs_file, sep="\t", index=False, header=False)

    return prepared_df, pairs_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--pairs_file", required=True)
    parser.add_argument("--hla_dir", required=True)
    parser.add_argument("--default_sample", default="")
    parser.add_argument("--input_source", default="custom_peptides")

    args = parser.parse_args()
    prepare_mixmhc2_inputs(
        input_file=args.input_file,
        output_file=args.output_file,
        pairs_file=args.pairs_file,
        hla_dir=args.hla_dir,
        default_sample=args.default_sample,
        input_source=args.input_source,
    )


if __name__ == "__main__":
    main()
