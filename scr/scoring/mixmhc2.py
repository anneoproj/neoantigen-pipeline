import re
import pandas as pd


MIXMHC2_COLUMNS = [
    "MixMHC2pred_Score",
    "MixMHC2pred_PercentileRank",
]


def _normalize(name):
    return re.sub(r"\W+", "", str(name).strip().lower())


def _pick_column(columns, candidates):
    normalized = {_normalize(column): column for column in columns}
    for candidate in candidates:
        for name in columns:
            if _normalize(name) == candidate:
                return name
        if candidate in normalized:
            return normalized[candidate]
        if candidate in normalized.values():
            return candidate
    return None


def parse_mixmhc2_output(path):
    try:
        df = pd.read_csv(path)
    except Exception:
        df = pd.read_csv(path, sep="\t")

    if df.empty:
        return pd.DataFrame(columns=["HLA", "Peptide", *MIXMHC2_COLUMNS])

    hla_column = _pick_column(df.columns, {"hla", "allele", "mhc", "ml"})
    peptide_column = _pick_column(df.columns, {"peptide", "sequence"})
    if not hla_column or not peptide_column:
        raise ValueError(f"MixMHC2 output missing HLA/Peptide columns: {path}")

    parsed = pd.DataFrame({
        "HLA": df[hla_column],
        "Peptide": df[peptide_column],
    })

    score_candidates = [
        "score", "mixmhc2score", "prediction",
        "predscore", "affinity", "binding_score",
    ]
    rank_candidates = [
        "percentilerank", "percentile", "rank",
        "percentile_rank", "rank_percentile", "rank_el",
        "percentilerank_el", "mixmhc2percentilerank",
    ]

    score_col = _pick_column(df.columns, score_candidates)
    rank_col = _pick_column(df.columns, rank_candidates)

    parsed["MixMHC2pred_Score"] = (
        pd.to_numeric(df[score_col], errors="coerce")
        if score_col else pd.NA
    )
    parsed["MixMHC2pred_PercentileRank"] = (
        pd.to_numeric(df[rank_col], errors="coerce")
        if rank_col else pd.NA
    )

    return parsed
