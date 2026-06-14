process MHCNUGGETS {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(peptides_csv)

    output:
    tuple val(sample_id), path("${sample_id}_mhcnuggets.csv")

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 - <<'PY'
    import pandas as pd

    from scoring.hla import normalize_hla, normalize_peptide

    df = pd.read_csv("${peptides_csv}")
    hla_col = next((name for name in df.columns if name.lower() == "hla"), None)
    peptide_col = next(
        (name for name in df.columns if name.lower() in {"peptide", "tumor_peptide", "mutant_seq"}),
        None,
    )
    if hla_col is None or peptide_col is None:
        raise ValueError("Input must contain HLA and Peptide columns.")

    output = pd.DataFrame({
        "HLA": df[hla_col].map(normalize_hla),
        "Peptide": df[peptide_col].map(normalize_peptide),
    })
    output["MHCnuggets_IC50"] = pd.NA
    output = output.drop_duplicates()
    output.to_csv("${sample_id}_mhcnuggets.csv", index=False)
    PY
    """
}

process PSSMHCPAN {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(peptides_csv)

    output:
    tuple val(sample_id), path("${sample_id}_pssmhcpan.csv")

    script:
    """
    PYTHONPATH=${projectDir}/scr python3 - <<'PY'
    import pandas as pd

    from scoring.hla import normalize_hla, normalize_peptide

    df = pd.read_csv("${peptides_csv}")
    hla_col = next((name for name in df.columns if name.lower() == "hla"), None)
    peptide_col = next(
        (name for name in df.columns if name.lower() in {"peptide", "tumor_peptide", "mutant_seq"}),
        None,
    )
    if hla_col is None or peptide_col is None:
        raise ValueError("Input must contain HLA and Peptide columns.")

    output = pd.DataFrame({
        "HLA": df[hla_col].map(normalize_hla),
        "Peptide": df[peptide_col].map(normalize_peptide),
    })
    output["PSSMHCpan_IC50"] = pd.NA
    output = output.drop_duplicates()
    output.to_csv("${sample_id}_pssmhcpan.csv", index=False)
    PY
    """
}
