import re


def normalize_hla(value):
    text = str(value).strip().upper()
    if not text or text == "NAN":
        return ""

    text = re.sub(r"^HLA[-_]?", "", text)
    return re.sub(r"[^A-Z0-9]", "", text)


def normalize_peptide(value):
    text = str(value).strip().upper()
    if not text or text == "NAN":
        return ""
    return text.replace("-", "")


def hla_to_mhcflurry(value):
    normalized = normalize_hla(value)
    match = re.fullmatch(r"([A-Z]+)(\d{2})(\d{2})", normalized)
    if not match:
        return str(value).strip()
    locus, first, second = match.groups()
    return f"HLA-{locus}*{first}:{second}"


def hla_to_mixmhcpred(value):
    return normalize_hla(value)


def hla_to_mixmhc2(value):
    return hla_to_mixmhcpred(value)
