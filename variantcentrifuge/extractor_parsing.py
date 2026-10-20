"""Parsing of SnpEff ANN and NMD subfields from raw extracted columns."""

import pandas as pd

# SnpEff ANN field positions (pipe-delimited)
# See: http://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files
ANN_FIELD_POSITIONS = {
    "Allele": 0,
    "Annotation": 1,  # EFFECT
    "Annotation_Impact": 2,  # IMPACT
    "Gene_Name": 3,  # GENE
    "Gene_ID": 4,
    "Feature_Type": 5,
    "Feature_ID": 6,  # FEATUREID
    "Transcript_BioType": 7,
    "Rank": 8,
    "HGVS.c": 9,  # HGVS_C
    "HGVS.p": 10,  # HGVS_P
    "cDNA.pos/cDNA.length": 11,
    "CDS.pos/CDS.length": 12,
    "AA.pos/AA.length": 13,  # Contains both AA_POS and AA_LEN
    "Distance": 14,
    "ERRORS/WARNINGS/INFO": 15,
}

# Map config field names to ANN positions
ANN_SUBFIELD_MAP = {
    "ANN[0].ALLELE": ("Allele", 0),
    "ANN[0].GENE": ("Gene_Name", 3),
    "ANN[0].FEATUREID": ("Feature_ID", 6),
    "ANN[0].EFFECT": ("Annotation", 1),
    "ANN[0].IMPACT": ("Annotation_Impact", 2),
    "ANN[0].HGVS_C": ("HGVS.c", 9),
    "ANN[0].HGVS_P": ("HGVS.p", 10),
    "ANN[0].AA_POS": ("AA.pos/AA.length", 13),  # Will need split
    "ANN[0].AA_LEN": ("AA.pos/AA.length", 13),  # Will need split
}

# NMD field positions (pipe-delimited SnpEff LOF/NMD format:
# Gene | ID | num_transcripts | percent_affected)
NMD_SUBFIELD_MAP = {
    "NMD[0].GENE": 0,
    "NMD[0].GENEID": 1,
    "NMD[0].NUMTR": 2,
    "NMD[0].PERC": 3,
}


def parse_ann_subfields(df: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    """
    Parse ANN subfields from raw ANN column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with raw "ANN" column containing pipe-delimited annotations
    fields : list[str]
        Original field list to determine which ANN subfields to extract

    Returns
    -------
    pd.DataFrame
        DataFrame with ANN subfields added and raw ANN column removed
    """
    if "ANN" not in df.columns:
        return df

    # Determine which ANN subfields are needed
    ann_fields_needed = [f for f in fields if f.startswith("ANN[0].")]

    if not ann_fields_needed:
        # No ANN subfields requested, just drop the raw column
        return df.drop(columns=["ANN"])

    # Select first annotation matching ALT (when biallelic ALT is present)
    if "ALT" in df.columns:
        first_ann = pd.Series(
            [
                next(
                    (e for e in str(ann).split(",") if e and e.split("|", 1)[0] == str(alt)),
                    "" if str(alt) and "," not in str(alt) else str(ann).split(",", 1)[0],
                )
                for ann, alt in zip(df["ANN"].fillna(""), df["ALT"].fillna(""), strict=False)
            ],
            index=df.index,
            dtype="object",
        )
    else:
        first_ann = df["ANN"].fillna("").str.split(",", n=1).str[0]

    first_ann = first_ann.replace({"": pd.NA, ".": pd.NA})

    # Split by pipe
    ann_split = first_ann.str.split("|", expand=True)

    # Extract each requested subfield
    for field in ann_fields_needed:
        if field in ANN_SUBFIELD_MAP:
            _, pos = ANN_SUBFIELD_MAP[field]

            # Special handling for AA_POS and AA_LEN (combined in one field)
            if field == "ANN[0].AA_POS":
                # Extract position part (before "/")
                if pos < ann_split.shape[1]:
                    s = ann_split[pos].fillna("").str.split("/").str[0]
                    df[field] = s.where(ann_split[pos].notna(), "NA")
                else:
                    df[field] = "NA"
            elif field == "ANN[0].AA_LEN":
                # Extract length part (after "/")
                if pos < ann_split.shape[1]:
                    s = ann_split[pos].fillna("").str.split("/").str[1].fillna("")
                    df[field] = s.where(ann_split[pos].notna(), "NA")
                else:
                    df[field] = "NA"
            else:
                # Normal extraction
                if pos < ann_split.shape[1]:
                    df[field] = ann_split[pos].fillna("NA")
                else:
                    df[field] = "NA"

    # Drop raw ANN column
    df = df.drop(columns=["ANN"])

    return df


def parse_nmd_subfields(df: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    """
    Parse NMD subfields from raw NMD column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with raw "NMD" column containing pipe-delimited data
    fields : list[str]
        Original field list to determine which NMD subfields to extract

    Returns
    -------
    pd.DataFrame
        DataFrame with NMD subfields added and raw NMD column removed
    """
    if "NMD" not in df.columns:
        return df

    # Determine which NMD subfields are needed
    nmd_fields_needed = [f for f in fields if f.startswith("NMD[0].")]

    if not nmd_fields_needed:
        return df.drop(columns=["NMD"])

    # Take first NMD annotation and remove SnpEff's surrounding parentheses:
    # NMD=(GENE|GENEID|NUMTR|PERC)
    first_nmd = df["NMD"].fillna("").str.split(",", n=1).str[0]
    first_nmd = first_nmd.str.strip().str.replace(r"^\((.*)\)$", r"\1", regex=True)

    # Split by pipe
    nmd_split = first_nmd.str.split("|", expand=True)

    # Extract each requested subfield
    for field in nmd_fields_needed:
        if field in NMD_SUBFIELD_MAP:
            pos = NMD_SUBFIELD_MAP[field]
            if pos < nmd_split.shape[1]:
                df[field] = nmd_split[pos].fillna("NA")
            else:
                df[field] = "NA"

    # Drop raw NMD column
    df = df.drop(columns=["NMD"])

    return df
