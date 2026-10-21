"""Parsing of SnpEff ANN, LOF and NMD subfields from raw extracted columns."""

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

# SnpEff LOF and NMD share one pipe-delimited format:
# (Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected)
LOF_NMD_TAGS = ("LOF", "NMD")
LOF_NMD_SUBFIELD_POSITIONS = {"GENE": 0, "GENEID": 1, "NUMTR": 2, "PERC": 3}
NMD_SUBFIELD_MAP = {f"NMD[0].{k}": v for k, v in LOF_NMD_SUBFIELD_POSITIONS.items()}


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
    Parse LOF[0].* and NMD[0].* subfields from the raw LOF and NMD columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with raw "LOF"/"NMD" columns containing pipe-delimited data
    fields : list[str]
        Original field list to determine which subfields to extract

    Returns
    -------
    pd.DataFrame
        DataFrame with subfields added and raw LOF/NMD columns removed
    """
    for tag in LOF_NMD_TAGS:
        if tag not in df.columns:
            continue
        # Take the first annotation and strip SnpEff's parentheses: (GENE|GENEID|NUMTR|PERC)
        first = df[tag].fillna("").str.split(",", n=1).str[0]
        first = first.str.strip().str.replace(r"^\((.*)\)$", r"\1", regex=True)
        split = first.str.split("|", expand=True)
        for field in fields:
            if field.startswith(f"{tag}[0]."):
                pos = LOF_NMD_SUBFIELD_POSITIONS.get(field.split(".", 1)[1])
                if pos is not None:
                    df[field] = split[pos].fillna("NA") if pos < split.shape[1] else "NA"
        if tag not in fields:  # a bare LOF/NMD field keeps its raw column
            df = df.drop(columns=[tag])
    return df
