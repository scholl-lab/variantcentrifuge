"""
Sample pseudonymization module for privacy-preserving data sharing.

Supports multiple naming schemas suitable for different publication contexts:
- Sequential: STUDY_001, STUDY_002, etc.
- Categorical: CASE_001, CONTROL_001, etc.
- Anonymous: A001, B002, etc.
- Custom: User-defined pattern
"""

import hashlib
import io
import json
import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class PseudonymSchema:
    """Base class for pseudonymization schemas."""

    def generate(self, sample_id: str, index: int, metadata: dict | None = None) -> str:
        """Generate a pseudonym for a sample."""
        raise NotImplementedError


class SequentialSchema(PseudonymSchema):
    """Sequential numbering schema: PREFIX_001, PREFIX_002, etc."""

    def __init__(self, prefix: str = "SAMPLE", padding: int = 3):
        """Initialize sequential schema.

        Parameters
        ----------
        prefix : str
            Prefix for generated IDs
        padding : int
            Number of digits to pad with zeros
        """
        self.prefix = prefix
        self.padding = padding

    def generate(self, sample_id: str, index: int, metadata: dict | None = None) -> str:
        """Generate sequential pseudonym."""
        return f"{self.prefix}_{str(index).zfill(self.padding)}"


class CategoricalSchema(PseudonymSchema):
    """Category-based schema: CASE_001, CONTROL_001, etc."""

    def __init__(self, category_field: str = "phenotype", padding: int = 3):
        """Initialize categorical schema.

        Parameters
        ----------
        category_field : str
            Metadata field to use for categorization
        padding : int
            Number of digits to pad with zeros
        """
        self.category_field = category_field
        self.padding = padding
        self._category_counters: dict[str, int] = {}

    def generate(self, sample_id: str, index: int, metadata: dict | None = None) -> str:
        """Generate category-based pseudonym."""
        category = "UNKNOWN"
        if metadata and self.category_field in metadata:
            category = str(metadata[self.category_field]).upper()
            # Sanitize category name
            category = re.sub(r"[^A-Z0-9]", "", category)[:10]

        if category not in self._category_counters:
            self._category_counters[category] = 0
        self._category_counters[category] += 1

        return f"{category}_{str(self._category_counters[category]).zfill(self.padding)}"


class AnonymousSchema(PseudonymSchema):
    """Anonymous letter-number schema: A001, B002, etc."""

    def __init__(self, use_hash: bool = False):
        """Initialize anonymous schema.

        Parameters
        ----------
        use_hash : bool
            If True, use MD5 hash for ID generation
        """
        self.use_hash = use_hash

    def generate(self, sample_id: str, index: int, metadata: dict | None = None) -> str:
        """Generate anonymous pseudonym."""
        if self.use_hash:
            # Generate deterministic hash-based ID
            hash_val = hashlib.md5(sample_id.encode()).hexdigest()[:6].upper()
            return f"ID{hash_val}"
        else:
            # Simple letter-number combination
            letter = chr(65 + ((index - 1) // 999))  # A, B, C, etc.
            number = ((index - 1) % 999) + 1
            return f"{letter}{str(number).zfill(3)}"


class CustomSchema(PseudonymSchema):
    """Custom pattern-based schema with placeholders."""

    def __init__(self, pattern: str):
        """Initialize with a pattern containing placeholders.

        Parameters
        ----------
        pattern : str
            Pattern with placeholders:
            - {prefix}: Custom prefix
            - {index}: Sequential index
            - {category}: Category from metadata
            - {hash}: First 6 chars of MD5 hash

        Example: "{prefix}_{category}_{index:03d}"
        """
        self.pattern = pattern

    def generate(self, sample_id: str, index: int, metadata: dict | None = None) -> str:
        """Generate custom pattern pseudonym."""
        replacements = {
            "index": index,
            "hash": hashlib.md5(sample_id.encode()).hexdigest()[:6].upper(),
        }

        if metadata:
            replacements.update(metadata)

        try:
            return self.pattern.format(**replacements)
        except KeyError as e:
            logger.warning(f"Missing placeholder in pattern: {e}")
            return f"SAMPLE_{index:03d}"


class SamplePseudonymizer:
    """Main class for sample pseudonymization."""

    def __init__(self, schema: PseudonymSchema, deterministic: bool = True):
        """Initialize pseudonymizer with a schema.

        Parameters
        ----------
        schema : PseudonymSchema
            The naming schema to use
        deterministic : bool
            If True, sort samples before assignment for reproducibility
        """
        self.schema = schema
        self.deterministic = deterministic
        self._mapping: dict[str, str] = {}
        self._reverse_mapping: dict[str, str] = {}

    def create_mapping(
        self, sample_list: list[str], metadata: dict[str, dict] | None = None
    ) -> dict[str, str]:
        """Create pseudonym mapping for samples.

        Parameters
        ----------
        sample_list : List[str]
            List of original sample IDs
        metadata : Optional[Dict[str, Dict]]
            Sample metadata for category-based schemas

        Returns
        -------
        Dict[str, str]
            Mapping from original to pseudonym IDs
        """
        if self.deterministic:
            sample_list = sorted(set(sample_list))

        self._mapping.clear()
        self._reverse_mapping.clear()

        for i, sample_id in enumerate(sample_list, 1):
            sample_meta = metadata.get(sample_id, {}) if metadata else {}
            pseudonym = self.schema.generate(sample_id, i, sample_meta)

            # Ensure uniqueness
            if pseudonym in self._reverse_mapping:
                logger.warning(f"Duplicate pseudonym {pseudonym}, adding suffix")
                suffix = 2
                while f"{pseudonym}_{suffix}" in self._reverse_mapping:
                    suffix += 1
                pseudonym = f"{pseudonym}_{suffix}"

            self._mapping[sample_id] = pseudonym
            self._reverse_mapping[pseudonym] = sample_id

        logger.info(f"Created pseudonyms for {len(self._mapping)} samples")
        return self._mapping

    def pseudonymize_gt_column(self, gt_value: str) -> str:
        """Replace sample IDs in a GT column value.

        Handles formats like:
        - "SAMPLE1(0/1);SAMPLE2(1/1)"
        - "SAMPLE1(0/1:50,50);SAMPLE2(./.)"
        """
        if not isinstance(gt_value, str) or not gt_value:
            return gt_value

        # Pattern to match SAMPLE_ID(genotype_info)
        pattern = re.compile(r"([^;()\s]+)\(([^)]+)\)")

        def replace_func(match):
            original_id = match.group(1)
            genotype_info = match.group(2)
            pseudonym = self._mapping.get(original_id, original_id)
            return f"{pseudonym}({genotype_info})"

        return pattern.sub(replace_func, gt_value)

    def pseudonymize_inheritance_column(self, inheritance_value: str) -> str:
        """Replace sample IDs in inheritance analysis columns.

        Handles formats like:
        - "child(0/1); father(0/1), partner:chr2:3000:C>T"
        - "proband1(0/1)"
        """
        if not isinstance(inheritance_value, str) or not inheritance_value:
            return inheritance_value

        # Apply mapping to each word that matches a known sample ID
        parts = []
        for part in inheritance_value.split():
            # Check if this part contains a sample ID
            for original_id, pseudonym in self._mapping.items():
                if part.startswith(f"{original_id}("):
                    part = part.replace(original_id, pseudonym, 1)
                    break
            parts.append(part)

        return " ".join(parts)

    def pseudonymize_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply pseudonymization to all relevant columns in a DataFrame."""
        df_copy = df.copy()

        # Handle GT column
        if "GT" in df_copy.columns:
            df_copy["GT"] = df_copy["GT"].apply(self.pseudonymize_gt_column)
            logger.debug("Pseudonymized GT column")

        # Handle inheritance-related columns
        inheritance_cols = ["Inheritance_Samples", "Inheritance_Description"]
        for col in inheritance_cols:
            if col in df_copy.columns:
                df_copy[col] = df_copy[col].apply(self.pseudonymize_inheritance_column)
                logger.debug(f"Pseudonymized {col} column")

        # Handle any column that might contain sample IDs (configurable)
        # This could be extended based on specific needs

        return df_copy

    def pseudonymize_ped_file(self, ped_path: str, output_path: str):
        """Create a pseudonymized version of a PED file."""
        ped_df = pd.read_csv(
            ped_path,
            sep="\t",
            header=None,
            names=["family", "individual", "father", "mother", "sex", "phenotype"],
        )

        # Apply mapping to all ID columns
        for col in ["individual", "father", "mother"]:
            ped_df[col] = ped_df[col].apply(lambda x: self._mapping.get(x, x) if x != "0" else "0")

        # Optionally pseudonymize family IDs too
        family_mapping = {}
        for i, family in enumerate(sorted(ped_df["family"].unique()), 1):
            family_mapping[family] = f"FAM{i:03d}"
        ped_df["family"] = ped_df["family"].map(family_mapping)

        # Apply compression if output path ends with .gz
        compression = "gzip" if str(output_path).endswith(".gz") else None
        ped_df.to_csv(output_path, sep="\t", header=False, index=False, compression=compression)
        logger.info(f"Pseudonymized PED file saved to {output_path}")

    def save_mapping(self, filepath: str, include_metadata: bool = True):
        """Save the pseudonymization mapping to a TSV file.

        Parameters
        ----------
        filepath : str
            Output file path
        include_metadata : bool
            Include schema information and timestamp
        """
        # Create mapping dataframe
        map_df = pd.DataFrame(list(self._mapping.items()), columns=["original_id", "pseudonym_id"])

        # Save primary mapping with compression if path ends with .gz
        compression = "gzip" if str(filepath).endswith(".gz") else None
        map_df.to_csv(filepath, sep="\t", index=False, compression=compression)

        # Save metadata if requested
        if include_metadata:
            metadata = {
                "schema": self.schema.__class__.__name__,
                "deterministic": self.deterministic,
                "creation_date": pd.Timestamp.now().isoformat(),
                "num_samples": len(self._mapping),
            }

            meta_path = Path(filepath).with_suffix(".meta.json")
            with open(meta_path, "w") as f:
                json.dump(metadata, f, indent=2)

        logger.info(f"Pseudonymization mapping saved to {filepath}")

    def load_mapping(self, filepath: str) -> dict[str, str]:
        """Load a previously saved mapping."""
        map_df = pd.read_csv(filepath, sep="\t", low_memory=False)
        self._mapping = dict(zip(map_df["original_id"], map_df["pseudonym_id"], strict=False))
        self._reverse_mapping = {v: k for k, v in self._mapping.items()}
        logger.info(f"Loaded {len(self._mapping)} sample mappings")
        return self._mapping


def create_pseudonymizer(schema_type: str, **kwargs) -> SamplePseudonymizer:
    """Create a pseudonymizer with the specified schema.

    Parameters
    ----------
    schema_type : str
        One of: 'sequential', 'categorical', 'anonymous', 'custom'
    **kwargs
        Schema-specific parameters

    Returns
    -------
    SamplePseudonymizer
        Configured pseudonymizer instance
    """
    schemas = {
        "sequential": SequentialSchema,
        "categorical": CategoricalSchema,
        "anonymous": AnonymousSchema,
        "custom": CustomSchema,
    }

    if schema_type not in schemas:
        raise ValueError(f"Unknown schema type: {schema_type}")

    schema_class = schemas[schema_type]
    schema = schema_class(**kwargs)

    return SamplePseudonymizer(schema)


def apply_pseudonymization(
    buffer: list[str],
    sample_list: list[str],
    cfg: dict,
    ped_data: pd.DataFrame | None = None,
) -> tuple[list[str], SamplePseudonymizer | None]:
    """Apply pseudonymization to the output buffer if configured.

    Parameters
    ----------
    buffer : List[str]
        Output buffer lines
    sample_list : List[str]
        List of sample IDs
    cfg : Dict
        Configuration dictionary
    ped_data : Optional[pd.DataFrame]
        Pedigree data for metadata extraction

    Returns
    -------
    Tuple[List[str], Optional[SamplePseudonymizer]]
        Modified buffer and pseudonymizer instance (if used)
    """
    if not cfg.get("pseudonymize"):
        return buffer, None

    logger.info("Applying sample pseudonymization...")

    # Create pseudonymizer based on schema
    schema_type = cfg.get("pseudonymize_schema", "sequential")

    if schema_type == "sequential":
        pseudonymizer = create_pseudonymizer(
            "sequential", prefix=cfg.get("pseudonymize_prefix", "SAMPLE")
        )
    elif schema_type == "categorical":
        # Collect metadata from PED if available
        metadata = {}
        if ped_data is not None:
            for _, row in ped_data.iterrows():
                metadata[row["individual"]] = {
                    "phenotype": "CASE" if row["phenotype"] == 2 else "CONTROL",
                    "sex": "M" if row["sex"] == 1 else "F",
                    "family": row["family"],
                }

        pseudonymizer = create_pseudonymizer(
            "categorical", category_field=cfg.get("pseudonymize_category_field", "phenotype")
        )
        pseudonymizer.create_mapping(sample_list, metadata)
    elif schema_type == "anonymous":
        pseudonymizer = create_pseudonymizer("anonymous")
        pseudonymizer.create_mapping(sample_list)
    elif schema_type == "custom":
        if not cfg.get("pseudonymize_pattern"):
            raise ValueError("--pseudonymize-pattern required for custom schema")
        pseudonymizer = create_pseudonymizer("custom", pattern=cfg["pseudonymize_pattern"])
        pseudonymizer.create_mapping(sample_list)
    else:
        raise ValueError(f"Unknown schema type: {schema_type}")

    # Create mapping if not already done
    if not pseudonymizer._mapping:
        pseudonymizer.create_mapping(sample_list)

    # Convert buffer to DataFrame
    if len(buffer) > 1:
        df = pd.read_csv(io.StringIO("\n".join(buffer)), sep="\t", dtype=str)

        # Apply pseudonymization
        df = pseudonymizer.pseudonymize_dataframe(df)

        # Convert back to buffer
        output = io.StringIO()
        df.to_csv(output, sep="\t", index=False, na_rep="")
        buffer = output.getvalue().strip().split("\n")

    return buffer, pseudonymizer
