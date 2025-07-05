#!/usr/bin/env python3
"""Example script demonstrating the pseudonymization functionality."""

import pandas as pd
from variantcentrifuge.pseudonymizer import create_pseudonymizer

# Example 1: Sequential Schema
print("=== Sequential Schema Example ===")
pseudonymizer = create_pseudonymizer("sequential", prefix="STUDY")
samples = ["Patient_A", "Patient_B", "Patient_C", "Control_1", "Control_2"]
mapping = pseudonymizer.create_mapping(samples)
print("Mapping:")
for orig, pseudo in mapping.items():
    print(f"  {orig} -> {pseudo}")

# Example 2: Categorical Schema
print("\n=== Categorical Schema Example ===")
pseudonymizer = create_pseudonymizer("categorical", category_field="phenotype")
metadata = {
    "Patient_A": {"phenotype": "case", "sex": "M"},
    "Patient_B": {"phenotype": "case", "sex": "F"},
    "Patient_C": {"phenotype": "case", "sex": "M"},
    "Control_1": {"phenotype": "control", "sex": "F"},
    "Control_2": {"phenotype": "control", "sex": "M"},
}
mapping = pseudonymizer.create_mapping(samples, metadata)
print("Mapping:")
for orig, pseudo in mapping.items():
    print(f"  {orig} -> {pseudo}")

# Example 3: Anonymous Schema
print("\n=== Anonymous Schema Example ===")
pseudonymizer = create_pseudonymizer("anonymous")
mapping = pseudonymizer.create_mapping(samples)
print("Mapping:")
for orig, pseudo in mapping.items():
    print(f"  {orig} -> {pseudo}")

# Example 4: Custom Schema
print("\n=== Custom Schema Example ===")
pattern = "{study}_{phenotype}_{index:03d}"
pseudonymizer = create_pseudonymizer("custom", pattern=pattern)
# Add study to metadata
for sample, meta in metadata.items():
    meta["study"] = "GENO2025"
mapping = pseudonymizer.create_mapping(samples, metadata)
print("Mapping:")
for orig, pseudo in mapping.items():
    print(f"  {orig} -> {pseudo}")

# Example 5: Pseudonymizing GT Column
print("\n=== Pseudonymizing Genotype Data ===")
df = pd.DataFrame(
    {
        "CHROM": ["chr1", "chr2", "chr3"],
        "POS": [1000, 2000, 3000],
        "GT": ["Patient_A(0/1);Control_1(0/0)", "Patient_B(1/1);Patient_C(0/1)", "Control_2(0/1)"],
    }
)
print("Original:")
print(df[["CHROM", "POS", "GT"]])

# Apply pseudonymization using sequential schema
pseudonymizer = create_pseudonymizer("sequential", prefix="ID")
pseudonymizer.create_mapping(samples)
df_pseudo = pseudonymizer.pseudonymize_dataframe(df)
print("\nPseudonymized:")
print(df_pseudo[["CHROM", "POS", "GT"]])

# Save mapping
print("\n=== Saving Mapping Table ===")
pseudonymizer.save_mapping("example_mapping.tsv")
print("Mapping saved to example_mapping.tsv")

# Load and display
mapping_df = pd.read_csv("example_mapping.tsv", sep="\t")
print("\nMapping table contents:")
print(mapping_df)
