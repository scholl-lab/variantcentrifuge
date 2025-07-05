"""Test suite for sample pseudonymization functionality."""

import json
import os
import tempfile

import pandas as pd
import pytest

from variantcentrifuge.pseudonymizer import (
    AnonymousSchema,
    CategoricalSchema,
    CustomSchema,
    SequentialSchema,
    apply_pseudonymization,
    create_pseudonymizer,
)


class TestPseudonymSchemas:
    """Test different pseudonymization schemas."""

    def test_sequential_schema(self):
        """Test sequential numbering schema."""
        schema = SequentialSchema(prefix="STUDY", padding=4)

        assert schema.generate("sample1", 1) == "STUDY_0001"
        assert schema.generate("sample2", 99) == "STUDY_0099"
        assert schema.generate("sample3", 1000) == "STUDY_1000"

    def test_sequential_schema_default(self):
        """Test sequential schema with defaults."""
        schema = SequentialSchema()

        assert schema.generate("s1", 1) == "SAMPLE_001"
        assert schema.generate("s2", 10) == "SAMPLE_010"

    def test_categorical_schema(self):
        """Test category-based schema."""
        schema = CategoricalSchema(category_field="phenotype")

        # Test with metadata
        meta1 = {"phenotype": "case"}
        meta2 = {"phenotype": "control"}

        assert schema.generate("s1", 1, meta1) == "CASE_001"
        assert schema.generate("s2", 2, meta1) == "CASE_002"
        assert schema.generate("s3", 3, meta2) == "CONTROL_001"

    def test_categorical_schema_sanitization(self):
        """Test category name sanitization."""
        schema = CategoricalSchema()

        # Test with special characters
        meta = {"phenotype": "case-control_test!@#"}
        result = schema.generate("s1", 1, meta)
        assert result == "CASECONTRO_001"  # Sanitized to first 10 alphanumeric chars

    def test_categorical_schema_unknown(self):
        """Test categorical schema without metadata."""
        schema = CategoricalSchema()

        # No metadata provided
        assert schema.generate("s1", 1) == "UNKNOWN_001"
        assert schema.generate("s2", 2, {}) == "UNKNOWN_002"

    def test_anonymous_schema(self):
        """Test anonymous schema."""
        schema = AnonymousSchema(use_hash=False)

        assert schema.generate("s1", 1) == "A001"
        assert schema.generate("s2", 999) == "A999"
        assert schema.generate("s3", 1000) == "B001"
        assert schema.generate("s4", 1998) == "B999"
        assert schema.generate("s5", 1999) == "C001"

    def test_anonymous_schema_hash(self):
        """Test hash-based anonymous schema."""
        schema = AnonymousSchema(use_hash=True)

        # Test deterministic hash generation
        id1 = schema.generate("sample1", 1)
        id2 = schema.generate("sample1", 2)  # Same input
        assert id1 == id2  # Should be deterministic
        assert id1.startswith("ID")
        assert len(id1) == 8  # ID + 6 hash chars

        # Different sample should give different hash
        id3 = schema.generate("sample2", 1)
        assert id3 != id1

    def test_custom_schema(self):
        """Test custom pattern schema."""
        schema = CustomSchema("{prefix}_{phenotype}_{index:04d}")

        meta = {"prefix": "PROJ", "phenotype": "CASE"}
        assert schema.generate("s1", 1, meta) == "PROJ_CASE_0001"

    def test_custom_schema_missing_placeholder(self):
        """Test custom schema with missing metadata."""
        schema = CustomSchema("{study}_{index:03d}")

        # Missing 'study' in metadata
        result = schema.generate("s1", 1, {"other": "value"})
        assert result == "SAMPLE_001"  # Fallback

    def test_custom_schema_hash_placeholder(self):
        """Test custom schema with hash placeholder."""
        schema = CustomSchema("ID_{hash}_{index:02d}")

        result = schema.generate("test_sample", 5)
        assert result.startswith("ID_")
        assert len(result.split("_")[1]) == 6  # Hash is 6 chars
        assert result.endswith("_05")


class TestSamplePseudonymizer:
    """Test the main pseudonymizer functionality."""

    def test_create_mapping_deterministic(self):
        """Test deterministic mapping creation."""
        pseudonymizer = create_pseudonymizer("sequential", prefix="TEST")

        samples = ["Charlie", "Alice", "Bob"]
        mapping = pseudonymizer.create_mapping(samples)

        # Should be sorted alphabetically
        assert mapping["Alice"] == "TEST_001"
        assert mapping["Bob"] == "TEST_002"
        assert mapping["Charlie"] == "TEST_003"

        # Test reproducibility
        mapping2 = pseudonymizer.create_mapping(samples)
        assert mapping == mapping2

    def test_create_mapping_non_deterministic(self):
        """Test non-deterministic mapping."""
        from variantcentrifuge.pseudonymizer import SamplePseudonymizer

        schema = SequentialSchema()
        pseudonymizer = SamplePseudonymizer(schema, deterministic=False)

        samples = ["Charlie", "Alice", "Bob"]
        mapping = pseudonymizer.create_mapping(samples)

        # Should preserve original order
        assert mapping["Charlie"] == "SAMPLE_001"
        assert mapping["Alice"] == "SAMPLE_002"
        assert mapping["Bob"] == "SAMPLE_003"

    def test_create_mapping_duplicates(self):
        """Test handling of duplicate sample IDs."""
        pseudonymizer = create_pseudonymizer("sequential")

        samples = ["S1", "S2", "S1", "S3", "S2"]
        mapping = pseudonymizer.create_mapping(samples)

        # Should handle duplicates
        assert len(mapping) == 3
        assert "S1" in mapping
        assert "S2" in mapping
        assert "S3" in mapping

    def test_pseudonymize_gt_column(self):
        """Test GT column pseudonymization."""
        pseudonymizer = create_pseudonymizer("sequential", prefix="S")
        pseudonymizer.create_mapping(["John", "Jane", "Jim"])

        # Test various GT formats
        # After sorting: Jane->S_001, Jim->S_002, John->S_003
        assert pseudonymizer.pseudonymize_gt_column("John(0/1)") == "S_003(0/1)"
        assert pseudonymizer.pseudonymize_gt_column("Jane(0/1);Jim(1/1)") == "S_001(0/1);S_002(1/1)"
        assert (
            pseudonymizer.pseudonymize_gt_column("John(0/1:50,50);Jane(./.)")
            == "S_003(0/1:50,50);S_001(./.)"
        )
        assert pseudonymizer.pseudonymize_gt_column("") == ""
        assert (
            pseudonymizer.pseudonymize_gt_column("Unknown(0/1)") == "Unknown(0/1)"
        )  # Unknown sample

    def test_pseudonymize_gt_column_edge_cases(self):
        """Test GT column edge cases."""
        pseudonymizer = create_pseudonymizer("sequential")
        pseudonymizer.create_mapping(["S1"])

        # Test None and non-string inputs
        assert pseudonymizer.pseudonymize_gt_column(None) is None
        assert pseudonymizer.pseudonymize_gt_column(123) == 123

    def test_pseudonymize_inheritance_column(self):
        """Test inheritance column pseudonymization."""
        pseudonymizer = create_pseudonymizer("anonymous")
        pseudonymizer.create_mapping(["child1", "father1", "mother1"])

        # Test inheritance samples format
        input_str = "child1(0/1); father1(0/1), partner:chr2:3000:C>T"
        expected = input_str.replace("child1", "A001").replace("father1", "A002")
        assert pseudonymizer.pseudonymize_inheritance_column(input_str) == expected

    def test_pseudonymize_inheritance_column_complex(self):
        """Test complex inheritance patterns."""
        pseudonymizer = create_pseudonymizer("sequential", prefix="P")
        pseudonymizer.create_mapping(["proband", "parent1", "parent2", "sibling"])

        # Test compound het pattern
        input_str = (
            "proband(0/1); parent1(0/1), partner:chr2:3000:C>T; parent2(0/1), partner:chr2:4000:G>A"
        )
        result = pseudonymizer.pseudonymize_inheritance_column(input_str)

        assert "P_003" in result  # proband
        assert "P_001" in result  # parent1
        assert "P_002" in result  # parent2
        assert "partner:chr2:3000:C>T" in result  # Variant info preserved

    def test_pseudonymize_dataframe(self):
        """Test full dataframe pseudonymization."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2"],
                "POS": [1000, 2000],
                "GT": ["Sample1(0/1);Sample2(1/1)", "Sample3(0/1)"],
                "Inheritance_Samples": ["Sample1(0/1)", "Sample2(1/1); Sample3(0/1)"],
                "Inheritance_Description": ["De novo in Sample1", "Inherited from Sample2"],
                "Other": ["data1", "data2"],
            }
        )

        pseudonymizer = create_pseudonymizer("sequential", prefix="ID")
        pseudonymizer.create_mapping(["Sample1", "Sample2", "Sample3"])

        result_df = pseudonymizer.pseudonymize_dataframe(df)

        # Check GT column
        assert "ID_001" in result_df["GT"].iloc[0]
        assert "ID_002" in result_df["GT"].iloc[0]
        assert "ID_003" in result_df["GT"].iloc[1]

        # Check inheritance columns
        assert "ID_001" in result_df["Inheritance_Samples"].iloc[0]
        assert "ID_002" in result_df["Inheritance_Samples"].iloc[1]

        # Check other columns unchanged
        assert result_df["CHROM"].equals(df["CHROM"])
        assert result_df["POS"].equals(df["POS"])
        assert result_df["Other"].equals(df["Other"])

    def test_pseudonymize_dataframe_missing_columns(self):
        """Test dataframe pseudonymization with missing columns."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [1000], "REF": ["A"], "ALT": ["T"]})

        pseudonymizer = create_pseudonymizer("sequential")
        pseudonymizer.create_mapping(["S1"])

        # Should handle missing columns gracefully
        result_df = pseudonymizer.pseudonymize_dataframe(df)
        assert result_df.equals(df)

    def test_pseudonymize_ped_file(self):
        """Test PED file pseudonymization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test PED
            ped_path = os.path.join(tmpdir, "test.ped")
            with open(ped_path, "w") as f:
                f.write("FAM1\tChild1\tFather1\tMother1\t1\t2\n")
                f.write("FAM1\tFather1\t0\t0\t1\t1\n")
                f.write("FAM1\tMother1\t0\t0\t2\t1\n")
                f.write("FAM2\tChild2\tFather2\tMother2\t2\t2\n")
                f.write("FAM2\tFather2\t0\t0\t1\t1\n")
                f.write("FAM2\tMother2\t0\t0\t2\t1\n")

            pseudonymizer = create_pseudonymizer("sequential", prefix="IND")
            pseudonymizer.create_mapping(
                ["Child1", "Child2", "Father1", "Father2", "Mother1", "Mother2"]
            )

            output_path = os.path.join(tmpdir, "pseudo.ped")
            pseudonymizer.pseudonymize_ped_file(ped_path, output_path)

            # Read and verify
            result_df = pd.read_csv(output_path, sep="\t", header=None)

            # Check family IDs are pseudonymized
            assert all(result_df[0].str.startswith("FAM"))
            assert "FAM001" in result_df[0].values
            assert "FAM002" in result_df[0].values

            # Check individual IDs
            assert "IND_001" in result_df[1].values  # Child1
            assert "IND_002" in result_df[1].values  # Child2
            assert "IND_003" in result_df[1].values  # Father1

            # Check relationships maintained
            child1_row = result_df[result_df[1] == "IND_001"].iloc[0]
            assert child1_row[2] == "IND_003"  # Father1
            assert child1_row[3] == "IND_005"  # Mother1

    def test_save_load_mapping(self):
        """Test saving and loading mappings."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create and save
            pseudonymizer = create_pseudonymizer("sequential", prefix="TEST")
            original_mapping = pseudonymizer.create_mapping(["A", "B", "C"])

            mapping_path = os.path.join(tmpdir, "mapping.tsv")
            pseudonymizer.save_mapping(mapping_path)

            # Verify files created
            assert os.path.exists(mapping_path)
            meta_path = os.path.join(tmpdir, "mapping.meta.json")
            assert os.path.exists(meta_path)

            # Verify metadata content
            with open(meta_path, "r") as f:
                metadata = json.load(f)
            assert metadata["schema"] == "SequentialSchema"
            assert metadata["deterministic"] is True
            assert metadata["num_samples"] == 3

            # Load in new instance
            new_pseudonymizer = create_pseudonymizer("sequential")
            loaded_mapping = new_pseudonymizer.load_mapping(mapping_path)

            assert loaded_mapping == original_mapping

    def test_save_mapping_no_metadata(self):
        """Test saving mapping without metadata."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pseudonymizer = create_pseudonymizer("sequential")
            pseudonymizer.create_mapping(["S1", "S2"])

            mapping_path = os.path.join(tmpdir, "mapping.tsv")
            pseudonymizer.save_mapping(mapping_path, include_metadata=False)

            # Only TSV should exist
            assert os.path.exists(mapping_path)
            meta_path = os.path.join(tmpdir, "mapping.meta.json")
            assert not os.path.exists(meta_path)

    def test_ensure_unique_pseudonyms(self):
        """Test that duplicate pseudonyms are handled."""
        from variantcentrifuge.pseudonymizer import SamplePseudonymizer

        # Create a schema that might generate duplicates
        class DuplicateSchema(SequentialSchema):
            def generate(self, sample_id, index, metadata=None):
                return "DUP_001"  # Always return same ID

        pseudonymizer = SamplePseudonymizer(DuplicateSchema())

        mapping = pseudonymizer.create_mapping(["S1", "S2", "S3"])

        # Should add suffixes to ensure uniqueness
        assert mapping["S1"] == "DUP_001"
        assert mapping["S2"] == "DUP_001_2"
        assert mapping["S3"] == "DUP_001_3"


class TestFactoryFunction:
    """Test the create_pseudonymizer factory function."""

    def test_create_sequential(self):
        """Test creating sequential pseudonymizer."""
        p = create_pseudonymizer("sequential", prefix="TEST", padding=5)
        assert isinstance(p.schema, SequentialSchema)
        assert p.schema.prefix == "TEST"
        assert p.schema.padding == 5

    def test_create_categorical(self):
        """Test creating categorical pseudonymizer."""
        p = create_pseudonymizer("categorical", category_field="group")
        assert isinstance(p.schema, CategoricalSchema)
        assert p.schema.category_field == "group"

    def test_create_anonymous(self):
        """Test creating anonymous pseudonymizer."""
        p = create_pseudonymizer("anonymous", use_hash=True)
        assert isinstance(p.schema, AnonymousSchema)
        assert p.schema.use_hash is True

    def test_create_custom(self):
        """Test creating custom pseudonymizer."""
        p = create_pseudonymizer("custom", pattern="{prefix}_{index}")
        assert isinstance(p.schema, CustomSchema)
        assert p.schema.pattern == "{prefix}_{index}"

    def test_create_invalid_schema(self):
        """Test creating with invalid schema type."""
        with pytest.raises(ValueError, match="Unknown schema type"):
            create_pseudonymizer("invalid_type")


class TestApplyPseudonymization:
    """Test the apply_pseudonymization integration function."""

    def test_apply_sequential(self):
        """Test applying sequential pseudonymization."""
        buffer = [
            "CHROM\tPOS\tGT",
            "chr1\t1000\tSample1(0/1);Sample2(1/1)",
            "chr2\t2000\tSample3(0/1)",
        ]

        cfg = {
            "pseudonymize": True,
            "pseudonymize_schema": "sequential",
            "pseudonymize_prefix": "TEST",
        }

        result_buffer, pseudonymizer = apply_pseudonymization(
            buffer, ["Sample1", "Sample2", "Sample3"], cfg
        )

        assert pseudonymizer is not None
        assert "TEST_001" in result_buffer[1]  # Sample1
        assert "TEST_002" in result_buffer[1]  # Sample2
        assert "TEST_003" in result_buffer[2]  # Sample3

    def test_apply_categorical_with_ped(self):
        """Test applying categorical pseudonymization with PED data."""
        buffer = ["CHROM\tPOS\tGT", "chr1\t1000\tCase1(0/1);Control1(1/1)"]

        ped_data = pd.DataFrame(
            {
                "family": ["F1", "F1"],
                "individual": ["Case1", "Control1"],
                "father": ["0", "0"],
                "mother": ["0", "0"],
                "sex": [1, 2],
                "phenotype": [2, 1],  # 2=case, 1=control
            }
        )

        cfg = {"pseudonymize": True, "pseudonymize_schema": "categorical"}

        result_buffer, pseudonymizer = apply_pseudonymization(
            buffer, ["Case1", "Control1"], cfg, ped_data
        )

        assert pseudonymizer is not None
        assert "CASE_001" in result_buffer[1]
        assert "CONTROL_001" in result_buffer[1]

    def test_apply_no_pseudonymization(self):
        """Test when pseudonymization is disabled."""
        buffer = ["CHROM\tPOS\tGT", "chr1\t1000\tSample1(0/1)"]
        cfg = {"pseudonymize": False}

        result_buffer, pseudonymizer = apply_pseudonymization(buffer, ["Sample1"], cfg)

        assert pseudonymizer is None
        assert result_buffer == buffer  # Unchanged

    def test_apply_custom_missing_pattern(self):
        """Test custom schema without pattern raises error."""
        buffer = ["CHROM\tPOS\tGT"]
        cfg = {"pseudonymize": True, "pseudonymize_schema": "custom"}

        with pytest.raises(ValueError, match="--pseudonymize-pattern required"):
            apply_pseudonymization(buffer, ["S1"], cfg)

    def test_apply_empty_buffer(self):
        """Test handling empty buffer."""
        buffer = []
        cfg = {"pseudonymize": True}

        result_buffer, pseudonymizer = apply_pseudonymization(buffer, ["S1"], cfg)

        assert result_buffer == []
        assert pseudonymizer is not None  # Still created

    def test_apply_header_only_buffer(self):
        """Test handling header-only buffer."""
        buffer = ["CHROM\tPOS\tGT"]
        cfg = {"pseudonymize": True}

        result_buffer, pseudonymizer = apply_pseudonymization(buffer, ["S1"], cfg)

        assert len(result_buffer) == 1
        assert result_buffer[0] == "CHROM\tPOS\tGT"
        assert pseudonymizer is not None


class TestIntegrationScenarios:
    """Test realistic usage scenarios."""

    def test_publication_ready_output(self):
        """Test creating publication-ready pseudonymized output."""
        # Simulate a cohort with cases and controls
        samples = [f"CASE_{i}" for i in range(1, 6)] + [f"CTRL_{i}" for i in range(1, 6)]

        metadata = {}
        for s in samples:
            metadata[s] = {
                "phenotype": "case" if s.startswith("CASE") else "control",
                "sex": "M" if int(s.split("_")[1]) % 2 == 1 else "F",
            }

        # Test categorical schema for publication
        pseudonymizer = create_pseudonymizer("categorical", category_field="phenotype")
        mapping = pseudonymizer.create_mapping(samples, metadata)

        # Verify case/control separation
        case_ids = [v for k, v in mapping.items() if k.startswith("CASE")]
        ctrl_ids = [v for k, v in mapping.items() if k.startswith("CTRL")]

        assert all(id.startswith("CASE_") for id in case_ids)
        assert all(id.startswith("CONTROL_") for id in ctrl_ids)
        assert len(case_ids) == 5
        assert len(ctrl_ids) == 5

    def test_family_study_pseudonymization(self):
        """Test pseudonymization for family studies."""
        # Simulate a trio study
        families = []
        for i in range(1, 4):
            families.extend([f"FAM{i}_child", f"FAM{i}_father", f"FAM{i}_mother"])

        # Use custom schema to maintain family structure
        pattern = "F{family}_{role}_{index:02d}"
        pseudonymizer = create_pseudonymizer("custom", pattern=pattern)

        # Add family metadata
        metadata = {}
        for sample in families:
            parts = sample.split("_")
            metadata[sample] = {"family": parts[0][-1], "role": parts[1].upper()}  # Family number

        mapping = pseudonymizer.create_mapping(families, metadata)

        # Verify family structure is preserved
        for original, pseudo in mapping.items():
            if "FAM1" in original:
                assert pseudo.startswith("F1_")
            if "child" in original:
                assert "_CHILD_" in pseudo
            if "father" in original:
                assert "_FATHER_" in pseudo

    def test_longitudinal_study_pseudonymization(self):
        """Test pseudonymization for longitudinal studies with timepoints."""
        # Simulate samples from multiple timepoints
        samples = []
        for patient in ["P001", "P002", "P003"]:
            for timepoint in ["T0", "T3", "T6", "T12"]:
                samples.append(f"{patient}_{timepoint}")

        # Use custom schema to preserve structure
        pattern = "SUBJ{patient_num:03d}_{timepoint}"
        pseudonymizer = create_pseudonymizer("custom", pattern=pattern)

        # Extract metadata
        metadata = {}
        for sample in samples:
            parts = sample.split("_")
            metadata[sample] = {"patient_num": int(parts[0][1:]), "timepoint": parts[1]}

        mapping = pseudonymizer.create_mapping(samples, metadata)

        # Verify structure is maintained
        assert mapping["P001_T0"] == "SUBJ001_T0"
        assert mapping["P001_T12"] == "SUBJ001_T12"
        assert mapping["P003_T6"] == "SUBJ003_T6"

    def test_multi_cohort_study(self):
        """Test pseudonymization for multi-cohort studies."""
        # Simulate samples from different cohorts
        cohorts = {
            "STUDY_A": ["SA_001", "SA_002", "SA_003"],
            "STUDY_B": ["SB_001", "SB_002"],
            "CLINIC_1": ["C1_P01", "C1_P02", "C1_P03", "C1_P04"],
        }

        all_samples = [s for samples in cohorts.values() for s in samples]

        # Use anonymous schema for complete de-identification
        pseudonymizer = create_pseudonymizer("anonymous", use_hash=True)
        mapping = pseudonymizer.create_mapping(all_samples)

        # Verify all samples get unique hash-based IDs
        assert len(set(mapping.values())) == len(all_samples)
        assert all(v.startswith("ID") and len(v) == 8 for v in mapping.values())


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
