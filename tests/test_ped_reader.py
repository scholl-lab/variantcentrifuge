"""Tests for PED file reader."""

import pytest
import tempfile
import os
from variantcentrifuge.ped_reader import read_pedigree, get_parents, is_affected, get_family_members


class TestPedReader:

    @pytest.fixture
    def valid_ped_content(self):
        """Standard trio PED file content."""
        return """FAM1 father 0 0 1 1
FAM1 mother 0 0 2 1
FAM1 child father mother 1 2"""

    @pytest.fixture
    def complex_ped_content(self):
        """Complex pedigree with multiple families."""
        return """FAM1 F1_father 0 0 1 1
FAM1 F1_mother 0 0 2 1
FAM1 F1_child1 F1_father F1_mother 1 2
FAM1 F1_child2 F1_father F1_mother 2 1
FAM2 F2_father 0 0 1 1
FAM2 F2_mother 0 0 2 1
FAM2 F2_child F2_father F2_mother 1 2"""

    @pytest.fixture
    def invalid_ped_content(self):
        """Invalid PED with wrong number of columns."""
        return """FAM1 father 0 0 1
FAM1 mother 0 0"""

    def test_read_valid_ped_file(self, valid_ped_content):
        """Test reading a valid PED file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(valid_ped_content)
            f.flush()

            try:
                pedigree = read_pedigree(f.name)

                assert len(pedigree) == 3
                assert "father" in pedigree
                assert "mother" in pedigree
                assert "child" in pedigree

                assert pedigree["father"]["family_id"] == "FAM1"
                assert pedigree["father"]["father_id"] == "0"
                assert pedigree["father"]["mother_id"] == "0"
                assert pedigree["father"]["sex"] == "1"
                assert pedigree["father"]["affected_status"] == "1"

                assert pedigree["child"]["father_id"] == "father"
                assert pedigree["child"]["mother_id"] == "mother"
                assert pedigree["child"]["affected_status"] == "2"

            finally:
                os.unlink(f.name)

    def test_read_complex_ped_file(self, complex_ped_content):
        """Test reading a complex PED file with multiple families."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(complex_ped_content)
            f.flush()

            try:
                pedigree = read_pedigree(f.name)

                assert len(pedigree) == 7

                # Check family 1
                fam1_members = [
                    sid for sid, info in pedigree.items() if info["family_id"] == "FAM1"
                ]
                assert len(fam1_members) == 4

                # Check family 2
                fam2_members = [
                    sid for sid, info in pedigree.items() if info["family_id"] == "FAM2"
                ]
                assert len(fam2_members) == 3

            finally:
                os.unlink(f.name)

    def test_read_invalid_ped_file(self, invalid_ped_content):
        """Test reading an invalid PED file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(invalid_ped_content)
            f.flush()

            try:
                with pytest.raises(ValueError, match="Failed to parse PED file"):
                    read_pedigree(f.name)
            finally:
                os.unlink(f.name)

    def test_read_empty_ped_file(self):
        """Test reading an empty PED file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.flush()

            try:
                with pytest.raises(ValueError, match="PED file is empty"):
                    read_pedigree(f.name)
            finally:
                os.unlink(f.name)

    def test_get_parents(self, valid_ped_content):
        """Test getting parent IDs."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(valid_ped_content)
            f.flush()

            try:
                pedigree = read_pedigree(f.name)

                # Test child with parents
                father_id, mother_id = get_parents("child", pedigree)
                assert father_id == "father"
                assert mother_id == "mother"

                # Test individual without parents
                father_id, mother_id = get_parents("father", pedigree)
                assert father_id is None
                assert mother_id is None

                # Test non-existent sample
                father_id, mother_id = get_parents("unknown", pedigree)
                assert father_id is None
                assert mother_id is None

            finally:
                os.unlink(f.name)

    def test_is_affected(self, valid_ped_content):
        """Test checking affected status."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(valid_ped_content)
            f.flush()

            try:
                pedigree = read_pedigree(f.name)

                assert is_affected("child", pedigree) is True
                assert is_affected("father", pedigree) is False
                assert is_affected("mother", pedigree) is False
                assert is_affected("unknown", pedigree) is False

            finally:
                os.unlink(f.name)

    def test_get_family_members(self, complex_ped_content):
        """Test getting family members."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(complex_ped_content)
            f.flush()

            try:
                pedigree = read_pedigree(f.name)

                # Test family 1
                fam1_members = get_family_members("F1_child1", pedigree)
                assert len(fam1_members) == 4
                assert "F1_father" in fam1_members
                assert "F1_mother" in fam1_members
                assert "F1_child1" in fam1_members
                assert "F1_child2" in fam1_members

                # Test family 2
                fam2_members = get_family_members("F2_child", pedigree)
                assert len(fam2_members) == 3
                assert "F2_father" in fam2_members
                assert "F2_mother" in fam2_members
                assert "F2_child" in fam2_members

                # Test non-existent sample
                members = get_family_members("unknown", pedigree)
                assert members == []

            finally:
                os.unlink(f.name)

    def test_ped_with_comments(self):
        """Test reading PED file with comment lines."""
        ped_content = """# This is a comment
# Another comment
FAM1 father 0 0 1 1
FAM1 mother 0 0 2 1
FAM1 child father mother 1 2"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(ped_content)
            f.flush()

            try:
                pedigree = read_pedigree(f.name)
                assert len(pedigree) == 3
                assert "child" in pedigree

            finally:
                os.unlink(f.name)

    def test_ped_with_missing_values(self):
        """Test reading PED file with missing values."""
        ped_content = """FAM1 father . . 1 1
FAM1 mother . . 2 1
FAM1 child father mother . ."""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(ped_content)
            f.flush()

            try:
                pedigree = read_pedigree(f.name)

                assert pedigree["father"]["father_id"] == "0"
                assert pedigree["father"]["mother_id"] == "0"
                assert pedigree["child"]["sex"] == "0"
                assert pedigree["child"]["affected_status"] == "0"

            finally:
                os.unlink(f.name)
