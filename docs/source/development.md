# Development Guide

This guide provides information for developers who want to contribute to VariantCentrifuge.

## Development Setup

### Prerequisites

- Python 3.7+
- Git
- External bioinformatics tools (bcftools, snpEff, SnpSift, bedtools)

### Setting Up Development Environment

1. **Clone the repository:**
   ```bash
   git clone https://github.com/scholl-lab/variantcentrifuge.git
   cd variantcentrifuge
   ```

2. **Create development environment:**
   ```bash
   # Using conda (recommended)
   mamba env create -f conda/environment.yml
   mamba activate annotation

   # Or using pip with virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r docs/requirements.txt
   ```

3. **Install in development mode:**
   ```bash
   pip install -e .
   ```

4. **Install pre-commit hooks:**
   ```bash
   pre-commit install
   ```

## Code Quality

VariantCentrifuge maintains high code quality standards using automated tools:

### Formatting and Linting

VariantCentrifuge uses a comprehensive linting setup to maintain code quality:

#### Tools Used

- **Black**: Code formatting with 100-character line length
- **isort**: Import statement organization (compatible with Black)
- **flake8**: Style checking and error detection with docstring requirements
- **pre-commit**: Automated quality checks before commits

#### Running Linting Tools

```bash
# Format code with Black (100 character line length)
black .

# Sort imports with isort (compatible with Black)
isort .

# Check code style with flake8
flake8 .

# Run all pre-commit hooks on all files
pre-commit run --all-files

# Install pre-commit hooks (run once after cloning)
pre-commit install
```

#### Automated Quality Assurance

Pre-commit hooks automatically run Black, isort, and flake8 on every commit to maintain code quality. The hooks will:

1. **Format code** automatically with Black
2. **Sort imports** according to isort configuration
3. **Check style** and fail commit if flake8 finds issues
4. **Require docstrings** for all functions and classes

#### Linting Configuration

Linting behavior is configured in:

- **pyproject.toml**: Black configuration (line length, target versions)
- **setup.cfg**: flake8 and isort configuration
- **.pre-commit-config.yaml**: Pre-commit hook definitions and versions

#### Common Linting Issues and Solutions

| Issue | Solution |
|-------|----------|
| Line too long | Black will auto-fix, or break long lines manually |
| Missing docstring | Add Google-style docstring to function/class |
| Import order | Run `isort .` to automatically fix |
| Unused imports | Remove unused imports or add `# noqa` comment |
| Trailing whitespace | Pre-commit will automatically remove |

#### Docstring Requirements

All functions and classes must have docstrings following Google style:

```python
def example_function(param1: str, param2: int = 10) -> bool:
    """
    Brief description of the function.

    Longer description if needed, explaining the function's purpose,
    algorithm, or important implementation details.

    Parameters
    ----------
    param1 : str
        Description of the first parameter.
    param2 : int, default=10
        Description of the second parameter with default value.

    Returns
    -------
    bool
        Description of the return value.

    Raises
    ------
    ValueError
        When param1 is empty or param2 is negative.

    Examples
    --------
    >>> example_function("test", 5)
    True
    """
    # Implementation here
    pass
```

### Configuration

Code quality settings are configured in:

- `pyproject.toml` - Black configuration
- `.pre-commit-config.yaml` - Pre-commit hook definitions
- `setup.cfg` or `pyproject.toml` - flake8 and isort settings

## Testing

### Running Tests

VariantCentrifuge has a comprehensive test suite with multiple categories and detailed coverage.

#### Basic Test Commands

```bash
# Run all tests with verbose output and colored results
pytest

# Run tests with maximum verbosity
pytest -v

# Run tests and stop on first failure
pytest -x

# Run tests in parallel (if pytest-xdist installed)
pytest -n auto
```

#### Test Categories

Tests are organized using pytest markers:

```bash
# Run only unit tests (fast, isolated tests)
pytest -m unit

# Run only integration tests (test component interactions)
pytest -m integration

# Run slow tests specifically (may involve large files or external tools)
pytest -m slow

# Run all tests except slow ones
pytest -m "not slow"

# Combine markers
pytest -m "unit or integration"
```

#### Test Coverage

```bash
# Run tests with coverage reporting
pytest --cov=variantcentrifuge

# Generate HTML coverage report
pytest --cov=variantcentrifuge --cov-report=html

# Generate XML coverage report (for CI)
pytest --cov=variantcentrifuge --cov-report=xml

# Show missing lines in coverage
pytest --cov=variantcentrifuge --cov-report=term-missing
```

#### Running Specific Tests

```bash
# Run specific test file
pytest tests/test_cli.py

# Run specific test class
pytest tests/test_cli.py::TestCLI

# Run specific test function
pytest tests/test_cli.py::TestCLI::test_basic_functionality

# Run tests matching a pattern
pytest -k "test_filter"

# Run tests with specific substring in name
pytest -k "vcf and not slow"
```

#### Test Output and Debugging

```bash
# Show print statements and logging output
pytest -s

# Drop into debugger on failures
pytest --pdb

# Drop into debugger on first failure
pytest --pdb -x

# Capture only failures (show output only for failed tests)
pytest --tb=short

# Show local variables on failures
pytest --tb=long
```

### Test Organization

Tests are organized by functionality in the `tests/` directory:

#### Test Files Structure

- **`test_cli.py`** - Command-line interface tests
- **`test_filters.py`** - Variant filtering tests  
- **`test_gene_lists.py`** - Gene list processing tests
- **`test_igv.py`** - IGV integration tests
- **`test_utils.py`** - Utility function tests
- **`conftest.py`** - Pytest configuration and shared fixtures
- **`pytest.ini`** - Pytest configuration file

#### Test Data

Test data is organized in subdirectories:

- **`tests/data/`** - Sample VCF files, configuration files, expected outputs
- **`tests/fixtures/`** - Pytest fixtures for common test setup
- **`tests/integration/`** - Integration test data and scenarios

#### Test Categories and Markers

| Marker | Purpose | Examples |
|--------|---------|----------|
| `@pytest.mark.unit` | Fast, isolated unit tests | Function input/output validation |
| `@pytest.mark.integration` | Component interaction tests | Pipeline workflow tests |
| `@pytest.mark.slow` | Tests that take significant time | Large file processing, external tools |
| `@pytest.mark.external_tools` | Tests requiring external tools | bcftools, snpEff integration tests |

#### Configuration Files

- **`pytest.ini`**: Main pytest configuration
  ```ini
  [tool:pytest]
  markers =
      unit: Unit tests (fast, isolated)
      integration: Integration tests
      slow: Slow tests that may take significant time
      external_tools: Tests requiring external tools
  testpaths = tests
  python_files = test_*.py
  python_classes = Test*
  python_functions = test_*
  ```

### Writing Tests

Follow these comprehensive guidelines when writing tests:

#### Test Design Principles

1. **Use descriptive test names** that explain what is being tested
2. **Follow the AAA pattern** (Arrange, Act, Assert)
3. **Use pytest fixtures** for setup and teardown
4. **Mock external dependencies** (file system, external tools)
5. **Test both success and failure cases**
6. **Keep tests independent** - each test should be able to run in isolation
7. **Use appropriate markers** to categorize tests

#### Test Naming Conventions

Use descriptive names that explain the scenario:

```python
# Good test names
def test_filter_variants_with_valid_quality_threshold():
def test_filter_variants_raises_error_with_invalid_expression():
def test_gene_normalization_handles_case_insensitive_input():
def test_vcf_extraction_preserves_header_order():

# Avoid generic names
def test_filter():  # Too vague
def test_success():  # Doesn't describe what succeeds
def test_error():   # Doesn't describe what causes error
```

#### Test Structure and Examples

**Basic Unit Test:**
```python
import pytest
from variantcentrifuge.filters import apply_filter

@pytest.mark.unit
def test_filter_variants_with_valid_quality_threshold():
    # Arrange
    vcf_file = "test_input.vcf"
    filter_expr = "QUAL >= 30"
    expected_variant_count = 5
    
    # Act
    result = apply_filter(vcf_file, filter_expr)
    
    # Assert
    assert result.returncode == 0
    assert "filtered variants" in result.output
    assert result.variant_count == expected_variant_count
```

**Integration Test with Fixtures:**
```python
@pytest.mark.integration
def test_full_pipeline_with_sample_data(sample_vcf, temp_output_dir):
    # Arrange
    config = {
        "gene_name": "BRCA1",
        "filters": ["rare", "coding"],
        "output_format": "tsv"
    }
    
    # Act
    result = run_pipeline(sample_vcf, config, temp_output_dir)
    
    # Assert
    assert result.success
    assert (temp_output_dir / "output.tsv").exists()
    assert result.variant_count > 0
```

**Error Testing:**
```python
@pytest.mark.unit
def test_filter_raises_error_with_invalid_expression():
    # Arrange
    vcf_file = "test_input.vcf"
    invalid_filter = "INVALID_FIELD >= 30"
    
    # Act & Assert
    with pytest.raises(ValueError, match="Invalid filter expression"):
        apply_filter(vcf_file, invalid_filter)
```

**Parametrized Tests:**
```python
@pytest.mark.unit
@pytest.mark.parametrize("input_gene,expected_output", [
    ("brca1", "BRCA1"),
    ("BRCA1", "BRCA1"),
    ("BrCa1", "BRCA1"),
    ("tp53", "TP53"),
])
def test_gene_name_normalization(input_gene, expected_output):
    # Act
    result = normalize_gene_name(input_gene)
    
    # Assert
    assert result == expected_output
```

#### Mocking External Dependencies

```python
from unittest.mock import Mock, patch

@pytest.mark.unit
@patch('variantcentrifuge.utils.subprocess.run')
def test_run_command_handles_tool_failure(mock_subprocess):
    # Arrange
    mock_subprocess.return_value.returncode = 1
    mock_subprocess.return_value.stderr = "Tool error"
    
    # Act & Assert
    with pytest.raises(subprocess.CalledProcessError):
        run_command(["failing_tool", "--option"])
```

#### Fixtures for Common Setup

```python
# In conftest.py
@pytest.fixture
def sample_vcf(tmp_path):
    """Create a sample VCF file for testing."""
    vcf_content = '''
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
CHR1	100	.	A	T	50	PASS	AC=1
'''
    vcf_file = tmp_path / "sample.vcf"
    vcf_file.write_text(vcf_content)
    return str(vcf_file)

@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir
```

#### Test Quality Checklist

Before submitting tests, verify:

- [ ] Test names clearly describe the scenario
- [ ] Tests follow AAA pattern
- [ ] External dependencies are mocked appropriately
- [ ] Both success and failure cases are covered
- [ ] Tests use appropriate markers (`@pytest.mark.unit`, etc.)
- [ ] Tests are independent and can run in any order
- [ ] Test data is properly cleaned up (use fixtures)
- [ ] Assertions are specific and meaningful
- [ ] Coverage includes edge cases and error conditions

## Documentation

### Building Documentation Locally

```bash
cd docs
pip install -r requirements.txt
sphinx-build -b html source build/html
```

Open `docs/build/html/index.html` in your browser.

### Documentation Structure

- `docs/source/` - Documentation source files
- `docs/source/api/` - Auto-generated API documentation
- `docs/source/guides/` - User guides and tutorials

### Writing Documentation

- Use **Markdown** for all documentation files
- Follow **Google-style docstrings** in Python code
- Include **code examples** in user-facing documentation
- Update **API documentation** when adding new modules or functions

## Architecture Overview

### Core Design Principles

1. **Modularity** - Each module has a single, well-defined responsibility
2. **Separation of Concerns** - Clear boundaries between data processing, analysis, and reporting
3. **Testability** - Code is structured to enable comprehensive testing
4. **Extensibility** - New functionality can be added without breaking existing code

### Key Components

```
variantcentrifuge/
├── cli.py              # Command-line interface
├── pipeline.py         # Main workflow orchestration
├── config.py           # Configuration management
├── filters.py          # Variant filtering
├── extractor.py        # Field extraction
├── analyze_variants.py # Statistical analysis
├── generate_*_report.py # Report generation
└── utils.py            # Common utilities
```

### Data Flow

1. **Input Validation** - CLI validates arguments and files
2. **Configuration Loading** - Load and merge config from file and CLI
3. **Gene Processing** - Convert genes to BED regions
4. **Variant Extraction** - Extract variants from VCF using external tools
5. **Filtering** - Apply SnpSift filters
6. **Field Extraction** - Extract specified fields
7. **Analysis** - Perform statistical analysis and gene burden testing
8. **Report Generation** - Create output files and reports

## Contributing

### Workflow

1. **Fork** the repository on GitHub
2. **Create a feature branch** from main: `git checkout -b feature/my-feature`
3. **Make changes** following the coding standards
4. **Add tests** for new functionality
5. **Update documentation** as needed
6. **Run tests and quality checks**
7. **Commit changes** with descriptive messages
8. **Push** to your fork and create a pull request

### Pull Request Guidelines

- **Describe the change** clearly in the PR description
- **Reference any issues** that the PR addresses
- **Include tests** for new functionality
- **Update documentation** for user-facing changes
- **Ensure CI passes** (tests, linting, documentation build)

### Commit Message Format

Use clear, concise commit messages:

```
type(scope): brief description

Longer explanation if needed

Fixes #123
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

## Release Process

### Version Management

Versions are managed in `variantcentrifuge/version.py` following semantic versioning (MAJOR.MINOR.PATCH).

### Creating a Release

1. **Update version** in `version.py`
2. **Update CHANGELOG.md** with release notes
3. **Create a release tag**: `git tag -a v0.5.0 -m "Release v0.5.0"`
4. **Push tag**: `git push origin v0.5.0`
5. **Create GitHub release** with release notes

## Debugging

### Common Development Issues

1. **Import errors** - Check that package is installed in development mode
2. **Test failures** - Ensure external tools are available in PATH
3. **Documentation build failures** - Check for syntax errors in docstrings
4. **Pre-commit failures** - Run tools manually to fix issues

### Debugging Tools

- **pdb** - Python debugger for interactive debugging
- **pytest --pdb** - Drop into debugger on test failures  
- **logging** - Use appropriate log levels for debugging
- **--keep-intermediates** - Retain intermediate files for inspection

## Getting Help

- **GitHub Issues** - Report bugs and request features
- **GitHub Discussions** - Ask questions and discuss development
- **Code Review** - Request feedback on complex changes
- **Documentation** - Check existing docs before asking questions