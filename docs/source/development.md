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

- **Black**: Code formatting with 100-character line length
- **isort**: Import statement organization  
- **flake8**: Style checking and error detection

Run these tools manually:

```bash
# Format code
black .

# Sort imports
isort .

# Check style
flake8 .

# Run all pre-commit hooks
pre-commit run --all-files
```

### Configuration

Code quality settings are configured in:

- `pyproject.toml` - Black configuration
- `.pre-commit-config.yaml` - Pre-commit hook definitions
- `setup.cfg` or `pyproject.toml` - flake8 and isort settings

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific test categories
pytest -m unit
pytest -m integration
pytest -m slow

# Run with coverage
pytest --cov=variantcentrifuge --cov-report=html

# Run specific test file
pytest tests/test_cli.py
```

### Test Organization

Tests are organized by functionality:

- `test_cli.py` - Command-line interface tests
- `test_filters.py` - Variant filtering tests
- `test_utils.py` - Utility function tests
- `test_igv_*.py` - IGV integration tests

### Writing Tests

Follow these guidelines when writing tests:

1. **Use descriptive test names** that explain what is being tested
2. **Follow the AAA pattern** (Arrange, Act, Assert)
3. **Use pytest fixtures** for setup and teardown
4. **Mock external dependencies** (file system, external tools)
5. **Test both success and failure cases**

Example test structure:

```python
def test_filter_variants_with_valid_expression():
    # Arrange
    vcf_file = "test_input.vcf"
    filter_expr = "QUAL >= 30"
    
    # Act
    result = apply_filter(vcf_file, filter_expr)
    
    # Assert
    assert result.returncode == 0
    assert "filtered variants" in result.output
```

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