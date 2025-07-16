# Contributing to VariantCentrifuge

We welcome contributions to VariantCentrifuge! This document provides guidelines for contributing to the project.

## Getting Started

### Prerequisites

Before contributing, ensure you have:

- Python 3.7+ installed
- Git configured with your GitHub account
- External bioinformatics tools (bcftools, snpEff, SnpSift, bedtools) installed
- Familiarity with the VariantCentrifuge codebase and documentation

### Development Setup

1. **Fork the repository** on GitHub
2. **Clone your fork locally:**
   ```bash
   git clone https://github.com/YOUR-USERNAME/variantcentrifuge.git
   cd variantcentrifuge
   ```
3. **Add upstream remote:**
   ```bash
   git remote add upstream https://github.com/scholl-lab/variantcentrifuge.git
   ```
4. **Set up development environment:**
   ```bash
   mamba env create -f conda/environment.yml
   mamba activate annotation
   pip install -e .
   pre-commit install
   ```

## Ways to Contribute

### Reporting Bugs

When reporting bugs, please include:

- **Clear description** of the problem
- **Steps to reproduce** the issue
- **Expected vs actual behavior**
- **Environment details** (OS, Python version, tool versions)
- **Log output** or error messages
- **Example data** if possible (anonymized)

Use the [bug report template](https://github.com/scholl-lab/variantcentrifuge/issues/new?template=bug_report.md) on GitHub.

### Suggesting Features

For feature requests, please provide:

- **Clear description** of the proposed feature
- **Use case** and motivation
- **Proposed implementation** (if you have ideas)
- **Examples** of how it would be used

Use the [feature request template](https://github.com/scholl-lab/variantcentrifuge/issues/new?template=feature_request.md) on GitHub.

### Contributing Code

#### Types of Contributions

- **Bug fixes** - Fix identified issues
- **New features** - Add functionality
- **Performance improvements** - Optimize existing code
- **Documentation** - Improve or expand documentation
- **Tests** - Add or improve test coverage
- **Refactoring** - Improve code structure without changing functionality

#### Development Workflow

1. **Create an issue** for significant changes (optional for small fixes)
2. **Create a feature branch:**
   ```bash
   git checkout -b feature/descriptive-name
   ```
3. **Make your changes** following the coding standards
4. **Add tests** for new functionality
5. **Update documentation** as needed
6. **Run the test suite:**
   ```bash
   pytest
   black .
   isort .
   flake8 .
   ```
7. **Commit your changes** with clear messages
8. **Push to your fork** and create a pull request

#### Code Standards

- **Follow PEP 8** with 100-character line length
- **Use type hints** for function signatures
- **Write docstrings** in Google format for all public functions
- **Add unit tests** for new functionality
- **Keep functions focused** and modular
- **Use meaningful variable names**
- **Handle errors gracefully** with appropriate logging

#### Testing Guidelines

- **Write tests for all new code**
- **Ensure tests pass** before submitting PR
- **Use pytest fixtures** for setup/teardown
- **Mock external dependencies** (files, external tools)
- **Test both success and failure cases**
- **Aim for high test coverage**

Example test:

```python
def test_extract_variants_success(mock_run_command):
    """Test successful variant extraction."""
    # Arrange
    mock_run_command.return_value = CompletedProcess(
        args=[], returncode=0, stdout="Success"
    )

    # Act
    result = extract_variants("input.vcf", "output.vcf", "chr1:1000-2000")

    # Assert
    assert result == "output.vcf"
    mock_run_command.assert_called_once()
```

#### Documentation Guidelines

- **Update user documentation** for new features
- **Add docstrings** to all new functions and classes
- **Include code examples** in docstrings where helpful
- **Update API documentation** automatically via docstrings
- **Test documentation builds** locally before submitting

Example docstring:

```python
def filter_variants(vcf_file: str, filter_expr: str, output_file: str) -> str:
    """Filter variants using SnpSift filter expression.

    Args:
        vcf_file: Path to input VCF file
        filter_expr: SnpSift filter expression
        output_file: Path to output filtered VCF file

    Returns:
        Path to the filtered output file

    Raises:
        FileNotFoundError: If input VCF file doesn't exist
        RuntimeError: If SnpSift command fails

    Example:
        >>> filter_variants("input.vcf", "QUAL >= 30", "filtered.vcf")
        "filtered.vcf"
    """
```

## Pull Request Process

### Before Submitting

- [ ] Code follows project style guidelines
- [ ] Tests pass locally (`pytest`)
- [ ] Code is properly formatted (`black .` and `isort .`)
- [ ] Linting passes (`flake8 .`)
- [ ] Documentation builds successfully
- [ ] Commit messages are clear and descriptive
- [ ] Changes are focused and atomic

### Pull Request Template

When creating a pull request, include:

- **Description** of changes made
- **Motivation** for the changes
- **Testing** performed
- **Breaking changes** (if any)
- **Related issues** (use "Fixes #123" format)

### Review Process

1. **Automated checks** run on all PRs (CI/CD)
2. **Code review** by maintainers
3. **Discussion** and iteration if needed
4. **Approval** and merge by maintainers

### After Submission

- **Respond to feedback** promptly
- **Make requested changes** in the same branch
- **Keep the PR up to date** with main branch if needed
- **Participate in discussions** constructively

## Community Guidelines

### Code of Conduct

We are committed to fostering an inclusive and welcoming community. Please:

- **Be respectful** and considerate
- **Welcome newcomers** and help them get started
- **Focus on what's best** for the community
- **Show empathy** towards other community members
- **Be collaborative** and constructive in discussions

### Communication

- **GitHub Issues** - Bug reports, feature requests, and project discussions
- **GitHub Discussions** - General questions and community discussions
- **Pull Request Comments** - Code-specific discussions
- **Email** - Contact maintainers for sensitive issues

### Recognition

Contributors are recognized through:

- **Git commit history** - Your contributions are permanently recorded
- **Changelog entries** - Significant contributions noted in releases
- **Contributors file** - Recognition of ongoing contributors
- **GitHub contributors graph** - Visual representation of contributions

## Getting Help

### For Contributors

- **Read the documentation** thoroughly before asking questions
- **Search existing issues** to see if your question has been answered
- **Check the development guide** for technical details
- **Ask specific questions** with context and examples

### For Maintainers

If you become a regular contributor, you may be invited to join the maintainer team. Maintainer responsibilities include:

- **Reviewing pull requests** thoroughly and constructively
- **Triaging issues** and labeling them appropriately
- **Helping newcomers** get started with contributing
- **Making decisions** about project direction and priorities
- **Maintaining code quality** and project standards

## Thank You!

Thank you for considering contributing to VariantCentrifuge! Your contributions help make this tool better for the entire genomics community.

Whether you're fixing a typo, adding a feature, or improving documentation, every contribution is valuable and appreciated.
