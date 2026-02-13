.PHONY: help install install-dev clean lint format format-check typecheck test test-unit test-integration test-fast ci-check fix

# Colors for output
CYAN := \033[0;36m
GREEN := \033[0;32m
YELLOW := \033[0;33m
NC := \033[0m

# Python
PYTHON := python

help: ## Show this help message
	@echo "$(CYAN)Available targets:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  $(GREEN)%-20s$(NC) %s\n", $$1, $$2}'

##@ Installation

install: ## Install the package
	uv pip install -e .

install-dev: ## Install with dev dependencies
	uv pip install -e ".[dev]"

##@ Cleaning

clean: ## Clean build artifacts and cache files
	@echo "$(YELLOW)Cleaning build artifacts...$(NC)"
	rm -rf build/ dist/ *.egg-info .pytest_cache .mypy_cache .ruff_cache htmlcov/ .coverage coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	@echo "$(GREEN)Cleaned$(NC)"

##@ Code Quality

lint: ## Run ruff linter (same as CI)
	@echo "$(CYAN)Running ruff linter...$(NC)"
	$(PYTHON) -m ruff check variantcentrifuge/ tests/
	@echo "$(GREEN)Lint passed$(NC)"

format: ## Format code with ruff
	@echo "$(CYAN)Formatting code with ruff...$(NC)"
	$(PYTHON) -m ruff format variantcentrifuge/ tests/
	$(PYTHON) -m ruff check --fix variantcentrifuge/ tests/
	@echo "$(GREEN)Code formatted$(NC)"

format-check: ## Check code formatting (same as CI)
	@echo "$(CYAN)Checking code formatting...$(NC)"
	$(PYTHON) -m ruff format --check --diff variantcentrifuge/ tests/
	@echo "$(GREEN)Format check passed$(NC)"

typecheck: ## Run mypy type checker (informational during gradual adoption)
	@echo "$(CYAN)Running mypy type checker...$(NC)"
	-$(PYTHON) -m mypy variantcentrifuge/
	@echo "$(YELLOW)Type check complete (non-blocking during gradual adoption)$(NC)"

##@ Testing

test: ## Run all tests with coverage (same as CI)
	@echo "$(CYAN)Running tests with coverage...$(NC)"
	$(PYTHON) -m pytest tests/ \
		--verbose \
		--cov=variantcentrifuge \
		--cov-report=xml \
		--cov-report=term-missing
	@echo "$(GREEN)Tests passed$(NC)"

test-unit: ## Run unit tests only
	$(PYTHON) -m pytest -m unit tests/ -v

test-integration: ## Run integration tests only
	$(PYTHON) -m pytest -m integration tests/ -v

test-fast: ## Run non-slow, non-integration tests
	$(PYTHON) -m pytest -m "not slow and not integration" tests/ --tb=short -q

##@ CI Verification

ci-check: ## Run ALL CI checks locally (lint, format-check, typecheck, test-fast) - matches GitHub Actions
	@echo "$(CYAN)Running Complete CI Check (mirroring GitHub Actions)$(NC)"
	@echo ""
	@echo "$(CYAN)[1/4] Linting with ruff...$(NC)"
	@$(MAKE) lint
	@echo ""
	@echo "$(CYAN)[2/4] Checking code format...$(NC)"
	@$(MAKE) format-check
	@echo ""
	@echo "$(CYAN)[3/4] Type checking with mypy...$(NC)"
	@$(MAKE) typecheck
	@echo ""
	@echo "$(CYAN)[4/4] Running tests...$(NC)"
	@$(MAKE) test-fast
	@echo ""
	@echo "$(GREEN)ALL CI CHECKS PASSED - Safe to push!$(NC)"

##@ Quick Commands

fix: format lint ## Auto-fix formatting and linting issues
	@echo "$(GREEN)Auto-fixes applied$(NC)"
