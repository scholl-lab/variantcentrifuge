.PHONY: help install install-dev clean lint lint-loc format format-check typecheck test test-unit test-integration test-fast ci-check fix

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

lint: ## Run ruff linter across the repository
	@echo "$(CYAN)Running ruff linter...$(NC)"
	$(PYTHON) -m ruff check .
	@echo "$(GREEN)Lint passed$(NC)"

lint-loc: ## Enforce the 600-line production Python module budget
	$(PYTHON) scripts/check_file_size.py variantcentrifuge --allowlist .loc-allowlist

format: ## Format code with ruff
	@echo "$(CYAN)Formatting code with ruff...$(NC)"
	$(PYTHON) -m ruff format .
	$(PYTHON) -m ruff check --fix .
	@echo "$(GREEN)Code formatted$(NC)"

format-check: ## Check code formatting across the repository
	@echo "$(CYAN)Checking code formatting...$(NC)"
	$(PYTHON) -m ruff format --check --diff .
	@echo "$(GREEN)Format check passed$(NC)"

typecheck: ## Run mypy type checker
	@echo "$(CYAN)Running mypy type checker...$(NC)"
	$(PYTHON) -m mypy variantcentrifuge/
	@echo "$(GREEN)Type check passed$(NC)"

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

ci-check: ## Run local lint, LOC, format, type, and fast-test gates
	@echo "$(CYAN)Running complete local CI check$(NC)"
	@echo ""
	@echo "$(CYAN)[1/5] Linting with ruff...$(NC)"
	@$(MAKE) lint
	@echo ""
	@echo "$(CYAN)[2/5] Checking production module line counts...$(NC)"
	@$(MAKE) lint-loc
	@echo ""
	@echo "$(CYAN)[3/5] Checking code format...$(NC)"
	@$(MAKE) format-check
	@echo ""
	@echo "$(CYAN)[4/5] Type checking with mypy...$(NC)"
	@$(MAKE) typecheck
	@echo ""
	@echo "$(CYAN)[5/5] Running tests...$(NC)"
	@$(MAKE) test-fast
	@echo ""
	@echo "$(GREEN)ALL LOCAL CI CHECKS PASSED$(NC)"

##@ Quick Commands

fix: format lint ## Auto-fix formatting and linting issues
	@echo "$(GREEN)Auto-fixes applied$(NC)"
