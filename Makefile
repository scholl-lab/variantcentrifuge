.PHONY: help install install-dev clean lint lint-loc format format-check typecheck check-version test test-unit test-integration test-fast test-e2e fixture-1000g ci-check fix build-ext

# Colors for output
CYAN := \033[0;36m
GREEN := \033[0;32m
YELLOW := \033[0;33m
NC := \033[0m

# Python
PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python)
# xdist workers; lower under a memory cap, e.g. make ci-check PYTEST_WORKERS=8
PYTEST_WORKERS ?= auto

help: ## Show this help message
	@echo "$(CYAN)Available targets:$(NC)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  $(GREEN)%-20s$(NC) %s\n", $$1, $$2}'

##@ Installation

install: ## Install the package
	uv sync --locked

install-dev: ## Install with dev dependencies
	uv sync --locked --extra dev

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

check-version: ## Verify version synchronization across metadata files
	@echo "$(CYAN)Checking version synchronization...$(NC)"
	$(PYTHON) scripts/check_version_sync.py
	@echo "$(GREEN)Version check passed$(NC)"

##@ Testing

DAVIES_EXT := variantcentrifuge/_qfc.abi3.so

build-ext: $(DAVIES_EXT) ## Build the Davies C extension in place (fresh worktrees lack it)

# The .so is untracked, so fresh worktrees silently fall back to Liu p-values without this.
$(DAVIES_EXT): variantcentrifuge/_davies_build.py variantcentrifuge/association/data/qfc.cpp
	$(PYTHON) setup.py -q build_ext --inplace

test: $(DAVIES_EXT) ## Run all tests with coverage (same as CI)
	@echo "$(CYAN)Running tests with coverage...$(NC)"
	$(PYTHON) -m pytest -n $(PYTEST_WORKERS) --dist worksteal tests/ \
		--verbose \
		--cov=variantcentrifuge \
		--cov-report=xml \
		--cov-report=term-missing
	@echo "$(GREEN)Tests passed$(NC)"

test-unit: $(DAVIES_EXT) ## Run unit tests only
	$(PYTHON) -m pytest -n $(PYTEST_WORKERS) --dist worksteal -m unit tests/ -v

test-integration: $(DAVIES_EXT) ## Run integration tests only
	$(PYTHON) -m pytest -n $(PYTEST_WORKERS) --dist worksteal -m integration tests/ -v

test-fast: $(DAVIES_EXT) ## Run non-slow, non-integration tests
	$(PYTHON) -m pytest -n $(PYTEST_WORKERS) --dist worksteal -m "not slow and not integration" tests/ --tb=short -q

# Real-data e2e: run inside the variantcentrifuge image (bcftools, snpEff, SnpSift).
E2E_IMAGE ?= ghcr.io/scholl-lab/variantcentrifuge:latest
SNPEFF_DATA ?= /mnt/nas-bioinformatics/bioinformatics/variantcentrifuge/resources/snpeff

fixture-1000g: ## Rebuild tests/fixtures/1000g_cohort from 1000 Genomes (needs SNPEFF_DATA with GRCh38.p14)
	docker run --rm --user $$(id -u):$$(id -g) -e HOME=/tmp -v $(CURDIR):/repo -w /repo \
		-v $(SNPEFF_DATA):/data/snpeff --entrypoint micromamba $(E2E_IMAGE) run -n base \
		bash scripts/build_1000g_fixture.sh

# Mounts only tests/ and pyproject.toml so the image's installed package (with its compiled
# Davies extension) is under test, not the uncompiled checkout. Rebuild the image to test changes.
test-e2e: ## Run real-data e2e tests (1000 Genomes cohort meta-analysis) in the image
	docker run --rm --user root -e PYTHONDONTWRITEBYTECODE=1 \
		-v $(CURDIR)/tests:/work/tests:ro -v $(CURDIR)/pyproject.toml:/work/pyproject.toml:ro -w /work \
		--entrypoint micromamba $(E2E_IMAGE) run -n base bash -c \
		'micromamba install -q -y -n base -c conda-forge pytest >/dev/null && \
		python -m pytest -p no:cacheprovider -o addopts= -m e2e tests/e2e -v'

##@ CI Verification

ci-check: ## Run local lint, LOC, format, type, and fast-test gates
	@echo "$(CYAN)Running complete local CI check$(NC)"
	@echo ""
	@echo "$(CYAN)[1/6] Linting with ruff...$(NC)"
	@$(MAKE) lint
	@echo ""
	@echo "$(CYAN)[2/6] Checking production module line counts...$(NC)"
	@$(MAKE) lint-loc
	@echo ""
	@echo "$(CYAN)[3/6] Checking code format...$(NC)"
	@$(MAKE) format-check
	@echo ""
	@echo "$(CYAN)[4/6] Type checking with mypy...$(NC)"
	@$(MAKE) typecheck
	@echo ""
	@echo "$(CYAN)[5/6] Checking version synchronization...$(NC)"
	@$(MAKE) check-version
	@echo ""
	@echo "$(CYAN)[6/6] Running tests...$(NC)"
	@$(MAKE) test-fast
	@echo ""
	@echo "$(GREEN)ALL LOCAL CI CHECKS PASSED$(NC)"

##@ Quick Commands

fix: format lint ## Auto-fix formatting and linting issues
	@echo "$(GREEN)Auto-fixes applied$(NC)"
