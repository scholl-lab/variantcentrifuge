# syntax=docker/dockerfile:1
# =============================================================================
# variantcentrifuge Docker image
# Multi-stage build using micromamba for bioconda tool installation
# =============================================================================

# ---------------------------------------------------------------------------
# Stage 1: Build — install conda environment and pip package
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:2.0-debian12-slim AS build

# Install conda environment from the Docker-specific env file first (layer cache)
COPY --chown=$MAMBA_USER:$MAMBA_USER conda/environment-docker.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Copy source and install the package (no-deps: all deps satisfied by conda)
# README.md is required by pyproject.toml metadata
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /tmp/src/
COPY --chown=$MAMBA_USER:$MAMBA_USER variantcentrifuge/ /tmp/src/variantcentrifuge/
# Activate conda env for RUN commands (micromamba convention)
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-deps --no-cache-dir /tmp/src

# ---------------------------------------------------------------------------
# Stage 2: Runtime — lean image with only the conda env and config data
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:2.0-debian12-slim AS runtime

LABEL org.opencontainers.image.title="variantcentrifuge" \
      org.opencontainers.image.description="Filter, extract, and analyze variants from VCF files" \
      org.opencontainers.image.source="https://github.com/scholl-lab/variantcentrifuge" \
      org.opencontainers.image.license="MIT" \
      org.opencontainers.image.authors="Bernt Popp <bernt.popp.md@gmail.com>"

# Copy the fully built conda environment from the build stage
COPY --from=build /opt/conda /opt/conda

# Include LICENSE for compliance
COPY --chown=$MAMBA_USER:$MAMBA_USER LICENSE /app/LICENSE

# Copy scoring models and stats configs into the image
COPY --chown=$MAMBA_USER:$MAMBA_USER scoring/ /app/scoring/
COPY --chown=$MAMBA_USER:$MAMBA_USER stats_configs/ /app/stats_configs/

# Working directory for user data (VCF files, output)
WORKDIR /data

# Health check — uses full path since HEALTHCHECK bypasses ENTRYPOINT
HEALTHCHECK --interval=60s --timeout=10s --start-period=5s --retries=3 \
    CMD ["/opt/conda/bin/variantcentrifuge", "--version"]

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "variantcentrifuge"]
CMD ["--help"]
