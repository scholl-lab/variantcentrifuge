# syntax=docker/dockerfile:1
# =============================================================================
# variantcentrifuge Docker image
# Build patched Java tools, install the conda environment, then copy only runtime
# artifacts into the production image.
# =============================================================================

# ---------------------------------------------------------------------------
# Stage 1: Java build — verify and package patched SnpEff and SnpSift
# ---------------------------------------------------------------------------
FROM maven:3.9.11-eclipse-temurin-11@sha256:c095e2421eaf3e5cb1573fd0474e68e17062866d454362349e17bbb75f44031e AS java-build

ADD --checksum=sha256:8633ecad8dbf06af19d5002818b65dc9e7b78419b94866dd4b1daed9d1e9d2b2 https://github.com/pcingola/SnpEff/archive/0c5e74f9b6ca6ed3db720177eb1f95b9d47d45f2.tar.gz /tmp/snpeff.tar.gz
ADD --checksum=sha256:74c37b85e74a390a27f122164fd06c5b56131e3fa1eca205636d3de4a8c94934 https://github.com/pcingola/SnpSift/archive/20978614457f14ec7a0c70539d5a7a2b7e754f60.tar.gz /tmp/snpsift.tar.gz

RUN mkdir -p /build/snpeff /build/snpsift /out && \
    tar -xzf /tmp/snpeff.tar.gz --strip-components=1 -C /build/snpeff && \
    tar -xzf /tmp/snpsift.tar.gz --strip-components=1 -C /build/snpsift

COPY docker/java/snpeff-pom.xml /build/snpeff/pom.xml
COPY docker/java/snpsift-pom.xml /build/snpsift/pom.xml
COPY docker/java/src/main/java/org/apache/commons/lang/ArrayUtils.java /build/snpeff/src/main/java/org/apache/commons/lang/ArrayUtils.java
COPY docker/java/src/test/java/org/apache/commons/lang/ArrayUtilsTest.java /build/snpeff/src/test/java/org/apache/commons/lang/ArrayUtilsTest.java
COPY docker/java/assert-runtime-dependencies.sh /usr/local/bin/assert-runtime-dependencies.sh

RUN cd /build/snpeff && \
    mvn -B verify && \
    set -- target/*-jar-with-dependencies.jar && \
    test "$#" -eq 1 && \
    test -f "$1" && \
    verified_sha256=$(sha256sum "$1" | cut -d ' ' -f 1) && \
    printf 'Verified SnpEff fat JAR SHA256: %s\n' "$verified_sha256" && \
    /usr/local/bin/assert-runtime-dependencies.sh /build/snpeff && \
    mvn -B install -DskipTests -Dassembly.skipAssembly=true && \
    set -- target/*-jar-with-dependencies.jar && \
    test "$#" -eq 1 && \
    test -f "$1" && \
    installed_sha256=$(sha256sum "$1" | cut -d ' ' -f 1) && \
    printf 'Installed SnpEff fat JAR SHA256: %s\n' "$installed_sha256" && \
    test "$verified_sha256" = "$installed_sha256" && \
    cp "$1" /out/snpEff.jar

RUN cd /build/snpsift && \
    mvn -B verify && \
    /usr/local/bin/assert-runtime-dependencies.sh /build/snpsift && \
    set -- target/*-jar-with-dependencies.jar && \
    test "$#" -eq 1 && \
    test -f "$1" && \
    cp "$1" /out/SnpSift.jar

# ---------------------------------------------------------------------------
# Stage 2: Conda build — install Python and native runtime dependencies
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:2.8.1-debian12-slim@sha256:c8198d53228ad7cfd7adcc0704e8837f9d1c9327fb363c6a7d3c5b1a51a4b561 AS conda-build

# Install conda environment from the Docker-specific env file first (layer cache)
COPY --chown=$MAMBA_USER:$MAMBA_USER conda/environment-docker.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# A builder-only C++ compiler is required for the optional Davies CFFI module.
# Only /opt/conda is copied into the runtime stage, so the compiler is excluded.
USER root
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -y --no-install-recommends install g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Copy source and install the package (no-deps: all deps satisfied by conda)
# README.md is required by pyproject.toml metadata
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md setup.py /tmp/src/
COPY --chown=$MAMBA_USER:$MAMBA_USER variantcentrifuge/ /tmp/src/variantcentrifuge/
# Activate conda env for RUN commands (micromamba convention)
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-deps --no-cache-dir /tmp/src

RUN /opt/conda/bin/python -m pip check && \
    /opt/conda/bin/python - <<'PY'
import io

import cffi
import numpy as np
import pandas as pd
import pyarrow as pa
import xlsxwriter

from variantcentrifuge import _qfc
from variantcentrifuge.association.backends.davies import (
    _try_load_davies,
    davies_pvalue,
)

table = pa.table({"variant": ["1-100-A-G", "1-200-C-T"], "score": [1.0, 2.0]})
frame = table.to_pandas()
assert frame.to_dict(orient="list") == {
    "variant": ["1-100-A-G", "1-200-C-T"],
    "score": [1.0, 2.0],
}

excel_buffer = io.BytesIO()
with pd.ExcelWriter(excel_buffer, engine="xlsxwriter") as writer:
    frame.to_excel(writer, sheet_name="Results", index=False)
assert excel_buffer.getvalue().startswith(b"PK")

assert cffi.__version__
assert xlsxwriter.__version__
assert _qfc.ffi is not None and _qfc.lib is not None
assert _try_load_davies()
pvalue, ifault = davies_pvalue(1.0, np.array([1.0]))
assert pvalue is not None and 0.0 <= pvalue <= 1.0
print(
    "Python runtime gate:",
    f"cffi={cffi.__version__}",
    f"pyarrow={pa.__version__}",
    f"xlsxwriter={xlsxwriter.__version__}",
    f"davies_pvalue={pvalue}",
    f"davies_ifault={ifault}",
)
PY

# Replace the Bioconda Java tool payloads with the verified patched builds.
COPY --from=java-build --chown=$MAMBA_USER:$MAMBA_USER /out/snpEff.jar /opt/conda/share/snpeff-5.2-3/snpEff.jar
COPY --from=java-build --chown=$MAMBA_USER:$MAMBA_USER /out/SnpSift.jar /opt/conda/share/snpsift-5.2-0/SnpSift.jar

RUN python - <<'PY'
import zipfile

jars = {
    "/opt/conda/share/snpeff-5.2-3/snpEff.jar": "org.snpeff.SnpEff",
    "/opt/conda/share/snpsift-5.2-0/SnpSift.jar": "org.snpsift.SnpSift",
}
for jar_path, expected_main_class in jars.items():
    with zipfile.ZipFile(jar_path) as jar:
        manifest_text = jar.read("META-INF/MANIFEST.MF").decode("utf-8")
    manifest = dict(
        line.split(": ", maxsplit=1)
        for line in manifest_text.splitlines()
        if ": " in line
    )
    assert manifest.get("Main-Class") == expected_main_class, (jar_path, manifest)
    assert manifest.get("Multi-Release") == "true", (jar_path, manifest)
PY

# Keep only the Java launcher and keytool.
RUN java_path=$(readlink -f "$(command -v java)") && \
    jvm_bin=$(dirname "$java_path") && \
    jvm_home=$(dirname "$jvm_bin") && \
    for link in /opt/conda/bin/*; do \
        test -L "$link" || continue; \
        link_target=$(readlink "$link"); \
        case "$link_target" in \
            ../lib/jvm/bin/*|/opt/conda/lib/jvm/bin/*) \
                tool=${link_target##*/}; \
                case "$tool" in \
                    java|keytool) ;; \
                    *) rm -f "$link" ;; \
                esac \
                ;; \
        esac; \
    done && \
    for tool_path in "$jvm_bin"/*; do \
        test -e "$tool_path" || continue; \
        tool=${tool_path##*/}; \
        case "$tool" in \
            java|keytool) ;; \
            *) rm -f "$jvm_bin/$tool" ;; \
        esac; \
    done && \
    rm -rf "$jvm_home/include" "$jvm_home/jmods" \
        "$jvm_home/man" "$jvm_home/demo" "$jvm_home/sample" && \
    rm -f "$jvm_home/src.zip" "$jvm_home/lib/src.zip" && \
    jvm_inventory=$(find "$jvm_bin" -mindepth 1 -maxdepth 1 -printf '%f\n' | sort) && \
    expected_jvm_inventory=$(printf '%s\n%s' java keytool) && \
    if test "$jvm_inventory" != "$expected_jvm_inventory"; then \
        printf 'Unexpected JVM runtime tool inventory:\n%s\n' "$jvm_inventory" >&2; \
        exit 1; \
    fi && \
    for link in /opt/conda/bin/*; do \
        test -L "$link" || continue; \
        link_target=$(readlink "$link"); \
        case "$link_target" in \
            ../lib/jvm/bin/*|/opt/conda/lib/jvm/bin/*) \
                tool=${link_target##*/}; \
                case "$tool" in \
                    java|keytool) ;; \
                    *) printf 'Unexpected conda JVM tool link: %s -> %s\n' \
                        "$link" "$link_target" >&2; exit 1 ;; \
                esac \
                ;; \
        esac; \
    done && \
    printf 'Final JVM bin inventory:\n%s\n' "$jvm_inventory" && \
    rm -rf /opt/conda/pkgs && \
    test ! -e /opt/conda/pkgs && \
    test -f /opt/conda/share/snpeff-5.2-3/snpEff.jar && \
    test -f /opt/conda/share/snpsift-5.2-0/SnpSift.jar && \
    jar_list=$(find /opt/conda -type f \( -name 'snpEff.jar' -o -name 'SnpSift.jar' \) -print | sort) && \
    expected_jars=$(printf '%s\n%s' \
        /opt/conda/share/snpeff-5.2-3/snpEff.jar \
        /opt/conda/share/snpsift-5.2-0/SnpSift.jar) && \
    test "$jar_list" = "$expected_jars" && \
    ! command -v javac

# Normalize modes in the builder so the runtime COPY remains the sole conda layer.
USER root
RUN chmod -R go-w /opt/conda && \
    test -z "$(find /opt/conda -xdev \
        \( ! -type l -a -perm /022 \) -print -quit)"

# ---------------------------------------------------------------------------
# Stage 3: Runtime — lean production image
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:2.8.1-debian12-slim@sha256:c8198d53228ad7cfd7adcc0704e8837f9d1c9327fb363c6a7d3c5b1a51a4b561 AS runtime

USER root
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -y --no-install-recommends dist-upgrade && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.title="variantcentrifuge" \
      org.opencontainers.image.description="Filter, extract, and analyze variants from VCF files" \
      org.opencontainers.image.source="https://github.com/scholl-lab/variantcentrifuge" \
      org.opencontainers.image.license="MIT" \
      org.opencontainers.image.authors="Bernt Popp <bernt.popp.md@gmail.com>"

# The runtime base pre-creates this mount point as 0777; normalize only the
# empty destination directory before copying the immutable environment into it.
RUN chmod go-w /opt/conda

# Copy the fully built and cleaned conda environment as root-owned runtime data.
COPY --from=conda-build --chown=root:root /opt/conda /opt/conda

# Include LICENSE for compliance
COPY --chown=0:0 LICENSE /app/LICENSE

# Copy scoring models and stats configs into the image
COPY --chown=0:0 scoring/ /app/scoring/
COPY --chown=0:0 stats_configs/ /app/stats_configs/

RUN chmod -R go-w /app && \
    mkdir -p /data && \
    chown $MAMBA_USER:$MAMBA_USER /data && \
    chmod 0750 /data && \
    invalid_runtime_path=$(find /opt/conda /app -xdev \
        \( ! -user root -o ! -group root -o \
            \( ! -type l -a -perm /022 \) \) -print -quit) && \
    if test -n "$invalid_runtime_path"; then \
        stat -c 'Unexpected mutable runtime path: %n owner=%U group=%G mode=%a' \
            "$invalid_runtime_path" >&2; \
        exit 1; \
    fi

USER $MAMBA_USER

# Working directory for user data (VCF files, output)
WORKDIR /data

# Health check — uses full path since HEALTHCHECK bypasses ENTRYPOINT
HEALTHCHECK --interval=60s --timeout=10s --start-period=5s --retries=3 \
    CMD ["/opt/conda/bin/variantcentrifuge", "--version"]

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "variantcentrifuge"]
CMD ["--help"]
