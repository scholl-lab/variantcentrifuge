# Code-Scanning Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task. Use
> `superpowers:test-driven-development` for every behavior change and
> `superpowers:verification-before-completion` before claiming success. Steps use checkbox
> (`- [ ]`) syntax for tracking.

**Goal:** Remediate every vendor-fixed vulnerability reported for the VariantCentrifuge
container, preserve complete visibility of vendor-unfixed findings, retain SnpEff/SnpSift 5.2
behavior, and gate pull requests before an image can be published or signed.

**Architecture:** Build security-patched SnpEff and SnpSift 5.2 fat JARs from checksum-verified,
commit-pinned sources in an isolated Maven stage. Install exact Bioconda wrapper/config builds in a
current digest-pinned Micromamba stage, replace only their vulnerable JAR payloads, remove caches,
and copy that environment into a Debian runtime that has applied all available package upgrades.
Exercise that exact production image with deterministic golden fixtures, then produce a complete
Trivy JSON audit and a separate fixed-only SARIF/gate before any registry push.

**Tech Stack:** Docker BuildKit, Micromamba 2.8.1 on Debian 12, conda-forge/Bioconda, Maven 3.9.11,
Eclipse Temurin Java 11, SnpEff/SnpSift 5.2, Bash, Python 3.12, pytest, GitHub Actions, Trivy 0.70.0,
GHCR, and Cosign.

**Global Constraints:**

- The authoritative design is
  `docs/superpowers/specs/2026-07-19-code-scanning-remediation-design.md`.
- The alert baseline is commit `335d971dfa813b8a3c3d126cf586ecce70c90bd5`: 203 open alerts,
  80 vendor-fixed and 123 vendor-unfixed, read on 2026-07-19.
- Never dismiss, PATCH, or otherwise mutate GitHub code-scanning alerts. Let a new SARIF analysis
  under `.github/workflows/docker.yml:scan` drive lifecycle changes.
- Never add `.trivyignore`, VEX exemptions, severity omissions, package/path omissions, or
  `continue-on-error` to make a fixed vulnerability pass.
- Scan `UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL` for both `os` and `library` vulnerabilities.
- Preserve exact SnpEff/SnpSift 5.2 CLI behavior; do not substitute a 5.4 upgrade.
- Before every commit step below, obtain explicit user permission. If permission has not been
  granted, leave the verified changes uncommitted and continue to the next task.
- Do not publish or sign during local implementation. Publication/signing is exercised only by the
  post-merge/tag workflow after its security gate passes.

---

## File Map

### Add

- `docker/java/snpeff-pom.xml` — full security overlay for the pinned SnpEff 5.2 source.
- `docker/java/snpsift-pom.xml` — full security overlay for the pinned SnpSift 5.2 source.
- `docker/java/assert-runtime-dependencies.sh` — fail a Maven build if a vulnerable/banned runtime
  coordinate survives dependency resolution.
- `docker/java/src/main/java/org/apache/commons/lang/ArrayUtils.java` — the one-method legacy binary
  compatibility class required by `mmtf-codec`.
- `docker/java/src/test/java/org/apache/commons/lang/ArrayUtilsTest.java` — compatibility contract.
- `scripts/test_container_image.sh` — production-image behavioral and filesystem contract.
- `scripts/summarize_trivy.py` — deterministic complete/actionable vulnerability accounting.
- `tests/unit/test_container_security_contract.py` — static policy tests for conda, Maven,
  Dockerfile, and workflow files.
- `tests/unit/test_summarize_trivy.py` — unit tests for Trivy JSON accounting.
- `tests/fixtures/container/snpeff/snpEff.config` — local test-genome configuration.
- `tests/fixtures/container/snpeff/data/testGenome/genes.gff` — 60-base coding transcript.
- `tests/fixtures/container/snpeff/data/testGenome/sequences.fa` — deterministic reference.
- `tests/fixtures/container/snpeff/input.vcf` — one missense-to-stop test variant.
- `tests/fixtures/container/snpeff/expected.vcf` — canonical normalized SnpEff output.
- `tests/fixtures/container/snpsift/input.vcf` — one passing and one rejected record.
- `tests/fixtures/container/snpsift/expected.vcf` — canonical normalized SnpSift output.
- `tests/fixtures/container/snpsift/assert_fail_closed.py` — invokes VariantCentrifuge's real
  malformed-filter guard in the container.

### Modify

- `conda/environment-docker.yml` — exact tool builds and current secure `pip` floor.
- `Dockerfile` — pinned Java builder, current pinned Micromamba stages, JAR replacement, cache
  removal, Debian upgrades, and explicit final non-root user.
- `.github/workflows/docker.yml` — build/load/test/scan before push, two reports, fixed-only gate,
  digest output, and signing dependency.
- `docs/source/installation.md` — explain the container's patched 5.2 tools and vulnerability-report
  policy.
- `docs/source/changelog.md` — record the security remediation and PR gate.

### Preserve Unchanged

- `.dockerignore` already includes the new `docker/` directory in build context and intentionally
  excludes `tests/` and `scripts/`; smoke fixtures/scripts run from the checked-out host after the
  build and are mounted into the container, so no ignore-rule change is needed.
- VariantCentrifuge Python pipeline modules remain unchanged.

---

## Task 1: Lock the Container Dependency Contract

**Files:**

- Add `tests/unit/test_container_security_contract.py`
- Modify `conda/environment-docker.yml`

**Interfaces:**

- `conda/environment-docker.yml` must contain exact scalar dependency entries
  `snpeff=5.2=hdfd78af_3`, `snpsift=5.2=hdfd78af_0`, and `pip>=26.1.2`.
- The test module reads repository files through `Path`; it must not import Docker, conda, or a YAML
  parser.

- [ ] Add the failing conda contract tests first:

```python
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_docker_environment_pins_release_era_java_tool_builds() -> None:
    environment = _text("conda/environment-docker.yml")
    assert "  - snpeff=5.2=hdfd78af_3\n" in environment
    assert "  - snpsift=5.2=hdfd78af_0\n" in environment


def test_docker_environment_has_secure_pip_floor() -> None:
    environment = _text("conda/environment-docker.yml")
    assert "  - pip>=26.1.2\n" in environment
```

- [ ] Run the focused test and confirm all three assertions fail against the current unqualified
  package entries:

```bash
pytest tests/unit/test_container_security_contract.py -q
```

Expected: `2 failed`.

- [ ] Change only the relevant conda entries:

```yaml
  - snpeff=5.2=hdfd78af_3
  - snpsift=5.2=hdfd78af_0
  - pip>=26.1.2
```

Keep the existing `pip:` subsection for `intervaltree` and `smart-open`; the top-level pin controls
the installer executable while the subsection controls packages installed through it.

- [ ] Run the focused test again.

Expected: `2 passed`.

- [ ] If explicit commit permission exists, commit this task:

```bash
git add conda/environment-docker.yml tests/unit/test_container_security_contract.py
git commit -m "test: lock secure container dependency versions"
```

Otherwise, leave the changes uncommitted.

---

## Task 2: Build Security-Patched SnpEff and SnpSift 5.2 JARs

**Files:**

- Add `docker/java/snpeff-pom.xml`
- Add `docker/java/snpsift-pom.xml`
- Add `docker/java/assert-runtime-dependencies.sh`
- Add `docker/java/src/main/java/org/apache/commons/lang/ArrayUtils.java`
- Add `docker/java/src/test/java/org/apache/commons/lang/ArrayUtilsTest.java`
- Modify `tests/unit/test_container_security_contract.py`

**Interfaces:**

- `ArrayUtils.toPrimitive(Integer[]) -> int[]` is binary-compatible with the only method referenced
  by `org.rcsb:mmtf-codec:1.0.10`.
- SnpEff artifact: `org.snpeff:SnpEff:5.2`, main class `org.snpeff.SnpEff`.
- SnpSift artifact: `org.snpsift:SnpSift:5.2`, main class `org.snpsift.SnpSift`.
- `assert-runtime-dependencies.sh <project-directory>` exits 0 only when every required fixed
  version is present and `commons-lang:commons-lang` is absent.

- [ ] Extend the static test module with failing XML assertions. Parse with
  `xml.etree.ElementTree`, using namespace `{"m": "http://maven.apache.org/POM/4.0.0"}`. Assert:

  - SnpEff and SnpSift are version `5.2` and have their expected assembly main classes;
  - all four JUnit/JUnit Platform dependencies have scope `test`;
  - SnpSift depends on locally installed `org.snpeff:SnpEff:5.2`, not 5.1;
  - Jackson Core/Databind are `2.22.1` and Annotations is `2.22`;
  - Gson is `2.14.0`, Commons IO `2.22.0`, Commons Compress `1.28.0`, and Log4j API,
    Core, and `log4j-slf4j-impl` are `2.25.4`;
  - the SnpEff `biojava-structure` dependency excludes `commons-lang:commons-lang`.

Use a coordinate helper so missing and duplicate dependencies fail clearly. Assert the version
property values separately, then assert the managed coordinates refer to the corresponding
property names:

```python
def _dependency_versions(path: str) -> dict[tuple[str, str], list[str]]:
    root = ElementTree.parse(ROOT / path).getroot()
    versions: dict[tuple[str, str], list[str]] = {}
    for dependency in root.findall(".//m:dependency", XML_NS):
        group = dependency.findtext("m:groupId", namespaces=XML_NS)
        artifact = dependency.findtext("m:artifactId", namespaces=XML_NS)
        version = dependency.findtext("m:version", namespaces=XML_NS)
        if group and artifact and version:
            versions.setdefault((group, artifact), []).append(version)
    return versions
```

- [ ] Run the focused tests and confirm failure because the overlay files do not exist.

Expected: the two conda tests pass and the new Java-policy tests fail with `FileNotFoundError`.

- [ ] Add the compatibility implementation exactly as follows:

```java
package org.apache.commons.lang;

/** Minimal binary compatibility for mmtf-codec 1.0.10. */
public final class ArrayUtils {
    private ArrayUtils() {
    }

    public static int[] toPrimitive(Integer[] array) {
        if (array == null) {
            return null;
        }

        int[] result = new int[array.length];
        for (int index = 0; index < array.length; index++) {
            result[index] = array[index].intValue();
        }
        return result;
    }
}
```

Do not add overloads, `ClassUtils`, reflection, class-name parsing, or any other Commons Lang code.

- [ ] Add the JUnit 5 compatibility tests:

```java
package org.apache.commons.lang;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

final class ArrayUtilsTest {
    @Test
    void nullInputReturnsNull() {
        assertNull(ArrayUtils.toPrimitive(null));
    }

    @Test
    void emptyArrayRemainsEmpty() {
        assertArrayEquals(new int[0], ArrayUtils.toPrimitive(new Integer[0]));
    }

    @Test
    void ordinaryAndNegativeValuesConvertExactly() {
        assertArrayEquals(new int[] {7, 0, -4},
                ArrayUtils.toPrimitive(new Integer[] {7, 0, -4}));
    }

    @Test
    void nullElementRetainsLegacyNullPointerBehavior() {
        assertThrows(NullPointerException.class,
                () -> ArrayUtils.toPrimitive(new Integer[] {1, null, 3}));
    }
}
```

- [ ] Create full repository-owned POM overlays from the POMs at the pinned source commits. Preserve
  every upstream compile dependency and repository required by the 5.2 source, then make these
  exact security changes:

```xml
<properties>
  <jackson.version>2.22.1</jackson.version>
  <jackson.annotations.version>2.22</jackson.annotations.version>
  <gson.version>2.14.0</gson.version>
  <commons.io.version>2.22.0</commons.io.version>
  <commons.compress.version>1.28.0</commons.compress.version>
  <log4j.version>2.25.4</log4j.version>
</properties>
```

Add fixed runtime coordinates under `dependencyManagement`, give all JUnit dependencies
`<scope>test</scope>`, and add the following exclusion to SnpEff's
`org.biojava:biojava-structure:6.0.4` dependency:

```xml
<exclusions>
  <exclusion>
    <groupId>commons-lang</groupId>
    <artifactId>commons-lang</artifactId>
  </exclusion>
</exclusions>
```

Set the SnpSift SnpEff dependency to `5.2`. Retain Java release 11 and ISO-8859-1 source/resource
encoding. Configure `maven-surefire-plugin:3.5.4` and `maven-assembly-plugin:3.7.1`, binding the
assembly `single` goal to `package`, so `mvn verify` both runs tests and emits one executable
`*-jar-with-dependencies.jar`.

- [ ] Add `assert-runtime-dependencies.sh`. It must run `mvn -B dependency:list` with runtime scope,
  save the resolved list, reject `commons-lang:commons-lang`, reject any resolved version of the
  seven secured families other than the exact versions below, and require each exact coordinate:

```text
com.fasterxml.jackson.core:jackson-core:2.22.1
com.fasterxml.jackson.core:jackson-databind:2.22.1
com.fasterxml.jackson.core:jackson-annotations:2.22
com.google.code.gson:gson:2.14.0
commons-io:commons-io:2.22.0
org.apache.commons:commons-compress:1.28.0
org.apache.logging.log4j:log4j-api:2.25.4
org.apache.logging.log4j:log4j-core:2.25.4
org.apache.logging.log4j:log4j-slf4j-impl:2.25.4
```

Use anchored Maven-coordinate matching; do not search binary JAR bytes or rely on filenames alone.
SnpSift may receive these through its local SnpEff 5.2 dependency, so apply the same checks to both
projects.

- [ ] Run static policy tests.

Expected: all conda and Java overlay tests pass.

- [ ] Prototype the builder independently of the production Dockerfile using the pinned Maven
  image. The source inputs and their required SHA-256 values are:

```text
SnpEff  0c5e74f9b6ca6ed3db720177eb1f95b9d47d45f2
        8633ecad8dbf06af19d5002818b65dc9e7b78419b94866dd4b1daed9d1e9d2b2
SnpSift 20978614457f14ec7a0c70539d5a7a2b7e754f60
        74c37b85e74a390a27f122164fd06c5b56131e3fa1eca205636d3de4a8c94934
```

Run `mvn -B verify`, the dependency assertion, and `mvn -B install -DskipTests` for SnpEff before
building SnpSift. Do not skip the first `verify`; SnpSift must resolve the locally installed patched
SnpEff artifact.

Expected: four compatibility tests plus upstream tests pass; both fat JARs are created; dependency
assertions report no banned coordinate.

- [ ] If an upstream test requires a downloadable genome or another non-hermetic external dataset,
  identify that test explicitly and split Maven execution into hermetic upstream unit tests plus the
  committed container golden tests in Task 4. Do not silently use `-DskipTests` for the builder.

- [ ] If explicit commit permission exists, commit this task:

```bash
git add docker/java tests/unit/test_container_security_contract.py
git commit -m "build: patch snpeff 5.2 runtime dependencies"
```

Otherwise, leave the changes uncommitted.

---

## Task 3: Harden the Production Docker Image

**Files:**

- Modify `Dockerfile`
- Modify `tests/unit/test_container_security_contract.py`

**Interfaces:**

- Stage names: `java-build`, `conda-build`, and `runtime`.
- Final image still exposes the existing entrypoint, command, health check, labels, `/app/LICENSE`,
  `/app/scoring`, `/app/stats_configs`, and `/data` workdir.
- Final user is `$MAMBA_USER`, never root.

- [ ] Add failing Dockerfile policy tests asserting:

  - the Java builder is
    `maven:3.9.11-eclipse-temurin-11@sha256:c095e2421eaf3e5cb1573fd0474e68e17062866d454362349e17bbb75f44031e`;
  - both Micromamba stages are
    `mambaorg/micromamba:2.8.1-debian12-slim@sha256:c8198d53228ad7cfd7adcc0704e8837f9d1c9327fb363c6a7d3c5b1a51a4b561`;
  - both remote archives use the BuildKit `ADD --checksum` form with the exact SHA-256 values from
    Task 2;
  - `apt-get update`, `dist-upgrade`, apt-index removal, and `/opt/conda/pkgs` removal are present;
  - both builder JARs replace the exact conda share paths;
  - JDK-only compiler/tooling removal and a negative `javac` assertion are present;
  - the last `USER` instruction is `USER $MAMBA_USER`;
  - the entrypoint, command, workdir, and health-check command are unchanged.

- [ ] Run the test and confirm the Dockerfile assertions fail against the two-stage 2.0 image.

- [ ] Replace the Dockerfile structure with these stages and ordering:

```dockerfile
FROM maven:3.9.11-eclipse-temurin-11@sha256:c095e2421eaf3e5cb1573fd0474e68e17062866d454362349e17bbb75f44031e AS java-build

ADD --checksum=sha256:8633ecad8dbf06af19d5002818b65dc9e7b78419b94866dd4b1daed9d1e9d2b2 \
    https://github.com/pcingola/SnpEff/archive/0c5e74f9b6ca6ed3db720177eb1f95b9d47d45f2.tar.gz /tmp/snpeff.tar.gz
ADD --checksum=sha256:74c37b85e74a390a27f122164fd06c5b56131e3fa1eca205636d3de4a8c94934 \
    https://github.com/pcingola/SnpSift/archive/20978614457f14ec7a0c70539d5a7a2b7e754f60.tar.gz /tmp/snpsift.tar.gz
```

Extract with `tar --strip-components=1` into separate directories. Copy the overlay POM,
compatibility source/test, and dependency assertion into SnpEff; run `mvn -B verify`, the assertion,
then `mvn -B install -DskipTests`. Copy the SnpSift POM, run `mvn -B verify`, and run the assertion.
Finally, copy each fat JAR to a stable `/out/snpEff.jar` or `/out/SnpSift.jar` path and verify there
is exactly one match before copying.

- [ ] Rename the current conda stage to `conda-build`, pin its base by digest, and keep the existing
  install order. After the Python package install:

```dockerfile
COPY --from=java-build --chown=$MAMBA_USER:$MAMBA_USER \
    /out/snpEff.jar /opt/conda/share/snpeff-5.2-3/snpEff.jar
COPY --from=java-build --chown=$MAMBA_USER:$MAMBA_USER \
    /out/SnpSift.jar /opt/conda/share/snpsift-5.2-0/SnpSift.jar
```

Use Python's standard-library `zipfile` module in a BuildKit heredoc to assert exactly these
manifest lines:

```text
Main-Class: org.snpeff.SnpEff
Main-Class: org.snpsift.SnpSift
```

Then remove `/opt/conda/pkgs` and fail if it still exists. Also fail if the two expected runtime
JAR paths are missing or if another `snpEff.jar`/`SnpSift.jar` exists under `/opt/conda`.

After manifest validation, remove JDK-only tooling inherited through Bioconda's `openjdk` package:
`javac`, `javadoc`, `javap`, `jar`, `jarsigner`, `jconsole`, `jdeps`, `jlink`, `jmod`, and `jshell`
from both `/opt/conda/bin` and the resolved JVM `bin` directory, plus the JVM `include`, `jmods`,
and `src.zip` paths. Fail the build if `command -v javac` still resolves. Keep `java`, `keytool`, and
the runtime module image needed by the SnpEff/SnpSift wrappers.

- [ ] Pin the runtime base by the same Micromamba digest. Before copying the conda environment:

```dockerfile
USER root
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y --no-install-recommends dist-upgrade \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
```

Copy `/opt/conda` from `conda-build`, preserve current application assets and metadata, and restore
`USER $MAMBA_USER` before `WORKDIR`, `HEALTHCHECK`, `ENTRYPOINT`, and `CMD`.

- [ ] Run static policy tests.

Expected: all tests in `test_container_security_contract.py` pass except workflow tests not yet
added.

- [ ] Build the production image, not a test-only target:

```bash
docker build --pull --tag variantcentrifuge:security .
```

Expected: checksummed source downloads succeed, Maven tests and dependency assertions pass, conda
installs exact builds, JAR manifest/cache assertions pass, and Docker emits
`naming to docker.io/library/variantcentrifuge:security`.

- [ ] Inspect the final image directly:

```bash
docker run --rm --entrypoint /usr/bin/id variantcentrifuge:security -u
docker run --rm --entrypoint /bin/sh variantcentrifuge:security -c \
  'test ! -e /opt/conda/pkgs && test "$(id -u)" -ne 0'
```

Expected: the first command prints a nonzero UID and the second exits 0.

- [ ] If explicit commit permission exists, commit this task:

```bash
git add Dockerfile tests/unit/test_container_security_contract.py
git commit -m "build: harden the production container"
```

Otherwise, leave the changes uncommitted.

---

## Task 4: Preserve the Container's Behavioral Contract

**Files:**

- Add `scripts/test_container_image.sh`
- Add `tests/fixtures/container/snpeff/snpEff.config`
- Add `tests/fixtures/container/snpeff/data/testGenome/genes.gff`
- Add `tests/fixtures/container/snpeff/data/testGenome/sequences.fa`
- Add `tests/fixtures/container/snpeff/input.vcf`
- Add `tests/fixtures/container/snpeff/expected.vcf`
- Add `tests/fixtures/container/snpsift/input.vcf`
- Add `tests/fixtures/container/snpsift/expected.vcf`
- Add `tests/fixtures/container/snpsift/assert_fail_closed.py`

**Interfaces:**

- `scripts/test_container_image.sh IMAGE_REF` exits 0 only if every behavioral and filesystem
  assertion passes.
- The script accepts exactly one positional image reference, creates a task-specific temporary
  directory with `mktemp -d`, and always removes it with `trap`.
- Golden comparisons remove only the path-bearing `##SnpEffCmd` and `##SnpSiftCmd` headers.
  Version/build headers, variant rows, and all other headers remain exact.

- [ ] Add the SnpEff fixture with this exact 60-base model:

```text
##gff-version 3
1	test	gene	1	60	.	+	.	ID=gene1;Name=GENE1
1	test	mRNA	1	60	.	+	.	ID=transcript1;Parent=gene1
1	test	exon	1	60	.	+	.	Parent=transcript1
1	test	CDS	1	60	.	+	0	ID=cds1;Parent=transcript1
```

The FASTA sequence is:

```text
>1
ATGGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAATAA
```

The input variant is `1:4 G>T`. The normalized golden row must be:

```text
1	4	.	G	T	.	PASS	ANN=T|stop_gained|HIGH|GENE1|gene1|transcript|transcript1|protein_coding|1/1|c.4G>T|p.Glu2*|4/60|4/60|2/19||
```

Build the local database with:

```bash
snpEff build -c snpEff.config -dataDir ./data -gff3 \
  -noCheckCds -noCheckProtein testGenome
```

Annotate with `snpEff -c snpEff.config -dataDir ./data -noStats testGenome input.vcf`.

- [ ] Add the SnpSift input with QUAL 50 and QUAL 10 records. Filter with
  `( QUAL >= 20 )`; the normalized golden output contains only the position 4 row, with QUAL
  serialized by SnpSift as `50.0`.

- [ ] Add `assert_fail_closed.py` using the production Python API:

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from variantcentrifuge.filters import apply_snpsift_filter

fixture = Path(__file__).with_name("input.vcf")
with TemporaryDirectory(prefix="vc-snpsift-") as directory:
    output = Path(directory) / "filtered.vcf.gz"
    try:
        apply_snpsift_filter(
            str(fixture),
            "(( QUAL >= 20 )",
            {"threads": 1},
            str(output),
        )
    except RuntimeError as error:
        if "SnpSift filter reported parser diagnostics" not in str(error):
            raise
    else:
        raise SystemExit("malformed SnpSift filter did not fail closed")

print("SnpSift parser diagnostics failed closed")
```

This expression is intentional: SnpSift 5.2 reports `missing ')' at '<EOF>'` but exits 0, so the
test proves VariantCentrifuge's stderr guard rather than merely checking a Java exit code.

- [ ] Implement the smoke script in fail-fast Bash (`set -euo pipefail`). It must check, in order:

  1. `variantcentrifuge --version` and `--help`;
  2. the literal health-check executable `/opt/conda/bin/variantcentrifuge --version`;
  3. `id -u` is nonzero;
  4. `bcftools --version` starts with 1.21 and `bedtools --version` contains 2.31.1;
  5. `snpEff -version` and `SnpSift -version` contain 5.2;
  6. Python `zipfile` reads the expected JAR manifest main classes;
  7. the valid SnpSift output matches its normalized golden VCF;
  8. `assert_fail_closed.py` prints the expected marker;
  9. the SnpEff output matches its normalized golden VCF;
  10. `/opt/conda/pkgs` is absent, the default UID is nonzero, exactly two expected tool JARs
      exist, `javac` is unavailable, and no Maven/compiler/source/cache path exists.

All arbitrary commands must set an explicit `docker run --entrypoint` value; do not accidentally
feed them to the production VariantCentrifuge entrypoint. Copy fixtures to the temporary directory and mount
that copy read/write so SnpEff database generation never changes committed files.

- [ ] First run the golden annotation/filter portions against the stock 5.2 images to prove the
  expected outputs describe upstream behavior:

```bash
docker pull quay.io/biocontainers/snpeff:5.2--hdfd78af_3
docker pull quay.io/biocontainers/snpsift:5.2--hdfd78af_0
```

Expected source image digests:

```text
snpeff  sha256:93ce8d13ec526514fdab68258c1b8ca9903d6bc6366e5560753b7e42fd789aca
snpsift sha256:599e98a362e32b122ebb519ff38ebd7e27533d0dd34d65b3e64e3a089c8f3770
```

- [ ] Run the complete contract against the patched production image:

```bash
bash scripts/test_container_image.sh variantcentrifuge:security
```

Expected final line: `container contract passed: variantcentrifuge:security`.

- [ ] Run the existing real-tool regression inside the image's installed environment if practical,
  or prove its two exact cases through the helper and valid golden test above:

```bash
pytest tests/integration/test_snpsift_filter_guards.py -q
```

Expected: `2 passed` in an environment with SnpSift, bgzip, and bcftools.

- [ ] If explicit commit permission exists, commit this task:

```bash
git add scripts/test_container_image.sh tests/fixtures/container
git commit -m "test: add production container behavior contract"
```

Otherwise, leave the changes uncommitted.

---

## Task 5: Make Trivy Evidence Deterministic and Auditable

**Files:**

- Add `scripts/summarize_trivy.py`
- Add `tests/unit/test_summarize_trivy.py`

**Interfaces:**

- `summarize_trivy.py COMPLETE_JSON ACTIONABLE_JSON` prints Markdown to stdout.
- A vulnerability is fixed/actionable only when `FixedVersion` is a nonempty string.
- Severity order is `UNKNOWN`, `LOW`, `MEDIUM`, `HIGH`, `CRITICAL`.
- Valid Trivy documents always exit 0 after printing counts; malformed input exits 2. The separate
  Trivy invocation with `--exit-code 1` remains the only authoritative workflow gate.

- [ ] Write failing tests using minimal Trivy documents with multiple `Results`, a result without a
  `Vulnerabilities` key, one empty `FixedVersion`, and one nonempty `FixedVersion`. Assert the
  summary contains:

```text
Complete findings: 3
Vendor-fixed findings: 1
Vendor-unfixed findings: 2
Actionable findings: 0
```

Add a second test asserting an actionable finding is counted under its severity while `main()`
still returns 0 for valid input. Test functions directly; use `capsys` only for CLI output.

- [ ] Run the focused tests and confirm they fail because the module is absent.

- [ ] Implement these pure interfaces before the CLI wrapper: `vulnerabilities(document)` returns
  the flattened vulnerability dictionaries; `has_fix(vulnerability)` returns whether
  `FixedVersion` is a nonempty string; `severity_counts(items)` returns all five severity buckets;
  `render_markdown(complete, actionable)` returns the report; and `main(argv=None)` reads both JSON
  paths, prints the report, and returns 0 for valid input. Give each interface the concrete typed
  arguments and returns described in this task's **Interfaces** section.

Handle missing/unknown severity as `UNKNOWN`; reject malformed top-level JSON with a clear error.
Do not infer fixed status from CVE text, installed version ordering, or severity.

- [ ] Run:

```bash
pytest tests/unit/test_summarize_trivy.py -q
ruff check scripts/summarize_trivy.py tests/unit/test_summarize_trivy.py
```

Expected: all focused tests pass and Ruff reports no errors.

- [ ] If explicit commit permission exists, commit this task:

```bash
git add scripts/summarize_trivy.py tests/unit/test_summarize_trivy.py
git commit -m "test: summarize complete and actionable trivy reports"
```

Otherwise, leave the changes uncommitted.

---

## Task 6: Gate Before Publishing in GitHub Actions

**Files:**

- Modify `.github/workflows/docker.yml`
- Modify `tests/unit/test_container_security_contract.py`

**Interfaces:**

- Job ID `scan` owns build, smoke test, both reports, gate, and conditional publication.
- Outputs: `image-digest` from the publish step and `image-tags` from Docker metadata.
- `sign` uses `needs: scan` and runs only outside pull requests.
- Local image reference: `variantcentrifuge:ci-${{ github.sha }}`.
- SARIF analysis category remains `.github/workflows/docker.yml:scan` because both workflow path and
  job ID are preserved.

- [ ] Add failing workflow-policy tests that read the YAML as text and assert:

  - there is no top-level `build` job and `scan` has no `needs: build` or pull-request skip;
  - PR path filters include `docker/**`, `scripts/test_container_image.sh`, and
    `tests/fixtures/container/**`;
  - build uses `load: true` and `push: false`;
  - smoke occurs before both Trivy reporting steps;
  - `trivy-complete.json` is generated without `ignore-unfixed` and uploaded with 90-day retention;
  - `trivy-actionable.sarif` uses `ignore-unfixed: true` and is uploaded;
  - the final Trivy invocation uses `ignore-unfixed: true`, all five severities, and
    `exit-code: "1"`;
  - every Trivy call pins `trivy-version: v0.70.0` and `scanners: vuln`;
  - login/push text occurs after the final gate;
  - `sign` contains `needs: scan`, not `needs: build`.

Use substring slices between unique step names so the complete-report test does not mistake the
actionable report's `ignore-unfixed` setting for its own.

- [ ] Run the test and confirm it fails against the current report-only, post-publish workflow.

- [ ] Replace the current build/scan split with one `scan` job. Keep metadata generation before the
  build, add `LOCAL_IMAGE: variantcentrifuge:ci-${{ github.sha }}` to workflow `env`, and build a
  single-platform local image with:

```yaml
- name: Build production image
  uses: docker/build-push-action@v6
  with:
    context: .
    load: true
    push: false
    tags: |
      ${{ env.LOCAL_IMAGE }}
      ${{ steps.meta.outputs.tags }}
    labels: ${{ steps.meta.outputs.labels }}
    cache-from: type=gha
    cache-to: type=gha,mode=max
```

- [ ] Immediately run:

```yaml
- name: Test production image
  run: bash scripts/test_container_image.sh "${LOCAL_IMAGE}"
```

- [ ] Generate and upload the complete report before any failing gate:

```yaml
- name: Generate complete vulnerability report
  uses: aquasecurity/trivy-action@v0.36.0
  with:
    image-ref: ${{ env.LOCAL_IMAGE }}
    format: json
    output: trivy-complete.json
    scanners: vuln
    vuln-type: os,library
    severity: UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL
    exit-code: "0"
    trivy-version: v0.70.0

- name: Upload complete vulnerability report
  if: always()
  uses: actions/upload-artifact@v7.0.1
  with:
    name: trivy-complete-${{ github.sha }}
    path: trivy-complete.json
    retention-days: 90
```

Do not add `ignore-unfixed` to this report.

- [ ] Generate the actionable SARIF with the same scanner/type/severity/version values plus
  `ignore-unfixed: true`, `format: sarif`, `output: trivy-actionable.sarif`, and `exit-code: "0"`.
  Upload it with `github/codeql-action/upload-sarif@v4.37.1`, `if: always()`, and
  `sarif_file: trivy-actionable.sarif`.

- [ ] Generate an actionable JSON report as well (`trivy-actionable.json`) with
  `ignore-unfixed: true` and `exit-code: "0"`, then append the summary to `$GITHUB_STEP_SUMMARY`:

```yaml
- name: Summarize vulnerability reports
  run: python scripts/summarize_trivy.py trivy-complete.json trivy-actionable.json \
    | tee -a "${GITHUB_STEP_SUMMARY}"
```

This reporting step stays green for any valid documents; the explicit final Trivy gate below is the
sole authoritative scanner exit status.

- [ ] Add the final gate after evidence generation/upload and before login:

```yaml
- name: Enforce zero vendor-fixed vulnerabilities
  uses: aquasecurity/trivy-action@v0.36.0
  with:
    image-ref: ${{ env.LOCAL_IMAGE }}
    format: table
    scanners: vuln
    vuln-type: os,library
    severity: UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL
    ignore-unfixed: true
    exit-code: "1"
    trivy-version: v0.70.0
```

- [ ] Only after the gate, and only when `github.event_name != 'pull_request'`, retain the existing
  three-attempt GHCR login. Push every newline-delimited metadata tag from
  `${{ steps.meta.outputs.tags }}` without rebuilding. Resolve the first pushed tag's registry
  digest from `docker image inspect`, validate it matches `^sha256:[0-9a-f]{64}$`, and write it to
  `$GITHUB_OUTPUT` from a step with ID `publish`.

- [ ] Configure job outputs:

```yaml
outputs:
  image-digest: ${{ steps.publish.outputs.digest }}
  image-tags: ${{ steps.meta.outputs.tags }}
```

Change the signing job to `needs: scan`; keep the non-PR condition, GHCR retry, keyless Cosign
installer, and `cosign sign --yes "${REGISTRY}/${IMAGE_NAME}@${DIGEST}"` behavior.

- [ ] Run static policy tests and validate workflow syntax:

```bash
pytest tests/unit/test_container_security_contract.py -q
docker run --rm -v "$PWD:/repo" -w /repo rhysd/actionlint:1.7.7
```

Expected: policy tests pass and actionlint exits 0 without diagnostics.

- [ ] Re-read the ordered workflow and prove these failure semantics manually:

  - a smoke or gate failure happens before login/push;
  - no PR step can publish or sign;
  - a failed `scan` prevents `sign` through `needs`;
  - complete JSON and actionable SARIF steps precede the final gate;
  - no second Docker build can create bytes different from the tested image.

- [ ] If explicit commit permission exists, commit this task:

```bash
git add .github/workflows/docker.yml tests/unit/test_container_security_contract.py
git commit -m "ci: gate container publication on trivy"
```

Otherwise, leave the changes uncommitted.

---

## Task 7: Document the Security and Risk Contract

**Files:**

- Modify `docs/source/installation.md`
- Modify `docs/source/changelog.md`

- [ ] Add a container-security subsection to installation documentation explaining:

  - the image retains SnpEff/SnpSift 5.2 wrappers/configs while rebuilding their JAR payloads from
    pinned 5.2 source commits with fixed runtime dependencies;
  - the complete JSON artifact includes fixed and unfixed findings and is retained 90 days;
  - code scanning and the merge gate contain vendor-fixed findings at every severity;
  - `--ignore-unfixed` partitions actionable reporting and is not a claim that unfixed findings do
    not exist;
  - no `.trivyignore` or alert dismissal is used for this remediation.

- [ ] Add an unreleased changelog item covering the pinned base refresh, patched Java JARs,
  behavioral smoke tests, pre-publication Trivy gate, complete audit artifact, and digest signing.

- [ ] Build documentation or at minimum run the repository's documented Sphinx link/format check.
  If the repo has no lightweight docs target, run:

```bash
ruff check tests/unit/test_container_security_contract.py tests/unit/test_summarize_trivy.py \
  scripts/summarize_trivy.py
```

Expected: no errors; visually verify Markdown/reStructuredText renders without malformed nesting.

- [ ] If explicit commit permission exists, commit this task:

```bash
git add docs/source/installation.md docs/source/changelog.md
git commit -m "docs: explain container vulnerability policy"
```

Otherwise, leave the changes uncommitted.

---

## Task 8: Run the Full Security Acceptance Gate

**Files:**

- Verify all changed files; only update a file when a failing check identifies a real issue.

- [ ] Run repository regression checks:

```bash
ruff check .
mypy variantcentrifuge
pytest -q
```

Expected: Ruff and mypy exit 0; the full pytest suite passes. Record any environment-based skip
count, but do not waive a failure as unrelated without reproducing it on baseline main.

- [ ] Rebuild from scratch so cached vulnerable layers cannot mask the result:

```bash
docker build --pull --no-cache --tag variantcentrifuge:security .
bash scripts/test_container_image.sh variantcentrifuge:security
```

Expected: fresh build and complete behavioral contract pass.

- [ ] Run the same pinned Trivy 0.70.0 policy locally. Use a task-specific temporary directory, not
  a repository output path:

```bash
security_tmp=$(mktemp -d /tmp/variantcentrifuge-trivy.XXXXXX)
docker run --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "${security_tmp}:/reports" \
  aquasec/trivy:0.70.0 image \
  --scanners vuln \
  --vuln-type os,library \
  --severity UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL \
  --format json \
  --output /reports/complete.json \
  variantcentrifuge:security
docker run --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "${security_tmp}:/reports" \
  aquasec/trivy:0.70.0 image \
  --scanners vuln \
  --vuln-type os,library \
  --severity UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL \
  --ignore-unfixed \
  --format json \
  --output /reports/actionable.json \
  --exit-code 1 \
  variantcentrifuge:security
python scripts/summarize_trivy.py \
  "${security_tmp}/complete.json" "${security_tmp}/actionable.json"
```

Expected: the actionable scan exits 0 because it contains zero vulnerabilities. The complete scan
may contain vendor-unfixed findings; record its total and severity distribution without calling
them remediated.

- [ ] If the actionable scan finds any implementation-time vendor-fixed vulnerability, do not add
  an exception. Identify its target/package/fixed version, add the minimal version/base/JAR change,
  add or tighten a regression assertion, and repeat Tasks 2-8 as applicable until actionable count
  is zero.

- [ ] Inspect the image history and filesystem:

```bash
docker history --no-trunc variantcentrifuge:security
docker run --rm --entrypoint /bin/sh variantcentrifuge:security -c \
  'test ! -e /opt/conda/pkgs && test ! -e /root/.m2 && test "$(id -u)" -ne 0'
```

Expected: no Maven/cache source layer is copied into runtime and the filesystem assertion exits 0.

- [ ] Re-query GitHub code scanning read-only immediately before opening the PR and save these facts
  in the PR description:

  - baseline commit and scan timestamp;
  - 203 baseline alerts / 80 fixed / 123 unfixed;
  - final complete count and severity distribution;
  - final actionable count `0`;
  - local image ID and, once main publishes it, registry digest;
  - workflow URL once CI runs;
  - explicit statement that no alert state was mutated.

- [ ] Review `git diff --check`, the full diff, and status:

```bash
git diff --check
git diff --stat
git status --short
```

Expected: no whitespace errors; only files in this plan plus the design/plan documents are changed.

- [ ] Request code review using `superpowers:requesting-code-review`, address technically valid
  findings with `superpowers:receiving-code-review`, and rerun the affected acceptance commands.

- [ ] If explicit commit permission exists and earlier task commits were intentionally deferred,
  create one reviewed implementation commit:

```bash
git add Dockerfile conda/environment-docker.yml docker/java \
  scripts/test_container_image.sh scripts/summarize_trivy.py \
  tests/fixtures/container tests/unit/test_container_security_contract.py \
  tests/unit/test_summarize_trivy.py .github/workflows/docker.yml \
  docs/source/installation.md docs/source/changelog.md \
  docs/superpowers/specs/2026-07-19-code-scanning-remediation-design.md \
  docs/superpowers/plans/2026-07-19-code-scanning-remediation.md
git commit -m "fix: remediate container code-scanning alerts"
```

Otherwise, leave all verified changes uncommitted for user review.

---

## Task 9: PR and Post-Merge Alert Closure Verification

This task requires explicit user authorization to publish a branch/open a PR and is not implied by
plan execution.

- [ ] After authorization, use `github:yeet` to push the reviewed branch and open one draft PR.
  Include the evidence from Task 8 and link the complete Trivy artifact/workflow run.

- [ ] Verify the PR `scan` job builds, runs smoke tests, uploads complete JSON and actionable SARIF,
  passes the fixed-only gate, and does not log in, push, or sign.

- [ ] After the PR merges, verify the main-branch `scan` job pushes the already tested local image,
  exports the immutable registry digest, and the `sign` job signs that digest.

- [ ] Read GitHub code-scanning alerts again without mutation. Confirm alerts 1 through 203 are
  closed/superseded by the new `.github/workflows/docker.yml:scan` analysis. If any remains open,
  compare its rule ID, package, location, commit, and analysis category against the new SARIF before
  deciding whether it represents a real remaining finding or stale analysis lifecycle state.

- [ ] Record the registry digest, main workflow URL, complete/unfixed counts, actionable count zero,
  and final dashboard closure status in the PR or follow-up release evidence.
