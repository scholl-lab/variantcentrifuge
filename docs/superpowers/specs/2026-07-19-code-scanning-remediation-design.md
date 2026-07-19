# Code-Scanning Remediation Design

## Status

- Date: 2026-07-19
- Repository: `scholl-lab/variantcentrifuge`
- Baseline commit: `335d971dfa813b8a3c3d126cf586ecce70c90bd5`
- Baseline source: GitHub code-scanning REST API, `state=open`, read on 2026-07-19
- Delivery boundary: one pull request

## Executive Summary

The GitHub security dashboard has 203 open Trivy alerts for the published container:

| Surface | Alerts | Vendor-fixed | Vendor-unfixed |
|---|---:|---:|---:|
| Debian runtime packages | 179 | 57 | 122 |
| SnpEff 5.2 fat JAR | 19 | 18 | 1 |
| SnpSift 5.2 fat JAR | 5 | 5 | 0 |
| **Total** | **203** | **80** | **123** |

The pull request will remediate every vendor-fixed finding, preserve complete visibility of
vendor-unfixed findings, and make the same standard a pull-request gate. It will not describe
vendor-unfixed findings as fixed or delete them through the GitHub API.

The remediation has three coordinated parts:

1. Rebuild SnpEff and SnpSift from pinned 5.2-era source commits with non-vulnerable runtime
   dependencies while preserving their CLI behavior.
2. Move the runtime to a current, digest-pinned Micromamba Debian image, apply available Debian
   security updates, and remove package/build caches from the final filesystem;
3. split CI evidence into a complete audit report and an actionable SARIF report, then fail every
   pull request that contains any vendor-fixed vulnerability at any severity.

## Evidence and Root Causes

### Dashboard snapshot

All 203 alerts were produced by Trivy 0.70.0 in
`.github/workflows/docker.yml:scan` against `refs/heads/main`. Severity distribution is:

| Severity | Count |
|---|---:|
| Critical | 3 |
| High | 40 |
| Medium | 79 |
| Low | 74 |
| Unknown | 7 |

There are 116 unique vulnerability identifiers across 47 package coordinates. Repeated CVEs on
different Debian binary packages are separate alerts and remain separately accounted for in the
appendix.

### Stale OS base

Both Docker stages use `mambaorg/micromamba:2.0-debian12-slim`. That tag resolves to a base built
in March 2025, while the current `2.8.1-debian12-slim` image was published in June 2026 and reports
Debian 12.14. The old runtime contains security-fixable versions of glibc, GnuPG, PAM, systemd,
libcap, dpkg, sed, and related packages.

Changing the base tag alone is insufficient. A same-version Trivy probe of the current base on
2026-07-19 still found 171 OS findings, including 14 with vendor fixes. The runtime must therefore
apply repository security upgrades during the image build and prove the result with a post-build
scan.

### Vulnerable Java fat JARs

The conda packages install fat JARs containing vulnerable transitive libraries:

- Jackson Core and Databind 2.10.5;
- Gson 2.8.7;
- Commons IO 2.11.0;
- Commons Lang 2.4;
- Commons Compress 1.19;
- Log4j Core 2.17.1.

The latest stock Bioconda SnpEff/SnpSift 5.4.0c images were probed with Trivy 0.70.0. They retain
the same vulnerable coordinates, so a version-only upgrade neither solves the alerts nor preserves
5.2 behavior.

The 5.2 dependency tree shows the origins:

| Vulnerable coordinate | Dependency path |
|---|---|
| Jackson Core/Databind | `biojava-structure -> mmtf-serialization -> jackson-dataformat-msgpack` |
| Gson | `biojava-structure -> ciftools-java-jdk8` |
| Commons Lang | `biojava-structure -> mmtf-codec` |
| Commons Compress | `htsjdk` |
| Log4j Core | `biojava-core` |
| Commons IO | direct SnpEff dependency |

`mmtf-codec` uses legacy Commons Lang only for
`ArrayUtils.toPrimitive(Integer[])`. The affected Commons Lang 2.x line has no fixed release;
Commons Lang 3.18.0 fixes the equivalent CVE but uses a different Java package namespace. A small,
source-owned compatibility implementation of that single conversion can remove the legacy library
without dropping MMTF functionality or retaining vulnerable `ClassUtils` code.

### CI cannot validate a security pull request

The Docker build loads the image for pull requests, but the `scan` job is explicitly skipped for
pull requests. On `main`, Trivy is report-only (`exit-code: "0"`). As a result:

- the pull request introducing a container vulnerability is not scanned;
- a scan with actionable findings remains green;
- vendor-fixed and vendor-unfixed findings are mixed in one dashboard stream;
- no complete machine-readable report is retained as a workflow artifact.

## Goals

1. Remediate all 80 currently vendor-fixed alerts and any additional vendor-fixed alerts found by
   the implementation-time Trivy database.
2. Preserve SnpEff and SnpSift 5.2 command-line behavior used by VariantCentrifuge.
3. Keep all vendor-unfixed findings visible in a complete JSON audit artifact.
4. Upload only actionable, vendor-fixed findings to GitHub code-scanning SARIF.
5. Fail pull requests and protected-branch builds when any actionable finding remains at
   `UNKNOWN`, `LOW`, `MEDIUM`, `HIGH`, or `CRITICAL` severity.
6. Remove vulnerable stock JARs, conda package caches, build tools, and package-manager indexes
   from the final runtime filesystem.
7. Keep the runtime non-root and retain the current entrypoint, health check, CLI, external-tool,
   scoring-model, and stats-config behavior.
8. Deliver the image changes, tests, workflow changes, documentation, and alert-accounting evidence
   in one pull request.

## Non-Goals

- Do not claim the final image has no known vulnerabilities. Vendor-unfixed findings remain known
  and auditable.
- Do not dismiss, close, or mutate alerts through the GitHub API. GitHub should close obsolete
  alerts from the next SARIF analysis.
- Do not hide vendor-fixed findings with `.trivyignore`, VEX, severity filters, or path filters.
- Do not upgrade SnpEff/SnpSift behavior to 5.4 in this security pull request.
- Do not remove SnpEff or SnpSift from the container.
- Do not rewrite VariantCentrifuge's Python pipeline.
- Do not promise a stable raw finding count: vulnerability databases and Debian advisories evolve.

## Definitions and Security Policy

### Actionable finding

An actionable finding is a Trivy vulnerability for which Trivy's selected vendor data provides a
non-empty fixed version. The actionable scan uses `--ignore-unfixed` and includes every severity.

### Vendor-unfixed finding

A vendor-unfixed finding has no fixed version in Trivy's selected vendor data. It is excluded from
SARIF so GitHub code scanning represents work this repository can act on, but it remains in the
complete JSON audit artifact.

`--ignore-unfixed` is a reporting partition, not a risk dismissal. The JSON artifact is mandatory,
uploaded even if a later gate fails, and retained for 90 days. A future scan automatically promotes
a finding into the actionable gate when vendor data gains a fixed version.

### Success criterion

The pull request succeeds only when:

```text
complete_scan = all severities + fixed + unfixed
actionable_scan = all severities + fixed only
count(actionable_scan.vulnerabilities) = 0
```

The complete scan may contain vendor-unfixed findings. The PR description must report their count
without calling them remediated.

## Chosen Design

### 1. Reproducible Java tool builder

Add a dedicated Maven build stage based on:

```text
maven:3.9.11-eclipse-temurin-11@sha256:c095e2421eaf3e5cb1573fd0474e68e17062866d454362349e17bbb75f44031e
```

Fetch and verify these source archives:

| Tool | Source commit | Archive SHA-256 |
|---|---|---|
| SnpEff | `0c5e74f9b6ca6ed3db720177eb1f95b9d47d45f2` (`v5.2`) | `8633ecad8dbf06af19d5002818b65dc9e7b78419b94866dd4b1daed9d1e9d2b2` |
| SnpSift | `20978614457f14ec7a0c70539d5a7a2b7e754f60` (5.2 release-era commit) | `74c37b85e74a390a27f122164fd06c5b56131e3fa1eca205636d3de4a8c94934` |

Repository-owned Maven overlays must:

- keep SnpEff/SnpSift source at those commits;
- mark JUnit and JUnit Platform dependencies as `test` so they are absent from runtime JARs;
- constrain Jackson Core and Databind to `2.22.1` and Jackson Annotations to its aligned `2.22`
  release (Annotations does not publish a `2.22.1` artifact);
- constrain Gson to `2.14.0`;
- constrain Commons IO to `2.22.0`;
- constrain Commons Compress to `1.28.0`;
- constrain Log4j API/Core/SLF4J implementation consistently to `2.25.4`;
- exclude `commons-lang:commons-lang` from `mmtf-codec`;
- add only the required behavior-compatible `Integer[]` to `int[]` conversion as
  `org.apache.commons.lang.ArrayUtils.toPrimitive(Integer[])`, satisfying the one legacy binary
  reference without copying the vulnerable Commons Lang implementation;
- build separate executable fat JARs with main classes `org.snpeff.SnpEff` and
  `org.snpsift.SnpSift`;
- run the upstream Maven tests before publishing either JAR to the next stage.

The compatibility class must not implement `ClassUtils`, accept class names, recurse, or copy
other Commons Lang behavior. Its tests preserve the original method contract: null input returns
null, a null element raises `NullPointerException`, and empty, negative, and ordinary arrays convert
exactly.

### 2. Hardened conda environment and runtime

Use the current Micromamba image pinned by digest in both conda-build and runtime stages:

```text
mambaorg/micromamba:2.8.1-debian12-slim@sha256:c8198d53228ad7cfd7adcc0704e8837f9d1c9327fb363c6a7d3c5b1a51a4b561
```

Keep exact conda package builds for wrapper/config compatibility:

```yaml
snpeff=5.2=hdfd78af_3
snpsift=5.2=hdfd78af_0
pip>=26.1.2
```

The explicit `pip` floor covers implementation-time Python-package findings already fixed by the
current conda-forge package. The final Trivy gate, rather than this single pin, remains authoritative
for any later package advisory.

After conda installation:

1. replace both stock runtime JARs with the verified builder outputs;
2. assert the JAR manifests contain the expected main classes;
3. remove the conda package cache rather than copying vulnerable cached stock JARs;
4. assert no stock JAR or file under `/opt/conda/pkgs` remains;
5. install VariantCentrifuge as today, without resolving a second dependency graph.

In the runtime stage, temporarily switch to root, run a non-interactive Debian package-index refresh
and distribution security upgrade, remove apt metadata, copy the prepared conda environment, then
switch back to `$MAMBA_USER`. No compiler, Maven cache, Java source, apt index, or conda package
cache may be present in the final filesystem.

The runtime must retain:

- `ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "variantcentrifuge"]`;
- `CMD ["--help"]`;
- the existing health check;
- working directories, labels, LICENSE, scoring models, and stats configs;
- non-root UID/GID behavior inherited from the Micromamba image.

### 3. Container behavioral contract

Add a deterministic container smoke-test script and small fixtures. It runs against the exact image
that Trivy scans and proves:

1. `variantcentrifuge --version` and `--help` succeed;
2. the health-check command succeeds;
3. the default user is not UID 0;
4. bcftools 1.21 and bedtools 2.31.1 remain available;
5. SnpEff and SnpSift report 5.2;
6. both JAR manifests have the correct main class;
7. SnpSift valid filtering produces the committed golden VCF;
8. SnpSift invalid syntax still triggers VariantCentrifuge's fail-closed diagnostic behavior;
9. SnpEff annotates a tiny local genome/VCF fixture to the committed canonical output;
10. no path under `/opt/conda/pkgs` and no vulnerable stock JAR remains.

Golden VCF comparison normalizes only path-bearing SnpEff/SnpSift command header lines. Version and
build headers, variant rows, ANN contents, filter results, exit codes, stdout, and stderr behavior
remain exact.

### 4. Build-before-publish, two-report CI model

The Docker workflow must scan the production image before merge or publication. Replace the
current build-before-scan job split with one job whose job ID remains `scan`. Keeping that ID
preserves the existing code-scanning analysis category
`.github/workflows/docker.yml:scan`, which is required for GitHub to close superseded alerts.

The `scan` job builds and loads one single-platform production image into the runner without
pushing it. Pull requests stop after the security gate. Pushes and tags log in and push the already
tested local image only after the gate passes, then expose its registry digest and tags as job
outputs. The signing job changes to `needs: scan` and signs that pushed digest. It cannot run after
a failed test or scan.

The job performs these ordered operations:

1. build and load the production image without publishing it;
2. run container behavioral tests;
3. generate `trivy-complete.json` with every severity and without `--ignore-unfixed`;
4. upload the complete JSON artifact with 90-day retention and `if: always()`;
5. generate `trivy-actionable.sarif` with every severity and `--ignore-unfixed`;
6. upload SARIF to GitHub code scanning with `if: always()`;
7. run a final actionable scan with every severity, `--ignore-unfixed`, and exit code 1;
8. on non-pull-request events only, authenticate and push every metadata-generated tag;
9. resolve and export the immutable registry digest for the signing job.

The evidence uploads precede the failing gate so a failure remains diagnosable. Trivy action and
upload actions remain version-pinned. No `.trivyignore` entry is part of this design. Pushing the
loaded image, rather than rebuilding after the gate, keeps tests, reports, publication, and signing
on one byte-identical image.

### 5. Alert lifecycle

The implementation does not PATCH alert state. After the fixed image is pushed from `main`, the
workflow uploads a new SARIF analysis under the same category
`.github/workflows/docker.yml:scan`. GitHub then determines which prior alerts are absent and closes
them according to normal code-scanning analysis semantics.

The PR completion report must record:

- the final image digest;
- complete and actionable finding counts;
- severity counts for remaining vendor-unfixed findings;
- the workflow run URL;
- whether all baseline actionable `(rule ID, package, location)` fingerprints are absent from
  actionable SARIF and, after the `main` analysis, whether GitHub closed alerts 1 through 203.

## Data and Control Flow

```text
pinned SnpEff source -----> security Maven overlay -----> patched snpEff.jar --+
                                                                              |
pinned SnpSift source ----> security Maven overlay -----> patched SnpSift.jar -+
                                                                              v
current pinned Micromamba -> conda environment -> replace JARs -> clean caches
                                                                              |
current pinned Micromamba -> Debian security upgrade -> copy clean env --------+
                                                                              v
                                                                  final runtime image
                                                                    |           |
                                                        behavioral tests       Trivy
                                                                                | \
                                                                  complete JSON  actionable SARIF+gate
                                                                                         |
                                                                            push passing image
                                                                                         |
                                                                                 resolve digest -> sign
```

The local image ID is the common identity used by tests and reports. After the gate, the pushed
registry digest becomes the common identity used by publication and signing; it must resolve back
to that tested local image.

## Error Handling and Failure Semantics

- A source checksum mismatch stops the Java build before extraction.
- Maven test, compilation, dependency resolution, or assembly failure stops the Docker build.
- A missing/duplicate JAR, wrong manifest main class, or residual stock JAR stops the Docker build.
- A Debian update failure stops the Docker build; it is never converted into a warning.
- A behavioral-contract mismatch stops the workflow before signing or release completion.
- Report upload uses `if: always()` so diagnostic evidence survives a later failure.
- Any fixed vulnerability at any severity makes the final scan return nonzero.
- Vendor-unfixed findings do not fail the actionable gate but must remain in the complete report.
- Authentication retry behavior for GHCR remains unchanged.
- A failed scan prevents publication and signing because push occurs after the gate and signing
  depends on the `scan` job.

## Testing Strategy

### Static repository tests

- Parse both security Maven overlays and assert every required version/exclusion/scope.
- Parse `conda/environment-docker.yml` and assert exact SnpEff/SnpSift package builds.
- Parse the Dockerfile and assert digest-pinned bases, non-root final user, cache removal, and expected
  stage names.
- Parse the Docker workflow and assert pull-request scanning, full audit upload, all-severity
  actionable scan, `ignore-unfixed`, and a nonzero gate.

These tests make policy regression fast and do not replace image-level tests.

### Java builder tests

- Run upstream SnpEff tests.
- Run upstream SnpSift tests against the locally installed patched SnpEff artifact.
- Run compatibility conversion unit tests.
- Emit `mvn dependency:tree` and fail if any banned coordinate/version remains.
- Inspect fat-JAR manifests and contents.

### Container tests

- Build the production target with no test-only Docker target substituted.
- Run the behavioral contract described above.
- Scan the production image with the same Trivy version/configuration as CI.

### Repository regression tests

The normal Python suite, Ruff, and mypy continue to run. The security PR does not change Python
runtime logic, so existing failures cannot be waived as unrelated without independent evidence.

## Acceptance Criteria and Authoritative Evidence

| Requirement | Evidence required |
|---|---|
| All vendor-fixed vulnerabilities remediated | Actionable Trivy JSON/SARIF contains zero vulnerabilities; gate exits 0 |
| Vendor-unfixed risk remains visible | Complete JSON artifact exists and includes unfixed findings |
| All 203 baseline alerts accounted for | Appendix maps every alert number 1-203 exactly once; final actionable report contains no baseline actionable fingerprint; post-merge dashboard marks baseline alerts closed |
| SnpEff/SnpSift behavior preserved | Golden annotation/filter container tests pass and both tools report 5.2 |
| Vulnerable Java dependencies removed | Dependency-policy test and Trivy JAR scan pass; no banned coordinates remain |
| OS fixes applied | Final-image Trivy scan has zero fixed OS findings |
| No scanner-path sleight of hand | Runtime filesystem assertion plus JAR-content scan; required tools still execute |
| Pull requests are gated | A PR workflow run builds, tests, scans, uploads evidence, and passes the nonzero gate |
| Published artifact is the tested artifact | Test, scan, sign, and push records reference the same image digest |
| Runtime stays non-root and functional | Container UID, health check, CLI, and external-tool tests pass |
| Existing application remains healthy | Python tests, Ruff, and mypy pass |

## Rollout and Recovery

1. Merge only after the pull-request image has zero actionable findings.
2. Let the `main` workflow build and scan locally, publish the passing image, and resolve its
   immutable digest before signing.
3. Confirm GitHub code scanning ingests the new same-category SARIF analysis.
4. Confirm baseline actionable alerts close without API mutation.
5. Retain the preceding signed image tag/digest for rollback.

If a patched JAR changes output, do not ship it and do not suppress its scan findings. Revert the
dependency substitution responsible, add a focused compatibility test, and choose a compatible
fixed release or source-level patch. If a new vendor-fixed CVE appears during implementation, the
all-severity actionable gate expands the PR scope automatically; the PR is not complete until that
finding is also removed.

## Alternatives Rejected

### Upgrade to stock SnpEff/SnpSift 5.4.0c

Rejected because a direct Trivy probe found the vulnerable Java coordinates still present. It also
introduces avoidable biological-output drift into a security PR.

### Suppress Java findings with VEX or `.trivyignore`

Rejected because 23 of the 24 current JAR alerts have fixed dependencies. Reporting policy is not
a substitute for upgrading fixable code.

### Remove SnpEff/SnpSift from the image

Rejected because it breaks the documented all-in-one container contract.

### Treat every stderr finding as fatal but leave CI report-only

Rejected because runtime behavior does not prevent vulnerable image publication, and a green
report-only workflow cannot enforce future remediation.

### Literal zero-finding dashboard without complete audit evidence

Rejected because it would hide vendor-unfixed risk. The chosen model gives GitHub an actionable
queue while retaining the unfiltered machine-readable record.

## Appendix A: Complete Baseline Alert Inventory

This table was generated from the GitHub code-scanning API on 2026-07-19. `Fixed` and `Unfixed`
refer to whether Trivy reported a non-empty fixed version. Each alert number from 1 through 203
appears exactly once.

| Surface | Package | Count | Severities | Fixed | Unfixed | GitHub alert numbers |
|---|---|---:|---|---:|---:|---|
| OS | `apt` | 1 | low | 0 | 1 | 1 |
| OS | `bash` | 1 | low | 0 | 1 | 2 |
| OS | `bsdutils` | 4 | low, medium | 0 | 4 | 3, 4, 139, 147 |
| OS | `coreutils` | 3 | low | 0 | 3 | 5, 6, 7 |
| OS | `dpkg` | 2 | low, medium | 2 | 0 | 8, 148 |
| OS | `gcc-12-base` | 2 | low | 1 | 1 | 9, 10 |
| OS | `gpgv` | 4 | high, low, medium | 1 | 3 | 11, 12, 13, 14 |
| OS | `libapt-pkg6.0` | 1 | low | 0 | 1 | 15 |
| OS | `libblkid1` | 4 | low, medium | 0 | 4 | 16, 17, 140, 149 |
| OS | `libc-bin` | 19 | high, low, medium | 8 | 11 | 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 150, 151, 152, 153, 154, 155, 156 |
| OS | `libc6` | 19 | high, low, medium | 8 | 11 | 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 157, 158, 159, 160, 161, 162, 163 |
| OS | `libcap2` | 2 | high, medium | 2 | 0 | 42, 164 |
| OS | `libgcc-s1` | 2 | low | 1 | 1 | 43, 44 |
| OS | `libgcrypt20` | 3 | high, low | 0 | 3 | 45, 46, 165 |
| OS | `libgnutls30` | 20 | critical, high, low, medium, unknown | 6 | 14 | 47, 48, 49, 50, 51, 52, 53, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178 |
| OS | `liblzma5` | 2 | high, medium | 1 | 1 | 54, 179 |
| OS | `libmount1` | 4 | low, medium | 0 | 4 | 55, 56, 141, 180 |
| OS | `libpam-modules` | 3 | high, medium | 2 | 1 | 57, 58, 59 |
| OS | `libpam-modules-bin` | 3 | high, medium | 2 | 1 | 60, 61, 62 |
| OS | `libpam-runtime` | 3 | high, medium | 2 | 1 | 63, 64, 65 |
| OS | `libpam0g` | 3 | high, medium | 2 | 1 | 66, 67, 68 |
| OS | `libsmartcols1` | 4 | low, medium | 0 | 4 | 69, 70, 142, 181 |
| OS | `libstdc++6` | 2 | low | 1 | 1 | 71, 72 |
| OS | `libsystemd0` | 10 | high, low, medium | 5 | 5 | 73, 74, 75, 76, 77, 182, 183, 184, 185, 186 |
| OS | `libtasn1-6` | 1 | medium | 0 | 1 | 78 |
| OS | `libtinfo6` | 3 | high, low, medium | 0 | 3 | 79, 80, 187 |
| OS | `libudev1` | 10 | high, low, medium | 5 | 5 | 81, 82, 83, 84, 85, 188, 189, 190, 191, 192 |
| OS | `libuuid1` | 4 | low, medium | 0 | 4 | 86, 87, 143, 193 |
| OS | `login` | 5 | low, medium | 2 | 3 | 88, 89, 90, 91, 92 |
| OS | `mount` | 4 | low, medium | 0 | 4 | 93, 94, 144, 194 |
| OS | `ncurses-base` | 3 | high, low, medium | 0 | 3 | 95, 96, 195 |
| OS | `ncurses-bin` | 3 | high, low, medium | 0 | 3 | 97, 98, 196 |
| OS | `passwd` | 5 | low, medium | 2 | 3 | 99, 100, 101, 102, 103 |
| OS | `perl-base` | 5 | high, low, medium | 3 | 2 | 104, 105, 106, 107, 108 |
| OS | `sed` | 1 | medium | 1 | 0 | 197 |
| OS | `sysvinit-utils` | 1 | low | 0 | 1 | 109 |
| OS | `tar` | 3 | low, medium | 0 | 3 | 110, 111, 198 |
| OS | `util-linux` | 4 | low, medium | 0 | 4 | 112, 113, 145, 199 |
| OS | `util-linux-extra` | 4 | low, medium | 0 | 4 | 114, 115, 146, 200 |
| OS | `zlib1g` | 2 | critical, medium | 0 | 2 | 116, 138 |
| SnpEff JAR | `com.fasterxml.jackson.core:jackson-core` | 3 | high, medium | 3 | 0 | 117, 118, 201 |
| SnpEff JAR | `com.fasterxml.jackson.core:jackson-databind` | 5 | high | 5 | 0 | 119, 120, 121, 122, 123 |
| SnpEff JAR | `com.google.code.gson:gson` | 1 | high | 1 | 0 | 124 |
| SnpEff JAR | `commons-io:commons-io` | 1 | high | 1 | 0 | 125 |
| SnpEff JAR | `commons-lang:commons-lang` | 1 | medium | 0 | 1 | 126 |
| SnpEff JAR | `org.apache.commons:commons-compress` | 5 | high, medium | 5 | 0 | 127, 129, 131, 133, 135 |
| SnpEff JAR | `org.apache.logging.log4j:log4j-core` | 3 | medium | 3 | 0 | 137, 202, 203 |
| SnpSift JAR | `org.apache.commons:commons-compress` | 5 | high, medium | 5 | 0 | 128, 130, 132, 134, 136 |

## Appendix B: Review Checklist

- [x] Every baseline alert number is represented.
- [x] Fixed and unfixed findings are not conflated.
- [x] The design changes vulnerable runtime content, not only scanner paths.
- [x] Pull-request and post-push scanning are both covered.
- [x] Functional and security regressions have authoritative acceptance evidence.
- [x] No alert dismissal API, `.trivyignore`, placeholder, or severity exception is required.
- [x] The work remains one coherent container-security pull request.
