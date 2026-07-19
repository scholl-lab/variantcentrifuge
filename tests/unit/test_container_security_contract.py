import os
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MAVEN_NAMESPACE = {"m": "http://maven.apache.org/POM/4.0.0"}

SECURED_PROPERTIES = {
    "jackson.version": "2.22.1",
    "jackson.annotations.version": "2.22",
    "gson.version": "2.14.0",
    "commons.io.version": "2.22.0",
    "commons.compress.version": "1.28.0",
    "log4j.version": "2.25.4",
}

MANAGED_RUNTIME_VERSIONS = {
    ("com.fasterxml.jackson.core", "jackson-core"): "${jackson.version}",
    ("com.fasterxml.jackson.core", "jackson-databind"): "${jackson.version}",
    ("com.fasterxml.jackson.core", "jackson-annotations"): "${jackson.annotations.version}",
    ("com.google.code.gson", "gson"): "${gson.version}",
    ("commons-io", "commons-io"): "${commons.io.version}",
    ("org.apache.commons", "commons-compress"): "${commons.compress.version}",
    ("org.apache.logging.log4j", "log4j-api"): "${log4j.version}",
    ("org.apache.logging.log4j", "log4j-core"): "${log4j.version}",
    ("org.apache.logging.log4j", "log4j-slf4j-impl"): "${log4j.version}",
}

JUNIT_COORDINATES = {
    ("org.junit.jupiter", "junit-jupiter-api"),
    ("org.junit.jupiter", "junit-jupiter-engine"),
    ("org.junit.platform", "junit-platform-suite-api"),
    ("org.junit.platform", "junit-platform-suite-engine"),
}

REQUIRED_RUNTIME_DEPENDENCIES = [
    "com.fasterxml.jackson.core:jackson-core:jar:2.22.1:runtime",
    "com.fasterxml.jackson.core:jackson-databind:jar:2.22.1:runtime",
    "com.fasterxml.jackson.core:jackson-annotations:jar:2.22:runtime",
    "com.google.code.gson:gson:jar:2.14.0:runtime",
    "commons-io:commons-io:jar:2.22.0:runtime",
    "org.apache.commons:commons-compress:jar:1.28.0:runtime",
    "org.apache.logging.log4j:log4j-api:jar:2.25.4:runtime",
    "org.apache.logging.log4j:log4j-core:jar:2.25.4:runtime",
    "org.apache.logging.log4j:log4j-slf4j-impl:jar:2.25.4:runtime",
]

# Derived from the pinned upstream POMs, with only the approved dependency,
# scope, and secured-version overlay changes applied.
EXPECTED_DIRECT_DEPENDENCIES = {
    "docker/java/snpeff-pom.xml": {
        ("org.apfloat", "apfloat", "1.10.1", "compile"),
        ("com.googlecode.charts4j", "charts4j", "1.3", "compile"),
        ("commons-cli", "commons-cli", "1.5.0", "compile"),
        ("org.junit.jupiter", "junit-jupiter-api", None, "test"),
        ("org.junit.jupiter", "junit-jupiter-engine", None, "test"),
        ("org.junit.platform", "junit-platform-suite-api", "1.8.2", "test"),
        ("org.junit.platform", "junit-platform-suite-engine", "1.8.2", "test"),
        ("net.sf.trove4j", "trove4j", "3.0.2", "compile"),
        ("org.freemarker", "freemarker", "2.3.31", "compile"),
        ("distlib", "distlib", "0.9.1", "compile"),
        ("org.apache.commons", "commons-math3", "3.6.1", "compile"),
        ("commons-io", "commons-io", "${commons.io.version}", "compile"),
        ("commons-codec", "commons-codec", "1.15", "compile"),
        ("org.biojava", "biojava-core", "6.0.4", "compile"),
        ("org.biojava", "biojava-structure", "6.0.4", "compile"),
        ("com.github.samtools", "htsjdk", "2.24.1", "compile"),
        ("javax.xml.bind", "jaxb-api", "2.3.1", "compile"),
        ("com.fasterxml.jackson.core", "jackson-core", "${jackson.version}", "compile"),
        (
            "com.fasterxml.jackson.core",
            "jackson-databind",
            "${jackson.version}",
            "compile",
        ),
        (
            "com.fasterxml.jackson.core",
            "jackson-annotations",
            "${jackson.annotations.version}",
            "compile",
        ),
        ("com.google.code.gson", "gson", "${gson.version}", "compile"),
        ("org.apache.commons", "commons-compress", "${commons.compress.version}", "compile"),
        ("org.apache.logging.log4j", "log4j-api", "${log4j.version}", "compile"),
        ("org.apache.logging.log4j", "log4j-core", "${log4j.version}", "compile"),
        (
            "org.apache.logging.log4j",
            "log4j-slf4j-impl",
            "${log4j.version}",
            "compile",
        ),
    },
    "docker/java/snpsift-pom.xml": {
        ("org.snpeff", "SnpEff", "5.2", "compile"),
        ("org.antlr", "antlr4", "4.9.3", "compile"),
        ("net.sf.trove4j", "trove4j", "3.0.2", "compile"),
        ("org.junit.jupiter", "junit-jupiter-api", None, "test"),
        ("org.junit.jupiter", "junit-jupiter-engine", None, "test"),
        ("org.junit.platform", "junit-platform-suite-api", "1.8.2", "test"),
        ("org.junit.platform", "junit-platform-suite-engine", "1.8.2", "test"),
        ("org.apache.commons", "commons-math3", "3.6.1", "compile"),
        ("com.github.samtools", "htsjdk", "2.24.1", "compile"),
        ("com.fasterxml.jackson.core", "jackson-core", "${jackson.version}", "compile"),
        (
            "com.fasterxml.jackson.core",
            "jackson-databind",
            "${jackson.version}",
            "compile",
        ),
        (
            "com.fasterxml.jackson.core",
            "jackson-annotations",
            "${jackson.annotations.version}",
            "compile",
        ),
        ("com.google.code.gson", "gson", "${gson.version}", "compile"),
        ("commons-io", "commons-io", "${commons.io.version}", "compile"),
        ("org.apache.commons", "commons-compress", "${commons.compress.version}", "compile"),
        ("org.apache.logging.log4j", "log4j-api", "${log4j.version}", "compile"),
        ("org.apache.logging.log4j", "log4j-core", "${log4j.version}", "compile"),
        (
            "org.apache.logging.log4j",
            "log4j-slf4j-impl",
            "${log4j.version}",
            "compile",
        ),
    },
}

EXPECTED_REPOSITORIES = {
    "docker/java/snpeff-pom.xml": {
        ("maven", "https://repo1.maven.org/maven2/"),
        ("central", "https://repo.maven.apache.org/maven2"),
        ("hadoop-bam", "https://hadoop-bam.sourceforge.net/maven/"),
        ("ncimvn-public", "https://ncimvn.nci.nih.gov/nexus/content/groups/public/"),
        ("typesafe", "https://repo.typesafe.com/typesafe/releases/"),
    },
    "docker/java/snpsift-pom.xml": {
        ("maven", "https://repo1.maven.org/maven2/"),
        ("ncimvn-public", "https://ncimvn.nci.nih.gov/nexus/content/groups/public/"),
    },
}


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def _workflow_step(workflow: str, name: str, next_name: str) -> str:
    """Return one workflow step, bounded by two unique step names."""
    start_marker = f"      - name: {name}\n"
    end_marker = f"      - name: {next_name}\n"
    assert workflow.count(start_marker) == 1
    assert workflow.count(end_marker) == 1
    return workflow.split(start_marker, maxsplit=1)[1].split(end_marker, maxsplit=1)[0]


def _docker_stage(dockerfile: str, stage: str) -> str:
    marker = f" AS {stage}\n"
    assert marker in dockerfile
    stage_body = dockerfile.split(marker, maxsplit=1)[1]
    return stage_body.split("\nFROM ", maxsplit=1)[0]


def _pom(relative: str) -> ET.Element:
    return ET.parse(ROOT / relative).getroot()


def _coordinate(dependency: ET.Element) -> tuple[str, str]:
    group_id = dependency.findtext("m:groupId", namespaces=MAVEN_NAMESPACE)
    artifact_id = dependency.findtext("m:artifactId", namespaces=MAVEN_NAMESPACE)
    assert group_id is not None
    assert artifact_id is not None
    return group_id, artifact_id


def _dependency_elements(pom: ET.Element, path: str) -> list[ET.Element]:
    return pom.findall(path, MAVEN_NAMESPACE)


def _dependency_versions(
    pom: ET.Element,
    path: str,
) -> dict[tuple[str, str], list[str | None]]:
    versions: dict[tuple[str, str], list[str | None]] = {}
    for dependency in _dependency_elements(pom, path):
        versions.setdefault(_coordinate(dependency), []).append(
            dependency.findtext("m:version", namespaces=MAVEN_NAMESPACE)
        )
    return versions


def _direct_dependency_contract(
    pom: ET.Element,
) -> set[tuple[str, str, str | None, str]]:
    return {
        (
            *_coordinate(dependency),
            dependency.findtext("m:version", namespaces=MAVEN_NAMESPACE),
            dependency.findtext("m:scope", "compile", MAVEN_NAMESPACE),
        )
        for dependency in _dependency_elements(pom, "m:dependencies/m:dependency")
    }


def _repository_contract(pom: ET.Element) -> set[tuple[str, str]]:
    repositories = set()
    for repository in pom.findall("m:repositories/m:repository", MAVEN_NAMESPACE):
        repository_id = repository.findtext("m:id", namespaces=MAVEN_NAMESPACE)
        repository_url = repository.findtext("m:url", namespaces=MAVEN_NAMESPACE)
        assert repository_id is not None
        assert repository_url is not None
        repositories.add((repository_id, repository_url))
    return repositories


def _plugin(pom: ET.Element, artifact_id: str) -> ET.Element:
    for plugin in pom.findall("m:build/m:plugins/m:plugin", MAVEN_NAMESPACE):
        if plugin.findtext("m:artifactId", namespaces=MAVEN_NAMESPACE) == artifact_id:
            return plugin
    raise AssertionError(f"missing Maven plugin: {artifact_id}")


def _run_dependency_assertion(
    tmp_path: Path,
    dependencies: list[str],
) -> subprocess.CompletedProcess[str]:
    return _run_dependency_assertion_output(
        tmp_path,
        [f"   {dependency} -- module fake.module [auto]" for dependency in dependencies],
    )


def _run_dependency_assertion_output(
    tmp_path: Path,
    dependency_lines: list[str],
) -> subprocess.CompletedProcess[str]:
    project = tmp_path / "project"
    project.mkdir()
    (project / "pom.xml").write_text("<project/>\n", encoding="utf-8")

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_maven = fake_bin / "mvn"
    fake_maven.write_text(
        """#!/usr/bin/env bash
set -euo pipefail
for argument in "$@"; do
    case "$argument" in
        -DoutputFile=*) output_file=${argument#-DoutputFile=} ;;
    esac
done
: "${output_file:?missing Maven dependency output path}"
printf '%s\n' "$FAKE_MAVEN_DEPENDENCIES" > "$output_file"
""",
        encoding="utf-8",
    )
    fake_maven.chmod(0o755)

    environment = {
        **os.environ,
        "PATH": f"{fake_bin}:{os.environ['PATH']}",
        "FAKE_MAVEN_DEPENDENCIES": "\n".join(dependency_lines),
    }
    return subprocess.run(
        [str(ROOT / "docker/java/assert-runtime-dependencies.sh"), str(project)],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )


def _contains_dependency_line(environment: str, expected: str) -> bool:
    return expected in environment.splitlines()


@pytest.mark.parametrize(
    "environment",
    [
        "#  - snpeff=5.2=hdfd78af_3\n",
        "    - snpeff=5.2=hdfd78af_3\n",
    ],
)
def test_dependency_matching_rejects_inactive_lines(environment: str) -> None:
    assert not _contains_dependency_line(
        environment,
        "  - snpeff=5.2=hdfd78af_3",
    )


def test_docker_environment_pins_release_era_java_tool_builds() -> None:
    environment = _text("conda/environment-docker.yml")
    assert _contains_dependency_line(environment, "  - snpeff=5.2=hdfd78af_3")
    assert _contains_dependency_line(environment, "  - snpsift=5.2=hdfd78af_0")


def test_docker_environment_has_secure_pip_floor() -> None:
    environment = _text("conda/environment-docker.yml")
    assert _contains_dependency_line(environment, "  - pip>=26.1.2")


def test_docker_environment_pins_compatible_smart_open() -> None:
    environment = _text("conda/environment-docker.yml")
    assert _contains_dependency_line(environment, "    - intervaltree")
    assert _contains_dependency_line(environment, "    - smart-open==7.7.1")


def test_docker_environment_installs_all_declared_binary_runtime_dependencies() -> None:
    environment = _text("conda/environment-docker.yml")
    for dependency in (
        "  - cffi=2.1.0",
        "  - pyarrow=25.0.0",
        "  - xlsxwriter=3.2.9",
    ):
        assert _contains_dependency_line(environment, dependency)


def test_dockerfile_pins_the_java_builder_image_by_digest() -> None:
    dockerfile = _text("Dockerfile")
    assert (
        "FROM maven:3.9.11-eclipse-temurin-11@sha256:"
        "c095e2421eaf3e5cb1573fd0474e68e17062866d454362349e17bbb75f44031e "
        "AS java-build"
    ) in dockerfile.splitlines()


def test_dockerfile_pins_both_micromamba_stages_by_digest() -> None:
    dockerfile = _text("Dockerfile")
    image = (
        "mambaorg/micromamba:2.8.1-debian12-slim@sha256:"
        "c8198d53228ad7cfd7adcc0704e8837f9d1c9327fb363c6a7d3c5b1a51a4b561"
    )
    assert f"FROM {image} AS conda-build" in dockerfile.splitlines()
    assert f"FROM {image} AS runtime" in dockerfile.splitlines()
    assert dockerfile.count(f"FROM {image} AS ") == 2


def test_dockerfile_fetches_pinned_java_sources_with_verified_checksums() -> None:
    dockerfile = _text("Dockerfile")
    assert (
        "ADD --checksum=sha256:"
        "8633ecad8dbf06af19d5002818b65dc9e7b78419b94866dd4b1daed9d1e9d2b2 "
        "https://github.com/pcingola/SnpEff/archive/"
        "0c5e74f9b6ca6ed3db720177eb1f95b9d47d45f2.tar.gz /tmp/snpeff.tar.gz"
    ) in dockerfile
    assert (
        "ADD --checksum=sha256:"
        "74c37b85e74a390a27f122164fd06c5b56131e3fa1eca205636d3de4a8c94934 "
        "https://github.com/pcingola/SnpSift/archive/"
        "20978614457f14ec7a0c70539d5a7a2b7e754f60.tar.gz /tmp/snpsift.tar.gz"
    ) in dockerfile


def test_dockerfile_upgrades_runtime_os_and_removes_package_caches() -> None:
    dockerfile = _text("Dockerfile")
    runtime = _docker_stage(dockerfile, "runtime")
    assert "apt-get update" in runtime
    assert "DEBIAN_FRONTEND=noninteractive" in runtime
    assert "apt-get -y --no-install-recommends dist-upgrade" in runtime
    assert "apt-get clean" in runtime
    assert "rm -rf /var/lib/apt/lists/*" in runtime

    conda_build = _docker_stage(dockerfile, "conda-build")
    assert "rm -rf /opt/conda/pkgs" in conda_build
    assert "test ! -e /opt/conda/pkgs" in conda_build


def test_dockerfile_replaces_the_exact_conda_java_tool_jars() -> None:
    dockerfile = _text("Dockerfile")
    assert (
        "COPY --from=java-build --chown=$MAMBA_USER:$MAMBA_USER /out/snpEff.jar "
        "/opt/conda/share/snpeff-5.2-3/snpEff.jar"
    ) in dockerfile
    assert (
        "COPY --from=java-build --chown=$MAMBA_USER:$MAMBA_USER /out/SnpSift.jar "
        "/opt/conda/share/snpsift-5.2-0/SnpSift.jar"
    ) in dockerfile


def test_dockerfile_routes_default_snpeff_databases_to_writable_data_storage() -> None:
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")
    runtime = _docker_stage(_text("Dockerfile"), "runtime")
    config_path = "/opt/conda/share/snpeff-5.2-3/snpEff.config"

    assert f"snpeff_config={config_path}" in conda_build
    assert "data_dir_pattern='^[[:space:]]*data[.]dir[[:space:]]*='" in conda_build
    assert "sed -i 's|^data.dir = \\./data/$|data.dir = /data/snpeff/|'" in conda_build
    assert conda_build.count('test "$(grep -Ec "$data_dir_pattern" "$snpeff_config")" -eq 1') == 2
    assert "grep -Fx 'data.dir = ./data/' \"$snpeff_config\"" in conda_build
    assert "grep -Fx 'data.dir = /data/snpeff/' \"$snpeff_config\"" in conda_build
    assert "mkdir -p /data/snpeff" in runtime
    assert "chown $MAMBA_USER:$MAMBA_USER /data /data/snpeff" in runtime
    assert "chmod 0750 /data /data/snpeff" in runtime


@pytest.mark.parametrize(
    "duplicate_assignment",
    [
        "data.dir=/shadow",
        " data.dir = /shadow",
        "\tdata.dir\t=\t/shadow",
    ],
)
def test_snpeff_data_assignment_guard_counts_whitespace_variants(
    tmp_path: Path, duplicate_assignment: str
) -> None:
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")
    pattern_match = re.search(r"data_dir_pattern='([^']+)'", conda_build)
    assert pattern_match is not None
    config = tmp_path / "snpEff.config"
    config.write_text(
        f"# data.dir = ignored\ndata.dir = ./data/\n{duplicate_assignment}\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        ["grep", "-Ec", pattern_match.group(1), str(config)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip() == "2"


def test_dockerfile_removes_jdk_only_tools_and_rejects_javac() -> None:
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")
    assert 'case "$tool" in' in conda_build
    assert "java|keytool" in conda_build
    assert "jspawnhelper" not in conda_build
    assert 'find "$jvm_bin" -mindepth 1 -maxdepth 1' in conda_build
    assert 'rm -f "$jvm_bin/$tool"' in conda_build
    assert "for link in /opt/conda/bin/*" in conda_build
    assert 'readlink "$link"' in conda_build
    assert 'rm -rf "$jvm_home/include" "$jvm_home/jmods"' in conda_build
    assert 'rm -f "$jvm_home/src.zip" "$jvm_home/lib/src.zip"' in conda_build
    assert "Unexpected JVM runtime tool" in conda_build
    assert "Unexpected conda JVM tool link" in conda_build
    assert "! command -v javac" in conda_build


def test_dockerfile_final_user_is_the_micromamba_user() -> None:
    dockerfile = _text("Dockerfile")
    user_lines = [line for line in dockerfile.splitlines() if line.startswith("USER ")]
    assert user_lines
    assert user_lines[-1] == "USER $MAMBA_USER"
    assert "USER root" in _docker_stage(dockerfile, "runtime")


def test_dockerfile_preserves_the_runtime_interface() -> None:
    dockerfile = _text("Dockerfile")
    assert "WORKDIR /data" in dockerfile
    assert 'CMD ["/opt/conda/bin/variantcentrifuge", "--version"]' in dockerfile
    assert 'ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "variantcentrifuge"]' in dockerfile
    assert 'CMD ["--help"]' in dockerfile
    assert "COPY --chown=0:0 LICENSE /app/LICENSE" in dockerfile
    assert "COPY --chown=0:0 scoring/ /app/scoring/" in dockerfile
    assert "COPY --chown=0:0 stats_configs/ /app/stats_configs/" in dockerfile


def test_dockerfile_build_gates_all_python_runtime_paths() -> None:
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")
    assert "COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md setup.py" in (conda_build)
    assert "/opt/conda/bin/python -m pip check" in conda_build
    for required_import in (
        "import cffi",
        "import pyarrow as pa",
        "import xlsxwriter",
    ):
        assert required_import in conda_build
    assert "pa.table" in conda_build
    assert ".to_pandas()" in conda_build
    assert "io.BytesIO()" in conda_build
    assert 'engine="xlsxwriter"' in conda_build
    assert "_try_load_davies" in conda_build
    assert "davies_pvalue" in conda_build


def test_production_package_and_image_gate_both_default_json_configs() -> None:
    pyproject = _text("pyproject.toml")
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")

    assert "[tool.setuptools.package-data]" in pyproject
    assert 'variantcentrifuge = ["config.json", "default_stats_config.json"]' in pyproject
    assert "from variantcentrifuge.config import load_config" in conda_build
    assert 'package_dir.is_relative_to(Path("/opt/conda"))' in conda_build
    assert 'not package_dir.is_relative_to(Path("/tmp/src"))' in conda_build
    assert "config = load_config()" in conda_build
    assert 'package_dir / "default_stats_config.json"' in conda_build


def test_builder_removes_conda_test_payloads_before_the_runtime_copy() -> None:
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")
    runtime = _docker_stage(_text("Dockerfile"), "runtime")
    cleanup = "rm -rf /opt/conda/etc/conda/test-files"
    normalization = "chmod -R go-w /opt/conda"

    assert cleanup in conda_build
    assert "test ! -e /opt/conda/etc/conda/test-files" in conda_build
    assert conda_build.index(cleanup) < conda_build.index(normalization)
    assert runtime.index("COPY --from=conda-build") < runtime.index("USER $MAMBA_USER")


def test_dockerfile_makes_runtime_payloads_immutable_to_the_service_user() -> None:
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")
    runtime = _docker_stage(_text("Dockerfile"), "runtime")
    conda_mode_normalization = "chmod -R go-w /opt/conda"
    assert conda_mode_normalization in conda_build
    assert conda_build.index("! command -v javac") < conda_build.index(conda_mode_normalization)
    runtime_copy = "COPY --from=conda-build --chown=root:root /opt/conda /opt/conda"
    runtime_mountpoint_normalization = "RUN chmod go-w /opt/conda"
    assert runtime_mountpoint_normalization in runtime
    assert runtime_copy in runtime
    assert runtime.index(runtime_mountpoint_normalization) < runtime.index(runtime_copy)
    after_runtime_copy = runtime.split(runtime_copy, maxsplit=1)[1].replace("\\\n", " ")
    assert not re.search(r"\b(?:chmod|chown)\s+-R\b[^\n]*\b/opt/conda\b", after_runtime_copy)
    assert "chmod -R go-w /app" in runtime
    assert "! -type l -a -perm /022" in runtime
    assert "chown $MAMBA_USER:$MAMBA_USER /data" in runtime
    assert "chmod 0750 /data" in runtime
    assert "USER $MAMBA_USER" in runtime


def test_dockerfile_validates_both_executable_multi_release_manifests() -> None:
    conda_build = _docker_stage(_text("Dockerfile"), "conda-build")
    assert "python - <<'PY'" in conda_build
    assert "zipfile.ZipFile" in conda_build
    assert '"/opt/conda/share/snpeff-5.2-3/snpEff.jar": "org.snpeff.SnpEff"' in conda_build
    assert '"/opt/conda/share/snpsift-5.2-0/SnpSift.jar": "org.snpsift.SnpSift"' in conda_build
    assert 'manifest.get("Main-Class") == expected_main_class' in conda_build
    assert 'manifest.get("Multi-Release") == "true"' in conda_build


def test_dockerfile_verifies_both_java_projects_and_runtime_dependencies() -> None:
    java_build = _docker_stage(_text("Dockerfile"), "java-build")
    assert java_build.count("mvn -B verify") == 2
    assert "mvn -B verify -DskipTests" not in java_build
    assert "mvn -B install -DskipTests -Dassembly.skipAssembly=true" in java_build
    assert "mvn -B install -DskipTests &&" not in java_build
    assert java_build.count("/usr/local/bin/assert-runtime-dependencies.sh /build/") == 2
    assert java_build.count("*-jar-with-dependencies.jar") == 3
    assert "verified_sha256=$(sha256sum" in java_build
    assert "installed_sha256=$(sha256sum" in java_build
    assert 'test "$verified_sha256" = "$installed_sha256"' in java_build
    assert java_build.index("mvn -B verify") < java_build.index(
        "mvn -B install -DskipTests -Dassembly.skipAssembly=true"
    )


@pytest.mark.parametrize(
    ("pom_path", "expected_group", "expected_artifact", "expected_main_class"),
    [
        ("docker/java/snpeff-pom.xml", "org.snpeff", "SnpEff", "org.snpeff.SnpEff"),
        ("docker/java/snpsift-pom.xml", "org.snpsift", "SnpSift", "org.snpsift.SnpSift"),
    ],
)
def test_patched_java_projects_build_executable_release_5_2_jars(
    pom_path: str,
    expected_group: str,
    expected_artifact: str,
    expected_main_class: str,
) -> None:
    pom = _pom(pom_path)
    assert pom.findtext("m:groupId", namespaces=MAVEN_NAMESPACE) == expected_group
    assert pom.findtext("m:artifactId", namespaces=MAVEN_NAMESPACE) == expected_artifact
    assert pom.findtext("m:version", namespaces=MAVEN_NAMESPACE) == "5.2"

    assembly = _plugin(pom, "maven-assembly-plugin")
    assert assembly.findtext("m:version", namespaces=MAVEN_NAMESPACE) == "3.7.1"
    assert (
        assembly.findtext(
            "m:configuration/m:archive/m:manifest/m:mainClass",
            namespaces=MAVEN_NAMESPACE,
        )
        == expected_main_class
    )
    assert (
        assembly.findtext("m:executions/m:execution/m:phase", namespaces=MAVEN_NAMESPACE)
        == "package"
    )
    assert (
        assembly.findtext("m:executions/m:execution/m:goals/m:goal", namespaces=MAVEN_NAMESPACE)
        == "single"
    )
    assert (
        assembly.findtext(
            "m:configuration/m:archive/m:manifestEntries/m:Multi-Release",
            namespaces=MAVEN_NAMESPACE,
        )
        == "true"
    )


@pytest.mark.parametrize(
    "pom_path",
    ["docker/java/snpeff-pom.xml", "docker/java/snpsift-pom.xml"],
)
def test_patched_java_projects_use_required_build_policy(pom_path: str) -> None:
    pom = _pom(pom_path)
    properties = pom.find("m:properties", MAVEN_NAMESPACE)
    assert properties is not None
    assert properties.findtext("m:maven.compiler.release", namespaces=MAVEN_NAMESPACE) == "11"
    assert (
        properties.findtext("m:project.build.sourceEncoding", namespaces=MAVEN_NAMESPACE)
        == "ISO-8859-1"
    )
    assert (
        properties.findtext("m:project.reporting.outputEncoding", namespaces=MAVEN_NAMESPACE)
        == "ISO-8859-1"
    )

    surefire = _plugin(pom, "maven-surefire-plugin")
    assert surefire.findtext("m:version", namespaces=MAVEN_NAMESPACE) == "3.5.4"


@pytest.mark.parametrize(
    "pom_path",
    ["docker/java/snpeff-pom.xml", "docker/java/snpsift-pom.xml"],
)
def test_patched_java_projects_preserve_stock_build_timestamp(pom_path: str) -> None:
    pom = _pom(pom_path)
    properties = pom.find("m:properties", MAVEN_NAMESPACE)
    assert properties is not None
    assert (
        properties.findtext("m:project.build.outputTimestamp", namespaces=MAVEN_NAMESPACE)
        == "2023-09-29T06:17:00Z"
    )


def test_snpeff_verify_runs_only_hermetic_upstream_and_compatibility_tests() -> None:
    pom = _pom("docker/java/snpeff-pom.xml")
    surefire = _plugin(pom, "maven-surefire-plugin")
    includes = {
        element.text
        for element in surefire.findall("m:configuration/m:includes/m:include", MAVEN_NAMESPACE)
    }
    assert includes == {
        "org/snpeff/snpEffect/testCases/unity/*.java",
        "org/apache/commons/lang/ArrayUtilsTest.java",
    }
    excludes = {
        element.text
        for element in surefire.findall("m:configuration/m:excludes/m:exclude", MAVEN_NAMESPACE)
    }
    assert excludes == {
        "org/snpeff/snpEffect/testCases/unity/TestCasesCytoBands.java",
        "org/snpeff/snpEffect/testCases/unity/TestCasesGenomicSequences.java",
    }


@pytest.mark.parametrize(
    "pom_path",
    ["docker/java/snpeff-pom.xml", "docker/java/snpsift-pom.xml"],
)
def test_patched_java_projects_keep_junit_dependencies_test_only(pom_path: str) -> None:
    pom = _pom(pom_path)
    dependencies = {
        _coordinate(dependency): dependency
        for dependency in _dependency_elements(pom, "m:dependencies/m:dependency")
    }
    assert dependencies.keys() >= JUNIT_COORDINATES
    for coordinate in JUNIT_COORDINATES:
        assert dependencies[coordinate].findtext("m:scope", namespaces=MAVEN_NAMESPACE) == "test"


def test_snpsift_uses_the_locally_patched_snpeff_release() -> None:
    pom = _pom("docker/java/snpsift-pom.xml")
    versions = _dependency_versions(pom, "m:dependencies/m:dependency")
    assert versions[("org.snpeff", "SnpEff")] == ["5.2"]


def test_snpsift_reclassifies_upstream_tests_as_test_sources() -> None:
    pom = _pom("docker/java/snpsift-pom.xml")
    assert (
        pom.findtext("m:build/m:testSourceDirectory", namespaces=MAVEN_NAMESPACE)
        == "${project.basedir}/src/main/java"
    )

    compiler = _plugin(pom, "maven-compiler-plugin")
    main_excludes = {
        element.text
        for element in compiler.findall("m:configuration/m:excludes/m:exclude", MAVEN_NAMESPACE)
    }
    assert main_excludes == {"org/snpsift/testCases/**"}

    test_includes = {
        element.text
        for element in compiler.findall(
            "m:configuration/m:testIncludes/m:testInclude", MAVEN_NAMESPACE
        )
    }
    assert test_includes == {"org/snpsift/testCases/**/*.java"}


def test_snpsift_verify_runs_only_archive_contained_upstream_tests() -> None:
    pom = _pom("docker/java/snpsift-pom.xml")
    surefire = _plugin(pom, "maven-surefire-plugin")
    includes = {
        element.text
        for element in surefire.findall("m:configuration/m:includes/m:include", MAVEN_NAMESPACE)
    }
    assert includes == {
        f"org/snpsift/testCases/unit/{test_class}.java"
        for test_class in {
            "TestCasesCaseControl",
            "TestCasesConcordance",
            "TestCasesExtractFields",
            "TestCasesFilter",
            "TestCasesFilterALL",
            "TestCasesFilterChrPos",
            "TestCasesFilterGt",
            "TestCasesGeneSets",
            "TestCasesGt",
            "TestCasesHwe",
            "TestCasesIndex",
            "TestCasesIntervals",
            "TestCasesLd",
            "TestCasesPrivate",
            "TestCasesSort",
            "TestCasesSplit",
            "TestCasesVarType",
            "TestCasesVcfCheck",
        }
    }


@pytest.mark.parametrize("pom_path", EXPECTED_DIRECT_DEPENDENCIES)
def test_patched_java_projects_lock_pinned_direct_dependency_graph(
    pom_path: str,
) -> None:
    assert _direct_dependency_contract(_pom(pom_path)) == EXPECTED_DIRECT_DEPENDENCIES[pom_path]


@pytest.mark.parametrize("pom_path", EXPECTED_REPOSITORIES)
def test_patched_java_projects_lock_pinned_repositories(pom_path: str) -> None:
    assert _repository_contract(_pom(pom_path)) == EXPECTED_REPOSITORIES[pom_path]


@pytest.mark.parametrize("pom_path", EXPECTED_REPOSITORIES)
def test_patched_java_projects_use_only_secure_repository_urls(pom_path: str) -> None:
    for _, repository_url in _repository_contract(_pom(pom_path)):
        assert repository_url.startswith("https://")


@pytest.mark.parametrize(
    "pom_path",
    ["docker/java/snpeff-pom.xml", "docker/java/snpsift-pom.xml"],
)
def test_patched_java_projects_pin_secured_runtime_families(pom_path: str) -> None:
    pom = _pom(pom_path)
    properties = pom.find("m:properties", MAVEN_NAMESPACE)
    assert properties is not None
    for property_name, expected_version in SECURED_PROPERTIES.items():
        assert (
            properties.findtext(f"m:{property_name}", namespaces=MAVEN_NAMESPACE)
            == expected_version
        )

    managed_versions = _dependency_versions(
        pom,
        "m:dependencyManagement/m:dependencies/m:dependency",
    )
    for coordinate, property_reference in MANAGED_RUNTIME_VERSIONS.items():
        assert managed_versions[coordinate] == [property_reference]

    direct_versions = _dependency_versions(pom, "m:dependencies/m:dependency")
    for coordinate, property_reference in MANAGED_RUNTIME_VERSIONS.items():
        assert direct_versions[coordinate] == [property_reference]


def test_snpeff_excludes_banned_commons_lang_from_biojava_structure() -> None:
    pom = _pom("docker/java/snpeff-pom.xml")
    dependencies = {
        _coordinate(dependency): dependency
        for dependency in _dependency_elements(pom, "m:dependencies/m:dependency")
    }
    biojava_structure = dependencies[("org.biojava", "biojava-structure")]
    exclusions = {
        _coordinate(exclusion)
        for exclusion in biojava_structure.findall("m:exclusions/m:exclusion", MAVEN_NAMESPACE)
    }
    assert ("commons-lang", "commons-lang") in exclusions


def test_runtime_dependency_assertion_accepts_only_the_fixed_runtime_set(
    tmp_path: Path,
) -> None:
    result = _run_dependency_assertion(
        tmp_path,
        [*REQUIRED_RUNTIME_DEPENDENCIES, "org.example:unrelated:jar:1.0:runtime"],
    )
    assert result.returncode == 0, result.stderr


def test_runtime_dependency_assertion_pins_dependency_plugin() -> None:
    script = _text("docker/java/assert-runtime-dependencies.sh")
    assert "org.apache.maven.plugins:maven-dependency-plugin:3.7.0:list" in script


def test_runtime_dependency_assertion_rejects_commons_lang(tmp_path: Path) -> None:
    result = _run_dependency_assertion(
        tmp_path,
        [*REQUIRED_RUNTIME_DEPENDENCIES, "commons-lang:commons-lang:jar:2.6:runtime"],
    )
    assert result.returncode != 0
    assert "commons-lang:commons-lang" in result.stderr


def test_runtime_dependency_assertion_rejects_other_secured_versions(tmp_path: Path) -> None:
    result = _run_dependency_assertion(
        tmp_path,
        [
            *REQUIRED_RUNTIME_DEPENDENCIES,
            "com.fasterxml.jackson.core:jackson-core:jar:2.22.0:runtime",
        ],
    )
    assert result.returncode != 0
    assert "jackson-core" in result.stderr


def test_runtime_dependency_assertion_requires_every_fixed_coordinate(tmp_path: Path) -> None:
    result = _run_dependency_assertion(tmp_path, REQUIRED_RUNTIME_DEPENDENCIES[:-1])
    assert result.returncode != 0
    assert "log4j-slf4j-impl" in result.stderr


@pytest.mark.parametrize(
    "malformed_dependency_line",
    [
        "[INFO] com.fasterxml.jackson.core:jackson-core:jar:2.22.1:runtime",
        "com.fasterxml.jackson.core:jackson-core:jar:2.22.1:runtime unexpected-trailer",
    ],
)
def test_runtime_dependency_assertion_rejects_malformed_coordinate_lines(
    tmp_path: Path,
    malformed_dependency_line: str,
) -> None:
    result = _run_dependency_assertion_output(
        tmp_path,
        [
            malformed_dependency_line,
            *[
                f"   {dependency} -- module fake.module [auto]"
                for dependency in REQUIRED_RUNTIME_DEPENDENCIES[1:]
            ],
        ],
    )
    assert result.returncode != 0
    assert "Required runtime dependency missing: com.fasterxml.jackson.core:jackson-core" in (
        result.stderr
    )


def test_container_image_contract_resolves_one_immutable_image_id() -> None:
    script = _text("scripts/test_container_image.sh")
    assert "image_id=$(docker image inspect --format '{{.Id}}' \"$image_ref\")" in script
    assert "[[ $image_id =~ ^sha256:[0-9a-f]{64}$ ]]" in script
    assert script.count('"$image_ref"') == 2
    assert script.count('"$image_id"') >= 10


def test_container_image_contract_keeps_golden_oracles_outside_writable_mount() -> None:
    script = _text("scripts/test_container_image.sh")
    for tool in ("snpeff", "snpsift"):
        assert f'"$fixture_source/{tool}/expected.vcf"' in script
        assert f'"$fixture_copy/{tool}/expected.vcf"' not in script
    assert 'cp -R "$fixture_source/." "$fixture_copy/"' not in script
    assert 'chmod -R u=rwX,go=rX "$fixture_copy"' in script
    assert 'chmod o+rwx "$fixture_copy/snpeff/data/testGenome"' in script


def test_container_image_contract_scans_broadly_for_build_executables() -> None:
    script = _text("scripts/test_container_image.sh")
    assert "find /opt /usr /bin /sbin /app -xdev \\" in script
    assert '\\( -type f -o -type l \\) -executable -printf "%p\\n"' in script
    assert (
        'compiler_name_pattern="(^|.*-)(gcc|g\\+\\+|cc|c\\+\\+|clang|clang\\+\\+|'
        'javac|mvn|mvnDebug|maven)(-?[0-9]+([.][0-9]+)*)?$"' in script
    )


def test_container_image_contract_checks_default_snpeff_data_storage() -> None:
    script = _text("scripts/test_container_image.sh")
    assert "checking default SnpEff database storage" in script
    assert 'grep -Fx "data.dir = /data/snpeff/"' in script
    assert "test -w /data/snpeff" in script
    assert "touch /data/snpeff/container-contract-write-test" in script
    assert "test ! -w /opt/conda/share/snpeff-5.2-3" in script
    assert "test ! -e /opt/conda/etc/conda/test-files" in script


def test_compose_and_installation_persist_the_default_snpeff_data_directory() -> None:
    compose = _text("docker-compose.yml")
    installation = _text("docs/source/installation.md")

    assert "snpeff_data:/data/snpeff" in compose
    assert "snpeff_data:/data/snpeff" in installation
    assert "/snpeff_data" not in compose
    assert "/snpeff_data" not in installation
    for documentation in (compose, installation):
        assert "first download" in documentation.lower()
        assert "fully populated" in documentation.lower()
        assert "snpeff_data:/data/snpeff:ro" in documentation


def test_scoring_integration_requires_every_invoked_bioinformatics_tool() -> None:
    scoring_integration = _text("tests/test_scoring_integration.py")
    required_executables = scoring_integration.split("_REQUIRED_EXECUTABLES = (", maxsplit=1)[
        1
    ].split(")", maxsplit=1)[0]
    assert '"sortBed"' in required_executables


def test_docker_workflow_preserves_scan_job_and_scans_pull_requests() -> None:
    workflow = _text(".github/workflows/docker.yml")
    assert "\n  build:\n" not in workflow
    scan_job = workflow.split("\n  scan:\n", maxsplit=1)[1].split("\n  sign:\n", maxsplit=1)[0]
    scan_header = scan_job.split("    steps:\n", maxsplit=1)[0]
    assert "needs: build" not in scan_job
    assert "if: github.event_name != 'pull_request'" not in scan_header
    assert "LOCAL_IMAGE: variantcentrifuge:ci-${{ github.sha }}" in workflow


def test_docker_workflow_pull_request_paths_cover_container_inputs() -> None:
    workflow = _text(".github/workflows/docker.yml")
    pull_request = workflow.split("  pull_request:\n", maxsplit=1)[1].split("\nenv:\n", maxsplit=1)[
        0
    ]
    configured_paths = {
        line.removeprefix('      - "').removesuffix('"')
        for line in pull_request.splitlines()
        if line.startswith('      - "')
    }
    assert configured_paths == {
        "Dockerfile",
        "docker/**",
        "docker-compose.yml",
        ".dockerignore",
        "conda/environment-docker.yml",
        "scripts/test_container_image.sh",
        "scripts/summarize_trivy.py",
        "tests/fixtures/container/**",
        "variantcentrifuge/**",
        "pyproject.toml",
        "setup.py",
        "README.md",
        "LICENSE",
        "scoring/**",
        "stats_configs/**",
        ".github/workflows/docker.yml",
    }


def test_docker_workflow_builds_and_tests_one_local_production_image() -> None:
    workflow = _text(".github/workflows/docker.yml")
    build = _workflow_step(workflow, "Build production image", "Test production image")
    assert workflow.count("uses: docker/build-push-action@v7") == 1
    assert "load: true" in build
    assert "push: false" in build
    assert "${{ env.LOCAL_IMAGE }}" in build
    assert "${{ steps.meta.outputs.tags }}" in build
    smoke = _workflow_step(
        workflow,
        "Test production image",
        "Generate complete vulnerability report",
    )
    assert 'run: bash scripts/test_container_image.sh "${LOCAL_IMAGE}"' in smoke
    assert workflow.index("- name: Test production image") < workflow.index(
        "- name: Generate complete vulnerability report"
    )
    assert workflow.index("- name: Test production image") < workflow.index(
        "- name: Generate actionable SARIF report"
    )


def test_docker_workflow_uses_node24_action_releases() -> None:
    workflow = _text(".github/workflows/docker.yml")
    uses = [
        line.strip().removeprefix("uses: ")
        for line in workflow.splitlines()
        if line.strip().startswith("uses: ")
    ]
    for action in (
        "actions/checkout@v6",
        "docker/setup-buildx-action@v4",
        "docker/metadata-action@v6",
        "docker/build-push-action@v7",
    ):
        family = action.rsplit("@", maxsplit=1)[0]
        assert [entry for entry in uses if entry.startswith(f"{family}@")] == [action]


def test_docker_workflow_records_the_built_image_identity_before_testing() -> None:
    workflow = _text(".github/workflows/docker.yml")
    assert (
        workflow.index("- name: Build production image")
        < workflow.index("- name: Record production image ID")
        < workflow.index("- name: Test production image")
    )
    record = _workflow_step(workflow, "Record production image ID", "Test production image")
    assert "id: production-image" in record
    assert "docker image inspect --format '{{.Id}}' \"${LOCAL_IMAGE}\"" in record
    assert "[[ $image_id =~ ^sha256:[0-9a-f]{64}$ ]]" in record
    assert 'echo "id=${image_id}" >> "${GITHUB_OUTPUT}"' in record


def test_docker_workflow_retains_complete_all_severity_audit() -> None:
    workflow = _text(".github/workflows/docker.yml")
    complete = _workflow_step(
        workflow,
        "Generate complete vulnerability report",
        "Upload complete vulnerability report",
    )
    for setting in (
        "uses: aquasecurity/trivy-action@v0.36.0",
        "image-ref: ${{ env.LOCAL_IMAGE }}",
        "format: json",
        "output: trivy-complete.json",
        "scanners: vuln",
        "vuln-type: os,library",
        "severity: UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL",
        'exit-code: "0"',
        "version: v0.70.0",
    ):
        assert setting in complete
    assert "ignore-unfixed" not in complete

    upload = _workflow_step(
        workflow,
        "Upload complete vulnerability report",
        "Generate actionable SARIF report",
    )
    assert "if: always()" in upload
    assert "uses: actions/upload-artifact@v7.0.1" in upload
    assert "path: trivy-complete.json" in upload
    assert "retention-days: 90" in upload
    assert "if-no-files-found: error" in upload


def test_docker_workflow_uploads_actionable_all_severity_sarif() -> None:
    workflow = _text(".github/workflows/docker.yml")
    sarif = _workflow_step(
        workflow,
        "Generate actionable SARIF report",
        "Upload actionable SARIF report",
    )
    for setting in (
        "uses: aquasecurity/trivy-action@v0.36.0",
        "image-ref: ${{ env.LOCAL_IMAGE }}",
        "format: sarif",
        "output: trivy-actionable.sarif",
        "scanners: vuln",
        "vuln-type: os,library",
        "severity: UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL",
        "ignore-unfixed: true",
        'exit-code: "0"',
        "version: v0.70.0",
    ):
        assert setting in sarif

    upload = _workflow_step(
        workflow,
        "Upload actionable SARIF report",
        "Generate actionable vulnerability report",
    )
    assert "if: always()" in upload
    assert "uses: github/codeql-action/upload-sarif@v4.37.1" in upload
    assert "sarif_file: trivy-actionable.sarif" in upload


def test_docker_workflow_summarizes_actionable_json_before_authoritative_gate() -> None:
    workflow = _text(".github/workflows/docker.yml")
    actionable = _workflow_step(
        workflow,
        "Generate actionable vulnerability report",
        "Summarize vulnerability reports",
    )
    for setting in (
        "format: json",
        "output: trivy-actionable.json",
        "ignore-unfixed: true",
        'exit-code: "0"',
    ):
        assert setting in actionable
    summary = _workflow_step(
        workflow,
        "Summarize vulnerability reports",
        "Enforce zero vendor-fixed vulnerabilities",
    )
    assert "python scripts/summarize_trivy.py trivy-complete.json trivy-actionable.json" in summary
    assert '| tee -a "${GITHUB_STEP_SUMMARY}"' in summary


def test_docker_workflow_final_gate_blocks_every_fixed_vulnerability() -> None:
    workflow = _text(".github/workflows/docker.yml")
    gate = _workflow_step(
        workflow,
        "Enforce zero vendor-fixed vulnerabilities",
        "Log in to GHCR for publication",
    )
    for setting in (
        "uses: aquasecurity/trivy-action@v0.36.0",
        "image-ref: ${{ env.LOCAL_IMAGE }}",
        "format: table",
        "scanners: vuln",
        "vuln-type: os,library",
        "severity: UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL",
        "ignore-unfixed: true",
        'exit-code: "1"',
        "version: v0.70.0",
    ):
        assert setting in gate


def test_docker_workflow_has_no_fail_open_security_steps() -> None:
    workflow = _text(".github/workflows/docker.yml")
    assert "continue-on-error" not in workflow

    smoke = _workflow_step(
        workflow,
        "Test production image",
        "Generate complete vulnerability report",
    )
    assert "if:" not in smoke
    assert "|| true" not in smoke

    gate = _workflow_step(
        workflow,
        "Enforce zero vendor-fixed vulnerabilities",
        "Log in to GHCR for publication",
    )
    assert "if:" not in gate
    assert "|| true" not in gate
    assert 'exit-code: "1"' in gate


def test_every_docker_workflow_trivy_call_has_the_same_scanner_scope() -> None:
    workflow = _text(".github/workflows/docker.yml")
    trivy_steps = workflow.split("uses: aquasecurity/trivy-action@v0.36.0")[1:]
    assert len(trivy_steps) == 4
    for step in trivy_steps:
        step = step.split("\n      - name:", maxsplit=1)[0]
        assert "scanners: vuln" in step
        assert "vuln-type: os,library" in step
        assert "severity: UNKNOWN,LOW,MEDIUM,HIGH,CRITICAL" in step
        assert "version: v0.70.0" in step
        assert "trivy-version:" not in step


def test_docker_workflow_publishes_only_after_the_gate_without_rebuilding() -> None:
    workflow = _text(".github/workflows/docker.yml")
    gate_index = workflow.index("- name: Enforce zero vendor-fixed vulnerabilities")
    login_index = workflow.index("- name: Log in to GHCR for publication")
    publish_index = workflow.index("- name: Push tested production image")
    assert gate_index < login_index < publish_index
    assert workflow.count("uses: docker/build-push-action@v7") == 1

    login = _workflow_step(
        workflow, "Log in to GHCR for publication", "Push tested production image"
    )
    assert "if: github.event_name != 'pull_request'" in login
    assert "for attempt in 1 2 3" in login
    publish = _workflow_step(workflow, "Push tested production image", "Install cosign")
    assert "id: publish" in publish
    assert "if: github.event_name != 'pull_request'" in publish
    assert "METADATA_TAGS: ${{ steps.meta.outputs.tags }}" in publish
    assert "PRODUCTION_IMAGE_ID: ${{ steps.production-image.outputs.id }}" in publish
    assert "local_image_id=$(docker image inspect --format '{{.Id}}' \"${LOCAL_IMAGE}\")" in publish
    assert '[[ $local_image_id == "$PRODUCTION_IMAGE_ID" ]]' in publish
    assert "tag_image_id=$(docker image inspect --format '{{.Id}}' \"$tag\")" in publish
    assert '[[ $tag_image_id == "$PRODUCTION_IMAGE_ID" ]]' in publish
    assert publish.index('[[ $tag_image_id == "$PRODUCTION_IMAGE_ID" ]]') < publish.index(
        'docker push "$tag"'
    )
    assert 'docker push "$tag"' in publish
    assert 'repository="${REGISTRY}/${IMAGE_NAME}"' in publish
    assert "${repository}@sha256:" in publish
    assert "[[ ${#matching_digests[@]} -eq 1 ]]" in publish
    assert "^sha256:[0-9a-f]{64}$" in publish
    assert 'echo "digest=${digest}" >> "${GITHUB_OUTPUT}"' in publish


def test_docker_workflow_exports_scan_outputs_and_signs_only_published_images() -> None:
    workflow = _text(".github/workflows/docker.yml")
    scan = workflow.split("\n  scan:\n", maxsplit=1)[1].split("\n  sign:\n", maxsplit=1)[0]
    assert "image-digest: ${{ steps.publish.outputs.digest }}" in scan
    assert "image-tags: ${{ steps.meta.outputs.tags }}" in scan

    sign = workflow.split("\n  sign:\n", maxsplit=1)[1]
    assert "needs: scan" in sign
    assert "needs: build" not in sign
    assert "if: github.event_name != 'pull_request'" in sign
    assert "DIGEST: ${{ needs.scan.outputs.image-digest }}" in sign
    assert 'cosign sign --yes "${{ env.REGISTRY }}/${{ env.IMAGE_NAME }}@${DIGEST}"' in sign


def test_docker_workflow_grants_oidc_only_to_the_signing_job() -> None:
    workflow = _text(".github/workflows/docker.yml")
    workflow_permissions = workflow.split("\njobs:\n", maxsplit=1)[0]
    scan = workflow.split("\n  scan:\n", maxsplit=1)[1].split("\n  sign:\n", maxsplit=1)[0]
    sign = workflow.split("\n  sign:\n", maxsplit=1)[1]

    assert "id-token:" not in workflow_permissions
    assert "id-token:" not in scan
    assert "permissions:\n      packages: write\n      id-token: write" in sign
    assert "security-events: write" not in sign
    assert "contents: write" not in sign
