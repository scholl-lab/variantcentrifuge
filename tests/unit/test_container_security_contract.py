import os
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
    assert assembly.findtext(
        "m:executions/m:execution/m:phase", namespaces=MAVEN_NAMESPACE
    ) == "package"
    assert assembly.findtext(
        "m:executions/m:execution/m:goals/m:goal", namespaces=MAVEN_NAMESPACE
    ) == "single"
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


def test_snpeff_verify_runs_only_hermetic_upstream_and_compatibility_tests() -> None:
    pom = _pom("docker/java/snpeff-pom.xml")
    surefire = _plugin(pom, "maven-surefire-plugin")
    includes = {
        element.text
        for element in surefire.findall(
            "m:configuration/m:includes/m:include", MAVEN_NAMESPACE
        )
    }
    assert includes == {
        "org/snpeff/snpEffect/testCases/unity/*.java",
        "org/apache/commons/lang/ArrayUtilsTest.java",
    }
    excludes = {
        element.text
        for element in surefire.findall(
            "m:configuration/m:excludes/m:exclude", MAVEN_NAMESPACE
        )
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
        assert (
            dependencies[coordinate].findtext("m:scope", namespaces=MAVEN_NAMESPACE)
            == "test"
        )


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
        for element in compiler.findall(
            "m:configuration/m:excludes/m:exclude", MAVEN_NAMESPACE
        )
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
        for element in surefire.findall(
            "m:configuration/m:includes/m:include", MAVEN_NAMESPACE
        )
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
        for exclusion in biojava_structure.findall(
            "m:exclusions/m:exclusion", MAVEN_NAMESPACE
        )
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
