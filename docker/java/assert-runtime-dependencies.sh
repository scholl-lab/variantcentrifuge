#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf 'Usage: %s <project-directory>\n' "$0" >&2
    exit 2
fi

project_directory=$1
if [[ ! -f "$project_directory/pom.xml" ]]; then
    printf 'Maven project does not contain pom.xml: %s\n' "$project_directory" >&2
    exit 2
fi

resolved_dependencies=$(mktemp)
trap 'rm -f "$resolved_dependencies"' EXIT

mvn -B -f "$project_directory/pom.xml" \
    org.apache.maven.plugins:maven-dependency-plugin:3.7.0:list \
    -DincludeScope=runtime \
    -DoutputFile="$resolved_dependencies" \
    -DappendOutput=false

declare -A required_versions=(
    [com.fasterxml.jackson.core:jackson-core]=2.22.1
    [com.fasterxml.jackson.core:jackson-databind]=2.22.1
    [com.fasterxml.jackson.core:jackson-annotations]=2.22
    [com.google.code.gson:gson]=2.14.0
    [commons-io:commons-io]=2.22.0
    [org.apache.commons:commons-compress]=1.28.0
    [org.apache.logging.log4j:log4j-api]=2.25.4
    [org.apache.logging.log4j:log4j-core]=2.25.4
    [org.apache.logging.log4j:log4j-slf4j-impl]=2.25.4
)
declare -A found_required=()

coordinate_pattern='^[[:space:]]*([^:[:space:]]+):([^:[:space:]]+):([^:[:space:]]+):([^:[:space:]]+):([^:[:space:]]+)([[:space:]]+--[[:space:]]+module[[:space:]].*)?[[:space:]]*$'
while IFS= read -r line; do
    if [[ $line =~ $coordinate_pattern ]]; then
        group_id=${BASH_REMATCH[1]}
        artifact_id=${BASH_REMATCH[2]}
        version=${BASH_REMATCH[4]}
        coordinate="$group_id:$artifact_id"

        if [[ $coordinate == "commons-lang:commons-lang" ]]; then
            printf 'Banned runtime dependency resolved: %s\n' "$coordinate" >&2
            exit 1
        fi

        if [[ -n ${required_versions[$coordinate]+present} ]]; then
            expected_version=${required_versions[$coordinate]}
            if [[ $version != "$expected_version" ]]; then
                printf 'Unexpected runtime version for %s: expected %s, found %s\n' \
                    "$coordinate" "$expected_version" "$version" >&2
                exit 1
            fi
            found_required[$coordinate]=1
        fi
    fi
done < "$resolved_dependencies"

for coordinate in "${!required_versions[@]}"; do
    if [[ -z ${found_required[$coordinate]+present} ]]; then
        printf 'Required runtime dependency missing: %s:%s\n' \
            "$coordinate" "${required_versions[$coordinate]}" >&2
        exit 1
    fi
done
