#!/usr/bin/env bash

set -euo pipefail

if (( $# != 1 )); then
    printf 'usage: %s IMAGE_REF\n' "${0##*/}" >&2
    exit 64
fi

image_ref=$1

fail() {
    printf 'container contract failure: %s\n' "$*" >&2
    exit 1
}

if ! image_id=$(docker image inspect --format '{{.Id}}' "$image_ref"); then
    fail 'IMAGE_REF did not resolve to an image ID'
fi
[[ $image_id =~ ^sha256:[0-9a-f]{64}$ ]] || \
    fail "resolved image ID is empty or invalid: '$image_id'"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd -- "$script_dir/.." && pwd)
fixture_source=$repository_root/tests/fixtures/container
work_dir=$(mktemp -d -t vc-container-contract-XXXXXX)
trap 'rm -rf -- "$work_dir"' EXIT

assert_contains() {
    local output=$1
    local expected=$2
    local description=$3
    if ! grep -Fq -- "$expected" <<<"$output"; then
        fail "$description"
    fi
}

fixture_copy=$work_dir/fixtures
mkdir -p "$fixture_copy/snpeff/data/testGenome" "$fixture_copy/snpsift"
cp "$fixture_source/snpeff/snpEff.config" "$fixture_source/snpeff/input.vcf" \
    "$fixture_copy/snpeff/"
cp "$fixture_source/snpeff/data/testGenome/genes.gff" \
    "$fixture_source/snpeff/data/testGenome/sequences.fa" \
    "$fixture_copy/snpeff/data/testGenome/"
cp "$fixture_source/snpsift/input.vcf" \
    "$fixture_source/snpsift/assert_fail_closed.py" "$fixture_copy/snpsift/"
# Inputs are read-only in the disposable mount. Only the local SnpEff database
# directory is writable by the image's arbitrary nonroot UID. Goldens are never
# copied into or mounted in the container.
chmod -R u=rwX,go=rX "$fixture_copy"
chmod o+rwx "$fixture_copy/snpeff/data/testGenome"

printf '%s\n' 'checking variantcentrifuge entrypoint behavior'
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" variantcentrifuge --version
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" variantcentrifuge --help >"$work_dir/variantcentrifuge-help.txt"
test -s "$work_dir/variantcentrifuge-help.txt"

printf '%s\n' 'checking configured health executable'
docker run --rm --entrypoint /opt/conda/bin/variantcentrifuge \
    "$image_id" --version

printf '%s\n' 'checking nonroot runtime user'
runtime_uid=$(docker run --rm --entrypoint /usr/bin/id "$image_id" -u)
if [[ ! $runtime_uid =~ ^[0-9]+$ ]] || (( runtime_uid == 0 )); then
    fail "runtime UID must be numeric and nonzero, got '$runtime_uid'"
fi

printf '%s\n' 'checking native tool versions'
bcftools_version=$(docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" /opt/conda/bin/bcftools --version)
if [[ $bcftools_version != 'bcftools 1.21'* ]]; then
    fail 'bcftools version does not start with 1.21'
fi
bedtools_version=$(docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" /opt/conda/bin/bedtools --version)
assert_contains "$bedtools_version" '2.31.1' 'bedtools version does not contain 2.31.1'

printf '%s\n' 'checking Java tool versions'
snpeff_version=$(docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" /opt/conda/bin/snpEff -version 2>&1)
assert_contains "$snpeff_version" '5.2' 'snpEff version does not contain 5.2'

printf '%s\n' 'checking default SnpEff database storage'
docker run --rm --user 0:0 --entrypoint /bin/sh "$image_id" -c '
set -eu
config=/opt/conda/share/snpeff-5.2-3/snpEff.config
test "$(grep -c "^data[.]dir = " "$config")" -eq 1
grep -Fx "data.dir = /data/snpeff/" "$config"
'
set +e
snpeff_default_data_output=$(docker run --rm \
    --entrypoint /usr/local/bin/_entrypoint.sh "$image_id" \
    /opt/conda/bin/snpEff genes2bed -v -nodownload -ud 1000 GRCh37.75 2>&1)
snpeff_default_data_status=$?
set -e
if (( snpeff_default_data_status == 0 )); then
    fail 'SnpEff default data resolution unexpectedly found a database in an empty /data'
fi
assert_contains "$snpeff_default_data_output" \
    "/data/snpeff/GRCh37.75/snpEffectPredictor.bin" \
    'SnpEff did not resolve its default database beneath /data/snpeff'

set +e
snpsift_version=$(docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" /opt/conda/bin/SnpSift -version 2>&1)
snpsift_version_status=$?
set -e
if (( snpsift_version_status != 1 )); then
    fail "SnpSift -version must preserve upstream status 1, got $snpsift_version_status"
fi
assert_contains "$snpsift_version" "Error: Unknown command '-version'" \
    'SnpSift -version missing upstream unknown-command diagnostic'
assert_contains "$snpsift_version" 'SnpSift -version -version' \
    'SnpSift -version missing upstream command diagnostic'
assert_contains "$snpsift_version" 'SnpSift version 5.2' \
    'SnpSift version output does not contain 5.2'

printf '%s\n' 'checking exact Java manifests'
docker run --rm --entrypoint /opt/conda/bin/python "$image_id" -c '
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
'

printf '%s\n' 'checking SnpSift golden filter behavior'
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    --mount "type=bind,src=$fixture_copy/snpsift,dst=/fixtures" \
    -w /fixtures "$image_id" /opt/conda/bin/SnpSift \
    filter '( QUAL >= 20 )' input.vcf >"$work_dir/snpsift-actual.vcf"
if ! diff -u \
    <(sed '/^##SnpSiftCmd=/d' "$fixture_source/snpsift/expected.vcf") \
    <(sed '/^##SnpSiftCmd=/d' "$work_dir/snpsift-actual.vcf"); then
    fail 'SnpSift output differs from the stock 5.2 golden'
fi

printf '%s\n' 'checking SnpSift parser diagnostics fail closed'
fail_closed_output=$(docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    --mount "type=bind,src=$fixture_copy/snpsift,dst=/fixtures" \
    -w /fixtures "$image_id" /opt/conda/bin/python \
    /fixtures/assert_fail_closed.py)
if [[ $fail_closed_output != 'SnpSift parser diagnostics failed closed' ]]; then
    fail "unexpected fail-closed marker: '$fail_closed_output'"
fi

printf '%s\n' 'checking SnpEff local database and annotation behavior'
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    --mount "type=bind,src=$fixture_copy/snpeff,dst=/fixtures" \
    -w /fixtures "$image_id" /opt/conda/bin/snpEff \
    build -c snpEff.config -dataDir ./data -gff3 \
    -noCheckCds -noCheckProtein testGenome
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    --mount "type=bind,src=$fixture_copy/snpeff,dst=/fixtures" \
    -w /fixtures "$image_id" /opt/conda/bin/snpEff \
    -c snpEff.config -dataDir ./data -noStats testGenome input.vcf \
    >"$work_dir/snpeff-actual.vcf"
if ! diff -u \
    <(sed '/^##SnpEffCmd=/d' "$fixture_source/snpeff/expected.vcf") \
    <(sed '/^##SnpEffCmd=/d' "$work_dir/snpeff-actual.vcf"); then
    fail 'SnpEff output differs from the stock 5.2 golden'
fi

printf '%s\n' 'checking runtime filesystem inventory and immutability'
docker run --rm --user 0:0 --entrypoint /bin/sh "$image_id" -c '
set -eu
test ! -e /opt/conda/pkgs

jar_inventory=$(find /opt/conda -type f \( -name snpEff.jar -o -name SnpSift.jar \) -print | sort)
expected_jars=$(printf "%s\n%s" \
    /opt/conda/share/snpeff-5.2-3/snpEff.jar \
    /opt/conda/share/snpsift-5.2-0/SnpSift.jar)
test "$jar_inventory" = "$expected_jars"

java_path=$(readlink -f /opt/conda/bin/java)
jvm_bin=$(dirname "$java_path")
jvm_inventory=$(find "$jvm_bin" -mindepth 1 -maxdepth 1 -printf "%f\n" | sort)
expected_jvm_inventory=$(printf "%s\n%s" java keytool)
test "$jvm_inventory" = "$expected_jvm_inventory"

conda_jvm_inventory=$(
    for link in /opt/conda/bin/*; do
        test -L "$link" || continue
        link_target=$(readlink "$link")
        case "$link_target" in
            ../lib/jvm/bin/*|/opt/conda/lib/jvm/bin/*) printf "%s\n" "${link##*/}" ;;
        esac
    done | sort
)
test "$conda_jvm_inventory" = "$expected_jvm_inventory"

for forbidden_path in \
    /root/.m2 /build /out /tmp/src \
    /tmp/snpeff.tar.gz /tmp/snpsift.tar.gz \
    "$jvm_bin/../include" "$jvm_bin/../jmods" \
    "$jvm_bin/../man" "$jvm_bin/../demo" "$jvm_bin/../sample" \
    "$jvm_bin/../src.zip" "$jvm_bin/../lib/src.zip"; do
    test ! -e "$forbidden_path"
done

compiler_name_pattern="(^|.*-)(gcc|g\+\+|cc|c\+\+|clang|clang\+\+|javac|mvn|mvnDebug|maven)(-?[0-9]+([.][0-9]+)*)?$"
compiler_inventory=$(
    find /opt /usr /bin /sbin /app -xdev \
        \( -type f -o -type l \) -executable -printf "%p\n" |
    while IFS= read -r tool_path; do
        tool_name=${tool_path##*/}
        if printf "%s\n" "$tool_name" | grep -Eq "$compiler_name_pattern"; then
            printf "%s\n" "$tool_path"
        fi
    done
)
test -z "$compiler_inventory"
test -z "$(find /var/lib/apt/lists -mindepth 1 -print -quit)"
test -z "$(find /var/cache/apt/archives -type f -name "*.deb" -print -quit)"

invalid_runtime_path=$(find /opt/conda /app -xdev \
    \( ! -user root -o ! -group root -o \
       \( ! -type l -a -perm /022 \) \) -print -quit)
test -z "$invalid_runtime_path"
'
docker run --rm --entrypoint /bin/sh "$image_id" -c '
set -eu
test "$(id -u)" -ne 0
test -w /data
test -w /data/snpeff
test ! -w /opt/conda/share/snpeff-5.2-3
touch /data/container-contract-write-test
touch /data/snpeff/container-contract-write-test
'

printf '%s\n' 'checking Python dependency and native runtime operations'
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" /opt/conda/bin/python -m pip check
docker run --rm --entrypoint /usr/local/bin/_entrypoint.sh \
    "$image_id" /opt/conda/bin/python -c '
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

table = pa.table({"variant": ["1-4-G-T"], "score": [1.0]})
frame = table.to_pandas()
assert frame.to_dict(orient="records") == [{"variant": "1-4-G-T", "score": 1.0}]

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
assert ifault == 0
'

printf '%s\n' 'checking image configuration and layer shape'
config_user=$(docker image inspect --format '{{.Config.User}}' "$image_id")
[[ -n $config_user && $config_user != 0 && $config_user != root ]] || \
    fail "configured image user must be nonroot, got '$config_user'"
[[ $(docker image inspect --format '{{json .Config.Entrypoint}}' "$image_id") \
    == '["/usr/local/bin/_entrypoint.sh","variantcentrifuge"]' ]] || \
    fail 'configured ENTRYPOINT differs from the production contract'
[[ $(docker image inspect --format '{{json .Config.Cmd}}' "$image_id") \
    == '["--help"]' ]] || fail 'configured CMD differs from the production contract'
[[ $(docker image inspect --format '{{.Config.WorkingDir}}' "$image_id") == /data ]] || \
    fail 'configured workdir differs from /data'
[[ $(docker image inspect --format '{{json .Config.Healthcheck.Test}}' "$image_id") \
    == '["CMD","/opt/conda/bin/variantcentrifuge","--version"]' ]] || \
    fail 'configured health command differs from the production contract'
[[ $(docker image inspect --format '{{.Config.Healthcheck.Interval}}' "$image_id") \
    == 1m0s ]] || fail 'configured health interval differs from 1m0s'
[[ $(docker image inspect --format '{{.Config.Healthcheck.Timeout}}' "$image_id") \
    == 10s ]] || fail 'configured health timeout differs from 10s'
[[ $(docker image inspect --format '{{.Config.Healthcheck.StartPeriod}}' "$image_id") \
    == 5s ]] || fail 'configured health start period differs from 5s'
[[ $(docker image inspect --format '{{.Config.Healthcheck.Retries}}' "$image_id") \
    == 3 ]] || fail 'configured health retries differs from 3'

image_size=$(docker image inspect --format '{{.Size}}' "$image_id")
if [[ ! $image_size =~ ^[0-9]+$ ]] || (( image_size >= 2000000000 )); then
    fail "image size must be below 2,000,000,000 bytes, got '$image_size'"
fi

# The conda runtime copy is expected to be the sole layer above 500 MB. The
# threshold is well above the 75 MB base layer; a second layer above it catches
# an accidental duplicate conda copy-up without relying on an exact layer size.
large_layer_count=$(
    docker image history --no-trunc --format '{{.Size}}' "$image_id" |
        python3 -c '
import re
import sys

factors = {"B": 1, "kB": 1_000, "MB": 1_000_000, "GB": 1_000_000_000, "TB": 1_000_000_000_000}
count = 0
for raw_line in sys.stdin:
    size = raw_line.strip()
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([kMGT]?B)", size)
    if match is None:
        raise SystemExit(f"unrecognized Docker history size: {size!r}")
    byte_count = float(match.group(1)) * factors[match.group(2)]
    count += byte_count > 500_000_000
print(count)
'
)
if [[ $large_layer_count != 1 ]]; then
    fail "expected exactly one Docker history layer above 500 MB, got $large_layer_count"
fi

printf 'container contract passed: %s\n' "$image_ref"
