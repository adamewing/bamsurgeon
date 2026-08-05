#!/usr/bin/env bash
#
# Install everything BAMSurgeon needs to run and to run its test suite.
# Idempotent: safe to re-run, skips anything already present.

set -euo pipefail

SUDO=""
if [ "$(id -u)" -ne 0 ]; then
    SUDO="sudo"
fi

log() { printf '\n=== %s ===\n' "$1"; }

have() { command -v "$1" >/dev/null 2>&1; }

APT_PKGS=(samtools bcftools bwa minimap2 tabix build-essential zlib1g-dev)

missing=()
for pkg in "${APT_PKGS[@]}"; do
    if ! dpkg -s "$pkg" >/dev/null 2>&1; then
        missing+=("$pkg")
    fi
done

if [ ${#missing[@]} -gt 0 ]; then
    log "apt: installing ${missing[*]}"
    $SUDO apt-get update -qq
    DEBIAN_FRONTEND=noninteractive $SUDO apt-get install -y -qq "${missing[@]}"
else
    log "apt: all packages already present"
fi

log "pip: pysam, pytest"
pip3 install --quiet --disable-pip-version-check pysam pytest

# There is no standalone wgsim package (the `wgsim` apt name is dwgsim, a
# different tool with a different CLI). On Ubuntu the samtools package ships
# /usr/bin/wgsim, so the apt step above usually covers it; the source build is
# the fallback for distros where it does not.
if have wgsim; then
    log "wgsim: already installed ($(command -v wgsim))"
else
    log "wgsim: building from source"
    tmp=$(mktemp -d)
    git clone --depth 1 https://github.com/lh3/wgsim.git "$tmp/wgsim"
    gcc -g -O2 -Wall -o "$tmp/wgsim/wgsim" "$tmp/wgsim/wgsim.c" -lz -lm
    $SUDO install -m 0755 "$tmp/wgsim/wgsim" /usr/local/bin/wgsim
    rm -rf "$tmp"
fi

# picard and the JRE it needed are gone: BAM to FASTQ conversion is done in
# pysam by bamsurgeon.common.bamtofastq(). --picardjar and
# $BAMSURGEON_PICARD_JAR are still accepted and ignored.

log "verifying"
python3 "$(dirname "$0")/check_dependencies.py"

printf '\nSetup complete.\n\n'
