#!/bin/bash
#
# Install BAMSurgeon's runtime and test dependencies so that a Claude Code on
# the web session can actually run the tools and their tests.
#
# None of samtools/bwa/pysam/picard ship in the base image, so without
# this every test in test/ fails at the first subprocess call.

set -euo pipefail

# Local checkouts are assumed to already have a working environment; running
# apt-get on someone's laptop is not our business.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

bash "$PROJECT_DIR/scripts/setup_env.sh"

# remap_bam() shells out to picard for every aligner that round-trips through
# FASTQ; both addsnv and addindel read this variable as their --picardjar default.
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo 'export BAMSURGEON_PICARD_JAR=/opt/picard/picard.jar' >> "$CLAUDE_ENV_FILE"
fi
