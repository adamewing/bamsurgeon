#!/bin/bash
#
# Install BAMSurgeon's runtime and test dependencies so that a Claude Code on
# the web session can actually run the tools and their tests.
#
# None of samtools/bwa/pysam ship in the base image, so without this every
# test in test/ fails at the first subprocess call.

set -euo pipefail

# Local checkouts are assumed to already have a working environment; running
# apt-get on someone's laptop is not our business.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

bash "$PROJECT_DIR/scripts/setup_env.sh"
