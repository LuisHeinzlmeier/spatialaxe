#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git diff --exit-code "${1:?provide a reference}" "${2:?provide a candidate}" -- . ':!AGENTS.md' ':!fork/'
