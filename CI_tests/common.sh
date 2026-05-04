#!/bin/bash
# Shared helpers for permutalib CI tests.
# Each ProcessExamples.sh script should source this file from its directory:
#     source ../common.sh
# and call ci_reset at the start and ci_write_ok at the end on success.
# The CI workflow checks for the presence of CI_CONCLUSION to decide whether
# the test passed.

CI_CONCLUSION_FILE="CI_CONCLUSION"

ci_reset() {
  rm -f "$CI_CONCLUSION_FILE"
}

ci_write_ok() {
  touch "$CI_CONCLUSION_FILE"
}

# Resolve the binary path against the src/ directory. Fail loudly if missing.
ci_require_binary() {
  local bin_path="$1"
  if [ ! -x "$bin_path" ]; then
    echo "Missing binary: $bin_path" >&2
    echo "Build it from src/ before running the CI tests" >&2
    exit 1
  fi
}
