#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
DATA=../../systematic_test/Groups_12perfect
ci_require_binary "$SRC/TestEnumerationElement"

if [ ! -f "$DATA" ]; then
  echo "Missing data file: $DATA" >&2
  exit 1
fi

echo "== TestEnumerationElement check (limit 1000) =="
"$SRC/TestEnumerationElement" "$DATA" 1000 check

ci_write_ok
echo "Normal case"
