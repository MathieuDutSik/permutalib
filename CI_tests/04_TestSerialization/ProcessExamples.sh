#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
DATA=../../systematic_test/Groups_12perfect
ci_require_binary "$SRC/TestSerialization"

if [ ! -f "$DATA" ]; then
  echo "Missing data file: $DATA" >&2
  exit 1
fi

echo "== TestSerialization on Groups_12perfect =="
"$SRC/TestSerialization" "$DATA"

ci_write_ok
echo "Normal case"
