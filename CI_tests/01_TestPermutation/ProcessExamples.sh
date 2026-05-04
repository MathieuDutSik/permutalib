#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestPermutation"

for n in 5 10 15; do
  echo "== TestPermutation n=$n =="
  "$SRC/TestPermutation" "$n"
done

ci_write_ok
echo "Normal case"
