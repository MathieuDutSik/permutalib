#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestSymmetricGroup"

for n in 6 8 10; do
  echo "== TestSymmetricGroup n=$n =="
  "$SRC/TestSymmetricGroup" "$n"
done

ci_write_ok
echo "Normal case"
