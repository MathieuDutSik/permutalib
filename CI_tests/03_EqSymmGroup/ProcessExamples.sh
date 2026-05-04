#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/EqSymmGroup"

for n in 5 8 10; do
  echo "== EqSymmGroup n=$n =="
  "$SRC/EqSymmGroup" "$n"
done

ci_write_ok
echo "Normal case"
