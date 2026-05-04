#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestComputeDoubleCosets"

echo "== TestComputeDoubleCosets exhaustive (G=S4, U=<(0 1 2 3)>, V=<(0 1)(2 3)>) =="
"$SRC/TestComputeDoubleCosets" G_U_V_S4 exhaustive

ci_write_ok
echo "Normal case"
