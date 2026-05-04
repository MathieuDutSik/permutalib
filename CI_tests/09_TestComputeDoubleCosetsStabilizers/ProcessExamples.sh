#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestComputeDoubleCosetsStabilizers"

echo "== TestComputeDoubleCosetsStabilizers (G=S4, U=<(0 1 2 3)>, V=<(0 1)(2 3)>) =="
"$SRC/TestComputeDoubleCosetsStabilizers" G_U_V_S4

ci_write_ok
echo "Normal case"
