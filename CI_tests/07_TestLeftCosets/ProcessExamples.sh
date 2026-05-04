#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestLeftCosets"

echo "== TestLeftCosets (H<G with G=S4, H=<(0 1)>) =="
"$SRC/TestLeftCosets" H_G_S4

ci_write_ok
echo "Normal case"
