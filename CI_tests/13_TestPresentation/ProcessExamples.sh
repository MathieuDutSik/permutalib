#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestPresentation"

# Builds the presentations of Presentation.h for the groups of CI test 12,
# from the deterministic and the randomized stabilizer chains. Checks that
# the relators hold, the index of each level by coset enumeration (which
# proves the order of the presented group) and, for the orders up to 50000,
# the order of the simplified presentation on the input generators by
# coset enumeration.
GROUPS_DIR=../12_TestStabChainRandomized
echo "== TestPresentation on the groups of CI test 12 =="
"$SRC/TestPresentation" "$GROUPS_DIR/ListGroups" "$GROUPS_DIR/Groups" 50000

ci_write_ok
echo "Normal case"
