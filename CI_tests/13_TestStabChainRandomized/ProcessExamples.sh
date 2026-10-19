#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestStabChainRandomized"

# Builds the stabilizer chain of each group of ListGroups (built by GAP, with
# their orders) deterministically and by StabChainRandomized with no order,
# the exact order and a wrong order (twice the exact one), and
# checks the orders, the equality of the groups, the chains, the canonical
# images and the set stabilizers. The binary is built with
# PERMUTALIB_BLOCKING_SANITY_CHECK, so a chain accepted at the given order
# is also completed and checked.
echo "== TestStabChainRandomized on ListGroups =="
"$SRC/TestStabChainRandomized" ListGroups Groups 5

ci_write_ok
echo "Normal case"
