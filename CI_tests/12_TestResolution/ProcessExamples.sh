#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestResolution"

# Computes free ZG-resolutions for a list of small groups and checks
# d o d = 0, the contracting homotopy identity and the integral homology
# against the values stored in GroupsHomology.
echo "== TestResolution check =="
"$SRC/TestResolution" GroupsHomology check

ci_write_ok
echo "Normal case"
