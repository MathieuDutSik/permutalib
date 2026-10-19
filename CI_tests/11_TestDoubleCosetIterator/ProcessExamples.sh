#!/bin/bash
set -euo pipefail

source ../common.sh
ci_reset

SRC=../../src
ci_require_binary "$SRC/TestDoubleCosetIterator"

# Checking method passed to TestDoubleCosetIterator. The comparison between the
# full list of double cosets and the one obtained from the iterators is always
# done. On top of that, METH_CHECK selects how the double cosets themselves are
# checked: "exhaustive" (the default, it enumerates the groups and is affordable
# on those entries), "fast_check_sizes", "fast_check_intersection" or "nothing"
# for only doing the comparison of the two lists.
METH_CHECK=${METH_CHECK:-exhaustive}

LOG_DIR=$(mktemp -d)
trap 'rm -rf "$LOG_DIR"' EXIT

echo "== TestDoubleCosetIterator on RawDC (meth_check=$METH_CHECK) =="
n_file=0
n_fail=0
list_fail=()
for eFile in RawDC/RawDoubleCoset_*; do
  n_file=$((n_file + 1))
  eLog="$LOG_DIR/$(basename "$eFile").log"
  if "$SRC/TestDoubleCosetIterator" "$eFile" "$METH_CHECK" > /dev/null 2> "$eLog"; then
    rm -f "$eLog"
  else
    n_fail=$((n_fail + 1))
    list_fail+=("$eFile")
    echo "FAIL $eFile"
    tail -n 20 "$eLog"
  fi
done

echo "n_file=$n_file n_fail=$n_fail"
if [ "$n_fail" -ne 0 ]; then
  echo "Failing entries:" >&2
  for eFile in "${list_fail[@]}"; do
    echo "  $eFile" >&2
  done
  exit 1
fi

ci_write_ok
echo "Normal case"
