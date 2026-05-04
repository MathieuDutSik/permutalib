# Permutalib CI tests

Each subdirectory holds one CI test that exercises a binary built from `src/`.
The structure follows the same convention as the `polyhedral_common` and
`basic_common_cpp` repositories:

* `common.sh` is sourced by every `ProcessExamples.sh` and provides
  `ci_reset` (delete `CI_CONCLUSION`) and `ci_write_ok` (touch `CI_CONCLUSION`).
* Each test directory contains a `ProcessExamples.sh` script. On success it
  writes a `CI_CONCLUSION` file inside its own directory; the absence of that
  file marks the test as failed.
* Test inputs that are required by the binary are stored next to the script
  in the same directory.

To run all tests locally:

```
./run_ci_tests.sh
```

This builds the required binaries via `make -C src` and then executes each
`ProcessExamples.sh` in lexicographic order. The GitHub workflow
`.github/workflows/ci_run_tests.yml` runs the same script on a daily
schedule and on manual dispatch.
