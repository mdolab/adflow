# Tests

## Dependencies
Running tests require the following packages, `testflo`, which can be installed via `pip install .[testing]` when installing ADflow, or separately via `pip`:
- `testflo`
- `parameterized`

## Testing
### Real tests
To run all real tests, use
```
testflo
```
To run all real tests in specific test file, use
```
testflo <path/to/test-file.py>
```
To run a specific test in a specific test file, use
```
testflo <path/to/test-file.py>:<test-class>.<test-function>
```
for example,
```
testflo reg_tests/test_functionals.py:TestFunctionals_1_euler_scalar_jst_tut_wing.test_jac_vec_prod_bwd
```
If you need help figuring out the full signature for a specific test, you can run
```
testflo --dryrun
```
which will list all real tests, and each line can be used directly to run that sepcific test.
You can also run
```
testflo <path/to/test-file.py> --dryrun
```
to do the same but for a single file.
Note that you can run `testflo` from any directory as long as the tests are contained in one of its sub-directories.

### Complex tests
Since the user may not have complex ADflow compiled, these are disabled by default.
To run complex tests, use
```
testflo -m "cmplx_test*"
```
Similarly, running
```
testflo -m "cmplx_test*" --dryrun
```
will show the list of complex tests.

## Training
ADflow tests are based on reference JSON files containing reference values.
To train all tests, i.e. generate these reference JSON files, use
```
testflo -m "train*"
```
To train all tests in a specific test file, use
```
testflo <path/to/test-file.py> -m "train*"`
```
To train specific regression test in specific testfile, use
```
testflo <path/to/test-file.py>:<test-class>.train
```
for example,
```
testflo reg_tests/test_functionals.py:TestFunctionals_1_euler_scalar_jst_tut_wing.train
```
To list all training functions, use
```
testflo -m "train*" --dryrun
```

## Important options
- To enable screen output: use `-s` argument
- To change number of concurrent tests to 1 for running tests: use `-n 1`
- To print message when tests pass use `-v`
- To have the tests announced before they are run (e.g. to debug hanging tests) use `--pre_announce`
- To list the tests/training functions without running anything, use `--dryrun`.
