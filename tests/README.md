# How to run the tests

```python run_all_tests.py```

This acts as an alias for invoking testflo from the terminal


# Dependencies
- testflo
- baseclasses


# Calling using Testflo directly examples

## Testing
- To test all regression tests: `testflo reg_tests`
- To test all regression tests with specific function signature: `testflo reg_tests -m "test_*"`
- To test all regression tests in specific testfile either works: `cd reg_tests; testflo test_DVGeometry.py -m "test_*"` or `testflo reg_tests/test_DVGeometry.py -m "test_*"` 
- To test specific regression test in a specific testfile either works: `testflo reg_tests/test_DVGeometry.py -m "test_1" or `testflo reg_tests/test_euler_tutorial_wing.py:TestFunctionals.test_jac_vec_prod_bwd`

## Training
- To train all regression tests with specific function signature: `testflo reg_tests -m "train_*"`
- To train all regression tests in specific testfile either works: `cd reg_tests; testflo test_DVGeometry.py -m "train_*"` or `testflo reg_tests/test_DVGeometry.py -m "train_*"` 
- To train specific regression test in specific testfile either works: `cd reg_tests; testflo test_DVGeometry.py -m "test_1` or `testflo reg_tests/test_DVGeometry.py -m "train_1"`

## important options 
- To enable screen output: use `-s` argument
- To change number of cores to 1 for running tests: use `-n 1`
- To print message when tests pass use `-v`
