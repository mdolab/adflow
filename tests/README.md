# Dependencies
- baseclasses
- testflo


# Instructions and examples
- To train all regression tests with specific function signature: `testflo reg_tests -m "train_*"`
- To train all regression tests in specific testfile either works: `cd reg_tests; testflo test_DVGeometry.py -m "train_*"` or `testflo reg_tests/test_DVGeometry.py -m "train_*"` 
- To train specific regression test in specific testfile either works: `cd reg_tests; testflo test_DVGeometry.py -m "test_1` or `testflo reg_tests/test_DVGeometry.py -m "train_1"`

- To test all regression tests: `testflo reg_tests`
- To test all regression tests with specific function signature: `testflo reg_tests -m "test_*"`
- To test all regression tests in specific testfile either works: `cd reg_tests; testflo test_DVGeometry.py -m "test_*"` or `testflo reg_tests/test_DVGeometry.py -m "test_*"` 
- To test specific regression test in a specific testfile either works: `cd reg_tests; testflo test_DVGeometry.py -m "test_1` or `testflo reg_tests/test_DVGeometry.py -m "test_1"`

- To enable screen output: use `-s` argument
- To change number of cores to 1 for running tests: use `-n 1`

## Other ways to run
The tests are also compatiable with general ```unittesting``` python framework and can be call individually using, 

```python unit_tests/test_basics.py```


,or can be called from within the in the vscode testing framework 

