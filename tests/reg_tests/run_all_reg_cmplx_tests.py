import os 

print('======================== running regression tests =====================')
os.system('testflo -n 2 ./ -v -s -m cmplx_test_* ')