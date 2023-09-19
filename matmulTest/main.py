# --- Python 3.9 ---
"""
@File    :   main.py
@Time    :   2023/07/19
@Desc    :   Should be a simple test of the f2py wrapper
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import time

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
# This is the f2py wrapped object but you have to access the module within it too
import test_f2py

# ==============================================================================
#                         MAIN DRIVER
# ==============================================================================
if __name__ == "__main__":
    np.random.seed(10)
    # Call the f2py wrapped matmult function
    n = 32


    a = np.zeros((n, n),order="F")
    b = np.zeros((n, n),order="F")
    a[:] = np.random.rand(n, n)[:]
    b[:] = np.random.rand(n, n)[:]

    # # Call the f2py wrapped matmult function
    cF = test_f2py.f2py_api.my_matmult_api(a, b, n)
    # print(np.isfortran(cF))
    # print(cF)
    print(np.mean(a),np.mean(b),np.mean(cF),np.mean(a@b))
    # # print(f"Output from my_matmult:\n{c}")

    # tic = time.perf_counter()
    # cPY = a @ b
    # toc = time.perf_counter()
    # print(f"Time taken for python matrix-matrix: {toc-tic} seconds")
    # res = np.amax(np.abs(cF - cPY))
    # print(res)

    # data = np.loadtxt("testoutput.txt")
    # a_test = np.zeros((n,n))
    # rows = data[:,3].astype(int)-1
    # cols = data[:,2].astype(int)-1
    # a_test[rows,cols] = data[:,0]
    # b_test = np.zeros((n,n))
    # b_test[rows,cols] = data[:,1]   
    # print(np.amax(np.abs(a_test-a)))
    # print(np.amax(np.abs(b_test-b)))