#!/bin/bash

cd input_files && ./get-input-files.sh
export PETSC_ARCH=complex-opt-\$COMPILERS-\$PETSCVERSION
testflo -v . -m "cmplx_test*" -n 1