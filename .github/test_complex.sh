#!/bin/bash
set -e
./input_files/get-input-files.sh
export PETSC_ARCH=$PETSC_ARCH_COMPLEX
testflo -v . -m "cmplx_test*" -n 1
