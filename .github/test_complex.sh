#!/bin/bash
set -e
./input_files/get-input-files.sh
testflo -v . -m "cmplx_test*" -n 1
