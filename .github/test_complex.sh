#!/bin/bash
set -e
./input_files/get-input-files.sh
export PETSC_ARCH=$PETSC_ARCH_COMPLEX
if [[ $AGENT_NAME == "Azure Pipelines"* ]]; then
    N_TEST=1
else
    N_TEST=2
fi
testflo -v -m "cmplx_test*" -n $N_TEST --coverage --coverpkg adflow
