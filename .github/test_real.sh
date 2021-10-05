#!/bin/bash
set -e
./input_files/get-input-files.sh
if [[ $AGENT_NAME == "Azure Pipelines"* ]]; then
    N_TEST=1
else
    N_TEST=2
fi

cd tests
testflo -v -n $N_TEST --coverage --coverpkg adflow
