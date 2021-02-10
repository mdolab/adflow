#!/bin/bash
set -e
./input_files/get-input-files.sh
testflo -v . -n 1
