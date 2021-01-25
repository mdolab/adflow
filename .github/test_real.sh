#!/bin/bash
set -e
cd input_files && ./get-input-files.sh
testflo -v . -n 1