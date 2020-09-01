#!/bin/bash
# this file will download the input files for ADflow regression tests
# and extract them to the right place.

wget http://umich.edu/~mdolaboratory/repo_files/ADflow/adflow_input_files.tar.gz
tar -xzf adflow_input_files.tar.gz -C ../
