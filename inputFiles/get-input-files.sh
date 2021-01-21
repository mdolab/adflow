#!/bin/bash
# this file will download the input files for ADflow regression tests
# and extract them to the right place.

DIR=$(dirname $0)
wget -O $DIR/input_files.tar.gz http://umich.edu/~mdolaboratory/repo_files/ADflow/adflow_input_files.tar.gz
tar -xzf $DIR/input_files.tar.gz -C $DIR/../
