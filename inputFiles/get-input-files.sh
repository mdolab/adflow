#!/bin/bash
# this file will download the input files for ADflow regression tests
# and extract them to the right place.

wget -O adflow_input_files.tar.gz https://umich.box.com/shared/static/vnra7aaepyo3elndakfyjgany15823qy
tar -xf adflow_input_files.tar.gz -C ../
