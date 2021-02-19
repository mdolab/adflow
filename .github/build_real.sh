#!/bin/bash
set -e
cp $CONFIG_FILE config/config.mk
make
pip install .
