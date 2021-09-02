#!/bin/bash
set -e
cp $CONFIG_FILE config/config.mk
make -f Makefile_CS PETSC_ARCH=$PETSC_ARCH_COMPLEX
pip install -e .
