#!/bin/bash
set -ev

make -f Makefile_CS PETSC_ARCH=complex-opt-\$COMPILERS-\$PETSCVERSION