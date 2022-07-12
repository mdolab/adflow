#!/bin/bash
set -e

cd src/adjoint
make -f Makefile_tapenade ad_forward ad_reverse ad_reverse_fast
