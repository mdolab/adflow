#!/bin/bash
set -e

cd src/adjoint
python run_tapenade.py -force_all
