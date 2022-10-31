#!/usr/bin/env python
# Standard Python modules
import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument(
    "name",
    type=str,
    help="Library name (example: libpackage.so). Note: This script must be run in the same dir as the library.",
)
args = parser.parse_args()

# Make sure the executing directory is always in the PATH before importing
sys.path.insert(0, os.getcwd())

# Only get the filename without the extension
name = os.path.splitext(args.name)[0]
print(f"Testing if module {name} can be imported...")

try:
    import_cmd = f"import {name}"
    exec(import_cmd)
except ImportError as e:
    print(f"Error: {e}")
    print(f"Error: library {args.name} was not imported correctly")
    sys.exit(1)

print(f"Module {name} was successfully imported")
