from __future__ import print_function
import sys

if __name__ == "__main__":
    import argparse

    # Add command line arguments
    parser = argparse.ArgumentParser(description="Test import of the Python shared library file")
    parser.add_argument("-n", "--name", type=str, help="Module name")

    # Parse the command line arguments
    args = parser.parse_args()
    name = args.name
    print("Testing if module %s can be imported..." % name)
    import_cmd = "import %s" % name
    try:
        exec(import_cmd)
    except:
        print("Error: %s.so was not imported correctly" % name)
        sys.exit(1)
    print("Module %s was successfully imported." % name)
