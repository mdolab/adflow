#!/usr/bin/python

"""
Text in the adflow.pyf cannot contain ANY upper case characters.
This script is intended to report any lines that contain UPPER
case characters. Text in comments and preprocessor tags are ignored
"""

import sys

def error(i,line):
    print("ERROR : No code related characters in adflow.pyf cannot contain upper case characters")
    print("ERROR : Uppercase found in adflow.pyf line {0:d}. String found is: {1:s} ".format(i,line))
    sys.exit(1)


# Define charachter set
ignoreChars = set("!#")

f = open("../f2py/adflow.pyf","r")
for i, line in enumerate(f):
    line = line.rstrip()
    if not line:
        # Line is empty
        continue
    elif line.islower():
        # Check if the line contains all lower case.
        continue
    elif any((c in ignoreChars) for c in line):
        # Check for comments and preprocessor tags

        # Comments are special since they might be trailing actual code, need to check further
        if "!" in line:
            tmp = line.split("!")[0]  # Get the text left of comment
            tmp = tmp.replace(" ","") # Clean all white spaces since they will show as not lower case later

            if not tmp:
                # Check if the rest is an empty string after the manipulation
                continue

            if not tmp.islower():
                # Check if the rest contains upper case characters
                error(i,line)
        # No upper case characters was found
        continue

    else:
         error(i,line)

# No errors in file exit with no errors
sys.exit(0)


#
