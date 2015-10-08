#!/usr/bin/python

"""
Text in the sumb.pyf cannot contain ANY upper case characters.
This script is intended to report any lines that contain UPPER 
case characters. Text in comments and preprocessor tags are ignored
"""

import sys

f = open("sumb.pyf","r")

ignoreChars = set("!#")

for i, line in enumerate(f):
    line = line.rstrip()
    if not line:
        # Line is empty
        continue
    elif any((c in ignoreChars) for c in line):
        # Check for comments and preprocessor tags
        continue
    elif line.islower():
        # Check if the line contains all lower case. If not we adjust
        continue
    else:
         print "ERROR : Uppercase found in sumb.pyf line {0:d}. String found is: {1:s} ".format(i,line)
         print "ERROR : No code related characters in sumb.pyf cannot contain upper case characters"
         sys.exit(1)

# No errors in file exit with no errors
sys.exit(0)
        

