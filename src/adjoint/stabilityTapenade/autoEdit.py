#! /usr/bin/env python
"""
autoEdit - A Python tool to automatically edit a set of files
           according to the specified user rules:

           1) Discard module files

           2) Rename module names in subroutines

Written by Andre C. Marta          Last updated: Apr 6, 2007

"""

# Import modules

import os
import string

# Specify file extension

EXT = '.f90'

# Specify directory containing the original source files
# and the output directory for edited files


#DIR_ORI = '/home/mader/UTIAS/pyMDO/pyHF/SUmbADjoint/trunk/src/adjoint/stabilityTapenade'
#DIR_MOD = '/home/mader/UTIAS/pyMDO/pyHF/SUmbADjoint/trunk/src/adjoint/stabilityOutput'
DIR_ORI = '/nfs/basalt/home/mdo/mader/svn/pyHF/SUmbADjoint/trunk/src/adjoint/stabilityTapenade'
DIR_MOD = '/nfs/basalt/home/mdo/mader/svn/pyHF/SUmbADjoint/trunk/src/adjoint/stabilityOutput'
# Specify line identifier

FILE_EXCL = 'MODULE'

LINE_ID = '  USE '

# Specify the string to be find and its replacement

STR_OLD = '_b'
STR_NEW = ''

# Some feedback

print "Directory of input source files  :", DIR_ORI
print "Directory of output source files :", DIR_MOD

# Clean the output directory

print "Cleaning output directory"
for f in os.listdir(DIR_MOD):
    if f.endswith(EXT):
        filePathComplete = os.path.join(DIR_MOD, f)
        os.remove(filePathComplete)

# Create a list with all the files to be parsed

for f in os.listdir(DIR_ORI):

    if f.endswith(EXT):

        # go to original directory
        os.chdir(DIR_ORI)
        # open original file in read mode
        file_object_ori = open(f,'r')
        print "\nParsing input file", file_object_ori.name

        # read to whole file to string and reposition the pointer
        # at the first byte for future reading
        all_src = file_object_ori.read()
        file_object_ori.seek(0)
        # check whether this is a fortran module file to discard
        if string.find(all_src, FILE_EXCL) > -1 :
            print " Module found -> discarded ! "
            continue

        # go to modified directory
        os.chdir(DIR_MOD)
        # open modified file in write mode
        file_object_mod = open(f,'w')

        # read the original file, line-by-line
        for line in file_object_ori:

            # parse original line for relevante identifier
            # and replace the string

            if line[0:len(LINE_ID)] == LINE_ID:
                line_mod = string.replace(line, STR_OLD, STR_NEW)
            else:
                line_mod = line

            # write the modified line to new file
            file_object_mod.write(line_mod)

        # close the files
        file_object_ori.close()
        file_object_mod.close()
        # success message
        print " Modified file saved", file_object_mod.name
