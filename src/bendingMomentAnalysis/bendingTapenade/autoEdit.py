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

EXT = '_b.f90' 

# Specify directory containing the original source files
# and the output directory for edited files
DIR_ORI = os.getcwd()+'/bendingTapenade'
DIR_MOD = DIR_ORI + '/../bendingOutput'
# Specify line identifier

# Specifiy the list of LINE ID's to find, what to replace and with what

LINE_ID = ['USE','USE']#,'CALL']
STR_OLD = ['_B' ,'_b'  ]# ,'_CB'    ]
STR_NEW = [''   ,''  ]# ,''       ]


STR_REPLACE_ALL = {'_CB':''}

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

        # go to modified directory
        os.chdir(DIR_MOD)
        # open modified file in write mode
        file_object_mod = open(f,'w')

        # read the original file, line-by-line
        nEdits = len(LINE_ID)
        for line in file_object_ori:

            # parse original line for relevante identifier
            # and replace the string
            line_mod = line.lstrip()
           
            for i in xrange(nEdits):
                if line_mod[0:len(LINE_ID[i])] == LINE_ID[i]:
                    line_mod = string.replace(line_mod, STR_OLD[i], STR_NEW[i])
                #end
            #end
            for key in STR_REPLACE_ALL:
                line_mod = string.replace(line_mod,key,STR_REPLACE_ALL[key])
            #end

            # write the modified line to new file
            file_object_mod.write('   '+ line_mod)

        # close the files
        file_object_ori.close()
        file_object_mod.close()
        # success message
        print " Modified file saved", file_object_mod.name




