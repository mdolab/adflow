#! /usr/bin/env python
"""
autoEdit - A Python tool to automatically edit a set of files
           according to the specified user rules:

           1) Discard module files

           2) Rename module names in subroutines

Written by Andre C. Marta          Last updated: Apr 6, 2007
"""

# Import modules
import os, sys
import string
import re
# Specify file extension
EXT = '_b.f90' 

DIR_ORI = sys.argv[1]
DIR_MOD = sys.argv[2]

# Specifiy the list of LINE ID's to find, what to replace and with what
patt_modules = re.compile(r'(\s*use\s*\w*)(_b)\s*')

del_patterns = [re.compile(r'(\s*call pushreal8)'), 
                re.compile(r'(\s*call popreal8)'),
                re.compile(r'(\s*call pushinteger4)'),
                re.compile(r'(\s*call popinteger4)')]
patt_pushcontrol1b = re.compile(r'(\s*call pushcontrol1b\()(.*)\)')
patt_popcontrol1b = re.compile(r'(\s*call popcontrol1b\()(.*)\)')
patt_comment = re.compile(r'\s*!.*')
print "Directory of input source files  :", DIR_ORI
print "Directory of output source files :", DIR_MOD

for f in os.listdir(DIR_ORI):
    if f.endswith(EXT):
        # open original file in read mode
        file_object_ori = open(DIR_ORI + '/' + f,'r')
        print "\nParsing input file", file_object_ori.name

        # read to whole file to string and reposition the pointer
        # at the first byte for future reading
        all_src = file_object_ori.read()
        file_object_ori.seek(0)

        # open modified file in write mode
        file_object_mod = open(os.path.join(DIR_MOD,f), 'w')

        # read the original file, line-by-line
        addedModule = False

        for line in file_object_ori:
            # Just deal with lower case string
            line = line.lower()

            
            # Replace modules
            m = patt_modules.match(line)
            if m:
                line = m.group(1) + '\n'

                if not addedModule:
                    line = '  use myPushPopLib\n'+line
                    addedModule = True

            # Replace _cb on calls
            if '_cb' in line:
                line = line.replace('_cb', '')
                
            # Delete patterns
            for p in del_patterns:
                if p.match(line):
                    line = ''
                    
            # Push control 1b's
            m = patt_pushcontrol1b.match(line)
            if m:
                num = m.group(2)
                line = 'myIntPtr = myIntPtr + 1\n myIntStack(myIntPtr) = %s\n'%num

            # Pop control 1b's
            m = patt_popcontrol1b.match(line)
            if m:
                num = m.group(2)
                line = '%s = myIntStack(myIntPtr)\n myIntPtr = myIntPtr - 1\n'%num
            
            # write the modified line to new file
            file_object_mod.write(line)

        # close the files
        file_object_ori.close()
        file_object_mod.close()

        # success message
        print " Modified file saved", file_object_mod.name


