#! /usr/bin/env python
"""
autoEdit - A Python tool to automatically edit a set of files
           according to the specified user rules:
G. Kenway
"""

# Import modules
import os, sys
import string
import re
# Specify file extension
EXT = '_d.f90'

DIR_ORI = sys.argv[1]
DIR_MOD = sys.argv[2]

# Specifiy the list of LINE ID's to find, what to replace and with what
patt_modules = re.compile(r'(\s*use\s*\w*)(_d)\s*')
patt_module = re.compile(r'\s*module\s\w*')
patt_module_start = re.compile('(\s*module\s)(\w*)(_d)\s*')
patt_module_end   = re.compile('(\s*end module\s)(\w*)(_d)\s*')
patt_subroutine = re.compile(r'\s*subroutine\s\w*')
patt_function = re.compile(r'\s*function\s\w*')

patt_subend = re.compile(r'\s*end\s*subroutine')
patt_funcend = re.compile(r'\s*end\s*function\n')

print("Directory of input source files  :", DIR_ORI)
print("Directory of output source files :", DIR_MOD)

useful_modules = ['bcroutines_d','turbbcroutines_d',
                  'utils_d', 'flowutils_d', 'walldistance_d', 'bcpointers_d',
                  'initializeflow_d', 'turbutils_d', 'sa_d', 'fluxes_d',
                  'solverutils_d', 'residuals_d', 'surfaceintegrations_d']
modSubToKeep = []

for f in os.listdir(DIR_ORI):
    if f.endswith(EXT):
        # open original file in read mode
        file_object_ori = open(os.path.join(DIR_ORI,f),'r')
        print("\nParsing input file", file_object_ori.name)

        # read to whole file to string and reposition the pointer
        # at the first byte for future reading
        all_src = file_object_ori.read()
        file_object_ori.seek(0)

        # First we want to dertmine if it is a 'useful' module or a
        # 'useless' module. A useful module is one that has
        # subroutines in it.
        isModule = False
        hasSubroutine = False
        for line in file_object_ori:
            line = line.lower()
            if patt_module.match(line):
                isModule = True
            if patt_subroutine.match(line):
                hasSubroutine = True
            if patt_function.match(line):
                hasSubroutine = True

        # If we have a module, close the input and cycle to next file.
        if isModule and not hasSubroutine:
            file_object_ori.close()
            continue
        elif isModule and hasSubroutine:
            f = f.replace('_d', '_d')

        # open modified file in write mode
        file_object_mod = open(os.path.join(DIR_MOD,f), 'w')

        # Go back to the beginning
        file_object_ori.seek(0)
        subActive = False

        for line in file_object_ori:
            # Just deal with lower case string
            line = line.lower()

            # Replace _cb on calls
            if '_cd' in line:
                line = line.replace('_cd', '')

            # Replace _d modules with normal -- except for the useful
            # ones.
            m = patt_modules.match(line)
            if m:
                found = False
                for m in useful_modules:
                    if m in line:
                        found = True
                if found:
                    line = line.replace('_d', '_d', 1)
                else:
                    line = line.replace('_d', '')

            # # See if we need to modify the line with changing the
            # # module names
            # m = patt_module_start.match(line)
            # if m:
            #     line = 'module %s_d2\n'%m.group(2)

            # m = patt_module_end.match(line)
            # if m:
            #     line = 'end module %s_d2\n'%m.group(2)

            file_object_mod.write(line)

        # close the files
        file_object_ori.close()
        file_object_mod.close()

        # success message
        print(" Modified file saved", file_object_mod.name)
