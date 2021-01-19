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
EXT = '_b.f90'

DIR_ORI = sys.argv[1]
DIR_MOD = sys.argv[2]

# Specifiy the list of LINE ID's to find, what to replace and with what
patt_modules = re.compile(r'(\s*use\s*\w*)(_b)\s*')
patt_module = re.compile(r'\s*module\s\w*')
patt_module_start = re.compile('(\s*module\s)(\w*)(_b)\s*')
patt_module_end   = re.compile('(\s*end module\s)(\w*)(_b)\s*')
del_patterns = [re.compile(r'(\s*call pushreal8)'),
                re.compile(r'(\s*call popreal8)'),
                re.compile(r'(\s*call pushinteger4)'),
                re.compile(r'(\s*call popinteger4)')]
patt_pushcontrol1b = re.compile(r'(\s*call pushcontrol1b\()(.*)\)')
patt_popcontrol1b = re.compile(r'(\s*call popcontrol1b\()(.*)\)')
patt_subroutine = re.compile(r'\s*subroutine\s\w*')
patt_comment = re.compile(r'\s*!.*')
print("Directory of input source files  :", DIR_ORI)
print("Directory of output source files :", DIR_MOD)

useful_modules = ['bcroutines_b','turbbcroutines_b',
                  'utils_b', 'flowutils_b', 'walldistance_b', 'bcpointers_b',
                  'initializeflow_b', 'turbutils_b', 'sa_b', 'fluxes_b',
                  'solverutils_b', 'residuals_b', 'surfaceintegrations_b']

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

        # If we have a module, close the input and cycle to next file.
        if isModule and not hasSubroutine:
            file_object_ori.close()
            continue
        elif isModule and hasSubroutine:
            # We need to change the name of the module file:
            f = f.replace('_b','_fast_b')

        # open modified file in write mode
        file_object_mod = open(os.path.join(DIR_MOD,f), 'w')

        # read the original file, line-by-line
        addedModule = False

        # Go back to the beginning
        file_object_ori.seek(0)

        for line in file_object_ori:
            # Just deal with lower case string
            line = line.lower()

            # # Replace modules
            # m = patt_modules.match(line)
            # if m:
            #     if not addedModule:
            #         file_object_mod.write('  use myPushPopLib\n')
            #         addedModule = True

            # Replace _cb on calls
            if '_cb' in line:
                line = line.replace('_cb', '')

            # Replace _b modules with normal -- except for the useful
            # ones.
            m = patt_modules.match(line)
            if m:
                found = False
                for m in useful_modules:
                    if m in line:
                        found = True
                if found:
                    line = line.replace('_b', '_fast_b', 1)
                else:
                    line = line.replace('_b', '')



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

            # See if we need to modify the line
            m = patt_module_start.match(line)
            if m:
                line = 'module %s_fast_b\n'%m.group(2)

            m = patt_module_end.match(line)
            if m:
                line = 'end module %s_fast_b\n'%m.group(2)

            # Delete patterns
            if not f.lower() == 'bcroutines_b.f90':
                for p in del_patterns:
                    if p.match(line):
                        line = ''


            # write the modified line to new file
            file_object_mod.write(line)

        # close the files
        file_object_ori.close()
        file_object_mod.close()

        # success message
        print(" Modified file saved", file_object_mod.name)


