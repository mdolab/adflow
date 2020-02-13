from __future__ import print_function
# This file contains two functions to help regression testing. The
# first is used to format float values with a specified absolute and
# relative tolerance. This information is used by the second function
# when it takes in two such formatted strings and decides if they are
# sufficiently close to be considered equal.
import os
import sys
REG_FILES_MATCH = 0
REG_FILES_DO_NOT_MATCH = 1
REG_ERROR = -1

def _reg_str_comp(str1, str2):
    '''Compare the float values in str1 and str2 and determine if they
    are equal. Returns True if they are the "same", False if different'''

    aux1 = str1.split()
    aux2 = str2.split()

    if not aux1[0] == aux2[0] == '@value':
        # This line does not need to be compared
        return True

    # Extract required tolerances and values
    rel_tol = float(aux1[2])
    abs_tol = float(aux1[3])
    val1 = float(aux1[1])
    val2 = float(aux2[1])

    rel_err = 0
    if val2 != 0:
        rel_err = abs((val1-val2)/val2)
    else:
        rel_err = abs((val1-val2)/(val2 + 1e-16))

    abs_err = abs(val1-val2)

    if abs_err < abs_tol or rel_err < rel_tol:
        return True
    else:
        return False


def reg_file_comp(ref_file, comp_file):
    '''Compare the reference file 'ref_file' with 'comp_file'. The
    order of these two files matter. The ref_file MUST be given
    first. Only values specified by reg_write() are compared.  All
    other lines are ignored. Floating point values are compared based
    on rel_tol and abs_tol'''

    all_ref_lines = []
    ref_values = []
    comp_values = []
    try:
        f = open(ref_file, 'r')
    except IOError:
        print('File %s was not found. Cannot do comparison.'% ref_file)
        return REG_ERROR
    for line in f.readlines():
        all_ref_lines.append(line)
        if line[0:6] == '@value':
            ref_values.append(line)

    f.close()

    try:
        f = open(comp_file, 'r')
    except IOError:
        print('File %s was not found. Cannot do comparison.'% comp_file)
        return REG_ERROR

    for line in f.readlines():
        if line[0:6] == '@value':
            comp_values.append(line)

    f.close()

    # Copy the comp_file to compe_file.orig
    os.system('cp %s %s.orig'% (comp_file, comp_file))

    # We must check that we have the same number of @value's to compare:
    if len(ref_values) != len(comp_values):
        print('Error: number of @value lines in file not the same!')
        return REG_FILES_DO_NOT_MATCH

    # Open the (new) comp_file:
    f = open(comp_file,'w')

    # Loop over all the ref_lines, for value lines, do the
    # comparison. If comparison is ok, write the ref line, otherwise
    # write orig line.

    j = 0
    res = REG_FILES_MATCH
    for i in range(len(all_ref_lines)):
        line = all_ref_lines[i]
        if line[0:6] == '@value':
            if _reg_str_comp(line, comp_values[j]) is False:
                f.write(comp_values[j])
                res = REG_FILES_DO_NOT_MATCH
            else:
                f.write(line)

            j += 1
        else:
            f.write(line)

    f.close()

    return res
