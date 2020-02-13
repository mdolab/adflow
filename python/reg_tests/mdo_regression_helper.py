from __future__ import print_function
# This file contains two functions to help regression testing. The
# first is used to format float values with a specified absolute and
# relative tolerance. This information is used by the second function
# when it takes in two such formatted strings and decides if they are
# sufficiently close to be considered equal.
import numpy, os
from mpi4py import MPI
import sys
REG_FILES_MATCH = 0
REG_FILES_DO_NOT_MATCH = 1
REG_ERROR = -1

def printHeader(testName):
    if MPI.COMM_WORLD.rank == 0:
        print('+' + '-'*78 + '+')
        print('| Test Name: ' + '%-66s'%testName + '|')
        print('+' + '-'*78 + '+')


def reg_write(values, rel_tol=1e-12, abs_tol=1e-12):
    '''Write values in special value format'''
    values = numpy.atleast_1d(values)
    values = values.flatten()
    for val in values:
        s = '@value %20.13e %g %g'% (val, rel_tol, abs_tol)
        print(s)

    return

def reg_par_write(values, rel_tol=1e-12, abs_tol=1e-12):
    """Write value(values) from parallel process in sorted order"""
    values = MPI.COMM_WORLD.gather(values)
    if MPI.COMM_WORLD.rank == 0:
        for i in xrange(len(values)):
            print ('Value(s) on processor: %d'%i)
            reg_write(values[i], rel_tol, abs_tol)

def reg_root_write(values, rel_tol=1e-12, abs_tol=1e-12):
    """Write values but only on the root proc"""
    if MPI.COMM_WORLD.rank == 0:
        reg_write(values, rel_tol, abs_tol)

def reg_par_write_sum(values, rel_tol=1e-12, abs_tol=1e-12):
    """Write the sum of sum of the values from all processors."""
    reducedSum = MPI.COMM_WORLD.reduce(numpy.sum(values))
    if MPI.COMM_WORLD.rank == 0:
        reg_write(reducedSum, rel_tol, abs_tol)

def reg_par_write_norm(values, rel_tol=1e-12, abs_tol=1e-12):
    """Write the sum of sum of the values from all processors."""
    reducedSum = MPI.COMM_WORLD.reduce(numpy.sum(values**2))
    if MPI.COMM_WORLD.rank == 0:
        reg_write(numpy.sqrt(reducedSum), rel_tol, abs_tol)

def reg_write_dict(d, rel_tol=1e-12, abs_tol=1e-12):
    """Write all values in a dictionary in sorted key order"""
    for key in sorted(d.keys()):
        print ('Dictionary Key: %s'%key)
        if isinstance(d[key],dict):
            reg_write_dict(d[key], rel_tol, abs_tol)
        elif type(d[key]) == bool:
            reg_write(int(d[key]), rel_tol, abs_tol)
        else:
            reg_write(d[key], rel_tol, abs_tol)

def reg_root_write_dict(d, rel_tol=1e-12, abs_tol=1e-12):
    """Only write from the root proc"""
    if MPI.COMM_WORLD.rank == 0:
        reg_write_dict(d, rel_tol, abs_tol)

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

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Single int write:')
        reg_write(1)

        print('Single float write:')
        reg_write(3.14159)

        print('List write:')
        reg_write([1.0, 3.5, 6.0], 1e-8, 1e-10)

        print('1D Numpy array write')
        vals = numpy.linspace(0, numpy.pi, 5)
        reg_write(vals, 1e-12, 1e-12)

        print('2D Numpy array write:')
        vals = numpy.linspace(0, 9.876, 4).reshape((2, 2))
        reg_write(vals)

        str1 = "@value    3.141592653589793 1e-12 1e-12"
        str2 = "@value    3.141592653589999 1e-12 1e-12"

        print('This comp should be True: ', _reg_str_comp(str1, str2))

        str1 = "@value    3.141592653589793 1e-12 1e-12"
        str2 = "@value    3.141592999999999 1e-12 1e-12"

        print('This comp should be False: ', _reg_str_comp(str1, str2))

    else:
        res = reg_file_comp(sys.argv[1], sys.argv[2])
        if res == 0:
            print ('Success!')
        elif res == 1:
            print ('Failure!')