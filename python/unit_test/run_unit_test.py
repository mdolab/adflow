''' This script is designed to test and verify the key functionality
in pySUmb. To execute the script run:

python run_unit_test.py [-mpiexec mpiexec] [-ref file] [-diff diff_program]

There are two optional argument: 

-mpiexec  which is used to specify a custom mpirun command. Default is
mpirun.

-ref      file to use as reference instead of the master reference. 

-diff     program to use to compare output. By default xxdiff is used


Several different examples are run and the output piped to a file
called output.out. A reference file, output.out.ref is included for
comparison if -ref is not specified. These two files are compared for
any differences. If any difference is found a comparison between the
two files are made to show where the differences are. 
'''


# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, argparse

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================

from mdo_import_helper import MPI

if __name__ == '__main__':    
    # First check if the user ran this in parallel:
    if MPI.COMM_WORLD.size > 1:
        mpiPrint('Error: run_unit_test.py must be run with only one processor')
        sys.exit(1)
    # end if

    # Get the optional commandline arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpiexec",default='mpirun',
                        help='Command to use for running mpi')
    parser.add_argument("--ref",default='output.out.ref',
                        help='File to use for comparison')
    parser.add_argument("--diff",default='xxdiff',
                        help='Command to run for displaying diff')
    args = parser.parse_args()

    mpirun = args.mpiexec
    ref_file = args.ref
    diff_cmd = args.diff

    f = open('output.out.ref','r') 
    g = open(ref_file,'r') 

    if len(f.readlines())!=len(g.readlines()): 
        print 'Local reference file different length than master reference\
 file. Please update local reference file.'
        sys.exit(1)
    # end if

    serial = 1
    try:
        import multiprocessing
        parallel = multiprocessing.cpu_count()
        parallel = min(4, parallel)
    except:
        parallel = 2
    # end if

    # Define the test matrix:
    test_matrix = [
        [serial  , 'bump', 'euler', 'NK', 'Steady', 'rev'],
        [serial  , 'bump', 'euler', 'NK', 'Steady', 'fwd'],
        [parallel, 'bump', 'euler', 'NK', 'Steady', 'rev'],
        [parallel, 'bump', 'euler', 'NK', 'Steady', 'fwd'],

        [serial  , 'wing', 'euler', 'NK', 'Steady', 'rev'],
        [serial  , 'wing', 'euler', 'NK', 'Steady', 'fwd'],
        [parallel, 'wing', 'euler', 'NK', 'Steady', 'rev'],
        [parallel, 'wing', 'euler', 'NK', 'Steady', 'fwd'],

        [serial  , 'wing', 'euler', 'NK', 'TS',     'rev'],
        [serial  , 'wing', 'euler', 'NK', 'TS',     'fwd'],
        [parallel, 'wing', 'euler', 'NK', 'TS',     'rev'],
        [parallel, 'wing', 'euler', 'NK', 'TS',     'fwd']
        ]

    f = open('output.out','w')

    for t in test_matrix:
        # Compose command 
        cmd = '%s -np %d python solve_script.py %s %s %s %s %s | tee tmp.out'%(
            mpirun, t[0], t[1], t[2], t[3], t[4], t[5])

        # Print some information about the case we're running:
        f.write('#'*80+'\n')
        f.write(' Running %s case with %d processors in %s mode with %s solver.\n'%(
            t[1], t[0], t[2], t[3]))
        f.write(' %s solution is performed. Adjoint is computed with %s mode.\n'%(
            t[4], t[5]))
        f.write('#'*80+'\n')

        # Run Command
        os.system(cmd)

        # Load in tmp and write to f
        g = open('tmp.out','r')
        for line in g:
            f.write(line)
        # end for
        
        g.close()
    # end for

    # Close full file
    f.close()

    # Remove tmp.out
    os.system('rm -fr tmp.out')

    # Comparison:
    try:
        f = open(ref_file,'r')
    except:
        print 'ERROR: %s is not found. This file should be \
included on the repository and not revomed. Revert file to restore'%ref_file
        sys.exit(1)
    # end try

    g = open('output.out','r')
    
    f_lines = f.readlines()
    g_lines = g.readlines()
    diff_files = False

    if len(f_lines) != len(g_lines):
        
        print '+'+'-'*78+'+'
        print ' Output File are not EXACTLY identical!'
        print '+'+'-'*78+'+'
        diff_files = True

    else:
        # Check each line 

        for i in xrange(len(f_lines)):
            if f_lines[i] != g_lines[i]:
                if 'time' in f_lines[i] or 'Time' in f_lines[i]:
                    # Lines that list times are ok...they can be different
                    pass
                else:
                    print '+'+'-'*78+'+'
                    print ' Output File are not EXACTLY identical!'
                    print '+'+'-'*78+'+'
                    diff_files = True
                    break
            # end if
        # end for
    # end if

    if diff_files:
        print '+'+'-'*78+'+'
        print ' Running Diff of files. Reference Output on left, current on right'
        print '+'+'-'*78+'+'
        cmd = '%s %s output.out &'%(diff_cmd,ref_file)
        os.system(cmd)
    else:
        print '+'+'-'*78+'+'
        print ' Output and reference file are identical! :)'
        print '+'+'-'*78+'+'
    # end if
        
    import datetime
    new_name = 'output.out_%s'%(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))
    os.system('cp output.out %s'%(new_name))
    
