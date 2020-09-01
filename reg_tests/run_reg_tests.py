# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, argparse, glob

# =============================================================================
# Extension modules
# =============================================================================
import mdo_regression_compare as reg

# define scripts to run:
module_name = 'adflow'

# Get the optional commandline arguments:
parser = argparse.ArgumentParser()
parser.add_argument("--mode",default='compare',choices=['train','compare'],
                    help='Train generates reference file. Compare runs test')

parser.add_argument("--procs",default=4, type=int,
                    help='The number of processors to use to run tests.')

parser.add_argument("--diff_cmd",default='xxdiff',
                    help='Command to run for displaying diff. Default: xxdiff')

parser.add_argument("--diff", action='store_true', default=False,
                    help='Display error diffs for each test.')

parser.add_argument("--mpiexec",default='mpirun',
                    help='Command to use for running mpi')

parser.add_argument('--test', metavar='test', type=int, nargs='+',
                    help='tests to run')

parser.add_argument("--solve", action='store_true', default=False,
                    help="Force solving on tests that use restart files.")
args = parser.parse_args()

# This is the total number of tests we have:
testDir = os.path.join('tests','test*.py')
testFiles = glob.glob(testDir)
nTestMax = len(testFiles)

if args.test is None:
    tests = range(1, nTestMax+1)
else:
    tests = args.test

solveStr = ""
if args.solve:
    solveStr = 'solve'

if args.mode == 'train':
    try:
        os.remove('%s_reg.ref'%(module_name))
    except OSError:
        pass

    # Run each script
    for iTest in tests:
        print('Running reference for test%d'%iTest)
        os.system('%s -np %d python tests/test%d.py %s > ref/%s_test%d_reg.ref 2>&1'%(
            args.mpiexec, args.procs, iTest, solveStr, module_name, iTest))

    # If we're training, we done (no comparison)
    sys.exit(0)
else:
    try:
        os.remove('%s_reg'%(module_name))
    except OSError:
        pass

    masterRes = 0

    os.system('rm -fr adflow_reg.ref adflow_reg adflow_reg.orig')
    os.system('touch adflow_reg.ref adflow_reg adflow_reg.orig')

    # Run each script
    for iTest in tests:
        os.system('%s -np %d python tests/test%d.py %s > %s_test%d_reg 2>&1'%(
            args.mpiexec, args.procs, iTest,  solveStr, module_name, iTest))

        refFile = 'ref/%s_test%d_reg.ref'%(module_name, iTest)
        curFile = '%s_test%d_reg'%(module_name, iTest)

        # Do the comparison (reference file must be first)
        res = reg.reg_file_comp(refFile, curFile)

        # Set the proper return codes for the script running this:
        if res == reg.REG_FILES_MATCH:
            print('%s test%d: Success!'%(module_name, iTest))
        elif res == reg.REG_FILES_DO_NOT_MATCH:
            print('%s test%d: Failure!'%(module_name, iTest))
            if  args.diff:
                os.system('%s %s %s'%(args.diff_cmd, refFile, curFile))
            masterRes += 1
        elif res == reg.REG_ERROR:
            print('%s test%d: Error in regression. Missing files.'%(module_name, iTest))
            masterRes += 1
        sys.stdout.flush()

        # Concentenate outputs for reference if it failed:
        if res == reg.REG_FILES_DO_NOT_MATCH or res == reg.REG_ERROR:
            os.system('cat %s >> adflow_reg.ref'%(refFile))
            os.system('cat %s >> adflow_reg'%(curFile))
            os.system('cat %s.orig >> adflow_reg.orig'%(curFile))

# Exit with code equal to the number of failures
sys.exit(masterRes)