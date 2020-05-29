import glob
import os 
from reg_test_utils import convertRegFileToJSONRegFile




# This is the total number of tests we have:
refDir = '../../python/reg_tests/ref'

refFiles = glob.glob(os.path.join(refDir, 'adflow_test*'))
outputDir = './reg_refs/ref'

for refFile in refFiles:

    basename = os.path.basename(refFile)
    basename = os.path.splitext(basename)[0]
    print(basename)
    outputFile = os.path.join(outputDir, basename)

    convertRegFileToJSONRegFile(refFile, outputFile + '.json')


# if args.ref is None:
#     refs = range(1, nrefMax+1)
# else:
#     refs = args.ref

# solveStr = ""
# if args.solve:
#     solveStr = 'solve'

# if args.mode == 'train':
#     try:
#         os.remove('%s_reg.ref'%(module_name))
#     except OSError:
#         pass

#     # Run each script
#     for iTest in tests:
#         print('Running reference for test%d'%iTest)
#         os.system('%s -np %d python tests/test%d.py %s > ref/%s_test%d_reg.ref 2>&1'%(
#             args.mpiexec, args.procs, iTest, solveStr, module_name, iTest))

#     # If we're training, we done (no comparison)
#     sys.exit(0)
# else: