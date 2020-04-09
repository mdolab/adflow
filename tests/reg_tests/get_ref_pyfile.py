from commonUtils import convertRegFileToPyRegFile
import glob
import os 

# This is the total number of refs we have:
refDir = '../../python/reg_tests/ref/'
refFiles = glob.glob(os.path.join(refDir,'*.ref'))

# Run each script
for refFile in refFiles:
    
    ref_id = int(refFile.split('_')[-2][4:])

    print('Running reference for ref%d'%ref_id)
    convertRegFileToPyRegFile(refFile, 'ref/ref%d.py'%(ref_id))
    