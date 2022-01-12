# ------------- VERY IMPORTANT ------------

# This script is necessary since f2py INSISTS on priting crap out when
# .f2py_f2cmap exists in the directory. Normally it get deleted, but
# if its still around a more naive approach will fail miseribly. Here
# we temporily reassign stdout such that when we import it, the output
# goes to stdout. Then we reassign stdout and simply puck off the
# include  directory.
import os
import sys
import numpy.f2py

tmp = sys.stdout
sys.stdout = sys.stderr

sys.stdout = tmp
print(os.path.dirname(os.path.abspath(numpy.f2py.__file__)))
