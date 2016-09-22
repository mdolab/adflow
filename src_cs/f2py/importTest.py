#! /usr/bin/env python

import sys

name = 'libadflow_cs'
print "Testing if module %s can be imported..." % name
import_cmd = "import %s" % name
try:
    exec import_cmd
except Exception, inst:
    print "Error: %s." % inst
    sys.exit(1)
# end try

print "Module %s was successfully imported." % name

