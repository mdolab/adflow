#! /usr/bin/env python

import sys

name = 'sumb'
print "Testing if module %s can be imported..." % name
try:
    import sumb
except Exception, inst:
    print "Error: %s." % inst
    sys.exit(1)
# end try

print "Module %s was successfully imported." % name
