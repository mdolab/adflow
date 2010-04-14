#! /usr/bin/env python

import os,sys,string

signature = sys.argv[1]
suffix = sys.argv[2]
#print 'signature',signature,suffix

# Determine the name of the module from the .pyf file
f = open(signature,'r')
lines = f.readlines()
f.close()
for i in range(len(lines)-1,-1,-1):
    if (lines[i][0:17] == 'end python module'):
        name = string.strip(lines[i][18:])
        break
    # end if
# end for
#print 'Name',name
# If necessary, create an _parallel.pyf signature file.
for i in range(len(lines)):
    if (string.count(lines[i],'python module') > 0):
        lines[i] = string.replace(lines[i],name,
                                  '%s%s'%(name,suffix))
signature = signature[:-4] + '%s.pyf'%(suffix)
name = '%s%s'%(name,suffix)
if (os.path.isfile(signature)):
    os.remove(signature)
#print 'sig.',signature,name,signature[:-4] + '%s.pyf'%(suffix)
f = open(signature,'w')
f.writelines(lines)
f.close()





