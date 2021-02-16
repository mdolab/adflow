#!/usr/bin/python

"""
The script will determine the number of processors that will result in the 
best loadbalancing between procs given a mesh and multigrid cycles. 

User can specify either specific number of procs that he wants to check
or a range of processors that he wants to check.

Note!
To get valid output the input needs to be valid as well. For example it 
does not make sense to check large number of procs and large number of
multigrid cycles for a small mesh as ADflow will crash on that.

Example usage:

#1 Get load balcance information for only two processor numbers
$ python numOfProcsLB.py fine_L1.cgns  32 64 --mgCycle 4w 

#2 Get load balcance information for only two processor numbers
python numOfProcsLB.py fine_L0.cgns  256 512 --sequence --mgCycle 4w

"""

from adflow import ADFLOW
import argparse
import numpy
import os

parser = argparse.ArgumentParser()
parser.add_argument("grid", type=str, help="Set the grid to be analyzed")
parser.add_argument("procs", type=int, nargs="+", help="Specify which procs you want to check the loadbalacning for. This can be a single number or list of numbers")
parser.add_argument("--sequence", action='store_true', help="If specified, the first two arguments from procs argument are used to create a range where all proc numbers in btween are tested.")
parser.add_argument("--mgCycle", type=str, default='2w', help="Set the level of multigrid intended to run on mesh [default: %(default)s]")
parser.add_argument("--outfile", type=str, help="Set the output filename for the analysis. Default: filename is generated from mesh name and other input.")
args = parser.parse_args()

# Parse input
if args.sequence:
    sidx = args.procs[0]
    eidx = args.procs[1]

    rge = range(sidx,eidx+1)
else:
    rge = args.procs


outName = args.outfile
if not args.outfile:
    # Construct a default 
    pth = os.path.dirname(args.grid)
    nme = os.path.basename(args.grid)
    outName = "lbAnalysis_{0}_{1}.log".format(os.path.splitext(nme)[0], args.mgCycle)
    outName = os.path.join(pth, outName)

# Define problem
aeroOptions = {'gridFile':args.grid,
                "mgcycle":args.mgCycle,
                "partitionOnly":True}

CFDSolver = ADFLOW(options=aeroOptions)

# Write results
loadImbalance =[]
faceImbalance =[]

fd = open(outName,"w")
fd.write(72*"=" + "\n")
fd.write("Load imbalance analysis for mesh {0}\n\n".format(args.grid))
fd.write("Multigrid cycle set to {0}\n\n".format(args.mgCycle))

for i in rge:   
    lIb, fIb = CFDSolver.checkPartitioning(i)
    loadImbalance.append(lIb)
    faceImbalance.append(fIb)
    fd.write("For {0:>5d} procs, Load imbalance : {1:>5.3f},  Face imbalance : {2:>5.3f}\n".format(i,lIb,fIb))

A = numpy.array((rge, loadImbalance, faceImbalance)).T

minLoad = A[A[:,1].argsort()]
minFace = A[A[:,2].argsort()]

fd.write(72*"-" + "\n")
fd.write("The folloing gives the smallest load imbalance\n")
numpy.savetxt(fd, minLoad, fmt="%5d %5.3f %5.3f")

fd.write(72*"-" + "\n")
fd.write("The folloing gives the smallest face imbalance\n")
numpy.savetxt(fd, minFace, fmt="%5d %5.3f %5.3f")

fd.close()    

print(90*"=" + "\n")
print("Loal balance analysis written to: \n" + outName)
print("\n" + 90*"=")
