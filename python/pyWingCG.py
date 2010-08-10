#!/usr/local/bin/python
'''
pyCG - A simple center of Gravity calculation class for use in stability optimizations.

Copyright (c) 2010 by Mr.C.A (Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 09/08/2010 11:00$


Developers:
-----------
- Mr. C.A.(Sandy) Mader (SM)

History
-------
	v. 1.0	- Original pyHF Framework Implementation (SM 2010)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''
# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import pdb
import time
from math import pi,cos,tan


# =============================================================================
# Extension modules
# =============================================================================
from mdo_import_helper import *
exec(import_modules('pyGeometry_liftingsurface_c'))

# =============================================================================
# Misc Definitions
# =============================================================================

## class DUMMYMAPPING(Mapping):

##     '''
##     Abstract Class for Surface to Surface mapping
##     '''

##     def __init__(self, *args, **kwargs):
		
## 		'''
## 		DUMMYMAPPING Class Initialization
		
## 		Documentation last updated:  July. 21, 2008 - C.A.(Sandy) Mader
                
## 		'''
##                 name = 'DUMMYMAPPING'
## 		category = 'Surface to surface mapping tool'
## 		def_opts = {
## 			}
## 		informs = {
## 			}

##                 #run the Generic surface defintion
## 		Mapping.__init__(self, name, category, def_opts, informs, *args, **kwargs)

##                 return


def calculateWingMAC(acg):
    '''
    Calculates the Wing Mean Aerodynamic Chord(MAC).

    acg - instance of pyACDT aircraft geometry class
    '''
    #determine the number of surfaces in geometry
    ncomp = len(acg)
    #loop over surfaces
    for i in xrange(ncomp):
        print 'ncomp',i,ncomp
        #deal only with the lifting surfaces
        print 'isisnstance',isinstance(acg[i],LiftingSurface)
        if isinstance(acg[i],LiftingSurface):
            #Deal only with the wing
            if acg[i].Name.lower()=='wing':
                #Determin the number of segments that make up the wing
                nseg=len(acg[i])
                SumCSquared = 0.0
                SumArea = 0.0
                SumChord = 0.0
                #loop over wing segments
                for j in xrange(nseg):
                   #Segments are linear therfore single trapezoid gives exact answer

                   #Copy parameters from geometry
                   yrLE = acg[i][j].yrLE
                   xrLE = acg[i][j].xrLE
                   Span = acg[i][j].Span
                   Dihedral = acg[i][j].Dihedral*(pi/180)
                   Area = acg[i][j].Area
                   Taper = acg[i][j].Taper
                   SweepLE = acg[i][j].SweepLE*(pi/180)

                   #Compute tip span location
                   ytLE = yrLE + Span*cos(Dihedral)
                   print 'ytle',ytLE
                   # Root Chord
                   #c_abs!!!!
                   C_root = (2.0*Area)/(abs(Span)*(1.0+Taper))
                   print 'croot',C_root
                   # Tip Chord
                   C_tip = Taper*C_root
                   print 'ctip',C_tip
                   #Integrate C**2 for this segment
                   I = ((ytLE-yrLE)/2.0*(C_root**2+C_tip**2))/Area
                   print 'I',I
                   #Sum up intgrations
                   SumCSquared = SumCSquared + I*Area
                   print 'area',Area,SumCSquared
                   SumArea= SumArea+Area
                   print 'Sumarea',SumArea
                   #Compute the chordwise tip offset of this segment
                   xtLE = xrLE + abs(Span)*tan(SweepLE)
                   
                   #MAC quarterchord contribution of this panel
                   xc4 = xrLE+(xtLE-xrLE)*((1.0+2.0*Taper)/(3.0*(1.0+Taper)))+I/4.0
                   
                   SumChord = SumChord+xc4*Area
                #endfor

                #Final MAC computation
                MAC = SumCSquared/SumArea

                #Quarter Chord location
                MACc4 = SumChord/SumArea
            #endif
        #endif
    #endfor
    
    return MAC,MACc4

def calculateWingCenterOfGravity(forwardSparPercent,rearSparPercent,MAC,MACc4):
    '''
    Empirical correlation used to estimate wing CG.
    '''

    #Store the quarter chord as a reference point
    refpoint = MACc4
    #Calculate the distance to the MAC leading edge
    leadingEdgeOffset = MAC*1.0/4.0
    MAC_LE= MACc4-leadingEdgeOffset

    #calculate spar locations
    forwardSparOffset= MAC*forwardSparPercent
    rearSparOffset= MAC*rearSparPercent

    x_ForwardSpar=MAC_LE+forwardSparOffset
    x_RearSpar=MAC_LE+rearSparOffset

    #Calculate the distance between spars
    deltaSpar=x_RearSpar-x_ForwardSpar

    #Take 55% of that value(reference for this?)
    deltaCG = deltaSpar*0.55

    #compute CG location
    x_CG = x_ForwardSpar+deltaCG 


    return x_CG

