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
#from cmath import pi,cos,tan,sin
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from mdo_import_helper import *
exec(import_modules('pyGeometry_liftingsurface_c'))
#exec(import_modules('pyGeometry_liftingsurface'))
#from pyGeometry_complex import catan,cabs,cmin,cmax
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
    for i in range(ncomp):
        #print 'ncomp',i,ncomp
        #deal only with the lifting surfaces
        #print 'isisnstance',isinstance(acg[i],LiftingSurface)
        if isinstance(acg[i],LiftingSurface):
            #Deal only with the wing
            if acg[i].Name.lower()=='wing':
                #Determin the number of segments that make up the wing
                nseg=len(acg[i])

                sumSpan = 0.0
                sumArea = 0.0
                for j in range(nseg):
                    Span = acg[i][j].Span
                    Area = acg[i][j].Area
                    sumSpan = sumSpan+Span
                    sumArea = sumArea+Area
                #endfor
                #print 'Span,Area',sumSpan,sumArea
                SumCSquared = 0.0
                
                SumChord = 0.0
                SumNumerator = 0.0
                SumDenomenator = 0.0
                sweepsum = 0.0
                xC4sum = 0.0
                #loop over wing segments
                for j in range(nseg):
                   #Segments are linear therfore single trapezoid gives exact answer

                   #Copy parameters from geometry
                   yrLE = acg[i][j].yrLE
                   xrLE = acg[i][j].xrLE
                   Span = acg[i][j].Span
                   Dihedral = acg[i][j].Dihedral*(numpy.pi/180)
                   Area = acg[i][j].Area
                   Taper = acg[i][j].Taper
                   SweepLE = acg[i][j].SweepLE*(numpy.pi/180)

                   #Compute tip span location
                   ytLE = yrLE + Span*numpy.cos(Dihedral)
                   #print 'ytle',ytLE
                   # Root Chord
                   #c_abs!!!!
                   C_root = (2.0*Area)/(numpy.abs(Span)*(1.0+Taper))
                   #print 'croot',C_root
                   # Tip Chord
                  # C_tip = Taper*C_root
                   #print 'ctip',C_tip

                   #original computation
                   #Integrate C**2 for this segment
                  # I = ((ytLE-yrLE)/2.0*(C_root**2+C_tip**2))/Area
                   #print 'I',I
                   #Sum up intgrations
                 #  SumCSquared = SumCSquared + I*Area
                  
                   #print 'Sumarea',SumArea
                   #Compute the chordwise tip offset of this segment
                 #  xtLE = xrLE + abs(Span)*tan(SweepLE)
                   
                   #MAC quarterchord contribution of this panel
                 #  xc4 = xrLE+(xtLE-xrLE)*((1.0+2.0*Taper)/(3.0*(1.0+Taper)))+I/4.0
                   
                 #  SumChord = SumChord+xc4*Area
                   #new xC4 computation
                   temp1 = (1+2*Taper)*((ytLE-yrLE)/sumSpan)\
                            *numpy.tan(SweepLE)

                   temp2 = 3*(1+Taper)*sweepsum

                   temp3 = ((ytLE-yrLE)/sumSpan)*C_root
                   #print 'temp',temp1,temp2,temp3
                   xC4sum = xC4sum+(temp3*(temp2+temp1))

                   
                   sweepsum = sweepsum+((ytLE-yrLE)/sumSpan)*numpy.tan(SweepLE)
                   

                   #new computation...
                   numerator = C_root**2*(1+Taper+Taper**2)*(ytLE-yrLE)/sumSpan
                   denomenator = C_root*(1+Taper)*(ytLE-yrLE)/sumSpan
                   SumNumerator= SumNumerator+numerator
                   SumDenomenator=SumDenomenator+denomenator
                #endfor

                #Final MAC computation
                #MAC = SumCSquared/sumArea
                MAC = (2.0/3.0)*(SumNumerator/SumDenomenator)
                #Quarter Chord location
                #MACc4 = SumChord/sumArea
                MACc4 = (((sumSpan*2)**2/(sumArea*2))/12)*xC4sum+MAC/4.0
            #endif
        #endif
    #endfor
 #   print 'MAC',MAC
 #   print 'xc4',MACc4
    return MAC,MACc4

def calculateWingCenterOfGravity(forwardSparPercent,rearSparPercent,CGPercentage,MAC,MACc4):
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

    #Take 55% of that value(reference for this?)Reference says anywhere between spars
    deltaCG = deltaSpar*CGPercentage

    #compute CG location
    x_CG = x_ForwardSpar+deltaCG 


    return x_CG

def calculateWingInertias(acg):
    '''
    Calculates the Wing mass moments of intertia.

    acg - instance of pyACDT aircraft geometry class

    note:
    x is chord direction
    y is span direction
    z is vertical direction
    '''


    Ix = 0.0
    Iy = 0.0
    Iz = 0.0
    #determine the number of surfaces in geometry
    ncomp = len(acg)
    #loop over surfaces
    for i in range(ncomp):
        #print 'ncomp',i,ncomp
        #deal only with the lifting surfaces
        #print 'isisnstance',isinstance(acg[i],LiftingSurface)
        if isinstance(acg[i],LiftingSurface):
            #Deal only with the wing
            if acg[i].Name.lower()=='wing':
                #Determin the number of segments that make up the wing
                nseg=len(acg[i])
                #loop over wing segments
                for j in range(nseg):
                   
                   #Copy parameters from geometry
                   tc_root = acg[i][j].root_Thickness #thickness to chord ratio...
                   tc_tip = acg[i][j].tip_Thickness
                   yrLE = acg[i][j].yrLE
                   xrLE = acg[i][j].xrLE
                   zrLE = acg[i][j].zrLE
                   Span = acg[i][j].Span
                   Dihedral = acg[i][j].Dihedral*(numpy.pi/180)
                   Area = acg[i][j].Area
                   Taper = acg[i][j].Taper
                   AR = Span**2/Area # Aspect Ratio
                   SweepLE = acg[i][j].SweepLE*(numpy.pi/180)
                   SweepTE = numpy.arctan(numpy.tan(SweepLE)-((4.0*(1.0/1.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)))
                   #print 'sweeps',SweepLE*(180.0/pi),SweepTE*(180.0/pi)
                   #Compute tip span location
                   ytLE = yrLE + Span*numpy.cos(Dihedral)
                   #print 'ytle',ytLE
                   # Root Chord
                   #c_abs!!!!
                   C_root = (2.0*Area)/(numpy.abs(Span)*(1.0+Taper))
                   #Root thickness
                   t_root = tc_root*C_root
                   #print 'croot',C_root,tc_root
                   # Tip Chord
                   C_tip = Taper*C_root
                   #tip thickness
                   t_tip = tc_tip*C_tip

                   #Weight
                   W = acg[i][j].Weight
                   #print 't,w',t_root,t_tip,W
                   #Volume???
                   V = Span*(t_root*(C_root+(Span/2.0)*(numpy.tan(SweepTE)-numpy.tan(SweepLE)))-((t_root-t_tip)*((C_root/2.0)+(Span/3.0)*(numpy.tan(SweepTE)-numpy.tan(SweepLE)))))
                   #print 'V',V
                   #I1x = (t_root-t_tip)*(C_root/4.0+(Span*tan(SweepTE)/5.0)-Span*tan(SweepLE)/5.0)#+(t_root*(C_root/3.0+(Span*tan(SweepTE)/4.0)-Span*tan(SweepLE)/4.0)))
                   #print 'I1x1',I1x
                   #I1x =(t_root*(C_root/3.0+(Span*tan(SweepTE)/4.0)-Span*tan(SweepLE)/4.0))
                   #print 'I1x2',I1x
                   I1x = (W*Span**3/V)*(((t_root-t_tip)*(C_root/4.0+(Span*numpy.tan(SweepTE)/5.0)-Span*numpy.tan(SweepLE)/5.0))+(t_root*(C_root/3.0+(Span*numpy.tan(SweepTE)/4.0)-Span*numpy.tan(SweepLE)/4.0)))
                   #I1x =27028033000
                   #print 'I1x',I1x
                   #Ok, except for I1x!!!!
                   I1y = ((W*Span)/V)*((t_root*((C_root**3/3.0)+Span*C_root*numpy.tan(SweepTE)*((C_root/2.0)+(Span*numpy.tan(SweepTE))/3.0)+(Span**3/12.0)*(numpy.tan(SweepTE)**3-numpy.tan(SweepLE)**3)))-((t_root-t_tip)*((C_root**3/6.0)+Span*C_root*numpy.tan(SweepTE)*((C_root/3.0)+(Span*numpy.tan(SweepTE)/4.0))+(Span**3/15.0)*(numpy.tan(SweepTE)**3-numpy.tan(SweepLE)**3))))
                   #print 'I1y',I1y
                   I1z = I1x+I1y
                   #print 'I1z',I1z
                   TI1y = (I1y*numpy.cos(Dihedral)+I1z*numpy.sin(Dihedral))
                   TI1z = I1y*numpy.sin(Dihedral)+I1z*numpy.cos(Dihedral)

                   I1y = TI1y
                   I1z = TI1z
                   #print 'TI1y',I1y
                   #print 'TI1z',I1z
                   a = numpy.array([C_root,Span*numpy.tan(SweepLE), Span*numpy.tan(SweepLE)+C_tip])
                   #print 'a',a
                   a.sort()
                   #print 'asort',a
                   Ca = a[0]#cmin(C_root,Span*tan(SweepLE), Span*tan(SweepLE)+C_tip)
                   #print 'ca',Ca
                   Cc = a[2]#cmax(C_root,Span*tan(SweepLE), Span*tan(SweepLE)+C_tip)
                   #print 'cc',Cc
                   Cb = a[1]#C_root+(Span*tan(SweepLE))+( Span*tan(SweepLE)+C_tip)-Ca-Cc
                   #print 'cb',Cb
                   K_o = 0.703 #(for a wing....)
                   acg[i][j].x_Centroid=((-Ca**2+Cb**2+Cc*Cb+Cc**2)/(3*(Cb+Cc-Ca)))*K_o**(1.0/2.0)
                   #print 'xs1', acg[i][j].x_Centroid
                   acg[i][j].y_Centroid=(Span**2/V)*((t_root*((C_root/2.0)+(Span/3.0)*(numpy.tan(SweepTE)-numpy.tan(SweepLE)))-(t_root-t_tip)*((C_root/3.0)+(Span/4.0)*(numpy.tan(SweepTE)-numpy.tan(SweepLE)))))
                   #print 'ys1', acg[i][j].y_Centroid
                   acg[i][j].z_Centroid=acg[i][j].y_Centroid*numpy.sin(Dihedral)

                   Xs = acg[i][j].x_Centroid
                   Xs4 = xrLE

                   Ys = acg[i][j].y_Centroid
                   Ys_dot = acg[i][j].y_Centroid*numpy.cos(Dihedral)
                   #print 'ydot',Ys_dot
                   Ysoff = yrLE

                   Zs1 = zrLE
                   Zs3 = acg[i][j].z_Centroid
                   #print 'zs3',Zs3
                   Zs4 = acg[i][j].z_Centroid
                   

                   Ix = Ix+I1x-W*(Ys_dot**2)-W*(Zs3**2)+W*(Ys_dot+Ysoff)**2 + W*(Zs3+Zs1)**2
                   Iy = Iy+I1y-W*(Xs**2)-W*(Zs3**2)+W*(Xs+Xs4)**2+W*(Zs3+Zs1)**2
                   Iz = Iz+I1z-W*(Xs**2+Ys_dot**2)+W*(Xs+Xs4)**2+W*(Ys_dot+Ysoff)**2
                   #Ixz = ...
                   
                #endfor
            #endif
        #endif
    #endfor

   
    return Ix,Iy,Iz

## def getAverageThickness(acg,surface,forwardSparPercent,rearSparPercent,npts):
##     percentchord = numpy.zeros([npts],numpy.float)
##     for i in range(npts):
##         percentchord[i] =forwardSparPercent+(rearSparPercent-forwardSparPercent)*( float(i)/float(npts-1))
##     #endfor
##     thickness =surface.getThickness(npts,percentchord)

##     ncomp = len(acg)
    
##     for i in range(ncomp):
##         if isinstance(acg[i],LiftingSurface):
##             #Deal only with the wing
##             if acg[i].Name.lower()=='wing':
##                 rootIndex = 0
##                 nseg=len(acg[i])
##                 for j in range(nseg):
                    
##                     if j ==0:
##                         tipIndex = rootIndex+acg[i][j].surface_SW_segments
##                     else:
##                         tipIndex = rootIndex+acg[i][j].surface_SW_segments-1
##                     #endif

##                     Span = acg[i][j].Span
##                     Area = acg[i][j].Area
##                     Taper = acg[i][j].Taper
                    
##                     #Chord Lengths
##                     C_root = (2.0*Area)/(numpy.abs(Span)*(1.0+Taper))
##                     C_tip = Taper*C_root
##                     #root thickness
##                     t_root = numpy.mean(thickness[i][rootIndex][:])
##                     #tip  thickness
##                     t_tip = numpy.mean(thickness[i][tipIndex][:])
                    
##                     tc_root = t_root/C_root
                    
##                     tc_tip = t_tip/C_tip
                    
##                     acg[i][j].root_Thickness  = tc_root#thickness to chord ratio...
##                     acg[i][j].tip_Thickness = tc_tip 
##                     rootIndex = tipIndex
##                 #endfor
##             #endif
##         #endif
##     #endfor

##     return


def getAverageThickness(acg,thickness):
    '''
    Take in the DVcontraints object and compute the updated
    thicknesses form the wing.
    '''
    ncomp = len(acg)
    
    for i in range(ncomp):
        if isinstance(acg[i],LiftingSurface):
            #Deal only with the wing
            if acg[i].Name.lower()=='wing':
                rootIndex = 0
                nseg=len(acg[i])
                for j in range(nseg):
                    
                   ##  if j ==0:
##                         tipIndex = rootIndex+acg[i][j].surface_SW_segments
##                     else:
##                         tipIndex = rootIndex+acg[i][j].surface_SW_segments-1
##                     #endif

                    Span = acg[i][j].Span
                    Area = acg[i][j].Area
                    Taper = acg[i][j].Taper
                    
                    #Chord Lengths
                    C_root = (2.0*Area)/(numpy.abs(Span)*(1.0+Taper))
                    C_tip = Taper*C_root
                    #root thickness
                    t_root = numpy.mean(thickness[j][:])
                    #tip  thickness
                    t_tip = numpy.mean(thickness[j+1][:])
                    
                    tc_root = t_root/C_root
                    
                    tc_tip = t_tip/C_tip
                    
                    acg[i][j].root_Thickness_act  = tc_root#thickness to chord ratio...
                    acg[i][j].tip_Thickness_act = tc_tip 
                    #rootIndex = tipIndex
                #endfor
            #endif
        #endif
    #endfor

    return

def calculateSegmentWeights(acg,Weight):
    '''
    Calculates the segment weight for each LiftingSegment
    '''

    ncomp = len(acg)
    
    for i in range(ncomp):
        if isinstance(acg[i],LiftingSurface):
            #Deal only with the wing
            if acg[i].Name.lower()=='wing':
                TotalVolume = 0.0
                nseg=len(acg[i])
                for j in range(nseg):

                    #thickness to chord ratios...
                    tc_root =acg[i][j].root_Thickness_act
                    tc_tip = acg[i][j].tip_Thickness_act

                    #Wing properties
                    Span = acg[i][j].Span
                    Area = acg[i][j].Area
                    Taper = acg[i][j].Taper
                    
                    #Chord Lengths
                    C_root = (2.0*Area)/(numpy.abs(Span)*(1.0+Taper))
                    C_tip = Taper*C_root

                    #Segment thicknesses
                    t_root = C_root* tc_root                    
                    t_tip = C_tip * tc_tip

                    #Approx Volume
                    V = Area*(t_root+t_tip)/2.0

                    acg[i][j].volumeWeight = V

                    TotalVolume = TotalVolume+V
                #endfor
            #endif
        #endif
    #endfor
    for i in range(ncomp):
        if isinstance(acg[i],LiftingSurface):
            #Deal only with the wing
            if acg[i].Name.lower()=='wing':
                nseg=len(acg[i])
                for j in range(nseg):
                    
                    acg[i][j].Weight = Weight*(acg[i][j].volumeWeight/TotalVolume)
                #endfor
            #endif
        #endif
    #endfor
                    
    return
