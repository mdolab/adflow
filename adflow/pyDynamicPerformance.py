#!/usr/local/bin/python
'''
pyDynamicPerformance - A set of routine to handle the computation of various dynamic handling qualities parameters.

Copyright (c) 2010 by Mr.C.A (Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 18/08/2010 11:00$


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
import numpy
#from cmath import pi,cos,sin

# =============================================================================
# Extension modules
# =============================================================================




def calculateThumbnailMethodConstraint(Wn,DampingRatio):
    '''
    Calculates the value of the thumnail plot constraint base on the natural
    frequncy and damping ratio.

    Wn - Natural frequency of the aircraft.

    This constraint is satisfied when z <= 0.  

    '''

    #Various function valus to allow the ellipse to coincide with
    #the satisfactory level on the thumbnail plot
    theta = 80 *numpy.pi/180.
    xcen = 0.76
    ycen = 3.01

    a = 2.2
    b=11
    c=0
    d=0.
    e=0.
    f=-1

    xp = DampingRatio
    yp = Wn

    #simple cordinate transformation
    x = (xp-xcen)*numpy.cos(theta)+(yp-ycen)*numpy.sin(theta)
    y = (yp-ycen)*numpy.cos(theta)-(xp-xcen)*numpy.sin(theta)

    #Equation of ellipse whose zero contour is roughly alligned with the
    #Satisfactory contour on the thumbnail plot.
    z = a*(x-d)*(x-d)+b*(y-e)*(y-e)+c*x*y+f
    

    return z

def calculateFrequencyAndDamping(Cmq,Clalpha,Cd,Cmalpha,Cmalphadot,mass,Iy,\
                                 rho,Area,U,c):
    '''
    Calculates the natural frequncy and damping for the aircraft using the
    second order 2 dof approximation from pg 309 of McRuer et al. 1973
    '''



    #Normalization of derivatives
    Mq = Cmq*rho*Area*U*c**2/(4*Iy)

    Zw = (-Clalpha-Cd)*rho*Area*U/(2*mass)

    Malpha = Cmalpha*rho*Area*U**2*c/(2*Iy)

    Malphadot = Cmalphadot*rho*Area*U*c**2/(4*Iy)
    #print 'mq',Mq,Zw,Malpha
    #Short period approximation
    Wsp = numpy.sqrt(Mq*Zw-Malpha)
    #print 'wsp',Wsp
    Damping = -(Zw+Mq+Malphadot)/(2*Wsp)
    #print 'Damping',Damping

    return Wsp,Damping

def calculateNAlpha(Clalpha,rho,Area,U,mass,g):

    '''
    calculate the g normalized lift derivative
    '''

    nalpha = rho*U**2*A*Clalpha/(2*mass*g)

    return nalpha


def calculateCAP(Wsp,nalpha):
    '''
    Compute the CAP parameter. Currently based on controls fixed
    interpretation.
    '''

    CAP = Wsp**2/nalpha

    return CAP
