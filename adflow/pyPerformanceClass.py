#!/usr/local/bin/python
'''
pyPerformanceClass - A class of routine to handle the computation of various dynamic handling qualities parameters.

Copyright (c) 2010 by Mr.C.A (Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 18/08/2010 11:00$


Developers:
-----------
- Mr. C.A.(Sandy) Mader (SM)

History
-------
	v. 1.0	- Original MACH Framework Implementation (SM 2010)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''
# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys,copy
import pdb
import time
import numpy 

# =============================================================================
# Extension modules
# =============================================================================


class PERFORMANCE(object):

    '''
    Basic Performance class to handle dynamic handling qualities constraints
    '''

    def __init__(self,*args,**kwargs):

        '''
        Performance class initialization

        '''

        name = 'Performance'
        category = 'Performance and Handling Qualities'
        def_opts = {}
        informs = { }

       
        return
    
    def calculateThumbnailMethodConstraint(self,Wn,DampingRatio):
        '''
        Calculates the value of the thumnail plot constraint base
        on the natural frequncy and damping ratio.

        Wn - Natural frequency of the aircraft.
        
        This constraint is satisfied when z <= 0.  
        
        '''
        #original function
       ##  #Various function valus to allow the ellipse to coincide with
##         #the satisfactory level on the thumbnail plot
##         theta = 80 *numpy.pi/180.
##         xcen = 0.76
##         ycen = 3.01
        
##         a = 2.2
##         b=11
##         c=0
##         d=0.
##         e=0.
##         f=-1.1
        
##         xp = DampingRatio
##         yp = Wn
        
##         #simple cordinate transformation
##         x = (xp-xcen)*numpy.cos(theta)+(yp-ycen)*numpy.sin(theta)
##         y = (yp-ycen)*numpy.cos(theta)-(xp-xcen)*numpy.sin(theta)
        
##         #Equation of ellipse whose zero contour is roughly alligned with the
##         #Satisfactory contour on the thumbnail plot.
##         z = a*(x-d)*(x-d)+b*(y-e)*(y-e)+c*x*y+f
        #Updated function
        #Various function values to allow the ellipse to coincide with
        #the satisfactory level on the thumbnail plot
        theta = 55 *numpy.pi/180.
        xcen = -0.4
        ycen = 3.01
        
        a = 0.3
        b=0.75
        c=-0.4
        d=0.
        e=0.
        f=-0.15
        
        xp = numpy.log(DampingRatio)
        yp = Wn
        
        #simple cordinate transformation
        x = (xp-xcen)*numpy.cos(theta)+(yp-ycen)*numpy.sin(theta)
        y = (yp-ycen)*numpy.cos(theta)-(xp-xcen)*numpy.sin(theta)
        
        #Equation of ellipse whose zero contour is roughly alligned with the
        #Satisfactory contour on the thumbnail plot.
        z = a*(x-d)*(x-d)+b*(y-e)*(y-e)+c*x*y+f
    
        
        return z

    def calculateFrequencyAndDamping(self,Cmq,Clalpha,Cd,Cmalpha,Cmalphadot,
                                     mass,Iy,rho,Area,U,c):
        '''
        Calculates the natural frequncy and damping for the aircraft using the
        second order 2 dof approximation from pg 309 of McRuer et al. 1973
        '''
        #Normalization of derivatives
        Mq = Cmq*rho*Area*U*c**2/(4*Iy)
        
        Zw = (-Clalpha-Cd)*rho*Area*U/(2*mass)
        
        Malpha = Cmalpha*rho*Area*U**2*c/(2*Iy)
        
        Malphadot = Cmalphadot*rho*Area*U*c**2/(4*Iy)
        
        #Short period approximation
        Wsp = numpy.sqrt(Mq*Zw-Malpha)
        
        Damping = -(Zw+Mq+Malphadot)/(2*Wsp)
        
        #print 'Frequerncy and Damping...:',Wsp,Damping
        return Wsp,Damping

    def calculateNAlpha(self,Clalpha,rho,Area,U,mass,g):

        '''
        calculate the g normalized lift derivative
        '''

        nalpha = rho*U**2*Area*Clalpha/(2*mass*g)
        
        return nalpha

    def computeStaticMargin(self,averagesol,geom,wbc,acg):
        '''
        compute the percent static margin and return
        '''
        #Compute the MAC
        [MAC,C4MAC] = wbc.calculateWingMAC(acg)
        SM = -averagesol['cmalpha']/averagesol['clalpha']

        percentSM = SM#/MAC
        #print 'SM',SM,MAC,-averagesol['cmalpha'],averagesol['clalpha']
        return percentSM

    
    def calculateCAP(self,Wsp,nalpha):
        '''
        Compute the CAP parameter. Currently based on controls fixed
        interpretation.
        '''
        
        CAP = Wsp**2/nalpha
        
        return CAP
    
    def thumbprintDriver(self,acg,wbc,geom,averagesol,rho,V,A,thickness):
        '''
        Driver to compute thumbprint constraint

        acg : pyACDT geometry definition
        wbc : weight and balance class instance
        '''

        #Get Aircraft weight
        W = wbc.estimateWeight(acg)
               
        #Compute the mass from the weight
        m = W/wbc.g

        #Compute the CG location as reference for the Iy calculation
        [MAC,MACc4] = wbc.calculateWingMAC(acg)

        xcg = wbc.calculateWingCenterOfGravity(geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4)

        #Compute the Wing Thickness (thickness constrain from pyGeo)
        #compute the average thickness of the wing for the Inertia Calculation
        wbc.getAverageThickness(acg,thickness)
        
        #Calculate the relative weights of the segments
        wbc.calculateSegmentWeights(acg,W)
    
        #Calculate Moment of Inertia
        [Ix,Iy,Iz]=wbc.calculateWingInertias(acg,xcg)

        #Calculate Frequency and Damping, use MAC as ref chord
        [Wn,DampingRatio]=self.calculateFrequencyAndDamping(averagesol['cmq'],averagesol['clalpha'],
                                            averagesol['cd0'],averagesol['cmalpha'],
                                            averagesol['cmalphadot'],m,Iy,
                                            rho,A,V,MAC)

        #Calculate Dynamic Stability constraint
        val = self.calculateThumbnailMethodConstraint(Wn,DampingRatio)
        

        return val

    def thumbprintDriverpyGeo(self,acg,wbc,geom,averagesol,rho,V,A,surface):
        '''
        run the routines to calculate the CAP values
        '''
        #Get Aircraft weight
        W = wbc.estimateWeight(acg)
               
        #Compute the mass from the weight
        m = W/wbc.g

        #Compute the CG location as reference for the Iy calculation
        [MAC,MACc4] = wbc.calculateWingMAC(acg)

        xcg = wbc.calculateWingCenterOfGravity(geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4)
        #print 'xcg Inertia',xcg,geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4
               
        #Calculate Moment of Inertia
        [Ix,Iz,Iy]=wbc.calculateWingInertiaspyGeo(surface,xcg)#acg,xcg)


        #Calculate the freqency and Damping for the aircraft, use MAC as ref chord
        [Wn,DampingRatio]=self.calculateFrequencyAndDamping(averagesol['cmq'],averagesol['clalpha'],
                                                            averagesol['cd0'],averagesol['cmalpha'],
                                                            averagesol['cmalphadot'],m,Iy,
                                                            rho,A,V,MAC)   
        #Calculate Dynamic Stability constraint
        val = self.calculateThumbnailMethodConstraint(Wn,DampingRatio)

       
        return val


    def CAPDriver(self,acg,wbc,geom,averagesol,rho,V,A,thickness):
        '''
        run the routines to calculate the CAP values
        '''
        #Get Aircraft weight
        W = wbc.estimateWeight(acg)
               
        #Compute the mass from the weight
        m = W/wbc.g

        #Compute the CG location as reference for the Iy calculation
        [MAC,MACc4] = wbc.calculateWingMAC(acg)

        xcg = wbc.calculateWingCenterOfGravity(geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4)
        #print 'xcg Inertia',xcg,geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4
        #Compute the Wing Thickness (thickness constrain from pyGeo)
        #compute the average thickness of the wing for the Inertia Calculation
        wbc.getAverageThickness(acg,thickness)

        #Calculate the relative weights of the segments
        wbc.calculateSegmentWeights(acg,W)
        
        #Calculate Moment of Inertia
        [Ix,Iy,Iz]=wbc.calculateWingInertias(acg,xcg)


        #Calculate the freqency and Damping for the aircraft, use MAC as ref chord
        [Wn,DampingRatio]=self.calculateFrequencyAndDamping(averagesol['cmq'],averagesol['clalpha'],
                                                            averagesol['cd0'],averagesol['cmalpha'],
                                                            averagesol['cmalphadot'],m,Iy/10.0,
                                                            rho,A,V,MAC)

        #Compute the change in g with alpha
        nalpha = self.calculateNAlpha(averagesol['clalpha'],rho,A,V,m,wbc.g)
        
        #Compute the Control Anticipation Parameter    
        CAP = self.calculateCAP(Wn,nalpha)

        #CAP = xcg#nalpha

        return CAP,DampingRatio

    def CAPDriverpyGeo(self,acg,wbc,geom,averagesol,rho,V,A,surface):
        '''
        run the routines to calculate the CAP values
        '''
        #Get Aircraft weight
        W = wbc.estimateWeight(acg)
               
        #Compute the mass from the weight
        m = W/wbc.g

        #Compute the CG location as reference for the Iy calculation
        [MAC,MACc4] = wbc.calculateWingMAC(acg)

        xcg = wbc.calculateWingCenterOfGravity(geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4)
        #print 'xcg Inertia',xcg,geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4
               
        #Calculate Moment of Inertia
        [Ix,Iz,Iy]=wbc.calculateWingInertiaspyGeo(surface,xcg)#acg,xcg)


        #Calculate the freqency and Damping for the aircraft, use MAC as ref chord
        [Wn,DampingRatio]=self.calculateFrequencyAndDamping(averagesol['cmq'],averagesol['clalpha'],
                                                            averagesol['cd0'],averagesol['cmalpha'],
                                                            averagesol['cmalphadot'],m,Iy,
                                                            rho,A,V,MAC)   
        #Compute the change in g with alpha
        nalpha = self.calculateNAlpha(averagesol['clalpha'],rho,A,V,m,wbc.g)
        
        #Compute the Control Anticipation Parameter    
        CAP = self.calculateCAP(Wn,nalpha)
        
        #CAP = xcg#nalpha

        return CAP,DampingRatio                        

    def computeStaticMarginDerivative(self,x,geo,dvFunc,geom,wbc,acg):
        '''
        compute the percent static margin derivative and return
        '''
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j
        SMderiv = {}
        for key in xw.keys():
            SMderiv[key] =[]
            for i in range(len(xw[key])):
                xw[key][i] = xref[key][i]+deltax
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                averagesol = dvFunc(xw)
                #Compute the MAC
                [MAC,C4MAC] = wbc.calculateWingMAC(acg)
                SM = -averagesol['cmalpha']/averagesol['clalpha']

                percentSM = SM/MAC
                SMderiv[key].append(numpy.imag(percentSM)/numpy.imag(deltax))
                xw = copy.deepcopy(xref)
            #end
        #end
        #print 'SM',SM,MAC,-averagesol['cmalpha'],averagesol['clalpha']
        return SMderiv

    def CAPDerivativeDriver(self,x,geo,con,acg,wbc,geom,dvFunc,rho,V,A):
        '''
        compute the derivative of the CAP function...
        '''
        #print 'dvlists',geo.DV_listGlobal,geo.DV_namesGlobal,'local', geo.DV_listLocal,geo.DV_namesLocal
        dCondx = con.getThicknessSensitivity(geo,'con')
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j

        CAPderiv = {}
        Dampderiv = {}
        #geo.complex = True
        for key in xw.keys():
            CAPderiv[key] =[]
            Dampderiv[key] =[]
            for i in range(len(xw[key])):
                #print 'key',key
                xw[key][i] = xref[key][i]+deltax
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                averagesol = dvFunc(xw)
                        
                con.setCoordinates(geo.update('con'))
                
                thick_con = con.getThicknessConstraints()
                thick_con_c = thick_con.astype('D')
                if key in geo.DV_namesGlobal.keys():
                    for j in range(len(thick_con)):
                        thick_con_c[j]+=dCondx[j,geo.DV_namesGlobal[key]]*deltax
                    #end
                
                #thick_con = numpy.real(thick_con)
                    #print 'thick',type(thick_con_c),thick_con_c,dCondx[:,geo.DV_namesGlobal[key]]*deltax
                #end
                CAP,Damp = self.CAPDriver(acg,wbc,geom,averagesol,
                                            rho,V,A,
                                            thick_con_c.reshape(con.thickConSizes[0]))          
   
                CAPderiv[key].append(numpy.imag(CAP)/numpy.imag(deltax))
                Dampderiv[key].append(numpy.imag(Damp)/numpy.imag(deltax))
                xw = copy.deepcopy(xref)
            #end
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')
        geo.update('con')

        return CAPderiv,Dampderiv

    def CAPDerivativeDriverpyGeo(self,x,geo,con,acg,wbc,geom,dvFunc,rho,V,A,surface,surfInst):
        '''
        compute the derivative of the CAP function...
        '''
        #print 'dvlists',geo.DV_listGlobal,geo.DV_namesGlobal,'local', geo.DV_listLocal,geo.DV_namesLocal
#        dCondx = con.getThicknessSensitivity(geo,'con')
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        wbc.setOption('dtype','D')
        geo.complex = True
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j

        CAPderiv = {}
        Dampderiv = {}
        #geo.complex = True
        for key in xw.keys():
            CAPderiv[key] =[]
            Dampderiv[key] =[]
            for i in range(len(xw[key])):
                #print 'key',key
                xw[key][i] = xref[key][i]+deltax
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                averagesol = dvFunc(xw)
                surface.setSurfaceCoordinates(surfInst,geo.update('X_coord'))#.reshape([wing.nSurf,nv,nu,3])
                surface.setCellCentroidCoordinates(surfInst,geo.update('Cent_coord'))
                #con.setCoordinates(geo.update('con'))
                
                #thick_con = con.getThicknessConstraints()
                #thick_con_c = thick_con.astype('D')
                #if key in geo.DV_namesGlobal.keys():
                #    for j in range(len(thick_con)):
                #        thick_con_c[j]+=dCondx[j,geo.DV_namesGlobal[key]]*deltax
                #    #end
               
                #thick_con = numpy.real(thick_con)
                    #print 'thick',type(thick_con_c),thick_con_c,dCondx[:,geo.DV_namesGlobal[key]]*deltax
                #end
                CAP,Damp = self.CAPDriverpyGeo(acg,wbc,geom,averagesol,
                                            rho,V,A,
                                            surfInst)              
          
   
                CAPderiv[key].append(numpy.imag(CAP)/numpy.imag(deltax))
                Dampderiv[key].append(numpy.imag(Damp)/numpy.imag(deltax))
                xw = copy.deepcopy(xref)
            #end
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')
        wbc.setOption('dtype','d')
        geo.complex = False
        geo.update('X_coord')
        geo.update('Cent_coord')

        return CAPderiv,Dampderiv

    def CAPDerivativeDriverpyGeoFD(self,x,geo,con,acg,wbc,geom,dvFunc,rho,V,A,surface,surfInst):
        '''
        compute the derivative of the CAP function...
        '''
        #print 'dvlists',geo.DV_listGlobal,geo.DV_namesGlobal,'local', geo.DV_listLocal,geo.DV_namesLocal
#        dCondx = con.getThicknessSensitivity(geo,'con')
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('d')
        #endfor
        wbc.setOption('dtype','d')
        geo.complex = False
        xref = copy.deepcopy(xw)
        deltax = 1.0e-7
        geo.setValues(xw,scaled=True)
        geo.update('wing')
        averagesol = dvFunc(xw)
        surface.setSurfaceCoordinates(surfInst,geo.update('X_coord'))
        surface.setCellCentroidCoordinates(surfInst,geo.update('Cent_coord'))
        CAPref,Dampref = self.CAPDriverpyGeo(acg,wbc,geom,averagesol,
                                             rho,V,A,
                                             surfInst)   
        CAPderiv = {}
        Dampderiv = {}
        #geo.complex = True
        for key in xw.keys():
            CAPderiv[key] =[]
            Dampderiv[key] =[]
            for i in range(len(xw[key])):
                #print 'key',key
                xw[key][i] = xref[key][i]+deltax
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                averagesol = dvFunc(xw)
                surface.setSurfaceCoordinates(surfInst,geo.update('X_coord'))
                surface.setCellCentroidCoordinates(surfInst,geo.update('Cent_coord'))
                #con.setCoordinates(geo.update('con'))
                
                #thick_con = con.getThicknessConstraints()
                #thick_con_c = thick_con.astype('D')
                #if key in geo.DV_namesGlobal.keys():
                #    for j in range(len(thick_con)):
                #        thick_con_c[j]+=dCondx[j,geo.DV_namesGlobal[key]]*deltax
                #    #end
               
                #thick_con = numpy.real(thick_con)
                    #print 'thick',type(thick_con_c),thick_con_c,dCondx[:,geo.DV_namesGlobal[key]]*deltax
                #end
                CAP,Damp = self.CAPDriverpyGeo(acg,wbc,geom,averagesol,
                                            rho,V,A,
                                            surfInst)              
          
   
                CAPderiv[key].append((CAP-CAPref)/(deltax))
                Dampderiv[key].append((Damp-Dampref)/(deltax))
                xw = copy.deepcopy(xref)
            #end
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')
        wbc.setOption('dtype','d')
        geo.complex = False
        geo.update('X_coord')
        geo.update('Cent_coord')

        return CAPderiv,Dampderiv

    def thumbprintDerivativeDriver(self,x,geo,con,acg,wbc,geom,dvFunc,rho,V,A):
        '''
        compute the derivative of the thumbprint function...
        '''
        #print 'dvlists',geo.DV_listGlobal,geo.DV_namesGlobal,'local', geo.DV_listLocal,geo.DV_namesLocal
        dCondx = con.getThicknessSensitivity(geo,'con')
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j

        TPderiv = {}
        
        #geo.complex = True
        for key in xw.keys():
            TPderiv[key] =[]
            for i in range(len(xw[key])):
                #print 'key',key
                xw[key][i] = xref[key][i]+deltax
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                averagesol = dvFunc(xw)
                        
                con.setCoordinates(geo.update('con'))
                
                thick_con = con.getThicknessConstraints()
                thick_con_c = thick_con.astype('D')
                if key in geo.DV_namesGlobal.keys():
                    for j in range(len(thick_con)):
                        thick_con_c[j]+=dCondx[j,geo.DV_namesGlobal[key]]*deltax
                    #end
                
                #thick_con = numpy.real(thick_con)
                    #print 'thick',type(thick_con_c),thick_con_c,dCondx[:,geo.DV_namesGlobal[key]]*deltax
                #end
                val = self.thumbprintDriver(acg,wbc,geom,averagesol,
                                            rho,V,A,
                                            thick_con_c.reshape(con.thickConSizes[0]))              
          
   
                TPderiv[key].append(numpy.imag(val)/numpy.imag(deltax))
                xw = copy.deepcopy(xref)
            #end
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')
        geo.update('con')

        return TPderiv
    
    def thumbprintDerivativeDriverpyGeo(self,x,geo,con,acg,wbc,geom,dvFunc,rho,V,A,surface,surfInst):
        '''
        compute the derivative of the thumbprint function...
        '''
        #print 'dvlists',geo.DV_listGlobal,geo.DV_namesGlobal,'local', geo.DV_listLocal,geo.DV_namesLocal
#        dCondx = con.getThicknessSensitivity(geo,'con')
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        wbc.setOption('dtype','D')
        geo.complex = True
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j

        TPderiv = {}
        
        #geo.complex = True
        for key in xw.keys():
            TPderiv[key] =[]
            for i in range(len(xw[key])):
                #print 'key',key
                xw[key][i] = xref[key][i]+deltax
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                averagesol = dvFunc(xw)
                surface.setSurfaceCoordinates(surfInst,geo.update('X_coord'))#.reshape([wing.nSurf,nv,nu,3])
                surface.setCellCentroidCoordinates(surfInst,geo.update('Cent_coord'))
                
                val  = self.thumbprintDriverpyGeo(acg,wbc,geom,averagesol,rho,V,A,surfInst)
                   
                TPderiv[key].append(numpy.imag(val)/numpy.imag(deltax))
                
                xw = copy.deepcopy(xref)
            #end
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')
        wbc.setOption('dtype','d')
        geo.complex = False
        geo.update('X_coord')
        geo.update('Cent_coord')

        return TPderiv
#==============================================================================
# PERFOMANCE Analysis Test
#==============================================================================
if __name__ == '__main__':
    
    # Test ADflow
    print('Testing ...')
    acp = PERFORMANCE()
    print(acp)
