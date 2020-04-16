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
import os, sys,copy
import pdb
import time
#from cmath import pi,cos,tan,sin
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from mdo_import_helper import *
exec(import_modules('pyBase_class'))
exec(import_modules('pyGeometry_liftingsurface_c'))
#exec(import_modules('pyGeometry_liftingsurface'))
from pyGeometry_complex import catan,cabs,cmin,cmax
# =============================================================================
# Misc Definitions
# =============================================================================

class WEIGHTANDBALANCE(Base):

    '''
    Implementation of basic preliminary level CG and inertia
    estimation methods.
    '''

    def __init__(self, *args, **kwargs):
		
        '''
        Weight and Balance Class Initialization
        
        '''
        name = 'WEIGHTANDBALANCE'
        category = 'CG and I estimation'
        def_opts = {'dtype':[str,'d'],
                    't':[float,0.001],
                    'rho':[float,2600],
                    'inertiaModifier':[float,1.5]
            }
        informs = {
            }

        #run the Generic surface defintion
        Base.__init__(self, name, category, def_opts, informs, *args, **kwargs)

        if 'g' in kwargs:
            self.g = kwargs['g']
        else:
            self.g = 9.80 #m/s**2
        #endif
        if 'Units' in kwargs:
            self.Units = kwargs['Units']
        else:
            self.Units ='metric'
        #end
        
        if 'AC_weight' in kwargs:
            self.AC_weight = kwargs['AC_weight']
        else:
            self.AC_weight = 436000 #Newtons. approx a global express.
        #endif
        return

    def estimateWeight(self,acg):
        '''
        Fill in a weight estimation function.
        For now just return the internal value
        acg - aircraft geometry defn.
        '''

        return self.AC_weight

    def calculateWingMAC(self,acg):
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
                    #Determine the number of segments that make up the wing
                    nseg=len(acg[i])

                    #Compute total area and span for this lifting surface
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
                    counter = 0
                    #loop over wing segments
                    for j in range(nseg):
                        #Segments are linear therfore single a
                        #trapezoid gives exact answer

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
                        C_root = (2.0*Area)/(cabs(Span)*(1.0+Taper))

                        #print 'croot',C_root
                        # Tip Chord
                        # C_tip = Taper*C_root
                        #print 'ctip',C_tip
                       
                        #new xC4 computation
                        temp1 = (1+2*Taper)*((ytLE-yrLE)/sumSpan)\
                                *numpy.tan(SweepLE)
                        #print 'temp1',(1+2*Taper),((ytLE-yrLE)/sumSpan),numpy.tan(SweepLE)
                        temp2 = 3*(1+Taper)*sweepsum
                        
                        temp3 = ((ytLE-yrLE)/sumSpan)*C_root
                        #print 'temp',temp1,temp2,temp3
                        xC4sum = xC4sum+(temp3*(temp2+temp1))
                        
                        
                        sweepsum = sweepsum+((ytLE-yrLE)/sumSpan)*numpy.tan(SweepLE)
                        #new computation...
                        numerator = C_root**2*(1+Taper+Taper**2)*(ytLE-yrLE)/sumSpan
                        #print 'numpart',C_root,C_root**2,(1+Taper+Taper**2),(ytLE-yrLE)/sumSpan
                        denomenator = C_root*(1+Taper)*(ytLE-yrLE)/sumSpan
                        SumNumerator= SumNumerator+numerator
                        SumDenomenator=SumDenomenator+denomenator
                        #print 'num',SumNumerator,numerator
                        #print 'den',SumDenomenator,denomenator
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
        #print 'MAC',MAC
        #print 'xc4',MACc4
        return MAC,MACc4

    def calculateWingCenterOfGravity(self,forwardSparPercent,rearSparPercent,CGPercentage,MAC,MACc4):
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

    def wingCGDriver(self,acg,geom):
        '''
        call routines to calculate wing CG
        '''
        #Calculate Wing CG
        [MAC,C4MAC] = self.calculateWingMAC(acg)
        ## print 'MAC',MAC,C4MAC

        xcg = self.calculateWingCenterOfGravity(geom.ForeSparPercent,
                                                    geom.RearSparPercent,
                                                    geom.CGPercent,MAC,C4MAC)
        #print 'xcg Inertia',xcg,geom.ForeSparPercent,geom.RearSparPercent,geom.CGPercent,MAC,MACc4
        return xcg


    def calculateWingInertias(self,acg,xcg):
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
                    #Determine the number of segments that make up the wing
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
                        #print 'span',Span
                        Dihedral = acg[i][j].Dihedral*(numpy.pi/180)
                        Area = acg[i][j].Area
                        Taper = acg[i][j].Taper
                        AR = Span**2/Area # Aspect Ratio
                        SweepLE = acg[i][j].SweepLE*(numpy.pi/180)
                        SweepTE = catan(numpy.tan(SweepLE)-((4.0*(1.0/1.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)))
                        #print 'sweeps',SweepLE*(180.0/numpy.pi),SweepTE*(180.0/numpy.pi)
                        #print 'SweepTE',numpy.tan(SweepLE),(4.0*(1.0/1.0)),(2.0*AR),((1.0-Taper)/(1.0+Taper)),numpy.tan(SweepLE)-((4.0*(1.0/1.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)),catan(numpy.tan(SweepLE)-((4.0*(1.0/1.0))/(2.0*AR))*((1.0-Taper)/(1.0+Taper)))
                        #Compute tip span location
                        ytLE = yrLE + Span*numpy.cos(Dihedral)
                        #print 'ytle',ytLE
                        # Root Chord
                        #c_abs!!!!
                        C_root = (2.0*Area)/(cabs(Span)*(1.0+Taper))
                        #Root thickness
                        t_root = tc_root*C_root
                        #print 'croot',C_root,tc_root
                        # Tip Chord
                        C_tip = Taper*C_root
                        #tip thickness
                        t_tip = tc_tip*C_tip

                        #Weight
                        if self.Units=='metric':
                            W = acg[i][j].Weight/self.g
                            #print 'mass',W
                        else:
                            W = acg[i][j].Weight
                        #end
                        #print 't,w',t_root,t_tip,W
                        #Volume???
                        V = Span*(t_root*(C_root+(Span/2.0)*(numpy.tan(SweepTE)-numpy.tan(SweepLE)))-((t_root-t_tip)*((C_root/2.0)+(Span/3.0)*(numpy.tan(SweepTE)-numpy.tan(SweepLE)))))
                        #print 'V',V
                        #I1x = (t_root-t_tip)*(C_root/4.0+(Span*tan(SweepTE)/5.0)-Span*tan(SweepLE)/5.0)#+(t_root*(C_root/3.0+(Span*tan(SweepTE)/4.0)-Span*tan(SweepLE)/4.0)))
                        #print 'I1x1',I1x
                        #I1x =(t_root*(C_root/3.0+(Span*tan(SweepTE)/4.0)-Span*tan(SweepLE)/4.0))
                        #print 'I1x2',I1x
                        I1x = (W*Span**3/(3*V))*(((t_root-t_tip)*(C_root/4.0+(Span*numpy.tan(SweepTE)/5.0)-Span*numpy.tan(SweepLE)/5.0))+(t_root*(C_root/3.0+(Span*numpy.tan(SweepTE)/4.0)-Span*numpy.tan(SweepLE)/4.0)))
                        #print 'Check',(W*Span**3/(3*V)),((t_root-t_tip)*(C_root/4.0+(Span*numpy.tan(SweepTE)/5.0)-Span*numpy.tan(SweepLE)/5.0)),(t_root*(C_root/3.0+(Span*numpy.tan(SweepTE)/4.0)-Span*numpy.tan(SweepLE)/4.0)),(((t_root-t_tip)*(C_root/4.0+(Span*numpy.tan(SweepTE)/5.0)-Span*numpy.tan(SweepLE)/5.0))+(t_root*(C_root/3.0+(Span*numpy.tan(SweepTE)/4.0)-Span*numpy.tan(SweepLE)/4.0))),(W*Span**3/(3*V))*(((t_root-t_tip)*(C_root/4.0+(Span*numpy.tan(SweepTE)/5.0)-Span*numpy.tan(SweepLE)/5.0))+(t_root*(C_root/3.0+(Span*numpy.tan(SweepTE)/4.0)-Span*numpy.tan(SweepLE)/4.0)))
                        #I1x1 =27028033000
                        #Note: changed(W*Span**3/(V)) to (W*Span**3/(3*V))
                        #to match published results
                        #print 'I1x',I1x#,I1x1,I1x/I1x1
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
                        #b = copy.deepcopy(a)
                        #print 'asort',b.sort(),b#a
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
                        Xs4 = xrLE-xcg
                        
                        Ys = acg[i][j].y_Centroid
                        Ys_dot = acg[i][j].y_Centroid*numpy.cos(Dihedral)
                        #print 'ydot',Ys_dot
                        Ysoff = yrLE
                        #print 'ysoff',Ysoff
                        Zs1 = zrLE
                        #print 'zrle',zrLE
                        Zs3 = acg[i][j].z_Centroid
                        #print 'zs3',Zs3
                        Zs4 = acg[i][j].z_Centroid
                        #print 'centroids',Xs,Ys,Zs3
                        
                        Ix = Ix+I1x-W*(Ys_dot**2)-W*(Zs3**2)+W*(Ys_dot+Ysoff)**2 + W*(Zs3+Zs1)**2
                        #print 'w',W,(Ys_dot),W*(Ys_dot**2)
                        #Ix = Ix-W*(Ys_dot**2)#-W*(Zs3**2)+W*(Ys_dot+Ysoff)**2 + W*(Zs3+Zs1)**2
                        #print 'ix',Ix,I1x
                        Iy = Iy+I1y-W*(Xs**2)-W*(Zs3**2)+W*(Xs+Xs4)**2+W*(Zs3+Zs1)**2
                        
                        Iz = Iz+I1z-W*(Xs**2+Ys_dot**2)+W*(Xs+Xs4)**2+W*(Ys_dot+Ysoff)**2
                        #Ixz = ...
                    
                    #endfor
                #endif
            #endif
        #endfor

   
        return Ix,Iy,Iz

    def calculateWingInertiaspyGeo(self,surface,xcg):
        '''
        Calculates the Wing mass moments of inertia from the pyGeo Surface
        
        wing - instance of pyGeo Surface
        
        note:
        x is chord direction
        y is vertical direction
        z is span direction
        t is thickness in meters
        rho is density of material in kg/m^3
        '''
        dtype = self.getOption('dtype')#'d'
        #a = time.time()
        X = surface.Xs
        Xc = surface.Xc
        nu = X.shape[2]#10#5#10
        nv = X.shape[1]#10#3#5
        nSurf = X.shape[0]
        #u = numpy.linspace(0,1,nu)
        #v = numpy.linspace(0,1,nv)
        t = self.getOption('t')#0.001 #thickness in meters
        rho = self.getOption('rho')#2700 #density of material in kg/m^3
        #[V,U] = numpy.meshgrid(u,v)
        #print 'u,v',u,v,'U',U,V
        v1 = numpy.zeros([3],dtype)
        v2 = numpy.zeros([3],dtype)
        s = numpy.zeros([nv-1,nu-1,3],dtype)
        r2 = numpy.zeros([nv-1,nu-1,3],dtype)
        inertia = numpy.zeros([nv-1,nu-1,3],dtype)
        sTot =numpy.zeros([nv-1,nu-1],dtype)
        Vol = numpy.zeros([nv-1,nu-1],dtype)
        mass = numpy.zeros([nv-1,nu-1],dtype)
        Area = numpy.zeros([3,nSurf],dtype)
        SurfaceInertia = numpy.zeros([3,nSurf],dtype)
        totalMass = 0
        totalInertia = numpy.zeros([3],dtype)
       #  Uc = numpy.zeros([nv-1,nu-1],dtype)
#         Vc = numpy.zeros([nv-1,nu-1],dtype)
#         for i in range(nv-1):
#             for j in range(nu-1):
#                 #print 'UV',V[i,j],V[i,j+1],U[i,j],U[i+1,j]
#                 Uc[i,j] = U[i,j]+0.5*(U[i+1,j]-U[i,j])
#                 #print 'Uc',Uc[i,j],i,j
#                 Vc[i,j] = V[i,j]+0.5*(V[i,j+1]-V[i,j])
#                 #print 'Vc',Vc[i,j]
#             #end
#         #end
        Xcg = [xcg,0.0,0.0]
        for i in range(nSurf):
            tempX =  X[i,:,:,:]
            tempXc = Xc[i,:,:,:]
            sz = tempX.shape
            for j in range(sz[0]-1):
                for k in range(sz[1]-1):
                    #print temp[j,k],j,k
                    v1[0] = tempX[j+1,k+1,0]-tempX[j,k,0]
                    v1[1] = tempX[j+1,k+1,1]-tempX[j,k,1]
                    v1[2] = tempX[j+1,k+1,2]-tempX[j,k,2]
                    #print 'v1',v1
                    v2[0] = tempX[j,k+1,0]-tempX[j+1,k,0]
                    v2[1] = tempX[j,k+1,1]-tempX[j+1,k,1]
                    v2[2] = tempX[j,k+1,2]-tempX[j+1,k,2]
                    s[j,k,0] = 0.5*(v1[1]*v2[2] - v1[2]*v2[1])
                    s[j,k,1] = 0.5*(v1[2]*v2[0] - v1[0]*v2[2])
                    s[j,k,2] = 0.5*(v1[0]*v2[1] - v1[1]*v2[0])
                    #s[j,k,:] = 0.5* numpy.cross(v1,v2)
                    sTot[j,k] = s[j,k,0]+s[j,k,1]+s[j,k,2]
                    #sTot[j,k] = numpy.sum(s[j,k,:])
                    Vol[j,k] = cabs(sTot[j,k])*t
                    mass[j,k] = Vol[j,k]*rho
                    totalMass = totalMass+mass[j,k]
                    dx = tempXc[j,k,0]-Xcg[0]
                    dy = tempXc[j,k,1]-Xcg[1]
                    dz = tempXc[j,k,2]-Xcg[2]
                    #print 'dist.',dx,dy,dz
                    r2[j,k,0] = (dy**2+dz**2)#aboutx axis
                    r2[j,k,1] = (dx**2+dz**2)# about y axis
                    r2[j,k,2] = (dy**2+dx**2)# about z axis
                    #print 'r2',r2[j,k,:]
                    inertia[j,k,0] = mass[j,k]*r2[j,k,0]#Ixx
                    inertia[j,k,1] = mass[j,k]*r2[j,k,1]#Iyy
                    inertia[j,k,2] = mass[j,k]*r2[j,k,2]#Izz
                    #inertia[j,k,:] = mass[j,k]*r2[j,k,:]
                    Area[:,i] = Area[:,i]+cabs(s[j,k,:])
                    SurfaceInertia[:,i] = SurfaceInertia[:,i]+inertia[j,k,:]
                #end
            #end
        #end     
        
        #print 'mass2',totalMass
        #print 'area2',Area
        #print SurfaceInertia
        self.totalMass = totalMass
        for i in range(nSurf):
            totalInertia = totalInertia+SurfaceInertia[:,i]
        #end
        temp = self.getOption('inertiaModifier')
        
        totalInertia = totalInertia*temp
        #print totalInertia
        #b = time.time()
        #print 'total time:',b-a   
        
        return totalInertia

    def computeRootBendingMoment(self,sol,ref,geom,liftIndex):
        #sum up the moments about the root elastic center to determine the effect of sweep on the moment
        momx = sol['mx']
        momy = sol['my']
        momz = sol['mz']
        fx = sol['fx']
        fy = sol['fy']
        fz = sol['fz']
        if liftIndex == 2:
            #z out wing sum momentx,momentz
            elasticMomentx = momx + fy*(geom.zRootec-ref.zref)-fz*(geom.yRootec-ref.yref)
            elasticMomentz = momz - fy*(geom.xRootec-ref.xref)+fx*(geom.yRootec-ref.yref)
                         
            BendingMoment = numpy.sqrt(elasticMomentx**2+elasticMomentz**2)

        elif liftIndex == 3:
            #y out wing sum momentx,momenty
            elasticMomentx = momx + fz*(geom.yRootec-ref.yref)+fy*(geom.zRootec-ref.zref)
            elasticMomenty = momy + fz*(geom.xRootec-ref.xref)+fx*(geom.zRootec-ref.zref)
                         
            BendingMoment = numpy.sqrt(elasticMomentx**2+elasticMomenty**2)
            
        return BendingMoment

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
##                     C_root = (2.0*Area)/(cabs(Span)*(1.0+Taper))
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


    def getAverageThickness(self,acg,thickness):
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
                        if (numpy.mod(len(thickness),nseg)==0):
                            #print 'thickness ok',numpy.mod(len(thickness),nseg)
                            multiple = len(thickness)/nseg
                            #print multiple
                        else:
                            print('Number of spanwise thicknesses must be a multiple of %d'%(nseg))
                        #endif
                        ##  if j ==0:
##                         tipIndex = rootIndex+acg[i][j].surface_SW_segments
##                     else:
##                         tipIndex = rootIndex+acg[i][j].surface_SW_segments-1
##                     #endif

                        Span = acg[i][j].Span
                        Area = acg[i][j].Area
                        Taper = acg[i][j].Taper
                        
                        #Chord Lengths
                        C_root = (2.0*Area)/(cabs(Span)*(1.0+Taper))
                        #print 'croot',C_root
                        C_tip = Taper*C_root
                        #root thickness
                        t_root = numpy.mean(thickness[j][:])
                        #print thickness[j][:]
                        #tip  thickness
                        t_tip = numpy.mean(thickness[j+multiple-1][:])
                        #print 't_tip',t_root,t_tip
                        tc_root = t_root/C_root
                    
                        tc_tip = t_tip/C_tip
                        
                        acg[i][j].root_Thickness_act  = tc_root#thickness to chord ratio...
                        acg[i][j].tip_Thickness_act =tc_tip
                        #print 'tc',tc_root,tc_tip
                        #rootIndex = tipIndex
                    #endfor
                #endif
            #endif
        #endfor
        
        return

    def calculateSegmentWeights(self,acg,Weight):
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
                        #print 'tc_act',tc_root,tc_tip
                        #Wing properties
                        Span = acg[i][j].Span
                        Area = acg[i][j].Area
                        Taper = acg[i][j].Taper
                        
                        #Chord Lengths
                        C_root = (2.0*Area)/(cabs(Span)*(1.0+Taper))
                        C_tip = Taper*C_root
                        
                        #Segment thicknesses
                        t_root = C_root* tc_root                    
                        t_tip = C_tip * tc_tip
                        #print 't_act',t_root,t_tip
                        #Approx Volume
                        V = Area*(t_root+t_tip)/2.0
                        #print 'V',V
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
                        #print 'weights',nseg,Weight,acg[i][j].volumeWeight,
                        acg[i][j].Weight = Weight*(acg[i][j].volumeWeight/TotalVolume)
                        #print 'acg weight',acg[i][j].Weight,acg[i][j].volumeWeight,TotalVolume,i,j
                    #endfor
                #endif
            #endif
        #endfor
                    
        return


    def calculateWingCenterOfGravityDerivatives(self,x,geo,acg,geom):
        '''
        Run the complex step derivatives required for the CG derivatives.
        '''
        #print 'DVs',x
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j
        CGderiv = {}
        #geo.complex = True
        for key in xw.keys():
            CGderiv[key] =[]
            for i in range(len(xw[key])):
                #if comm.rank == 0: print 'test',xw[key][i],xref[key][i],deltax,key,i
                xw[key][i] = xref[key][i]+deltax
                #if comm.rank == 0:print 'xw',xw
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                
                [MAC,C4MAC] = self.calculateWingMAC(acg)
                #print 'incgderiv',geom.CGPercent
                xcg = self.calculateWingCenterOfGravity(geom.ForeSparPercent,
                                                        geom.RearSparPercent,
                                                        geom.CGPercent,MAC,C4MAC)
                #if comm.rank == 0:print 'cgderiv',CGderiv,xcg
                CGderiv[key].append(numpy.imag(xcg)/numpy.imag(deltax))
                xw = copy.deepcopy(xref)
                #if comm.rank == 0:print 'xw2',xw
            #end
            
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')

        return CGderiv

    def computeWingMACDerivatives(self,x,geo,acg,geom):
        '''
        Run the complex step derivatives required for the CG derivatives.
        '''
        #print 'DVs',x
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j
        MACDeriv = {}
        MACc4Deriv = {}
        #geo.complex = True
        for key in xw.keys():
            MACDeriv[key] =[]
            MACc4Deriv[key] =[]
            for i in range(len(xw[key])):
                #if comm.rank == 0: print 'test',xw[key][i],xref[key][i],deltax,key,i
                xw[key][i] = xref[key][i]+deltax
                #if comm.rank == 0:print 'xw',xw
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                
                [MAC,C4MAC] = self.calculateWingMAC(acg)
                
                #if comm.rank == 0:print 'cgderiv',CGderiv,xcg
                MACDeriv[key].append(numpy.imag(MAC)/numpy.imag(deltax))
                MACc4Deriv[key].append(numpy.imag(C4MAC)/numpy.imag(deltax))
                xw = copy.deepcopy(xref)
                #if comm.rank == 0:print 'xw2',xw
            #end
            
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')

        return MACDeriv,MACc4Deriv

    def calculateBendingMomentDerivatives(self,x,geo,dvFunc,ref,geom,flow):
        '''
        Run the complex step derivatives required for the CG derivatives.
        '''
        #print 'DVs',x
        xw = copy.deepcopy(x)
        for i in xw.keys():
            xw[i] = numpy.atleast_1d(xw[i]).astype('D')
        #endfor
        xref = copy.deepcopy(xw)
        deltax = 1.0e-40j
        BMderiv = {}
        #geo.complex = True
        for key in xw.keys():
            BMderiv[key] =[]
            for i in range(len(xw[key])):
                #if comm.rank == 0: print 'test',xw[key][i],xref[key][i],deltax,key,i
                xw[key][i] = xref[key][i]+deltax
                #if comm.rank == 0:print 'xw',xw
                geo.setValues(xw,scaled=True)
                geo.update('wing')
                averagesol = dvFunc(xw)
                ref.xref = xw['xref']
                
                BM = self.computeRootBendingMoment(averagesol,ref,geom,flow.liftIndex)
               
                #if comm.rank == 0:print 'cgderiv',CGderiv,xcg
                BMderiv[key].append(numpy.imag(BM)/numpy.imag(deltax))
                xw = copy.deepcopy(xref)
                #if comm.rank == 0:print 'xw2',xw
            #end
            
        #endfor
        #geo.complex = False
        geo.setValues(x,scaled=True)
        geo.update('wing')

        return BMderiv

    def _on_setOption(self, name, value):

        '''
        Get Optimizer Option Value (Optimizer Specific Routine)

        Documentation last updated:  May. 21, 2008 - Ruben E. Perez
        '''

        pass

    def _on_getOption(self, name):

        '''
        Get Optimizer Option Value (Optimizer Specific Routine)

        Documentation last updated:  May. 21, 2008 - Ruben E. Perez
        '''

        pass

#==============================================================================
# weight and balance Analysis Test
#==============================================================================
if __name__ == '__main__':
    
    # Test ADflow
    print('Testing ...')
    wbc = WEIGHTANDBALANCE()
    print(wbc)
