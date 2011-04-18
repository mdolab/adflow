#!/usr/local/bin/python
'''
Aerodynamics

Holds the Aerodynamic Analysis Classes (base and inherited) 
used in the aircraft conceptual design implementation framework.

Copyright (c) 2007 by Dr. Ruben E. Perez
All rights reserved. Not to be used for commercial purposes.
Revision: 1.1   $Date: 06/10/2007 21:00$


History
-------
	v. 1.0 - Initial Class Creation
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, math, string
import pdb

# =============================================================================
# External Python modules
# =============================================================================
import numpy
from pylab import *
import matplotlib.axes3d

# =============================================================================
# Extension modules
# =============================================================================


# =============================================================================
# 
# =============================================================================
def getRe(Mach,Atmos,lref):
	
	'''
	Calculate Reynolds Number
	
	Inputs:
	-------
	Mach  - Flight Mach Number
	Atmos - Atmospheric Properties at given Altitude
		.rho - Density 
		.a   - Speed of Sound
		.Te  - Temperature
		.mu  - Viscocity
	lref  - Reference Length
		
	Outputs:
	--------
	Re   - Reynolds Number
	'''
	
	# Inputs
	rho = Atmos.Density
	a = Atmos.Sound_Speed
	mu = Atmos.Dynamic_Viscosity
	
	# Reynolds Number
	Re = (Mach*a)*lref*rho/mu
	
	return Re
	


# =============================================================================
# 
# =============================================================================
def getCf(Mach,Atmos,lref,xT,k):
	
	'''
	Calculate Friction Coefficient
	
	Inputs:
	-------
	Mach  - Flight Mach Number
	Atmos - Atmospheric Properties at given Altitude
		.rho - Density 
		.a   - Speed of Sound
		.Te  - Temperature
		.mu  - Viscocity
	lref  - Reference Length
	xt    - Point of Transition (x/l)
	k     - Skin Roughness Value (k=0 No Skin Roughness Calculation)
	
	Outputs:
	--------
	Cf   - Friction Coefficient
	Re   - Reynolds Number
	'''
	
	# Constants
	TwTaw = 1
	
	# Reynolds Number
	Re = getRe(Mach,Atmos,lref)
	
	# Skin Roughness Effect (Reynolds Cutoff Method - DATCOM 4.1.5.1-26)
	if (k != 0):
		if (Mach < 0.9):
			Recut = 38.21*(lref/k)**1.053
		else:
			Recut = 44.62*(lref/k)**1.053*Mach**1.16
		#end
		Re = min([Re,Recut])
	#end
	
	# Re Transition
	Rex = Re*xt    
	if (Rex <= 0):  
		Rex = 0.0001
	#end
	
	# theta = (0.671*xt)/math.sqrt(Rex)
	# xeff = (27.78*theta*Re**0.2)**1.25
	# Rext = Re*((1-xt)+xeff)
	
	# Turbulent and Laminar Friction Coefficients
	Cfturb = getCfturb(Re,Mach,TwTaw,Te)
	Cfxturb = getCfturb(Rex,Mach,TwTaw,Te)
	Cfxlam = getCflam(Rex,Mach,TwTaw,Te)
	
	# Cfturb = getCfturb(Rext,Mach,TwTaw,Te)
	# Cfstart = getCfturb(Re*xeff,Mach,TwTaw,Te)
	# Cflam = getCflam(Rex,Mach,TwTaw,Te)
	
	# Partial Laminar Flow Run (Schlicting's Formulation)
	Cf = Cfturb - xt*(Cfxturb-Cfxlam)
	# Cf = (Cflam*xt + Cfturb*((1-xt)+xeff)) - Cfstart*xeff

	return Cf, Re


	def getCflam(Rex,Mach,TwTaw,Te):
		
		'''
		Calculate Laminar Skin Friction Coefficient
		
		Inputs:
		-------
		Rex   - Reynolds Number
		Mach  - Mach Number
		TwTaw - Tw/Taw, the wall temperature ratio
		
		Outputs:
		--------
		Cf    - Laminar Friction Coefficient
		'''
		
		# Constants
		Pr = 0.72
		r = Pr**(1/2)
		Gamma = 1.4
		Ks = 198.594
		
		# Incompressible Laminar Friction Coefficient
		Cfincomp = 1.328/math.sqrt(Rex)
		
		# Compressibility Effects
		TwTe = TwTaw*(1.0 + r*((Gamma-1.0)/2.0)*Mach**2.0)
		TrTe = 1.0 + 0.035*Mach**2.0 + 0.45*(TwTe-1.0)
		# TrTe = 1.0 + 0.039*Mach**2.0 + 0.50*(TwTe-1.0)
		RrRe = ((Te+Ks)/(TrTe*Te+Ks))*(1/TrTe)**1.5		# Sutherland's Law
		Cf = Cfincomp*(1.0/(TrTe*(RrRe)^0.2))
		# Cf = Cfincomp*(1.0 + 0.104*Mach**2)**-0.773
		# Cf = Cfincomp*(1.0 + 0.144*Mach**2)**-0.650
		
		# Cstar = ((1.0 + Ks/Te)/(TrTe+Ks/Te))*(TrTe)**0.5
		# Cf = 2.0*(0.664*math.sqrt(Cstar)/math.sqrt(Rex))
		
		return Cf
		
		
	def getCfturb(Rex,Mach,TwTaw,Te):
		
		'''
		Calculate Turbulent Skin Friction Coefficient
		
		Inputs:
		-------
		Rex   - Reynolds Number
		Mach  - Mach Number
		TwTaw - Tw/Taw, the wall temperature ratio
		
		Outputs:
		--------
		Cf    - Laminar Friction Coefficient
		
		References:
		-----------
		- Van Driest II Method, NASA TN D-6945
		'''
		
		# Constants
		Gamma = 1.4
		r = 0.88
		
		# Van Driest II Method, NASA TN D-6945
		m = ((Gamma-1.0)/2.0)*Mach**2
		TawTe = 1.0 + r*m
		F = TwTaw*TawTe
		Tw = F*Te
		A = (r*m/F)**(1.0/2.0)
		B = (1.0 + r*m - F)/F
		Alpha = (2.0*A**2.0 - B)/((4.0*A**2.0 + B**2.0)^(1.0/2.0))
		Beta = B/((4.0*A**2.0 + B**2.0)**(1.0/2.0))
		if (Mach > 0.1): 
			Fc = (r*m)/((math.asin(Alpha) + math.asin(Beta))**2)
		elif (Mach <= 0.1):
			Fc = ((1.0 + math.sqrt(F))/2.0)**2.0
		#end
		Ftheta = math.sqrt(1.0/F)*((1.0 + 122.0/Tw*10**(-5.0/Tw))/(1.0 + 122.0/Te*10**(-5.0/Te)))
		Fx = Ftheta/Fc
		RexBar = Fx*Rex
		Cfb[0] = 0.074/RexBar**0.2
		for i in xrange(200):
			
			Cfb[i+1] = Cfb[i]*(1.0 + (0.242 - math.sqrt(Cfb[i])*math.log10(RexBar*Cfb[i]))/(0.121 + math.sqrt(Cfb[i])/math.log(10.0)))
			
			# Convergence Criteria
			if (abs(Cfb[i+1]-Cfb[i]) < 1e-9):
				Cf = Cfb[i+1]/Fc
				break
			elif (i == 200):
				print 'Cfb did not converge!'
				Cf = 0.455/(math.log10(Rex))**2.58   
			#end
			
		#end
		
		return Cf
	


#def plotCf():
#	
#
#% Variables Effects on Turbulent Skin Friction
#% i = 1;
#% for Mach = 0.01:0.01:4
#%     [Cf] = ACAero_Cf(0,Mach,10,0.05,3.33e-5);
#%     Cfi(i,1) = Cf;
#%     Cfi(i,2) = Mach;
#%     i=i+1;
#% end
#% plot(Cfi(:,2),Cfi(:,1));
#% title('Effect of Mach Number on Turbulent Skin Friction');
#% xlabel('Mach Number');
#% ylabel('Cf');
#% set(gcf,'color',[1 1 1]);



	
