.. _tutorial:

Tutorial
========

The following shows how to get started with SUmb by running the mdo_tutorial
wing problem. First, we show the complete program listing and then go through 
each statement line by line::

  # This is a template that should be used for setting up
  # RANS analysis scripts
  
  # ======================================================================
  #         Import modules
  # ======================================================================
  import numpy
  from mpi4py import MPI
  from baseclasses import *
  from sumb import SUMB
  
  # ======================================================================
  #         Input Information -- Modify accordingly!
  # ======================================================================
  outputDirectory = './'
  gridFile = '../INPUT/rans_grid.cgns'
  alpha = 2.0
  mach = 0.78
  areaRef = 45.5
  chordRef = 3.25
  MGCycle = '3w'
  altitude = 10000
  name = 'fc'
  
  aeroOptions = {
  # Common Parameters
  'gridFile':gridFile,
  'outputDirectory':outputDirectory,
  
  # Physics Parameters
  'equationType':'RANS',
    
  # Common Parameters
  'CFL':1.5,
  'CFLCoarse':1.25,
  'MGCycle':MGCycle, 
  'MGStartLevel':-1,
  'nCyclesCoarse':250,
  'nCycles':1000,
  'monitorvariables':['resrho','cl','cd'],
  'useNKSolver':False,
    
  # Convergence Parameters
  'L2Convergence':1e-6,
  'L2ConvergenceCoarse':1e-2,
  }

  # Aerodynamic problem description
  ap = AeroProblem(name=name, alpha=alpha, mach=mach, altitude=altitude,
  areaRef=areaRef, chordRef=chordRef,
  evalFuncs=['cl','cd'])
  # Create solver
  CFDSolver = SUMB(options=aeroOptions)
  
  # Solve and evaluate functions
  funcs = {}
  CFDSolver(ap)
  CFDSolver.evalFunctions(ap, funcs)
  
  # Print the evaluatd functions
  if MPI.COMM_WORLD.rank == 0:
  print funcs
  
Start by importing the ``sumb``, ``baseclasses``, and ``mpi4py``.::

  # ======================================================================
  #         Import modules
  # ======================================================================
  import numpy
  from mpi4py import MPI
  from baseclasses import *
  from sumb import SUMB

Then, we define inputs as well as options for SUmb. The grid files is in
CGNS format. We usually use ICEM CFD to generate grids. Here, we also define
a few flow conditions and reference values to be put into ``AeroProblem`` later.::

  # ======================================================================
  #         Input Information -- Modify accordingly!
  # ======================================================================
  outputDirectory = './'
  gridFile = '../INPUT/rans_grid.cgns'
  alpha = 2.0
  mach = 0.78
  areaRef = 45.5
  chordRef = 3.25
  MGCycle = '3w'
  altitude = 10000
  name = 'fc'
  
  aeroOptions = {
  # Common Parameters
  'gridFile':gridFile,
  'outputDirectory':outputDirectory,
  
  # Physics Parameters
  'equationType':'RANS',
    
  # Common Parameters
  'CFL':1.5,
  'CFLCoarse':1.25,
  'MGCycle':MGCycle, 
  'MGStartLevel':-1,
  'nCyclesCoarse':250,
  'nCycles':1000,
  'monitorvariables':['resrho','cl','cd'],
  'useNKSolver':False,
    
  # Convergence Parameters
  'L2Convergence':1e-6,
  'L2ConvergenceCoarse':1e-2,
  }

Now, this is the actually solution part. We start by defining the ``AeroProblem``, 
which is import from ``baseclasses``. We specify flow condtions and reference values
into the ``aeroProblem``. We also tell the solver which solution values that
we are interested in. In this case, we use the keyword ``evalFuncs``. ::

  # Aerodynamic problem description
  ap = AeroProblem(name=name, alpha=alpha, mach=mach, altitude=altitude,
  areaRef=areaRef, chordRef=chordRef,
  evalFuncs=['cl','cd'])

Then, we create the SUmb instant. We also provide SUmb all the options that we 
just specified above. ::

  # Create solver
  CFDSolver = SUMB(options=aeroOptions)
  
Now, we solve the CFD problem. ``CFDSolver(ap)`` is the command that actually
solve the CFD. You can see print out from SUmb of each iteration here. This 
example will take just a couple minutes. ``CFDSolver.evalFunctions()`` return
the function of interests we specified in ``AeroProblem``.::

  # Solve and evaluate functions
  funcs = {}
  CFDSolver(ap)
  CFDSolver.evalFunctions(ap, funcs)
  
Finally, we print out the value of `cd` and `cl`. We only print on the 
root processor. ::
 
  # Print the evaluatd functions
  if MPI.COMM_WORLD.rank == 0:
  print funcs

