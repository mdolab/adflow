.. SUmb documentation master file, created by
   sphinx-quickstart on Tue Mar 25 15:53:41 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====
SUmb
====

SUmb -- Stanford University MultiBlcok -- is structured multi-block 3D CFD solver. 

=====================================    =========   ========================      =======================================================================================
Parameter                                  Type       Default                                   Description
=====================================    =========   ========================      =======================================================================================
 `gridFile`                              `str` *     `None`                        This is the grid file to use. It must be a multi-block-structured CGNS file with
                                                                                   all block-to-block information and boundary condition information contained in the
		                       				                
`restartFile`                            `str` *     `restart.cgns`                Use this CGNS file to restart the solution from a previous solution. This file must
                                                                                   have been a volume solution file that was written by ``sumb``.  
                                                                                                                                                                         
`solRestart`                             `bool` *    `False`                       Flag to use specified restart file.
                                                                                                                                                                         
`storeRindLayer`                         `bool`      `True`                        Flag to use write halo or "rind cell" information into cgns files. This is required
                                                                                   to have some postprocessors compute contour lines correctly (Tecplot).                
		                       				                
`writeSymmetry`                          `bool`      `True`                        Flag to write symmetry plane in the output surface solution file.  
		                       				                
`writeFarfield`                          `bool`      `True`                        Flag to write farfield surfaces in output surface solution file. 
								                
`writeSurfaceSolution`                   `bool`      `True`                        Flag to write surface solution (automatically) after the end of each solution call. 
								                
`writeVolumeSolution`                    `bool`      `True`                        Flag to write volume solution(automatically) after the end of each solution call. 
								                
`solutionPrecision`                      `str`       `single`                      Must be `single` or `double`. Write solution data in surface and volume files in 
                                                                                   desired precision. If you are not planning on using restart files, it is best to use
                                                                                   single precision since the output files will be half the size and for visualization
                                                                                   purposes, this is more than sufficient. 
								                
`gridPrecision`                          `str`       `double`                      Must be `single` or `double`. Double precision is preferred especially for RANS grids
                                                                                   since lower precision can adversely affect the cells in the boundary layer. 
								                
`isoSurface`                             `dict`      `{}`                          Dictionary specifying the type and values to be used for isosurfaces. 
                                                                                   Any of the "volumeVariables" may be used. An example of the format is as 
								                   follows: 'isoSurface':{'Vx':-0.001, 'shock':1.0}
								                   This will place a isosurface at (essentially) 0 x-velocity and a iso-surface 
								                   at the shock sensor value of 1. (used to visualize shock region). 
								                
`viscousSurfaceVelocities`               `bool`      `True`                        Only applicable for RANS/laminar computations. Write surface velocities as the value
                                                                                   on the first cell off the wall. This means that the velocities on the surface will
								                   **not** be zero as enforced by the boundary conditions. The reason for this option
								                   is it lets post-processing software compute oil-flow surface stream-line patterns. 
								                
`discretization`                         `str`       `central plus`                Set the discretrization method. `central plus scalar dissipation` uses a central 
                                                     `scalar dissipation`          finite volume scheme with JST scalar dissipation. Other options are `central plus 
								                   matrix dissipation` (which generally does not work properly) and `upwind`.  The default
									       	   is generally recommended for both robustness and speed, at the cost of numerical 
									      	   accurary. 
									      
`coarseDiscretization`                   `str`       `central plus`                Set the discretrization method for the coarse grid. Generally should be the same as
                                                     `scalar dissipation`          the `discretization` option. 
									      
`limiter`                                `str`       `valalbeta`                   Type of flux limiter to use for `Upwind` scheme. Possible values are `vanalbeda` 
                                                                                   `minmod` and `nolimiter`. 
									      
`smoother`                               `str`       `runge kutta`                 Type of `solver` or smoother to use with Runge Kutta method. The default and most 
                                                                                   tested is `runge kutta`. Other possible values are `dadi` which uses a diagonalized
									      	   diagonal-dominant alternating-direction implicit (D3ADI) scheme. This is typically 
									      	   faster that the Runge Kutta method but may be less robust. 
									      
`equationType`                           `str` *     `euler`                       The type of equations to solve. Possible values are `euler`, `laminar NS`, or `RANS`. 
									      
`equationMode`                           `str` *     `steady`                      The temporal mode of the equations to solve. Possible values are `steady`, `unsteady`
                                                                                   or `time spectral`. The `steady` and `time spectral` methods have been extensively 
									      	   tested from the Python interface. The unsteady method has not been extensively 
									      	   tested from Python.
									      
`flowType`                               `str` *     `external`                    Type of flow simulation. Must be `internal` or `external`. Only external flow 
                                                                                   simulations have been tested with the Python interface. 
									      
`turbulenceModel`                        `str` *     `sa`                          Select the turbulence model to use. Possible values are `sa` for the Spalart Allmaras 
                                                                                   model, `sae` for the Sparart Allmaras-Edward model, `k omega wilcox`, 
									      	   `k omega modified`, `ktau`, `menter sst` and `v2f`. For external aerodynamic flow 
									      	   applications, `sa` is recommended as this is currently the only turbulence model
									      	   that has been differentiated. 
									      
`turbulenceOrder`                        `str` *     `first order`                 The numerical order of accuracy of the turbulence model. Possible values are 
                                                                                   `first order` or `second order`. Generally `first order` is recommended as the
									      	   adjoints systems are much easier to solve with the first order discretrization. 
									      
`useWallFunctions`                       `bool` *    `False`                       Flag specifying if wall functions are to be used. This is generally not recommended
                                                                                   since they give (potentially very) poor drag estimates. Furthermore, the required 
									      	   routines are differentiated so wall function simulations cannot be used for 
									      	   optimization
									      
`useApproximateWallDistance`             `bool` *    `True`                        Flag to use a `cheap` wall distance calculation. When this is true, the exact wall
                                                                                   distances are computed during initialization and the parametric location of the 
									      	   closest wall point is stored for each cell. After the geometry deforms (say during
									      	   an optimization) the spatial search algorithm is not run, but the distance between
									      	   the (new) parametric location and the (new) grid cell center is computed and taken
									      	   as the wall distance. This is substantially faster and permits efficient wall-distance
									      	   updates for use in aerostructural analysis. 
										
`wallTreatment`                          `str`       `linear pressure`             Specifies how the boundary conditions are implemented. Generally the default value
                                                     `extrapolation`               need not be changed. Other values include `constant pressure extrapolation`,
						                                   `quadratic pressure extrapolation` and `normal momentum`. Only `linear pressure 
									      	   extrapolation` and `constant pressure` extrapolation are know to work with the 
									      	   adjoint method. 
									      
`dissipationScalingExponent`             `float`     0.67                          Exponent factor to use in JST dissipation scheme. This value typically will not need
                                                                                   to be changed from its default value. The value of 2/3 is the theoretical best value
									      	   for this value assuming an orthogonal 3 dimensional grid. 
									      
`vis4`                                   `float`     0.0156                        Coefficient of the fourth order dissipation used in the scalar and matrix JST 
                                                                                   dissipation scheme. The default value is generally recommended if a converged solution
									      	   can be obtained. It may be raised slightly in the range of 0.02-0.025 which may help
									      	   achieve better convergence properties at the expense of numerical accuracy. 
									      
`vis2`                                   `float`     0.25                          Coefficient of the second order dissipation used in the scalar and matrix JST
                                                                                   dissipation schemes. This dissipation is only turned on at shocks, and thus may be 
									      	   set to 0.0 if the user knows a simulation will be entirely subsonic. 
									      
`vis2Coarse`                             `float`     0.50                          Set a difference vis2 for the coarse grid. This is typically larger than vis2. The
                                                                                   default value of 0.5 is generally sufficient for most cases. 
									      
`restrictionRelaxation`                  `float`     0.80                          The relaxation factor for the restriction operation in multigrid. Value must be 
                                                                                   between 0 and 1.0. A value of 1.0 will not perform any relaxation. On some problem
									      	   this may be faster, while slower on others. The default value of 0.80 appears to work
									      	   well for a wide variety of cases. 
									      
									      
`liftIndex`                              `int` *     None                          Specify the coordinate index that will be considered the 'lift' direction. 
                                                                                   If not supplied, this parameter will be determined automatically if there are 
									      	   symmetry planes present in the grid. Otherwise, it must be supplied. The applicable 
									      	   values are 2 for the y-axis as the lift direction and 3 for the z-axis as the lift
									      	   direction. 
									      
`nCycles`                                `int`       500                           Maximum Number of "iterations" to run. For the Runge Kutta solver this refers to the 
                                                                                   number of multigrid cycles to run on the fine grid. When the NK solver is used, it refers to the
									      	   total number of multi-grid cycles **plus** the number of function evaluations. Each
									      	   function evaluation corresponds roughly to single residual evaluation. 
									      
`nCyclesCoarse`                          `int`       500                           Maximum number of iterations to run on the coarse grid when performing a full-multigrid
                                                                                   start-up procedure. 
									      
`nSubIterTurb`                           `int`       1                             The number of **additional** iterations of the turbulent ADI solver to run. Only 
                                                                                   meaningful for RANS simulations. Certain RANS simulations may benefit from a slight
									      	   increase of this parameter to 2 or 3 which will lower the overall solution time. 
`CFL`                                    `float`     1.5                           The Courant–Friedrichs–Lewy (CFL) number to use for the Runge-Kutta simulations. This
                                                                                   is the main parameter that determines the overall speed and robustness of RK simulations.
									      	   Lower CFL numbers give more robust solutions but are slower. The default parameter of 
									      	   1.5 is a good place to start. Usually some experimentation is required to determine
									      	   the maximum CFL for a particular simulation. 
									      
`CFLCoarse`                              `float`     1.o                           The CFL number to use on the coarse grids of the multigrid simulations. It is often
                                                                                   desirable to have this number somewhat lower than the CFL number of the fine grid. 
									      
`mcCycle`                                `str` *     `3w`                          The type of multigrid cycle to use. The dimensions of the grid must be such that the 
                                                                                   requested multigrid level is possible. To run a single grid simulation (no multigrid)
									      	   use `sg`. To run 3 multigrid levels with a 'w' cycle use `3w`. To use a 'v' cycle use
									      	   `3v` etc. 
									      
`mgStartLevel`                           `int` *     -1                            Specifiy the starting grid level. This is ued to perform a "full multigrid startup"
                                                                                   procedure. This can lead to significally reduced simulation times since a good starting
									      	   point can be obtained from approximate solutions on the coarser grids. A -1 indicated
									      	   that the coarsest grid level should be used. For RANS simulations, it is often not
									      	   possible to start on the coarest grid, espeically if the coarse grid has very few 
									      	   cells. 
									      
`resAveraging`                           `str`       `alternateResAveraging`       Only perform residual averaging on every second stage of the RK procedure. This
                                                                                   save computation, but has very little impact on the convergence properties. 
									      
`smoothParameter`                        `float`     1.5                           Parameter used in residual smoothing. This value will typically not need to be
                                                                                   changed from the default. 
									      
`cflLimit`                               `float`     1.5                           The maximum CFL that could be run withiout residual smoothing. If the actual CFL
                                                                                   is lower than the CFLLimit, not smoothing will be applied, regardless of the `resAveraging`
									      	   option 
									      
`timeItegrationScheme`                   `str`       `bdf`                         The type of time integration scheme to use for unsteady analysis. Only the `bdf` option
                                                                                   is currently known to work. 
									      
`timeAccurary`                           `int`       2                             Order of accuracy of the time integration scheme. Valid values are 1, 2, or 3. 
									      
`nTimeStepsFine`                         `int`       100                           Number of time steps to run in an unsteady simulation. Note that MGStart level
                                                                                   should be 1 for a unsteady simulation.
									      
`deltaT`                                 `float`     0.01                          Time step to use for unsteady simulation.
									      
									      
`timeIntervals`                          `int`       1                             The number of "spectral instances" to use for a time spectral simulation. This 
                                                                                   option is only meaningful when `equationMode` is `time spectral`. 
									      
`alphaMode`                              `bool`      False                         Use a specified alpha motion for the Time spectral analysis. 
									      
`betaaMode`                              `bool`      False                         Use a specified beta motion for the Time spectral analysis.  Untested.
									      
`machMode`                               `bool`      False                         Use a specified mach number motion for the Time spectral analysis. Untested
									      
`pmode`   `                              `bool`      False                         Use a specified p-motion (rolling) motion for the Time spectral analysis. Untested.
									      
`qmode`                                  `bool`      False                         Use a specified q-motion (pictch) motion for the Time spectral analysis. 
									      
`rmode`                                  `bool`      False                         Use a specified r-motino (yaw) motion for the Time spectral analysis. Untested
									      
`altitudeMode`                           `bool`      False                         Use a specified h-variation  motion for the Time spectral analysis. Untested
									      
`windAxis`                               `bool`      False                         Not sure?
									      
`TSStability`                            `bool`      Flag                          Flag to compute time spectral stability information from a timespectral CFD solution
									      
`l2Convergence`                          `float`     1e-6                          This specifies the desired convergence factor. For the RK solver, this is taken 
                                                                                   relative initial residual on the **fine** grid. Since this prolonged solution
									      	   may be a fairly good starting point, the **actual** convergence relative to a 
									      	   free stream residual may be 1 to 2 orders magnitudes lower. For the NK solver, this
									      	   option also determines the convergence, but the reference is taken as free-stream 
									      	   residual. 
									      
`l2ConvergenceRel`                       `float`     1e-16                         This option is typically **only** used when SUmb is used in conjunction with an
                                                                                   aerostructural solver. This specifies the relative tolerance in relation to the
									      	   current starting point. 
									      
`l2ConvergneceCoarse`                    `float`     1e-2                          The convergence factor to perform on the coarse grids during multi-grid startup. 
                                                                                   Most of the benefits of the start-up procedure is obtained after converging 
									      	   between 2 and 3 orders of magnitude so this options is typically 1e-2 to 1e-3. 
									      
`maxL2DeviationFactor`                   `float`     1.0                           If the solver runs out of iterations, the maximum factor the residual can be 
                                                                                   above the target residual (as determined by l2Convergence) and still be considered
									      	   "converged". 
									      
`minIterationNum`                        `int`       10                            This option ensures that a minmum number of iterations are performed when using the 
                                                                                   RK solver. This can be useful when only changing the angle of attack; A small 
									      	   change in the anlge attack is not sufficient to increase the residual and the 
									      	   solver may stop prematurely before the peturbation is actually solved. 

`useNKSolver`                            `bool`      False                         Flag to turn on the Newton--Krylov solver. If this flag is `False`, the remainder of the
                                                                                   of the options that begin with `nk` will have no effect. The Newton solver only works 
										   with the Euler and Laminar NS equations, in either steady or time-spectral modes. 

`NKLinearSolver`                         `str`       `gmres`                       Type of PETSc KSP solver to use for the solution of the linear systems that arise
										   from Newton's method. For practically all cases, GMRES will perform the best. `TFQMR` --
										   Transpose-Free quasi minimial residual may also be used in certain situations which 
										   will use less memory that GMRES.

`NKSwtichTol`                            `float`     1e-2                          The relative tolerance to converge before the switch is made to the Newton solution
										   technique. This must be low enough that most of the difficult transients have beeen 
										   passed. If the NK solver stalls, this value can be set to a lower value which will
										   run the RK solver longer before switching. 

`NKSubSpaceSize`                         `int`       60                            The size of the GMRES subspace for the NK solver. For difficult problems, convergence
                                                                                   may be improved by increasing this value at the expense of more memory. 

`NKLinearSolveTol`                       `float`     1e-1                          The inital tolerance to solve the linear system resulting from the Newton approximation. 
                                                                                   This value is only used for the first solution; thereafter the forcing tolerance
										   is updated dynamically using the Einstat-Walker forcing criteria. 

`NKPC`                                   `str`       `additive schwartz`           The type of (global) preconditioner to use for the linearized system. The default
                                                                                   is recommended unless memory is a issue. In that case, `block jacobi` can be used
										   which is less efficient but, has a lower memory footprint.

`NKASMOverlap`                           `int`       1                             The number of overlap levels in the ASM preconditions. More overlap levels result in a
                                                                                   stronger preconditioner, at the expense of more expensive iterations and more memory. 
										   Typicaly values range from 1 for easy problems up to 2 or 3 for more difficult ones. 

`NKPCILUFill`                            `int`       1                             The number of levels of fill to use on the local (subdomain) Incomplete LU (ILU) factorization
                                                                                   Typical values are 1 for easy cases and up to 3 for more difficult cases. More levels
										   of fill result in a stronger precondtioner which will result in fewer (linear) 
										   iterations, but individual iterations will be more costly and consume more memory. 

`NKLocalOrdering`                        `str`       `rcm`                         The type of reordering algorithm to use on the local subdomains. For practically all
                                                                                   cases Reverse Cuthill McKee performs the best. 

`NKJacobianLag`                          `int`       10                            The option determines the frequency at which the precondition is reformed. In other words
                                                                                   the jacobian used for form the precondition is "lagged" behind the actual solution by
										   10 iterations. For simple problems, it may be possible to increase the jacobian lag
										   to such a high value that the precondition is never reformed at all during a solution. 
										   For more difficult cases, a lower value may help convergence. A lower value will 
										   result in more (preconditioner) jacobian assemblies that are fairly costly in SUmb. 

`RKReset`                                 `bool`      `False`                      ???

`NKReset`                                `int`       5                             ???

`applyPCSubSpaceSize`                    `int`       10                            This option is only used when SUmb is used in an aero-structural analysis. This parameter
                                                                                   determines the subspace **and** the total number of iterations to run when SUmb is only
										   being used to precondition residuals via the globalNKPreCon() function. 

`NKOuterPreConIts`                       `int`       1                             Number of times to apply the global (NKPC option) precondition. More iterations may help
                                                                                   converge the linear system faster. Typical values are from 1 to 3. 

`NKInnerPreConIts`                       `int`       1                             Number of time to apply the local precondition. More iterations may help converge the 
                                                                                   linear system faster. This should be left at 1, unless a very difficult problem is 
										   encountered. 

`blockSplitting`                         `bool` *    True                          Flag determining if the block may be split to obtain better load 
=====================================    =========   ========================      =======================================================================================


	    ..
	       'blocksplitting':[bool, True],*
	       'loadimbalance':[float, 0.1], *
	       'loadbalanceiter':[int, 10], *
	       'partitiononly':[bool, False], *
	       'metricconversion':[float, 1.0],
	       'autosolveretry':[bool, False],
	       'autoadjointretry':[bool, False], *
	       'storehistory':[bool, False], ????????
	       'numbersolutions':[bool, True],
	       'printiterations':[bool, True],
	       'printtiming':[bool, True],
	       'setmonitor':[bool, True],        
	       'monitorvariables':[list, ['cpu','resrho','cl', 'cd']], (CHECK)
	       'surfacevariables':[list, ['cp','vx', 'vy','vz', 'mach']], (chECK)
	       'volumevariables':[list, ['resrho']], (check)
	       'forcesastractions':[bool, True], 
	       'adjointl2convergence':[float, 1e-6],
	       'adjointl2convergencerel':[float, 1e-16],
	       'adjointl2convergenceabs':[float, 1e-16],
	       'adjointdivtol':[float, 1e5],
	       'approxpc': [bool, True],
	       'adpc': [bool, False],
	       'viscpc':[bool,False],
	       'usediagtspc':[bool, True],
	       'restartadjoint':[bool, True],
	       'adjointsolver': [str, 'gmres'],
	       'adjointmaxiter': [int, 500],
	       'adjointsubspacesize' : [int, 100],
	       'adjointmonitorstep': [int, 10],
	       'dissipationlumpingparameter':[float, 6.0],
	       'preconditionerside': [str, 'right'],
	       'matrixordering': [str, 'rcm'],
	       'globalpreconditioner': [str, 'additive schwartz'],
	       'localpreconditioner' : [str, 'ilu'],
	       'ilufill': [int, 1],
	       'asmoverlap' : [int, 1],
	       'innerpreconits':[int,1],
	       'outerpreconits':[int,3],
	       'usereversemodead':[bool, False],
	       'applyadjointpcsubspacesize':[int, 20],
	       'frozenturbulence':[bool, True], **** check *****
	       'firstrun':[bool, True],
	       'verifystate':[bool, True],
	       'verifyspatial':[bool, True],
	       'verifyextra':[bool, True],


Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

