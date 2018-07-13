.. _adflow_solvers:

Solvers
=======

This section contains some practical information about the available solver algorithms in ADflow, along with best practices.
Here, we assume that users started experimenting with ADflow, following the :ref:`tutorial <adflow_tutorial>` section.
The additional options we mention here can be added to the baseline runscripts users have to improve performance.

ADflow is capable of switching between solver algorithms during a solution procedure.
This process is controlled with the *relative convergence* metric, which is the ratio of current residual norm and the initial residual norm.
A 2 orders of magnitude relative convergence in this sense refers to converging the simulations such that the L2 norm of the current residual vector is one hundredth of the initial norm.
Technically, all three solver algorithms can be used in a single simulation, however, users will typically use either the multigrid or approximate Newton--Krylov (ANK) for the startup, combined with the Newton--Krylov (NK) for the final stages of convergence.

If all three solvers are used, ADflow uses the multigrid method until the relative convergence reaches ``'ankswitchtol'``.
Then the ANK method is used to reach ``'nkswitchtol'``, and after this relative convergence is reached, the solver switches to the NK solver.
Besides setting the *switch tolerances* for the ANK and NK solvers, users need to set the options ``'useanksolver'``, and ``'usenksolver'`` respectively to ``True`` in their runscripts.

With RANS simulations, the startup stage can be defined as the initial 4-6 orders of magnitude of convergence, while the terminal stage is the rest of the convergence until ``'L2Convergence'`` is reached.
For a RANS simulation, reaching steady state is equivalent to reducing the residual norm by 6-12 orders of magnitude.
The final relative convergence users prescribe is case dependent, however tighter tolerances (i.e. ``'L2Convergence': 1e-12``) will yield more accurate solutions and gradients, at the cost of higher computational effort.

Multigrid
---------

Multigrid is the baseline solver algorithm available in ADflow.
This method uses either the Runge--Kutta (RK), or the Diagonalized Diagonally-Dominant Alternating Direction Implicit (D3ADI) algorithm as the smoother.
The smoother algorithms only update the flow variables (i.e. density, momentum components, energy), and a Diagonalized Alternating Direction Implicit (DADI) method is used to update the turbulence model with RANS simulations.

Within a full multigrid solution process, the solution is initialized on the coarsest grid level, and the smoothers are used to reach the prescribed ``'L2ConvergenceCoarse'`` relative convergence is reached.
Once this target is reached, the solver moves to the next finer level, and repeats this process until the finest grid level is reached.
When the solver reaches to the finest grid level, it performs the prescribed multigrid cycle, that is set with the option ``'MGCycle'``.

Performance of the multigrid algorithm depends on a few critical options.
First of all, as the number of coarse grid levels are increased, the solver performance will generally improve.
However for practical cases, this is typically around 3-5, and coarser levels than these might not be available.
The number of grid levels is represented by the leading number with the ``'MGCycle'`` option (e.g. ``'3w'`` means performing a *w* multigrid cycle with 3 grid levels).

The second critical parameter is the CFL numbers selected.
CFL number is representative of the time-step size used with the smoother methods.
Larger CFL numbers translate into larger time steps, which in turn yield faster convergence to steady state.
However, both of the smoother algorithms are *conditionally stable*, meaning that there is a limit CFL number that can be used before the solution goes unstable.
A more stable solver can be achieved by reducing the CFL number, at the cost of higher computational effort, as the solver will take more iterations to reach the steady state.

The user can prescribe two different CFL numbers within ADflow.
The option ``'cfl'`` sets the CFL number used for operations on the finest grid level, while ``'cflCoarse'`` is used to set the CFL number used with the coarser grid levels.
The RK smoother is a fully explicit method, and the default values of ``'cfl': 1.7`` and ``'cflCoarse': 1.0`` should be reasonable.
When the D3ADI method is used by setting ``'smoother': 'DADI'``, the solver will benefit from a higher ``'cfl'`` value.
Even though the D3ADI smoother is an implicit method, the performance does not generally improve with higher CFL numbers, as the underlying factorization breaks down as the CFL number is increased.
Therefore, for practical cases, the CFL number used with the D3ADI implementation in ADflow is around 3-5.

Multigrid startup works well with multiblock meshes, as coarser levels of the grid is usually available.
However with overset meshes, coarser levels might not be available, and we recommend using the ANK solver for startup instead.
Furthermore, if the flowfield contains separated regions, the multigrid solver will likely stall, and not be able to complete the startup stage.
With these cases, we again recommend using the ANK solver.

Approximate Newton--Krylov
--------------------------

The approximate Newton--Krylov (ANK) solver is one of the fully implicit methods within ADflow.
It was developed based on the need for an efficient startup method for overset meshes, and it can also be used to converge difficult cases with heavy separation.
To use the ANK solver, users need to set ``'useanksolver': True`` in their scripts.
The default switch tolerance for the ANK solver (i.e. ``'ankswitchtol'``) is set to ``1.0``, therefore ADflow will start the simulations with the ANK solver if the solver is enabled.
However, if the coarser levels of the grid is available, and the ``'MGCycle'`` option is not set to ``'sg'`` (which stands for *single grid*), the solver might start with a full multigrid startup, and the ANK solver will only be activated once the finest grid level is reached.
This behavior can also be controlled with the ``'MGStartLevel'`` option, which controls the grid level that is used to initialize the simulations.
The default value for this option is set to ``-1``, which means that the solver will start with the coarsest grid level available.
See :ref:`ADflow options <adflow_options>` for a more detailed explanation of these parameters.
The ANK method can also be used after the multigrid method.
For example, if the user wants to obtain 2 orders of magnitude convergence with the multigrid, then switch to the ANK solver, the option ``'ankswitchtol'`` can be set to ``1e-2``, and the solver will switch to the ANK solver after this relative convergence value is reached.

The ANK solver has a large number of tunable parameters, and a few different modes.
The default set of tunings works well for a typical aeronautical application with transonic conditions, but there is always room for improvement.
For a quick and dirty solution, simply include ``'useanksolver': True`` option in the runscripts.
Furthermore, we recommend setting a relatively high turbulence sub-iteration number for RANS simulations, i.e. set ``'nsubiterturb'`` to somewhere between ``3`` and ``7``.
However, there will be cases where the default set of parameters will not yield good performance, and some tuning might be required.
To be able to understand why the solver is not working well, or even failing, the users need to understand the underlying algorithm to some level.
In the following subsections, we describe the main components within the ANK solver algorithm, and lay down the best practices for troubleshooting.


Solver Algorithm
~~~~~~~~~~~~~~~~

Matrix-Free Operations
~~~~~~~~~~~~~~~~~~~~~~

Turbulence Coupling
~~~~~~~~~~~~~~~~~~~

Interpreting the Output
~~~~~~~~~~~~~~~~~~~~~~~

Troubleshooting
~~~~~~~~~~~~~~~



Newton--Krylov
--------------

Solver Algorithm
~~~~~~~~~~~~~~~~

The Eisenstat--Walker Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Linear Solver Performance
~~~~~~~~~~~~~~~~~~~~~~~~~

Troubleshooting
~~~~~~~~~~~~~~~
