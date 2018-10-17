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

The ANK solver uses the backward Euler method to iteratively solve the system of equations.
As the method is fully implicit, the solver is unconditionally stable, meaning that there is not a stability limit for the CFL number.
In the context of backward Euler scheme, the CFL number used in the solver has very important outcomes.
Relatively small CFL numbers (around 1-10) translate into small time steps, and the solver yields very favorable stability properties.
However, convergence with low CFL numbers stagnate quickly after the first few orders of magnitude relative convergence.
At high CFL numbers (around 100-1e5), the solver approaches the Newton's method, and yields very favorable nonlinear convergence properties.
However, if the current state is far away from the solution, high CFL numbers might cause the solver to stagnate.
Keeping these outcomes in mind, we start the solver algorithm with a relatively small CFL number (``'ankcfl0': 5.0``), and adaptively ramp the CFL number to higher values as the simulation converges to steady state.
This process is called pseudo-transient continuation.
This continuation process is what enables the ANK solver to be used as a startup strategy, while yielding favorable convergence rates later on in the simulation.

During each nonlinear iteration, we determine the updates to the state variables by *inexactly* solving a large linear system.
In this context, inexactly solving the linear system refers to reducing the linear residual norm by a factor of ``'anklinearsolvetol'``, and the default value we use is ``0.05``.
We avoid exactly solving the linear system because it is often more beneficial to take a large number of cheaper nonlinear iterations, because the startup of RANS simulations often require tens or sometimes hundreds of nonlinear iterations.
Prescribing a higher linear solution tolerance on the other hand might destabilize the solver, because in this case the solution vector is too far from the actual solution.
The linear systems we solve with the ANK solver contains two components: the Jacobian and a time-stepping term.
This is simply a result of using the backward Euler method.
The time stepping term is calculated by selecting a global CFL number that is applied to each cell in the domain, while the Jacobian matrix contains the partial derivatives of residuals with respect to states in each cell.

The linear solver algorithm we use is called generalized minimal residual method (GMRES).
We use a matrix-based preconditioner to improve the convergence rate of the GMRES algorithm, while the actualy linear system we solve is never formed; instead we use a matrix-free approach.
This is enabled by the fact that the GMRES algorithm only requires matrix-vector products, rather than the full matrix itself.
The preconditioner we use is based on a *first order* Jacobian, that contains some approximations that are useful in reducing the memory requirements while storing the matrix.
Compared to the 33-cell stencil of the second order accurate scheme, this first order accurate Jacobian only requires a 7-cell stencil.
The matrix-free operations themselves also contain some approximations.
The resulting Jacobian we obtain with the matrix-free operations is somewhere between the full second, and first order accurate schemes.
We call the solver *approximate* Newton--Krylov because of this reason; the matrix-free operations contain some approximations.

The use of matrix-free operations enable the approximate Jacobian used in the solution process to be up to date on every nonlinear iteration.
However, keeping an up to date preconditioner for every nonlinear iteration is not a practical approach.
This is because forming and factorizing the preconditioner matrix is very expensive.
To alleviate this cost, we lag the preconditioner between nonlinear iterations.
The lagging process is done automatically, and user input is not required.
However, users should be aware that the CFL number for the ANK solver is only updated during iterations where we update the preconditioners.
These iterations are marked with a ``*`` leading the iteration type output.

After inexactly solving the linear system, we obtain an update vector for the state variables, however we usually do not take the full update vector.
Instead, we relax the update by a certain factor between 0 and 1.
This is similar to taking full or partial steps within an optimization.
A full step is equivalent to taking the update vector as it is, and a limited step is equivalent to multiplying each entry of the update vector with a relaxation parameter.

To determine this relaxation parameter, we first check the total physical change that the update yields.
We limit the step size such that the density and energy of each cell do not change more than ``'ankphysicallstol'`` fraction of the original value.
The default value is ``0.2``, which translates to limiting the physical change in these variables to 20% of the original value.
For the turbulence model, we follow a similar approach, however we only check the updates that reduce the value of the turbulence model variable, and limit the change to 99%, by using ``'ankphysicallstolturb': 0.99``.
We refer to to this process as the physicality check.

After the physicality check, we go into a backtracking line search, where the goal is to find a step size that yields a reduction in the unsteady residual norm.
This unsteady residual norm is different from the steady residual norm printed in the output.
As a result, the steady (or total) residual norm can actually increase, while the unsteady residual norm decreases.
This backtracking search starts with the step size calculated with the physicality check, and then traces this step back until it finds a step size that gives a reduction in the unsteady residual norm.
After the backtracking line search, the solver multiplies the update vector with the step size and updates the state vector.
We repeat this process until the simulation converges, or we reach ``'nkswitchtol'`` relative convergence.

Matrix-Free Operations
~~~~~~~~~~~~~~~~~~~~~~

The use of matrix-free operations for the actual linear system gives us the flexibility to be able to modify the Jacobian formulation on the go without any increased memory cost.
The default matrix-free operations contain some approximations compared to the exact residual routines.
However, users can switch to an exact Jacobian during runtime to improve nonlinear convergence.
This is achieved by using the ``'anksecondordswitchtol'`` option.
This option prescribes a relative convergence limit, above which the solver uses the default approximate Jacobian.
However, this relative convergence value is reached, the solver switches to using an exact Jacobian formulation for the matrix-free operations.
For example, setting ``'anksecondordswitchtol': 1e-2`` would cause the solver to use the approximate formulation for the initial 2 orders of magnitude convergence, and then switch to the exact formulation for the rest of convergence.
Note that this modification only changes how the implicit system is handled, and does not alter the baseline residual formulations.
Therefore the only effect will be in nonlinear convergence rate, and cost of each nonlinear iteration.

The approximate Jacobian is designed to have better conditioning properties, i.e. it is easier to solve numerically.
However, these approximations reduce the accuracy of the update vector, and nonlinear convergence rate suffers from this.
On the other hand, the exact Jacobian would be more difficult to solve compared to the approximate one, however the update vector obtained this way is expected to yield better nonlinear convergence.
As a result, the tradeoff is between cost of each nonlinear iteration, and rate of nonlinear convergence.

The second order switch is set to ``1e-16`` by default, meaning it is disabled.
However it can be manually set to improve performance.
For many practical cases, default approximate Jacobian is faster in the first 3-4 orders of magnitude convergence.
The solver can benefit from switching to second order formulation after this point.
However, the users should keep an eye on the linear residual during each nonlinear iteration.
As the second order Jacobian creates a more difficult linear system, the linear solver might fail and this might cause the solver to go unstable.
In cases where the prescribed linear solution tolerance cannot be met (e.g. linear residual above the ``0.05`` default value), the users are better off with just using the approximate formulation and disabling the second order switch.
The optimal switching point is case dependent, and users are encouraged to experiment with it.
Finally, the solver will print an ``S`` before the ``ANK`` identifier to state that it is using the second order Jacobian formulation.

.. _turbulence_coupling:

Turbulence Coupling
~~~~~~~~~~~~~~~~~~~

The turbulence models used with RANS equations can be notoriously difficult to converge.
To prevent issues in convergence, we solve the turbulence model separately from the flow variables.
In this context, the flow variables refer to the density, momentum components and energy, while the turbulence variable is typically the SA model working variable.
A decoupled algorithm updates the flow variables with the algorithms described above, and after updating the flow variables in each nonlinear iteration, we perform sub-iterations for the turbulence model before moving onto the next iteration.
Best way to diagnose if the turbulence model is causing problems with convergence is to print the turbulence residual norm with the output.
This can be done with including ``'resturb'`` in the list passed with the option ``'monitorvariables'``.

In the decoupled mode, the ANK solver has two turbulence solvers available.
The first one is the DADI based solver, which we refer to as ``turbDADI``.
This solver algorithm calculates the update vector for the turbulence model by using the diagonalized alternating direction implicit algorithm.
We typically use this algorithm with a large number of sub iterations.
The option ``'nsubiterturb'`` can be used to set the number of sub-iterations for the turbulence model to be performed after each flow update.
We typically recommend a value between 3 and 7, however more difficult cases might require up to 10 sub-iterations for the turbulence model.

The second turbulence solver available is called ``turbKSP``.
This is practically an isolated ANK algorithm just for the turbulence model.
We use the exact same options, and algorithms with the ANK solver, however the ``turbKSP`` solver has its own matrices, and linear system.
After each flow update, we repeat the similar ANK process for the turbulence model, and compute an update vector just for the turbulence model variable.
To use the ``turbKSP`` solver with ANK, users can set ``'ankuseturbdadi': False``, which implies that the solver will use the ``turbKSP`` solver instead of the default ``turbDADI`` solver.
To print useful debugging information about the turbulence solver, users can include ``'ankturbkspdebug': True`` in their runscripts.
The solver will then print some diagnostics about the turbulence solver in each iteration, such as linear convergence, step size, number of iterations GMRES algorithm takes etc.
This output will not look pretty but can be very useful while debugging.
We recommend using only 1 sub-iteration for this solver, as it is much more expensive, but more powerful at the same time compared to the ``turbDADI`` solver.
To use larger number of sub-iterations, users can set ``'anknsubiterturb'`` option to any integer larger than 1.

For smaller cases (<1M Cells) with multiblock meshes, we recommend using the turbDADI solver.
However, the performance of the ``turbDADI`` solver will deteriorate with overset meshes, therefore users can get better performance by switching to the ``turbKSP`` solver with more realistic cases (>1M Cells) with overset meshes.

Instead of the decoupled mode, the ANK solver is also capable of coupling the turbulence model to the flow variables.
We call this coupled ANK, and iterations in this mode is denoted with a ``C`` character before the ``ANK`` identifier.
The coupled mode can be beneficial because the solver now considers the coupling between the turbulence model and the flow variables, and this mode is expected to yield better convergence during the final stages of startup.
However, running in decoupled mode for the initial 4-5 orders of magnitude convergence is almost always going to yield better performance.
To start with the decoupled algorithm, and switch to the coupled algorithm, users can set a target relative convergence value with the ``'ankcoupledswitchtol'`` option.
Similar to the second order switch, the default for coupled switch tolerance is set to ``1e-16``.
To enable the coupled solver, users can pick a relative convergence value, e.g. setting ``'ankcoupledswitchtol': 1e-4`` will cause the solver to switch to the coupled formulation after 4 orders of magnitude of relative convergence is reached.
In this mode, the turbulence model and the flow variables are updated together with the ANK algorithm described above, and no sub-iterations for the turbulence model is performed.

There are two important aspects of converging the SA turbulence model.
First of all, with almost every case, the turbulence model residual norm will drop a few orders of magnitude within the initial 1-2 orders of magnitude relative convergence.
After this, the turbulence model residual will start increasing, until about 4 orders of magnitude relative convergence.
After this *hill*, the turbulence model residuals usually goes down monotonically.
Users should be aware that if coupled switch happens before the turbulence model goes over the hill, the solver might stall, or yield very bad performance.
It is usually better to use a decoupled method before this hill, and a coupled method after.

Second important aspect is related to the scaling of the turbulence model residuals.
The flow variables in ADflow are normalized with respect to the free stream reference values.
For example, density and velocity values of 1 represent the values that would be obtained in the free stream.
This is done to prevent precision loss with numerical algorithms, and with this normalization, the residuals of density, momentum and energy are typically around similar orders of magnitude.
However, performing the same normalization for the turbulence model is not straightforward, as the turbulence model variable can span a few orders of magnitude even after normalization.
As a result, the turbulence model variable is not normalized in the same way, and the turbulence residual norm is usually 4-5 orders of magnitude lower than the flow variables' residual norms.
To prevent numerical difficulties with the coupled ANK (and the NK) solver, we scale the turbulence model residual norm by ''1e4``.
This scaling has implications with the coupled solvers.
The coupled solvers will only yield good performance when the printed turbulence model residual norm is about 3-5 orders of magnitude lower than the flow variable residual norms (e.g. density).
When the difference is between 3-5 orders of mangitude, the scaling works as expected, and we can successfully solve the coupled linear systems.
However, if the difference is much larger, or smaller than ``1e4``, the scaling will be off, and the solver will encounter difficulties while solving the coupled linear systems.
In cases where the scaling differs greatly, users can manually set the turbulence scaling constant by setting the ``'turbresscale'`` option in their runscripts.
However, we recommend not modifying this variable, and using the decoupled ANK mode further, as the solver will eventually achieve this scaling where the default parameters work as expected.

.. _interpreting_output:

Interpreting the Output
~~~~~~~~~~~~~~~~~~~~~~~

ADflow prints a number of useful metrics for every nonlinear iteration within a simulation.
Understanding what these mean can be critical, especially with the ANK and NK solvers.
Here, we describe each relevant output with the ANK solver.

* ``Iter Tot``: The cumulative number of linear iterations, and residual evaluations. With the ANK solver, this number is calculated with the total number of GMRES iterations with the linear system for the flow variables in the decoupled mode (or the coupled linear system in the coupled mode), plus the residual evaluations required for the line search algorithm during each nonlinear iteration.

* ``Iter Type``: Solver type used. With the ANK solver, last three characters will always read ``ANK``. The leading characters determine what exactly the solver is doing. A ``*`` indicates that the solver updated the preconditioner during that nonlinear iteration. The additional characters ``C`` and ``S`` stand for coupled and second order modes respectively.

* ``CFL``: CFL number used for this nonlinear iteration. This parameter is only updated with the iterations where we update the preconditioner.

* ``Step``: The relaxation factor used. A step of ``1.0`` means the full update is taken, and any number less than this means that the update was relaxed using that factor.

* ``Lin Res``: The relative convergence achieved with the linear solver. The default linear convergence desired is ``0.05``. However, we limit the number of GMRES iterations for the sake of computational cost, and if the solver runs out of iterations, this number will go above the default value. This means that the linear solution failed, however as long as we solve the linear system to some degree, we can still use the update.

Furthermore, users can print some useful information if they are using the ``turbKSP`` solver with the decoupled ANK solver.
To enable this output, users can use the option ``'ankturbkspdebug': True``.
When enabled, the solver will print information related to the turbulence solver between nonlinear iteration outputs.
The turbulence information is printed first, then ADflow prints the default output.
So the turbulence output, and the following default output belong to the same nonlinear iteration.
The turbulence output will print ``LIN RES, ITER, INITRES, REASON, STEP``, followed by 5 numbers.
The numbers correspond to these variables at each ``turbKSP`` iteration.

* ``LIN RES``: Relative convergence of the linear solver. Note this is only the convergence of the linear system for the turbulence.

* ``ITER``: Number of iterations the linear solver took to reach the prescribed tolerance.

* ``INITRES``: Initial norm of the linear residual. This is only useful for developers.

* ``REASON``: The reason for terminating the linear solver. ``2`` means the desired relative convergence is reached, ``-3`` means the solver ran out of iterations. See `KSPConvergedReason <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPConvergedReason.html>`_ for more details.

* ``STEP``: Relaxation factor used for the update. Similar to the value printed with the default output. Note that this is only the relaxation for the turbulence update, and has nothing to do with the flow update.

Troubleshooting
~~~~~~~~~~~~~~~

The ANK solver is tuned for a typical aeronautical application with transonic conditions.
However, for many cases the solver performance can be improved.
Besides performance improvements, the solver might fail for a range of critical cases, and some troubleshooting might be required.
In this section, we talk about the common failure modes and how to fix them.
Before reading here, users should be familiar with the content presented in :ref:`interpreting_output`, as this will be the main source of information for our decisions.

Step size going to zero
***********************

CFL limit

low ma flows, incompressible etc

turbulence not converging

linear solver not converging



Newton--Krylov
--------------

The Newton--Krylov (NK) solver is the solver we recommend for using the final stages of convergence.
It yields the best nonlinear convergence if the initial guess is close to the *basin of attraction*.
With well behaving cases, it is typical to see the NK solver drop the residual norm by 2--3 orders of magnitude in one nonlinear iteration.
However, if the NK solver is used when the state is away from the solution, the solver will either stall or yield bad performance.
To use the NK solver, users can include ``'usenksolver': True`` in their runscrips.
We also recommend prescribing the relative convergence when ADflow will switch to the NK solver by setting ``'nkswitchtol'`` option.
This is a case dependent parameter, for RANS simulations, a relative convergence of ``1e-4`` would be a good case scenario.
For very difficult cases, this switch can be reduced down to ``1e-8`` levels to achieve reasonable performance with the NK solver.
Typically, a setting around ``1e-5`` and ``1e-6`` will yield good results.
The NK solver does not have as many tunable parameters as the ANK solver, however performance can still be improved.
Most of these parameters are related to the linear solver used with the NK solver, and we will describe these options under :ref:`linear_solver_performance`.
However, we first give a quick description of the solver algorithm, along with the important aspects that enable high nonlinear convergence rates.

Solver Algorithm
~~~~~~~~~~~~~~~~

The NK solver solves the nonlinear system of governing equations by simply using the Newton's method.
This involves solving a large linear system at each iteration to calculate the update to the state vector.
To solve this linear system, we use the GMRES algorithm, which is a Krylov subspace based solver; hence the name Newton--*Krylov*.
All state variables are handled in a coupled way, and we use the default scaling described in the :ref:`turbulence_coupling` section.
The method is equivalent to using Euler's method with an infinite time step, and as a result, we do not have a time step in our linear sytems; the implicit component is only composed of the Jacobian.
We stil use a matrix-based preconditioner based on an approximate Jacobian, however the main driver for the linear solver is the exact matrix-free residual operations.
As a result, we always solve for the exact Jacobian.
After solving for the update, we use a cubic line search by default to guaranee a reduction in the total residual norm.
A number of line search algorithms are available and can be specified with the option ``'nkls'``.
We recommend the default ``'cubic'`` line search, however setting this option to ``'non monotone'`` can help by relaxing the criteria to achieve a decrease in the residual norm, and users can even select ``'none'`` to default the solver to take the full step at each iteration.

Selecting the Linear Solver Tolerance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the most important aspects of the NK implementation in ADflow is how the linear solver tolerance is selected.
To do this, we use a method called Eisenstat--Walker (EW) algorithm (`Eisenstat and Walker, 1996 <http://epubs.siam.org/doi/10.1137/0917003>`_).
In simple words, the main idea is to avoid *over-solving* the linear system at each nonlinear iteration.
The Newton update is calculated using a linearization about the current state.
Therefore, if the state is away from the solution, the nonlinear convergence obtained with the update will be limited.
However, as the state approaches the solution, the Newton update can yield a few orders of magnitude of relative convergence in one iteration.
Besides these, the linear solutions are very expensive, and tighter linear convergence tolerances require large computational efforts due to the size of the problems.
Therefore, we face a trade-off: over-solving the linear system will yield better convergence, however the linear solution will be expensive.
Under-solving the linear system on the other hand will require lower computational effort, however the update vector will not be as accurate, and nonlinear convergence will suffer.

Keeping these outcomes in mind, Eisenstat and Walker developed an algorithm that monitors the linear and nonlinear convergence rates between iterations, and picks the optimal linear solution tolerance for the next nonlinear iteration.
When linear solver performance is not a limiting factor, the algorithm picks large linear solver tolerances that yield fast but inaccurate updates when the state is far from the solution, i.e. nonlinear convergence is not satisfactory.
However, as the nonlinear convergence rates improve, the algorithm picks a lower linear solver tolerance, in turn yielding a more expensive iteration, but a more accurate one at the same time.

The practical outcomes for users is as follows:
The solver will start with the default ``'nklinearsolvetol': 0.3``.
This is an option in ADflow, however if users want to prescribe a constant linear solution tolerance for each iteration, besides setting this option, they need to disable the EW algorithm by setting ``'nkuseew': False``.
If the default ``'nkuseew': True`` is preserved, the solver will only solve the linear system in the first NK iteration to 0.3 relative convergence.
After this step, the solver will monitor nonlinear convergence, and determine the linear solver tolerance for the next nonlinear iteration.
Users can monitor the linear convergence by reading the number under ``Lin Res`` that is printed with the default ADflow output.
The solver picking a lower linear solver tolerance means that the nonlinear convergence was satisfactory, and the performance can be improved with a tighter linear convergence.
This is the desired behavior, and the NK solver will gradually lower the linear solver tolerance.
As the linear solver tolerance gets lower values, each iteration will take more time, but nonlinear convergence between nonlinear iterations will improve.
When this happens, it is typical to see 2--3 orders of magnitude relative nonlinear convergence at each nonlinear iteration.

On the other hand, if the nonlinear convergence is not satisfactory, the solver will pick a larger linear solver tolerance, to avoid over-solving the linear system.
In this case, the solver will pick a linear convergence higher than the previous iteration, for example, the second iteration after the first 0.3 linear convergence will have a higher linear convergence tolerance.
This means that the state is not close to the solution, and the solver prefers to take more low-cost iterations, rather than taking fewer but more expensive ones.
This behavior is usually okay, as the solver will eventually start picking lower linear solver tolerances.
If the solver does not pick a lower linear solution tolerance after a handful of iterations, it is usually better to lower the ``'nkswitchtol'`` and try again.
The ANK solver will handle these *transients* better, and switching to the NK solver later on will help avoiding these issues.

There is a hard coded upper limit for the linear solver tolerance, which is set to ``0.8``.
This means that if the solver is consistently solving the linear system to 0.8 relative convergence, the state is far away from the solution, and the users should try again with a lower ``'nkswitchtol'``.

.. _linear_solver_performance:

Linear Solver Performance
~~~~~~~~~~~~~~~~~~~~~~~~~

All the scenarios described in the previous subsection assumes that the linear solver performance is not a limiting factor, i.e. the prescribed linear solution tolerance is reached on every nonlinear iteration.
However, this is usually not the outcome with difficult cases.
Especially with large overset meshes, the default linear solver might fail to meet the prescribed linear solver tolerances.
This can happen due to a weak or outdated preconditioner, accompanied with the solver running out of GMRES iterations.
To prevent these, users can tweak a number of ADflow options to obtain a stronger linear solver for the NK solver.

As stated in the previous subsection, there is a hard coded upper limit on linear solver tolerance, which is 0.8.
If the ``Lin Res`` outputs go above this value, it means that the linear solver is failing to meet the tolerances.
On the other hand, the linear solver can fail as the EW algorithm picks lower linear tolerances.
This case is usually okay, however users can monitor the health of the linear solver by observing the change under ``Iter Tot`` output.
This output prints the cumulative number of linear iterations.
If the change in ``Iter Tot`` between nonlinear iterations are larger than the specified linear iteration limit for the NK solver, i.e. ``'nksubspacesize'``, the solver linear solver is failing to reach the prescribed tolerance.

To obtain a stronger linear solver, there are a number of options.
Each option either increases memory requirements, CPU usage (more operations), or both.
Most likely, the adjoint solver will be the bottleneck in terms of memory usage, and users can read the :ref:`adflow_performance` section to get some estimates.
As a result, users will have bit of room to improve the linear solver used with the NK solver, as the default memory usage will be less than the adjoint solver.

Talk about options.

Troubleshooting
~~~~~~~~~~~~~~~

linear solver failing
linear solver not picking a tighter tolerance
step size going to zero (stall)
