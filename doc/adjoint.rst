.. _adflow_adjoint:

Adjoint overview
!!!!!!!!!!!!!!!!

This section is based on a session given by Gaetan Kenway and Charles Mader on 03-22-2017.

Derivative seed manipulation
============================

The workhorse routines for derivative seed manipulation are the master routines defined in ``src/adjoint/masterRoutines.F90``. These subroutines are ``master``, ``master_d``, and ``master_b``.

MASTER subroutine
-----------------
Simple diagram of showing input and outputs to the master routine. Note that these inputs are implicit to the master routine::

             +--------------------+
       w --->|                    |---> R (dw)
             |                    |
       x --->|                    |---> forces (or tractions)
             |       MASTER       |
    x_bc --->|                    |---> funcs
             |                    |
 x_extra --->|                    |
             +--------------------+
                 ^
                 |
              fam_list


The inputs to the master subroutine are all variables that matter for the flow solution. They are:

* ``w``: the flow states
* ``x``: the volume node coordinates
* ``x_bc`` : the coordinates of nodes in boundary conditions such as pressure and temperature faces;
* ``x_extra`` : aeroProblem variables such as alpha, Mach, reference point coordinates, and free-stream pressure, temperature, and density;
* ``fam_list`` : list of families to subdivide the surface integration of output quantities.

All these input variable, with the exception of fam_list, are active variables in the sense that they should be differentiated.  Every function get a familyList associated with it. The default one is used except this is defined specifically.

The outputs of the master function are the expected results from a flow solver. They are:

* ``R`` (dw in Fortran) : Flow solution residuals in every cell;
* ``forces`` : Nodal forces or tractions (forces/area), depending on the ``forcesAsTractions`` flag. This includes the zipper mesh nodes.
* ``funcs`` : ADflow functions requested by the user. They can be integrated quantities such as 'lift', 'cd', and 'mass flow', and the can be defined in arbitrary integration planes. User-defined functions are not included here since they are handled in the Python side.

The master function should perform all relevant operation to receive the solver inputs and compute the corresponding outputs. In other words, it is a collection of functions that will "reproduce" the solver residual and function calculation. Note the call order can be different than in analysis code. There are other calls as well that can be considered as preprocessing calls.
This code is not used during the actual flow solution, but it serves as a reference to develop the differentiated code, which should be done by hand.

In order to perform this task, the master routine should practically touch all ADflow modules. This is why it has such an extensive "use" list. If someone implements new modules to ADflow that changes the residual computation, they **should** also be added to this master subroutine, and its differentiated versions.

Preprocessing routines (such as ones that compute face areas and cell volumes) are included in the master routine since they will affect the linearization.

MASTER_B subroutine
-------------------

This is the reverse AD version of the master routine. Note that these differentiated versions are done manually for two main reasons:

* Tapenade cannot handle pointers such as the ones in ``setPointers`` and communication calls (MPI and PETSc).
* We can get a more efficient code by avoiding local calls to the non-linear code to push/pop values during the reverse AD run, since we already stored most of the necessary information during the forward pass.

.. important::
       The linearization of the low speed preconditioner is not implemented because it needs the residuals before the preconditioner update. However, we do not store this during the forward pass. All other residual operations are incremental (they just add a value to the residual), thus the reverse code does not need any intermediate residual values to accumulate seeds.

.. note::
       ``master_b`` is highly optimized to call only the necessary functions to get only the requested derivatives.

MASTER_STATE_B subroutine
-------------------------

An automatically differentiated code in reverse mode  (e.g. using Tapenade) can be up to 4x-6x slower than the original code due to the additional amount of data, the push/pop calls, and also the local forward pass calls to recompute intermediate states. The adjoint solution subroutine calls the reverse code several times (to compute dRdw vector product), therefore it needs to be as fast as possible.

This subroutine is an optimized version of master_b to compute dRdw products assuming frozen spatial variables (nodal volume coordinates, or x). This is the "fast" reverse code used during the solution of the adjoint systems. This subroutine also calls specialized "fast" versions of the differentiated modules, which also use the frozen spatial variables assumption (so they have fewer active variables during the differentiation).

Python never sees the fast mode. Only the Fortran adjoint solver sees it.

Python interface for derivatives
================================

The main Python methods for derivative seed manipulation in ADflow are ``computeJacobianVectorProductFwd`` and ``computeJacobianVectorProductBwd``, both located in ``pyADflow.py``. These methods of the ADflow object eventually call the master routines defined in Fortran.

ComputeJacobianVectorProductBwd
-------------------------------

This function can receive and return multiple subsets of seeds.

The :math:`\bar{X}_v` seeds are typically only used for derivative verification. The :math:`\bar{X}_s` seeds are only available if we include the mesh object in the solver (``ADFLOW.setMesh``). The :math:`\bar{X}_{Dv}` seeds need both the mesh object and the DVGeo object (``ADFLOW.setDVGeo``), also remember that :math:`\bar{X}_{Dv}` includes :math:`\bar{X}_{DvAero}` (design variables coming from the AeroProblem definition).

All ADflow's native cost functions (such as ``cl``, ``cd``, ``lift``, ``drag``, ...) are applied at all design families at first, but we can specify family subsets to compute these functions (if you want, for instance, a drag breakdown for wing and fuselage). One of the first procedures in the ``ComputeJacobianVectorProductBwd`` is the assembly of ``funcsBar``, which is a matrix that contains derivative seeds for all functions (columns) and every family subset (rows). This ensures that we will have the smallest number of reverse passes to get the derivatives.

The user can give a custom function that uses the intrinsic ADflow functions. ADflow will use complex-step to get the Jacobian of this user-defined function. This may reduce the number of adjoints. For instance, instead of solving one adjoint for cl and other for cd to get the L/D sensitivity, we can solve an adjoint for cl/cd directly::

    def userFunc(funcs):
        funcs['LD'] = funcs['cl']/funcs['cd']

The boundary condition definition and sensitivities are stored under ``AeroProblem`` to facilitate the use of different boundary conditions for multipoint cases.

The reverse master (``master_b``) routine is wrapped by a thin wrapper which is the ``computeMatrixFreeProductBwd`` subroutine defined in ``adjointAPI.F90``. This wrapper allocates memory etc. before the derivatives are calculated.

ComputeJacobianVectorProductFwd
-------------------------------

This function is the 'twin' of ``ComputeJacobianVectorProductBwd`` and can be used to calculate the derivatves using the tangent AD code (i.e. forward code) instead of the usual reverse (adjoint) AD code.
This is useful since computeJacobianVectorProductFwd() calls the function master_d() which is the hand-written recipe that calls all the differentiated code.
This means one can check master_d() against master_b() and verify that they produce consistent gradients.

Development tips:
.................

There is no current method that calculates the sensitivities with the forward method in ``pyADflow.py``, but a previous partial implementation to verify rotating frame derivatives can be found in `the early commits of this PR <https://github.com/mdolab/adflow/pull/121>`_.
The function should be identical to ``evalFunctionsSens()`` up to the point where the vector-product function is called:

.. code-block:: python

      for keys_, val_ in zip(self.DVGeo.getValues().keys(), self.DVGeo.getValues().values()):
                  residualDeriv_, funcDeriv_ = self.computeJacobianVectorProductFwd(
                        xDvDot={keys_: val_}, funcDeriv=True, residualDeriv=True
                  )

Note that the forward mode indeed requires a loop over all the design variables (only the DVGeo ones are reported in the example above).

Once ``ComputeJacobianVectorProductFwd()`` is called, we have all terms in the total derivative equation - i.e. the functions and residual partials, and the reverse mode seeds from ``self.getAdjoint()``.
However, before we construct the total derivative we must remember to call ``MPI_SUM`` on psi and residualDeriv with something like:

.. code-block:: python
            
            psi_dot_rDeriv_ = numpy.dot(psi, residualDeriv_)
            AllReduced_psi_dot_rDeriv = self.comm.allreduce(psi_dot_rDeriv_, op=MPI.SUM)

We do *NOT* ``MPI_SUM`` the funcDeriv since that has already been done in the hand-differentiated code from ``surfaceIntegrations.F90``in the subroutine called ``getSolution_d()``.
Now, we can finally compute the total sensitivity for every function of interest with:

.. code-block:: python

            TotalDeriv_ = funcDeriv_[myKey] + AllReduced_psi_dot_rDeriv

Note that the final derivatives must be processed with ``_processAeroDerivatives()`` into a usable dictionary.

Notes on dot-product tests
--------------------------
This example below uses the subroutine ``slipvelocitiesfineleve_block_d()`` which can be found in the file ``src/adjoint/outputForward/adfjointextra_d.f90``.
First we must insert a new entry in the ``adflow.pyf`` file.
To this end we cd to the directory for ``slipvelocitiesfineleve_block_d()``:

.. prompt:: bash

      $ cd /src/adjoint/outputForward/

and then we use f2py to generate the python interface:

.. prompt:: bash

      $ f2py -h TEMPORARY_FILE.pyf adjointextra_d.f90

Now we go into ``TEMPORARY_FILE.pyf`` and find the lines relevant for the routine of interest.
In this case we need the following lines:

.. code-block:: Fortran

      module adjointextra_d
        subroutine slipvelocitiesfinelevel_block_d(useoldcoor,t,sps) ! in adjointextra_d.f90:adjointextra_d
          use blockpointers
          use inputmotion
    
      #         A  L O T  O F   L I N E S  .....

          logical intent(in) :: useoldcoor
          real(kind=realtype) dimension(*),intent(in) :: t
          integer(kind=inttype) intent(in) :: sps
        end subroutine slipvelocitiesfinelevel_block_d
      end module adjointextra_d

Now we insert these in the ``/src/f2py/adflow.pyf file``.
Using the above defined function here in ``pyADflow.py`` we can now (thx to the changes we made in the ``f2py/pyADflow.pyf``) finally setup and call the dot test from the python layer.
