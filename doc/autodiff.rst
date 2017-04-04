.. adflow_autodiff

Differentiation overview
!!!!!!!!!!!!!!!!!!!!!!!!


Automatic Differentiation (AD)
==============================

When we run Tapenade in reverse mode, we explicitly ask to use "d" to flag derivative seeds, instead of the usual "b" identifier. We do this so that we can use the same variable definition for both forward and reverse modes. Note, however, that the subroutine names still have the "_b" for their reverse counterparts.

We should call Tapenade in a single command to generate all AD code. Having a single head will ensure that we get all dependencies correct.

We tell Tapenade to treat almost variables as input and output at the same time so that it does not zero the seeds after the subroutine calls. This saves time since zeroing vectors is a costly operation. However, we need to make sure that all seeds are set to zero at the beginning of the reverse AD pass. This is done in ``src/adjoint/adjointUtils.F90`` at the subroutine zeroADSeeds. This is done in ``adjointAPI.F90``.

If developers make important changes in the current ADflow modules, then the differentiated code should be updated as well. Use the following steps to run the automatic differentiation tools:

* Open a terminal at the src/adjoint folder
* Run the following commands::

    $ make -f Makefile_tapenade ad_forward
    $ make -f Makefile_tapenade ad_reverse
    $ make -f Makefile_tapenade ad_reverse_fast

Manual Differentiation
======================

The generation of the reverse code in ADflow is not fully automatic. Some subroutines need to be manually differentiated for two main reasons:

* Tapenade cannot handle pointers and communication calls (MPI and PETSc).
* We can get a more efficient code by avoiding local calls to the non-linear code to push/pop values during the reverse AD run, since we already stored most of the necessary information during the forward pass.

Examples of manually differentiated subroutines are master, ``haloExchange``, ``surfaceIntegration``, and ``getSolution``.

If developers change currently existing base modules, there is no need for manual differentiation: run the automatic differentiation command as mentioned above.

If new modules are added, or if the residual computation procedure drastically changes, then manual differentiation will be necessary (specially in ``master``, ``master_d``, and ``master_b``).

.. note::
       If routine is differenced by hand it should be labeled clearly at the top of that routine that it is manually differentiated. These routines however mostly call code that has been differentiated by Tapenade.

Differentiation of communication calls
---------------------------------------

Tapenade does not differentiate functions with MPI and PETSc calls. Therefore, we need to manually differentiate them, that is, for every MPI call in the non-linear (original) code, we need to add a corresponding MPI or PETSc call to transfer seeds. Here are some examples of the reverse counterparts of some communication procedures.

* MPI Sends become MPI Receives, and vice-versa.
* MPI Reduces become MPI Broadcasts, and vice-versa.
* All PETSc forward scatters with insert flag become PETSc reverse scatters with accumulate flag.
