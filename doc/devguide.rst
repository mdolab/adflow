.. _adflow_devguide:

Developers Guide
================

This guide is intended to be used by developers of ADflow. Higher level descriptions of the code itself can be found here.

Adjoint
-------

.. toctree::
    :maxdepth: 2

    adjoint


Differentiation
---------------

.. toctree::
    :maxdepth: 2

    autodiff




Extra Stuff
-----------

* Triangles of zipper mesh live only in the root processor.
* ``surfaceIntegration.F90``: Takes care of forces and flow-through integrations. User-defined surfaces can only be used for flow-through integrations
* We only assemble the full Jacobian for preconditioner. The adjoint is matrix-free. In the future we need a matrix-free preconditioner to avoid memory limitations.
* Overset interpolation is in the ``wOversetGeneric`` subroutine, located in ``src/utils/haloExchange.F90``.
