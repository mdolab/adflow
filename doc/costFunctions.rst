.. _adflow_costFunctions:

Cost Functions
==============

ADflow works in three dimensions (:math:`x, y, z`) and application of angles such as ``alpha`` and ``beta`` change the inflow angle, not the mesh itself.
For example a high lift airfoil at :math:`\alpha=20^{\circ}` has the same mesh as that same airfoil at :math:`\alpha=0^{\circ}`.
Another example is a moth T-foil at a :math:`\beta=5^{\circ}` leeway angle versus :math:`\beta=0^{\circ}`, the meshes are the same.

The relevant subroutines for the evaluation of cost functions are located in ``surfaceIntegrations.F90``, ``zipperIntegrations.F90``.
The primary cost functions that are computed in the fortran code are defined in the ``adflowCostFunctions`` dictionary in ``pyADflow.py``.
The following list describes each primary cost function and specifies its units.

.. costfunctionslist:: adflow.pyADflow.ADFLOW


References
==========

.. bibliography::
    :style: unsrt
