.. _adflow_costFunctions:

Cost Functions
==============

ADflow works in three dimensions (:math:`x, y, z`) and application of angles such as ``alpha`` and ``beta`` change the inflow angles, not the mesh itself.
For example a high lift airfoil at :math:`\alpha=20^{\circ}` has the same mesh as that same airfoil at :math:`\alpha=0^{\circ}`.
Another example is a moth T-foil at a :math:`\beta=5^{\circ}` leeway angle versus :math:`\beta=0^{\circ}`, the meshes are the same.

The relevant source code file for the evaluation of cost functions is ``surfaceIntegrations.F90``.
Cost functions are added to the ``adflowCostFunctions`` dictionary in ``pyADflow.py``.
The following list is in the same order as how it currently is in the source code.

.. costfunctionslist:: adflow.pyADflow.ADFLOW


References
==========

.. bibliography::
    :style: unsrt
