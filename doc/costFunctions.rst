.. _adflow_costFunctions:

Cost Functions
==============

ADflow works in three dimensions (x,y,z) and application of angles such as alpha and beta change the inflow angles, not the mesh itself.
So for example a high lift airfoil at 20 degrees has the same mesh as that same airfoil at 0 degrees.
Another example is a moth T-foil at a 10 degree beta (leeway) angle versus zero.

The relevant source code file for the evaluation of cost functions is ``surfaceIntegrations.F90``.
Cost functions are added to the ``adflowCostFunctions`` dictionary in ``pyADflow.py``.
The following list is in the same order as how it currently is in the source code.

.. costfunctionslist:: adflow.pyADflow.ADFLOW


References
==========

.. bibliography::
    :style: unsrt
