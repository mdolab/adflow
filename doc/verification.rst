.. _verification:

Verification and Validation
===========================
This guide is intended to be used to document major ADflow cases for validation and verification purposes.


Overset and Multiblock Meshes
-----------------------------
The following paper :cite:`Coder2017a` compares ADflow with OVERFLOW and elsA for the CRM and a Wing-Body-Nacelle-Pylon configuration. 
Both overset and multiblock structured meshes were used in the comparison. The cases are run at transonic mach numbers and Reynolds Numbers of ``Re = 5e6``. Figure 8 of in this paper details drag grid convergence of each of the flow solvers.

Hydrofoils
----------
The following paper :cite:`Garg2017a` couples hydrodynamic and structural analysis of hydrofoils and compares the results with experimental data.
The CFD analysis uses a low-speed preconditioner to solve hydrodynamic problems that have nearly incompressible flow. A NACA 0009 hydrofoil is analyzed at a Reynolds Number of ``Re = 1e6`` and compared to experimental results.
Figure 3 of this paper shows the comparison between the predicted force and moment coefficients as well as predicted tip deflection of ADflow with the experimental results. 
Figure 4 of this paper shows a drag convergence study and the results approach experimental values.

Wind Turbines
-------------
Analysis of wind turbines is conducted in this paper :cite:`Madsen2019a`. Appendix A compares ADflow to EllipSys3D at different flow conditions.
Figure 21 shows the Thrust and Torque calculations from ADflow and EllipSys3D. These results show ADflow consistently overshoots the EllipSys3D results at all wind speeds.

ONERA M6
--------
The MACH-Aero tutorial contains a comparison of ADflow with other flow solvers for the ONERA M6 wing :ref:`here <mach-aero:overset_analysis>`. 
The results section at the bottom of this page shows force and moment coefficient convergence of ADflow and other solvers. 


.. bibliography::
    :style: unsrt
