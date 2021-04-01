.. _Verification:

Verification and Validation
===========================
This guide is inteded to be used to document major ``ADflow`` cases for validation and verification purposes.


Overset and Multiblock Meshes
-----------------------------
The `following paper <https://arc.aiaa.org/doi/pdf/10.2514/1.C034486>`_ compares ``ADflow`` with OVERFLOW and elsA for the CRM and a Wing-Body-Nacelle-Pylon configuration. 
Both overset and multiblock structured meshes were used in the comparison. The cases are run at transonic mach numbers and Reynolds Numbers of ``Re = 5e6``. Figure 8 of in this paper details drag grid convergence of each of the flow solvers.

.. code-block:: bibtex

    @incollection{Coder2017a,
        author = {Coder, James G. and Pulliam, Thomas H. and Hue, David and Kenway, Gaetan K. W. and Sclafani, Anthony J.},
        booktitle = {AIAA SciTech Forum},
        doi = {10.2514/6.2017-0960},
        month = {January},
        publisher = {American Institute of Aeronautics and Astronautics},
        title = {Contributions to the 6th {AIAA} {CFD} {Drag Prediction Workshop} Using Structured Grid Methods},
        year = {2017}
    }

Hydrofoils
----------
The following `paper <https://www.sciencedirect.com/science/article/abs/pii/S088997461630144X?via%3Dihub>`_ couples hydrodynamic and structural analysis of hydrofoils and compares the resutls with experimental data.
The CFD analysis uses a low-speed preconditioner to solve hydrodynamic problems that have nearly incompressible flow. A NACA 0009 hydrofoil is analyzed at a Reynolds Number of ``Re = 1e6`` and compared to experimental results.
Figure 3 of this paper shows the comparison between the predicted force and moment coefficients as well as predicted tip deflection of ``ADflow`` with the experimental results. 
Figure 4 of this paper shows a drag convergence study and the results approach experimental values.

.. code-block:: bibtex

    @article{Garg2017a,
        author = {Nitin Garg and Gaetan K. W. Kenway and Joaquim R. R. A. Martins and Yin Lu Young},
        doi = {10.1016/j.jfluidstructs.2017.02.001},
        journal = {Journal of Fluids and Structures},
        month = {May},
        pages = {15--39},
        title = {High-fidelity Multipoint Hydrostructural Optimization of a {3-D} Hydrofoil},
        volume = {71},
        year = {2017}
    }


Wind Turbines
-------------
Analysis of wind turbines is conducted in `this paper <http://www.umich.edu/~mdolaboratory/pdf/Madsen2019a.pdf>`_. Appendix A compares ``ADflow`` to EllipSys3D at different flow conditions.
Figure 21 shows the Thrust and Torque calculations from ``ADflow`` and EllipSys3D. These results show ``ADflow`` consistently overshoots the EllipSys3D results at all wind speeds.

.. code-block:: bibtex

    @article{Madsen2019a,
        author = {Mads H. Aa. Madsen and Frederik Zahle and Niels N. S{\o}rensen and Joaquim R. R. A. Martins},
        doi = {10.5194/wes-4-163-2019},
        journal = {Wind Energy Science},
        keywords = {ank},
        month = {April},
        pages = {163--192},
        title = {Multipoint high-fidelity {CFD}-based aerodynamic shape optimization of a 10 {MW} wind turbine},
        volume = {4},
        year = {2019}
    }

ONERA M6
--------
The MACH-Aero tutorial contains a comparison of ADflow with other flow solvers for the ONERA M6 wing `here <https://github.com/mdolab/MACH-Aero/blob/2a17dc93f0d2ddfa779ee78e1105a56a71026a62/machAeroTutorials/overset_analysis.rst#id5>`_. 
The results section at the bottom of this page shows force and moment coefficient convergence of ``ADflow`` and other solvers. 
