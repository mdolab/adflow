!
!      ******************************************************************
!      *                                                                *
!      * File:          createTemplate.f90                              *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 12-12-2002                                      *
!      * Last modified: 11-28-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine createTemplate
!
!      ******************************************************************
!      *                                                                *
!      * create_template creates a template parameter file.             *
!      * This file contains all possible input parameters and must be   *
!      * adjusted by the user.                                          *
!      *                                                                *
!      ******************************************************************
!
       use allInputParam
       implicit none
!
!      local variables
!
       integer, parameter :: writeUnit = 35

       integer :: ios

       character(len=2*maxStringLen) :: string
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Open the file for writing and checks if it was succesful. If not
       ! terminate the program without writing the file.

       open(unit=writeUnit, file=paramFile, status="new", &
            action="write", iostat=ios)
       if(ios /= 0) then

         write(string,"(3a)") "Template parameter file ", &
                              trim(paramFile),            &
                              " could not be opened for writing"
         call terminate("createTemplate", string)

       endif

       ! Write the header of the file

       write(writeUnit,"(a)") "=========================================&
                              &======================================"
       write(writeUnit,"(a)") "This is an automatically created &
                              &template parameter file for SUmb."
       write(writeUnit,"(a)") "It contains all the possible parameters &
                              &in several sections."
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "This file must be adjusted by the user."
       write(writeUnit,"(a)") "Comments are indicated by #. All info &
                              &following a # is ignored."
       write(writeUnit,"(a)") "The sequence of the parameters is &
                              &arbitrary and keywords are case &
                              &insensitive."
       write(writeUnit,"(a)") "When a keyword occurs multiple times, &
                              &the last value is taken."
       write(writeUnit,"(a)") "=========================================&
                              &======================================"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the IO parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     IO Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "                       File format: CGNS"
       write(writeUnit,"(a)") "               # Other possibility: &
                              &PLOT3D"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "                         Grid file: &
                              &MISSING FILE NAME"
       write(writeUnit,"(a)") "          PLOT3D Connectivity file: &
                              &MISSING FILE NAME"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "                      Restart file: &
                              &MISSING FILE NAME"
       write(writeUnit,"(a)") "                           Restart: no"
       write(writeUnit,"(a)") "       Check nondimensionalization: yes"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "                     New grid file: &
                              &NewGrid.cgns"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "                     Solution file: &
                              &SolSUmb.cgns"
       write(writeUnit,"(a)") "             Surface solution file: &
                              &SolSUmb_surface.cgns"
       write(writeUnit,"(a)") "      Rind layer in solution files: no"
       write(writeUnit,"(a)") "        Write coordinates in meter: no"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "        Automatic parameter update: yes"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "                 Cp curve fit file: &
                              &MISSING FILE NAME"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "       Write precision grid: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "       # Default is executable precision"
       write(writeUnit,"(a)") "            # Possibilities: single"
       write(writeUnit,"(a)") "            #              : double"
        write(writeUnit,"(a)") "  Write precision solution: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "       # Default is executable precision"
       write(writeUnit,"(a)") "            # Possibilities: single"
       write(writeUnit,"(a)") "            #              : double"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "Store convergence inner iterations: no"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the physics parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Physics Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "                  Equations: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "            # Possibilities: Euler"
       write(writeUnit,"(a)") "            #              : Laminar NS"
       write(writeUnit,"(a)") "            #              : RANS"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                       Mode: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "            # Possibilities: Steady"
       write(writeUnit,"(a)") "            #              : Unsteady"
       write(writeUnit,"(a)") "            #              : Time spectral"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                  Flow type: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "            # Possibilities: Internal flow"
       write(writeUnit,"(a)") "            #              : External flow"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                   Cp model: Constant"
       write(writeUnit,"(a)") "        # Other possibility: &
                              &Temperature curve fits"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "           Turbulence model: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "            # Possibilities: Baldwin Lomax"
       write(writeUnit,"(a)") "            #              : &
                              &Spalart Allmaras"
       write(writeUnit,"(a)") "            #              : &
                              &Spalart Allmaras Edwards"
       write(writeUnit,"(a)") "            #              : KOmega Wilcox"
       write(writeUnit,"(a)") "            #              : &
                              &KOmega Modified"
       write(writeUnit,"(a)") "            #              : KTau"
       write(writeUnit,"(a)") "            #              : Menter SST"
       write(writeUnit,"(a)") "            #              : V2F"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "     V2F version (n1 or n6): 1"
       write(writeUnit,"(a)") "        # Other possibility: 6"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "       V2F with upper bound: yes"
       write(writeUnit,"(a)") "        # Other possibility: no"
       write(writeUnit,"(a)")


       write(writeUnit,"(a)") " Turbulence production term: Strain"
       write(writeUnit,"(a)") "      # Other possibilities: Vorticity"
       write(writeUnit,"(a)") "            #              : Kato-Launder"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                Use wall functions: no"
       write(writeUnit,"(a)") "Offset from wall in wall functions: 0.0"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "   Constant specific heat ratio: 1.4    #Air"
       write(writeUnit,"(a)") "   Gas constant (J/(kg K))     : 287.87 #Air"
       write(writeUnit,"(a)") "   Prandtl number              : 0.72   #Air"
       write(writeUnit,"(a)") "   Turbulent Prandtl number    : 0.90"
       write(writeUnit,"(a)") "   Max ratio k-prod/dest       : 20.0"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Free Stream Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "                             Mach: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "            Mach for coefficients: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "                          # Default is Mach"
       write(writeUnit,"(a)") "                         Reynolds: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "       Reynolds length (in meter): 1.0"
       write(writeUnit,"(a)") "   Free stream velocity direction: &
                              &1.0 0.0 0.0"
       write(writeUnit,"(a)") "                   Lift direction: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "     # Default is normal to free stream &
                                & without y-component"
       write(writeUnit,"(a)") "   Free stream temperature (in K): 288.15"
       write(writeUnit,"(a)") " Free stream eddy viscosity ratio: 0.01"
       write(writeUnit,"(a)") "  Free stream turbulent intensity: 0.001"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Reference State"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "            Reference pressure (in Pa): &
                              &101325.0"
       write(writeUnit,"(a)") "         Reference density (in kg/m^3): 1.25"
       write(writeUnit,"(a)") "          Reference temperature (in K): &
                              &273.15"
       write(writeUnit,"(a)") " Conversion factor grid units to meter: 1.0"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Geometrical Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "           Reference surface: 1.0"
       write(writeUnit,"(a)") "            Reference length: 1.0"
       write(writeUnit,"(a)") "    Moment reference point x: 0.0"
       write(writeUnit,"(a)") "    Moment reference point y: 0.0"
       write(writeUnit,"(a)") "    Moment reference point z: 0.0"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the discretization
       ! parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Fine Grid Discretization Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "       Discretization scheme: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "             # Possibilities: &
                              &Central plus scalar dissipation"
       write(writeUnit,"(a)") "             #              : &
                              &Central plus matrix dissipation"
       write(writeUnit,"(a)") "             #              : &
                              &Central plus CUSP dissipation"
       write(writeUnit,"(a)") "             #              : Upwind"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "   Order turbulent equations: First order"
       write(writeUnit,"(a)") "         # Other possibility: Second order"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "              Riemann solver: Roe"
       write(writeUnit,"(a)") "       # Other possibilities: Van Leer"
       write(writeUnit,"(a)") "       #                    : AUSMDV"

       write(writeUnit,"(a)") "   Limiter                  : No limiter"
       write(writeUnit,"(a)") "       # Other possibilities: First order"
       write(writeUnit,"(a)") "       #                    : Van Albeda"
       write(writeUnit,"(a)") "       #                    : Minmod"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "              Preconditioner: &
                              &No preconditioner"
       write(writeUnit,"(a)") "       # Other possibilities: Turkel"
       write(writeUnit,"(a)") "       #                    : Choi Merkle"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "     Wall boundary treatment: Normal momentum"
       write(writeUnit,"(a)") "       # Other possibilities: &
                              &Constant pressure"
       write(writeUnit,"(a)") "       #                    : &
                              &Linear extrapolation pressure"
       write(writeUnit,"(a)") "       #                    : &
                              &Quadratic extrapolation pressure"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "  Outflow boundary treatment: &
                              &Constant extrapolation"
       write(writeUnit,"(a)") "        # Other possibility: &
                              &Linear extrapolation"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "non-matching block to block treatment: &
                              &NonConservative"
       write(writeUnit,"(a)") "                  # Other possibility: &
                              &Conservative"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                           Vis2: 0.5"
       write(writeUnit,"(a)") "                           Vis4: &
                              &0.015625  # 1/64"
       write(writeUnit,"(a)") "Directional dissipation scaling: yes"
       write(writeUnit,"(a)") "   Exponent dissipation scaling: 0.0"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "   Total enthalpy scaling inlet: no"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "      Kappa interpolation value: 0.33333"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "              Vortex correction: no"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the unsteady
       ! parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Unsteady Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"

       write(writeUnit,"(a)") "     Time integration scheme: BDF"
       write(writeUnit,"(a)") "       # Other possibilities: explicit &
                              &Runge-Kutta"
       write(writeUnit,"(a)") "       # Other possibilities: implicit &
                              &Runge-Kutta"

       write(writeUnit,"(a)") "      Time accuracy unsteady: Second"
       write(writeUnit,"(a)") "       # Other possibilities: &
                              &First to Fifth"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "Number of unsteady time steps coarse grid: &
                              &-1  # Means same as on fine grid"
       write(writeUnit,"(a)") "  Number of unsteady time steps fine grid: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "              Unsteady time step (in sec): &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "      Update wall distance unsteady mode: &
                              &yes"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the time spectral
       ! parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Time Spectral Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "             Number time intervals spectral: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "            Write file for unsteady restart: &
                              &no"
       write(writeUnit,"(a)") "    Time step (in sec) for unsteady restart: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "       Write unsteady volume solution files: &
                              &no"
       write(writeUnit,"(a)") "      Write unsteady surface solution files: &
                              &no"
       write(writeUnit,"(a)") "          Number of unsteady solution files: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the iteration
       ! parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Iteration Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "                                   Smoother: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "                            # Possibilities: &
                              &Runge Kutta"
       write(writeUnit,"(a)") "                            #              : &
                              &Nonlinear LUSGS"
       write(writeUnit,"(a)") "                            #              : &
                              &Nonlinear LUSGS Line"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "               Number of Runge Kutta stages: &
                              &5"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "              Treatment turbulent equations: &
                              &Segregated"
       write(writeUnit,"(a)") "                        # Other possibility: &
                              &Coupled"
       write(writeUnit,"(a)") "    Number additional turbulence iterations: &
                              &0"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                         Turbulent smoother: &
                              &ADI"
       write(writeUnit,"(a)") "                        # Other possibility: &
                              &GMRES"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                        Update bleeds every: &
                              &50"
       write(writeUnit,"(a)") "Relaxation factor bleed boundary conditions: &
                              &0.1"

       write(writeUnit,"(a)") "                         Residual averaging: &
                              &all stages"
       write(writeUnit,"(a)") "                      # Other possibilities: &
                              &no"
       write(writeUnit,"(a)") "                      #                    : &
                              &alternate stages"
       write(writeUnit,"(a)") "     Residual averaging smoothing parameter: &
                              &1.5"

       write(writeUnit,"(a)") "                 Number of multigrid cycles: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "   Number of single grid startup iterations: &
                              &0"
       write(writeUnit,"(a)") "                                 Save every: &
                              &0"
       write(writeUnit,"(a)") "                         Save surface every: &
                              &0"
       write(writeUnit,"(a)") "                                 CFL number: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "                       Turbulent relaxation: &
                              &MISSING PARAMETER"
       write(writeUnit,"(a)") "                            # Possibilities: &
                              &Explicit"
       write(writeUnit,"(a)") "                            #              : &
                              &Implicit"
       write(writeUnit,"(a)") "                     Alpha turbulent DD-ADI: &
                              &0.8"
       write(writeUnit,"(a)") "                      Beta turbulent DD-ADI: &
                              &-1  # Same as alpha"
       write(writeUnit,"(a)") "           Relative L2 norm for convergence: &
                              &1.e-6"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Multigrid Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "      Number of multigrid cycles coarse grid: &
                              & -1 &
                                          & # Means same as on fine grid"
       write(writeUnit,"(a)") "                      CFL number coarse grid: &
                              &-1.0  # Means same as on fine grid"
       write(writeUnit,"(a)") "Relative L2 norm for convergence coarse grid: &
                              &1.e-2"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "        Discretization scheme coarse grid: &
                             &MISSING PARAMETER  # Default fine grid scheme"
       write(writeUnit,"(a)") "               Riemann solver coarse grid: &
                             &MISSING PARAMETER  # Default fine grid solver"
       write(writeUnit,"(a)") "                         Vis2 coarse grid: &
                              &0.5"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "        Freeze turbulent source terms in MG: &
                              &yes"
       write(writeUnit,"(a)") "                        # Other possibility: &
                              &no"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") " Treatment boundary multigrid corrections: &
                              &Zero Dirichlet"
       write(writeUnit,"(a)") "            Restriction relaxation factor: 1.0"
       write(writeUnit,"(a)") "                    Multigrid start level: &
                              &MISSING PARAMETER  # Default is coarsest &
                              &MG level"
       write(writeUnit,"(a)") "                 Multigrid cycle strategy: sg"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the overset parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Overset Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "      Input overset donors are guesses: no"
       write(writeUnit,"(a)") "Average restricted residual for blanks: no"
       write(writeUnit,"(a)") "            Overset interpolation type: &
                              &TriLinear"
       write(writeUnit,"(a)") "Overset interpolation type coarse grid: &
                              &TriLinear"
       write(writeUnit,"(a)") "              Allowable donor quality: 1.0"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the coupler parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Coupler Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"

       write(writeUnit,"(a)") "                            Code name: &
                              &SUmb"
       write(writeUnit,"(a)") "                 Get coarse-level sol: &
                              &no"
       write(writeUnit,"(a)") "              Mach for initialization: &
                              &0.5"
       write(writeUnit,"(a)") "          Pressure for initialization: &
                              &101325.0"
       write(writeUnit,"(a)") "           Density for initialization: &
                              &1.2"
       write(writeUnit,"(a)") "Velocity direction for initialization: &
                              &1.0 0.0 0.0"
       write(writeUnit,"(a)")

       ! Write the keywords and default values for the parallel, i.e.
       ! load balance parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Load balancing Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "        Allowable load imbalance: 0.1"
       write(writeUnit,"(a)") "   Split blocks for load balance: yes"
       write(writeUnit,"(a)")

       ! Write the visualization parameters.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Visualization Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "          pV3 visualization only: no"
       write(writeUnit,"(a)")

       ! Write the keywords and example values for the grid motion.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Grid motion Parameters"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"

       write(writeUnit,"(a)") "     Rotation point body (x,y,z): 0.0 0.0 0.0"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "    Degree polynomial x-rotation: 0"
       write(writeUnit,"(a)") "    Degree polynomial y-rotation: 0"
       write(writeUnit,"(a)") "    Degree polynomial z-rotation: 1"
       write(writeUnit,"(a)")
       
       write(writeUnit,"(a)") "Polynomial coefficients x-rotation: 0.0"
       write(writeUnit,"(a)") "Polynomial coefficients y-rotation: 0.0"
       write(writeUnit,"(a)") "Polynomial coefficients z-rotation: 0.0 1.e-3"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "       Degree fourier x-rotation: 0"
       write(writeUnit,"(a)") "       Degree fourier y-rotation: 0"
       write(writeUnit,"(a)") "       Degree fourier z-rotation: 1"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "        Omega fourier x-rotation: 0.25"
       write(writeUnit,"(a)") "        Omega fourier y-rotation: 0.32"
       write(writeUnit,"(a)") "        Omega fourier z-rotation: 0.41"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "Fourier cosine coefficients x-rotation: 0.0"
       write(writeUnit,"(a)") "Fourier cosine coefficients y-rotation: 0.0"
       write(writeUnit,"(a)") "Fourier cosine coefficients z-rotation: 0.0 0.0"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "Fourier sine coefficients z-rotation: 0.1"
       write(writeUnit,"(a)")

       ! Write the monitor, surface output and volume output variables.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Monitoring and output variables"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "                Monitoring variables: &
                              &resrho_cl_cd"
       write(writeUnit,"(a)") " Monitor massflow sliding interfaces: no"
       write(writeUnit,"(a)") "            Surface output variables: &
                              &rho_cp_vx_vy_vz_mach"
       write(writeUnit,"(a)") "             Volume output variables: &
                              &ptloss_resrho"
       write(writeUnit,"(a)")

       ! The section to overwrite the rotation info for the families.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Family rotation info "
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "                               Rotation &
                              &center  Rotation rate (rad/s)"
       write(writeUnit,"(a)") "Rotating family <family_name1> : 0.0 0.0 0.0 &
                              &   1.e+4 0.e+0 0.e+0"
       write(writeUnit,"(a)") "Rotating family <family_name2> : 0.0 0.0 0.0 &
                              &   1.e+3 0.e+0 0.e+0"
       write(writeUnit,"(a)") "Etc."
       write(writeUnit,"(a)")

       ! The section to overwrite the boundary condition data sets
       ! for the families.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Family boundary condition data sets "
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)")
       write(writeUnit,"(a)") "      # Subsonic outflow boundary with &
                                  &varying pressure in radial direction"
       write(writeUnit,"(a)") "Boundary family <family_name1> : CoordinateR &
                              & Pressure"
       write(writeUnit,"(a)") "            Npoints = 3"
       write(writeUnit,"(a)") "                                     0.2     &
                              &   85360  #Pa"
       write(writeUnit,"(a)") "                                     0.6     &
                              &   88900"
       write(writeUnit,"(a)") "                                     1.0     &
                              &   92730"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "      # Subsonic inflow boundary with a &
                                  &constant state"
       write(writeUnit,"(a)") "Boundary family <family_name2> :  &
                              &PressureStagnation TemperatureStagnation &
                              &VelocityAngleX VelocityAngleY VelocityAngleZ"
       write(writeUnit,"(a)") "            Npoints = 1"
       write(writeUnit,"(a)") "                                  &
                              &    112850                  300.0        &
                              &      0.0       1.570796327    1.570796327"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") "      # Isothermal wall boundary with a &
                              &variable temperature"
       write(writeUnit,"(a)") "Boundary family <family_name3> : CoordinateX &
                              &CoordinateY CoordinateZ Temperature"
       write(writeUnit,"(a)") "            Npoints = 2  2  # 2D structured"
       write(writeUnit,"(a)") "                                     0.0     &
                              &    0.0         0.0         312.5  # (1,1)"
       write(writeUnit,"(a)") "                                     0.5     &
                              &    0.0         0.1         315.0  # (2,1)"
       write(writeUnit,"(a)") "                                     0.0     &
                              &    0.5         0.0         310.5  # (1,2)"
       write(writeUnit,"(a)") "                                     0.5     &
                              &    0.5         0.1         313.5  # (2,2)"
       write(writeUnit,"(a)")

       write(writeUnit,"(a)") " # Whether or not to monitor the mass flow&
                              & for certain families."
       write(writeUnit,"(a)") "Boundary family <family_name4> : &
                              & monitor mass flow: yes   #no"
       write(writeUnit,"(a)")

       ! Write the header for updates, such that it is clear where the
       ! updated section of the parameter file starts.

       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)") "     Updates"
       write(writeUnit,"(a)") "-----------------------------------------&
                              &--------------------------------------"
       write(writeUnit,"(a)")

       ! close the file

       close(unit=writeUnit)

       end subroutine createTemplate
