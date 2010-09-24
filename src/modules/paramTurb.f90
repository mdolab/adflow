!
!      ******************************************************************
!      *                                                                *
!      * File:          paramTurb.f90                                  *
!      * Author:        Edwin van der Weide, Georgi Kalitzin            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module paramTurb
!
!      ******************************************************************
!      *                                                                *
!      * Module that contains the constants for the turbulence models   *
!      * as well as some global variables/parameters for the turbulent  *
!      * routines.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Spalart-Allmaras constants.                                    *
!      *                                                                *
!      ******************************************************************
!
       real(kind=realType), parameter :: rsaK   = 0.41_realType
       real(kind=realType), parameter :: rsaCb1 = 0.1355_realType
       real(kind=realType), parameter :: rsaCb2 = 0.622_realType
       real(kind=realType), parameter :: rsaCb3 = 0.66666666667_realType
       real(kind=realType), parameter :: rsaCv1 = 7.1_realType
       real(kind=realType), parameter :: rsaCw1 = rsaCb1/(rsaK**2) &
                                                  + (1.+rsaCb2)/rsaCb3
       real(kind=realType), parameter :: rsaCw2 = 0.3_realType
       real(kind=realType), parameter :: rsaCw3 = 2.0_realType
       real(kind=realType), parameter :: rsaCt1 = 1.0_realType
       real(kind=realType), parameter :: rsaCt2 = 2.0_realType
       real(kind=realType), parameter :: rsaCt3 = 1.2_realType
       real(kind=realType), parameter :: rsaCt4 = 0.5_realType
!
!      ******************************************************************
!      *                                                                *
!      * K-omega constants.                                             *
!      *                                                                *
!      ******************************************************************
!
       real(kind=realType), parameter :: rkwK     = 0.41_realType
       real(kind=realType), parameter :: rkwSigk1 = 0.5_realType
       real(kind=realType), parameter :: rkwSigw1 = 0.5_realType
       real(kind=realType), parameter :: rkwSigd1 = 0.5_realType
       real(kind=realType), parameter :: rkwBeta1 = 0.0750_realType
       real(kind=realType), parameter :: rkwBetas = 0.09_realType
!
!      ******************************************************************
!      *                                                                *
!      * K-omega SST constants.                                         *
!      *                                                                *
!      ******************************************************************
!
       real(kind=realType), parameter :: rSSTK     = 0.41_realType
       real(kind=realType), parameter :: rSSTA1    = 0.31_realType
       real(kind=realType), parameter :: rSSTBetas = 0.09_realType

       real(kind=realType), parameter :: rSSTSigk1 = 0.85_realType
       real(kind=realType), parameter :: rSSTSigw1 = 0.5_realType
       real(kind=realType), parameter :: rSSTBeta1 = 0.0750_realType

       real(kind=realType), parameter :: rSSTSigk2 = 1.0_realType
       real(kind=realType), parameter :: rSSTSigw2 = 0.856_realType
       real(kind=realType), parameter :: rSSTBeta2 = 0.0828_realType
!
!      ******************************************************************
!      *                                                                *
!      * K-tau constants.                                               *
!      *                                                                *
!      ******************************************************************
!
       real(kind=realType), parameter :: rktK     = 0.41_realType
       real(kind=realType), parameter :: rktSigk1 = 0.5_realType
       real(kind=realType), parameter :: rktSigt1 = 0.5_realType
       real(kind=realType), parameter :: rktSigd1 = 0.5_realType
       real(kind=realType), parameter :: rktBeta1 = 0.0750_realType
       real(kind=realType), parameter :: rktBetas = 0.09_realType
!
!      ******************************************************************
!      *                                                                *
!      * V2-f constants.                                                *
!      *                                                                *
!      ******************************************************************
!
       real(kind=realType), parameter :: rvfC1    = 1.4_realType
       real(kind=realType), parameter :: rvfC2    = 0.3_realType
       real(kind=realType), parameter :: rvfBeta  = 1.9_realType
       real(kind=realType), parameter :: rvfSigk1 = 1.0_realType
       real(kind=realType), parameter :: rvfSige1 = 0.7692307692_realType
       real(kind=realType), parameter :: rvfSigv1 = 1.00_realType
       real(kind=realType), parameter :: rvfCn    = 70.0_realType

       real(kind=realType), parameter :: rvfN1Cmu = 0.190_realType
       real(kind=realType), parameter :: rvfN1A   = 1.300_realType
       real(kind=realType), parameter :: rvfN1B   = 0.250_realType
       real(kind=realType), parameter :: rvfN1Cl  = 0.300_realType
       real(kind=realType), parameter :: rvfN6Cmu = 0.220_realType
       real(kind=realType), parameter :: rvfN6A   = 1.400_realType
       real(kind=realType), parameter :: rvfN6B   = 0.045_realType
       real(kind=realType), parameter :: rvfN6Cl  = 0.230_realType

       real(kind=realType) :: rvfLimitK, rvfLimitE, rvfCl
       real(kind=realType) :: rvfCmu
!
!      ******************************************************************
!      *                                                                *
!      * Variables to store the parameters for the wall functions fits. *
!      * As these variables depend on the turbulence model they are set *
!      * during runtime. Allocatables are used, because the number of   *
!      * fits could be different for the different models.              *
!      * The curve is divided in a number of intervals and is           *
!      * constructed such that both the function and the derivatives    *
!      * are continuous. Consequently cubic polynomials are used.       *
!      *                                                                *
!      ******************************************************************
!
       ! nFit:               Number of intervals of the curve.
       ! ypT(0:nFit):        y+ values at the interval boundaries.
       ! reT(0:nFit):        Reynolds number at the interval
       !                     boundaries, where the Reynolds number is
       !                     defined with the local velocity and the
       !                     wall distance.
       ! up0(nFit):          Coefficient 0 in the fit for the
       !                     nondimensional tangential velocity as a
       !                     function of the Reynolds number.
       ! up1(nFit):          Idem for coefficient 1.
       ! up2(nFit):          Idem for coefficient 2.
       ! up3(nFit):          Idem for coefficient 3.
       ! tup0(nFit,nt1:nt2): Coefficient 0 in the fit for the
       !                     nondimensional turbulence variables as a
       !                     function of y+.
       ! tup1(nFit,nt1:nt2): Idem for coefficient 1.
       ! tup2(nFit,nt1:nt2): Idem for coefficient 2.
       ! tup3(nFit,nt1:nt2): Idem for coefficient 3.
       ! tuLogFit(nt1:nt2):  Whether or not the logarithm of the variable
       !                     has been fitted.

       integer(kind=intType) :: nFit

       real(kind=realType), dimension(:), allocatable :: ypT, reT
       real(kind=realType), dimension(:), allocatable :: up0, up1
       real(kind=realType), dimension(:), allocatable :: up2, up3

       real(kind=realType), dimension(:,:), allocatable :: tup0, tup1
       real(kind=realType), dimension(:,:), allocatable :: tup2, tup3

       logical, dimension(:), allocatable :: tuLogFit

       end module paramTurb
