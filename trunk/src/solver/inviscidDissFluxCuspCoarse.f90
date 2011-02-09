!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidDissFluxCuspCoarse.f90                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-25-2003                                      *
!      * Last modified: 08-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidDissFluxCuspCoarse
!
!      ******************************************************************
!      *                                                                *
!      * inviscidDissFluxCuspCoarse computes the artificial dissipation *
!      * for the coarse grid, i.e. first order, CUSP scheme for a given *
!      * block. Therefore it is assumed that the pointers in            *
!      * blockPointers already point to the correct block.              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use inputDiscretization
       use iteration
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call terminate("inviscidDissFluxCuspCoarse", &
                      "not implemented yet")

       end subroutine inviscidDissFluxCuspCoarse
