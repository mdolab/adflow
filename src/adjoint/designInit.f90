!
!     ******************************************************************
!     *                                                                *
!     * File:          designInit.f90                                  *
!     * Authors:       C.A.(Sandy) Mader                               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine designInit
!
!
      use ADjointVars     
      use flowVarRefState 
      use inputTimeSpectral 
      implicit none
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      ! Determine the number of extra variables. This includes:
      ! - angle of attack
      ! - side-slip angle
      ! - Mach Number,mach NumberGrid
      ! - Rotational Rate (3)
      ! - Rotational Center (3)
      ! - Moment Reference Point (3)
      ! - will eventually include Flight Dynamic Derivatives.

      nDesignExtra = 13

      end subroutine designInit
