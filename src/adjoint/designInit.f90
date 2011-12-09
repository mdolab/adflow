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
!     ******************************************************************
!     *                                                                *
!     * designInit determines the number of design variables,          *
!     *   allocates the arrays that store the cost function names,     *
!     *   values and gradients and fills in the function names. It     *
!     *   also allocates the arrays that store the design variables,   *
!     *   their names, lower and upper bounds, filling their names.    *
!     *                                                                *
!     * The variables are split in the spatial variables (x) and       *
!     *   and extra variables, such as angle of attack.                *
!     *                                                                *
!     * The extra design variables are ordered as:                     *
!     *  - angle of attack                                             *
!     *  - side-slip angle                                             *
!     *  - ...                                                         *
!     *                                                                *
!     * The routines that compute dR/da (setupGradientMatrixExtra)     *
!     * and dJ/da (setupGradientRHSExtra) have to be consistent with   *
!     * this ordering...                                               *
!     *                                                                *
!     ******************************************************************
!
      use ADjointVars     ! nDesign, nDesignExtra, nDesignSpatial
                          ! nNodesGlobal, nDesignDipoles
      use communication   ! myID
      use flowVarRefState ! magnetic
      use inputTimeSpectral !nTimeIntervalsSpectral
      implicit none
!
!     Local variables.
!
      character(len=maxStringLen) :: varName, tmp
      integer(kind=intType)       :: n, nn,ierr
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

      nDesignExtra = 13!7
  
      ! Determine the total number of design variables.
      !   + nDesignSpatial (= 3 * nNodesGlobal)
      !   + nDesignExtra

      nDesignSpatial = (3 * nNodesGlobal)*nTimeIntervalsSpectral
      nDesign        = nDesignExtra + nDesignSpatial


      end subroutine designInit
