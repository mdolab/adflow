!
!     ******************************************************************
!     *                                                                *
!     * File:          getSolution.f90                                *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 23-07-2008                                      *
!     * Last modified: 23-07-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine getSolution
!
!     ******************************************************************
!     *                                                                *
!     * designExport compiles all the design data - functions and      *
!     *   their gradients, and design variable values - to export to   *
!     *   an optimizer.                                                *
!     *                                                                *
!     ******************************************************************
!
      use ADjointVars     ! functionValue, xDesignVar, nDesignDipoles
      use flowVarRefState ! magnetic
      use inputPhysics    ! velDirFreestream

      implicit none
!
!     Local variables.
!
      integer(kind=intType) :: level, sps, n, nn
      real(kind=realType)   :: CL, CD, CFx, CFy, CFz, CMx, CMy, CMz
      real(kind=realType)   :: alpha, beta
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      ! Set the relevant grid level and time instance.

      level = 1 ! finest grid level
      sps   = 1 ! first time instance
!
!     ******************************************************************
!     *                                                                *
!     * Function mapping.                                              *
!     *                                                                *
!     ******************************************************************
!
      ! Function values
      !print *,'calling computeAeroCoef'

      call computeAeroCoef(CL,CD, CFx, CFy, CFz,CMx,CMy,CMz,level,sps)

      functionValue(costFuncLiftCoef) = CL
      functionValue(costFuncDragCoef) = CD
      functionValue(costFuncForceXCoef) = CFx
      functionValue(costFuncForceYCoef) = CFy
      functionValue(costFuncForceZCoef) = CFz
      functionValue(costFuncMomXCoef) = CMx
      functionValue(costFuncMomYCoef) = CMy
      functionValue(costFuncMomZCoef) = CMz

    end subroutine getSolution
