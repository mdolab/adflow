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
      use inputTSStabDeriv !TSStability

      implicit none
!
!     Local variables.
!
      integer(kind=intType) :: level, sps, n, nn
      real(kind=realType)   :: CL, CD, CFx, CFy, CFz, CMx, CMy, CMz
      real(kind=realType)   :: alpha, beta

      real(kind=realType)::dcldp,dcldpdot,dcmzdp,dcmzdpdot         
      real(kind=realType)::dcldq,dcldqdot,dcmzdq,dcmzdqdot
      real(kind=realType)::dcldr,dcldrdot,dcmzdr,dcmzdrdot
      real(kind=realType)::dcldalpha,dcldalphadot,dcddalpha,dcmzdalpha,dcmzdalphadot
      real(kind=realType)::dcldMach,dcldMachdot,dcmzdMach,dcmzdMachdot
      real(kind=realType)::cl0,cl0dot,cD0,cmz0,cmz0dot
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
      
      if(TSStability)then
         
         call computeTSDerivatives(cl0,cd0,cmz0,dcldalpha,dcddalpha,&
           dcmzdalpha,dcmzdalphadot,dcmzdq)
         functionValue(costFuncCmzAlpha)     = dcmzdalpha
         functionValue( costFuncCm0)         = cmz0
         functionValue( costFuncClAlpha)     = dcldalpha
         functionValue( costFuncCl0  )       = cl0
         functionValue( costFuncCdAlpha )    = dcmzdalpha
         functionValue( costFuncCd0 )        = cd0
         functionValue( costFuncCmzAlphaDot) = dcmzdalphadot
         functionValue( costFuncCmzq)         = dcmzdq
      end if
    end subroutine getSolution
