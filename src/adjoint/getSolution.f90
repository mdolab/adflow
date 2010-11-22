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
      subroutine getSolution(sps)
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
      use communication !myID
      implicit none

!
!     Subroutine Variables
!
      integer(kind=intType):: sps
!
!     Local variables.
!
      integer(kind=intType) :: level, n, nn
      real(kind=realType)   :: CL, CD, CFx, CFy, CFz, CMx, CMy, CMz
      real(kind=realType)   :: alpha, beta

      real(kind=realType)::dcldp,dcldpdot,dcmzdp,dcmzdpdot         
      real(kind=realType)::dcldq,dcldqdot,dcddq,dcddqdot,dcmzdq,dcmzdqdot
      real(kind=realType)::dcldr,dcldrdot,dcmzdr,dcmzdrdot
      real(kind=realType)::dcldalpha,dcldalphadot,dcddalpha,dcddalphadot,dcmzdalpha,dcmzdalphadot
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
      !if (myid==0) print *,'getting solution',sps
      level = 1 ! finest grid level
      !sps   = 1 ! first time instance
!
!     ******************************************************************
!     *                                                                *
!     * Function mapping.                                              *
!     *                                                                *
!     ******************************************************************
!
      ! Function values
      !print *,'calling computeAeroCoef'
      functionValue(:) = 0.0
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
           dcmzdalpha,dcldalphadot,dcddalphadot,dcmzdalphadot,dcldq,&
           dcddq,dcmzdq,dcldqdot,dcddqdot,dcmzdqdot)

         functionValue( costFuncCmzAlpha)    = dcmzdalpha
         functionValue( costFuncCm0)         = cmz0
         functionValue( costFuncCmzAlphadot) = dcmzdalphadot
         functionValue( costFuncClAlpha)     = dcldalpha
         functionValue( costFuncCl0  )       = cl0
         functionValue( costFuncClAlphadot)  = dcldalphadot
         functionValue( costFuncCdAlpha )    = dcddalpha
         functionValue( costFuncCd0 )        = cd0
         functionValue( costFuncCdAlphaDot)  = dcddalphadot
         functionValue( costFuncCmzq)        = dcmzdq
         functionValue( costFuncCmzqdot)     = dcmzdqdot
         functionValue( costFuncClq)         = dcldq
         functionValue( costFuncClqdot)      = dcldqdot

      end if
    end subroutine getSolution
