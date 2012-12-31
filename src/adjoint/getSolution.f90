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

  use costFunctions
  use inputTSStabDeriv !TSStability
  use inputTimeSpectral !nTimeIntervalsSpectral
  use communication
  use inputIteration
  implicit none

  ! Input Variables
  integer(kind=intType) :: sps

  !  Local variables.
  real(kind=realType)   :: alpha, beta
  real(kind=realType), dimension(8) :: dcdp, dcdpdot, dcdq, dcdqdot, dcdr, dcdrdot
  real(kind=realType), dimension(8) :: dcdalpha, dcdalphadot, dcdbeta
  real(kind=realType), dimension(8) :: dcdbetadot, dcdMach, dcdMachdot
  real(kind=realType), dimension(8) :: Coef0,Coef0dot
  real(kind=realType), dimension(nCostFunction)::globalCFVals
  real(kind=realType), dimension(:),allocatable :: localVal,globalVal
  real(kind=realType)::bendingMoment,bendingSum
  real(kind=realType) :: value1,value2
  integer(kind=intType) :: ierr, i

  ! Function values
  if (.not. allocated(functionValue)) then
     allocate(functionValue(nCostFunction))
  end if

  functionValue(:) = 0.0
  bendingSum = 0.0
  do i =1,nTimeIntervalsSpectral
     call computeAeroCoef(globalCFVals,i)
     
     call computeRootBendingMoment(globalCFVals,bendingMoment)
     bendingsum = bendingsum+bendingMoment
     if(myid==0 .and. printIterations)then
        print *,'Bending Coefficient',bendingMoment,i
     end if
  end do
  if(myid==0 .and. printIterations)then
     print *,'bending average',bendingSum/nTimeIntervalsSpectral
  end if
  functionValue(costFuncBendingCoef)=bendingSum/nTimeIntervalsSpectral

  call computeAeroCoef(globalCFVals,sps)

  functionValue(costFuncLift) = globalCFVals(costFuncLift) 
  functionValue(costFuncDrag) = globalCFVals(costFuncDrag) 
  functionValue(costFuncLiftCoef) = globalCFVals(costFuncLiftCoef) 
  functionValue(costFuncDragCoef) = globalCFVals(costFuncDragCoef) 
  functionValue(costFuncForceX) = globalCFVals(costFuncForceX) 
  functionValue(costFuncForceY) = globalCFVals(costFuncForceY) 
  functionValue(costFuncForceZ) = globalCFVals(costFuncForceZ) 
  functionValue(costFuncForceXCoef) = globalCFVals(costFuncForceXCoef) 
  functionValue(costFuncForceYCoef) = globalCFVals(costFuncForceYCoef) 
  functionValue(costFuncForceZCoef) = globalCFVals(costFuncForceZCoef) 
  functionValue(costFuncMomX) = globalCFVals(costFuncMomX) 
  functionValue(costFuncMomY) = globalCFVals(costFuncMomY) 
  functionValue(costFuncMomZ) = globalCFVals(costFuncMomZ) 
  functionValue(costFuncMomXCoef) = globalCFVals(costFuncMomXCoef) 
  functionValue(costFuncMomYCoef) = globalCFVals(costFuncMomYCoef)
  functionValue(costFuncMomZCoef) = globalCFVals(costFuncMomZCoef)

  if(TSStability)then
     call computeTSDerivatives(coef0,dcdalpha,dcdalphadot,dcdq,dcdqdot)

     functionValue( costFuncCl0  )       = coef0(1)
     functionValue( costFuncCd0 )        = coef0(2)
     functionValue( costFuncCFy0 )       = coef0(4)
     functionValue( costFuncCm0 )        = coef0(8)

     functionValue( costFuncClAlpha)     = dcdalpha(1)
     functionValue( costFuncCdAlpha)     = dcdalpha(2)
     functionValue( costFuncCFyAlpha)    = dcdalpha(4)
     functionValue( costFuncCmzAlpha)    = dcdalpha(8)

     functionValue( costFuncClAlphaDot)     = dcdalphadot(1)
     functionValue( costFuncCdAlphaDot)     = dcdalphadot(2)
     functionValue( costFuncCFyAlphaDot)    = dcdalphadot(4)
     functionValue( costFuncCmzAlphaDot)    = dcdalphadot(8)
    
     functionValue( costFuncClq)         = dcdq(1)
     functionValue( costFuncCdq)         = dcdq(2)
     functionValue( costFuncCfyq)        = dcdq(4)
     functionValue( costFuncCmzq)        = dcdq(8)

     functionValue( costFuncClqDot)         = dcdqdot(1)
     functionValue( costFuncCdqDot)         = dcdqdot(2)
     functionValue( costFuncCfyqDot)        = dcdqdot(4)
     functionValue( costFuncCmzqDot)        = dcdqdot(8)
  end if
end subroutine getSolution
