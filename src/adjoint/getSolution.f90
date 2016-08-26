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
  use inputTSStabDeriv
  use inputPhysics
  use inputTimeSpectral 
  use communication
  use inputIteration
  use utils, only : computeTSDerivatives, computeRootBendingMoment
  implicit none

  ! Input Variables
  integer(kind=intType) :: sps

  !  Local variables.
  real(kind=realType)   :: alpha, beta
  real(kind=realType), dimension(8) ::  dcdq, dcdqdot
  real(kind=realType), dimension(8) :: dcdalpha, dcdalphadot
  real(kind=realType), dimension(8) :: Coef0
  real(kind=realType), dimension(nCostFunction)::globalCFVals
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment
  real(kind=realType)::bendingMoment,bendingSum, cf(3), cm(3)
  integer(kind=intType) :: i, liftIndex

  !Begin execution
  !determine the liftIndex from the flow and liftdirection
  call getDirAngle(velDirFreestream, LiftDirection,&
       liftIndex, alpha, beta)

  funcValues(:) = 0.0

  bendingSum = 0.0
  do i =1,nTimeIntervalsSpectral
     call computeAeroCoef(globalCFVals,i)

     force(1, i) = globalCFVals(costFuncForceX)
     force(2, i) = globalCFVals(costFuncForceY)
     force(3, i) = globalCFVals(costFuncForceZ)
     moment(1, i) = globalCFVals(costFuncMomX)
     moment(2, i) = globalCFVals(costFuncMomY)
     moment(3, i) = globalCFVals(costFuncMomZ)

     cf = (/globalCFVals(costFuncForceXCoef), &
          globalCFVals(costFuncForceYCoef), &
          globalCFVals(costFuncForceZCoef)/)
     
     cm = (/globalCFVals(costFuncMomXCoef), &
          globalCFVals(costFuncMomYCoef), &
          globalCFVals(costFuncMomZCoef)/)
     call computeRootBendingMoment(cf, cm, liftIndex, bendingMoment)
     bendingsum = bendingsum+bendingMoment
  end do
  
  call computeAeroCoef(globalCFVals,sps)
  funcValues(costFuncBendingCoef)=bendingSum/nTimeIntervalsSpectral

  funcValues(costFuncLift) = globalCFVals(costFuncLift) 
  funcValues(costFuncDrag) = globalCFVals(costFuncDrag) 
  funcValues(costFuncLiftCoef) = globalCFVals(costFuncLiftCoef) 
  funcValues(costFuncDragCoef) = globalCFVals(costFuncDragCoef) 
  funcValues(costFuncForceX) = globalCFVals(costFuncForceX) 
  funcValues(costFuncForceY) = globalCFVals(costFuncForceY) 
  funcValues(costFuncForceZ) = globalCFVals(costFuncForceZ) 
  funcValues(costFuncForceXCoef) = globalCFVals(costFuncForceXCoef) 
  funcValues(costFuncForceYCoef) = globalCFVals(costFuncForceYCoef) 
  funcValues(costFuncForceZCoef) = globalCFVals(costFuncForceZCoef) 
  funcValues(costFuncMomX) = globalCFVals(costFuncMomX) 
  funcValues(costFuncMomY) = globalCFVals(costFuncMomY) 
  funcValues(costFuncMomZ) = globalCFVals(costFuncMomZ) 
  funcValues(costFuncMomXCoef) = globalCFVals(costFuncMomXCoef) 
  funcValues(costFuncMomYCoef) = globalCFVals(costFuncMomYCoef)
  funcValues(costFuncMomZCoef) = globalCFVals(costFuncMomZCoef)
  funcValues(costFuncSepSensor) = globalCFVals(costFuncSepSensor)
  funcValues(costFuncSepSensorAvgX) = globalCFVals(costFuncSepSensorAvgX)
  funcValues(costFuncSepSensorAvgY) = globalCFVals(costFuncSepSensorAvgY)
  funcValues(costFuncSepSensorAvgZ) = globalCFVals(costFuncSepSensorAvgZ)

  funcValues(costFuncCavitation) = globalCFVals(costFuncCavitation)

  if(TSStability)then

     call computeTSDerivatives(force, moment, liftIndex, coef0, dcdalpha, &
          dcdalphadot, dcdq, dcdqdot)

     funcValues( costFuncCl0  )       = coef0(1)
     funcValues( costFuncCd0 )        = coef0(2)
     funcValues( costFuncCFy0 )       = coef0(4)
     funcValues( costFuncCm0 )        = coef0(8)

     funcValues( costFuncClAlpha)     = dcdalpha(1)
     funcValues( costFuncCdAlpha)     = dcdalpha(2)
     funcValues( costFuncCFyAlpha)    = dcdalpha(4)
     funcValues( costFuncCmzAlpha)    = dcdalpha(8)

     funcValues( costFuncClAlphaDot)     = dcdalphadot(1)
     funcValues( costFuncCdAlphaDot)     = dcdalphadot(2)
     funcValues( costFuncCFyAlphaDot)    = dcdalphadot(4)
     funcValues( costFuncCmzAlphaDot)    = dcdalphadot(8)
    
     funcValues( costFuncClq)         = dcdq(1)
     funcValues( costFuncCdq)         = dcdq(2)
     funcValues( costFuncCfyq)        = dcdq(4)
     funcValues( costFuncCmzq)        = dcdq(8)

     funcValues( costFuncClqDot)         = dcdqdot(1)
     funcValues( costFuncCdqDot)         = dcdqdot(2)
     funcValues( costFuncCfyqDot)        = dcdqdot(4)
     funcValues( costFuncCmzqDot)        = dcdqdot(8)
  end if
end subroutine getSolution
