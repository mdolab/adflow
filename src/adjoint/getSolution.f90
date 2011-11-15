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
  implicit none

  integer(kind=intType) :: sps,ierr,i
  !
  !     Local variables.
  !
  real(kind=realType)   :: alpha, beta
  real(kind=realType),dimension(8)::dcdp,dcdpdot,dcdq,dcdqdot,dcdr,dcdrdot
  real(kind=realType),dimension(8)::dcdalpha,dcdalphadot,dcdbeta,dcdbetadot,&
       dcdMach,dcdMachdot
  real(kind=realType),dimension(8)::Coef0,Coef0dot

  real(kind=realType), dimension(nCostFunction)::globalCFVals
  real(kind=realType),dimension(:),allocatable :: localVal,globalVal
  real(kind=realType)::bendingMoment,bendingSum

  real(kind=realType) :: value1,value2
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
     !if(myid==0)then
     !   print *,'Bending Coefficient',bendingMoment,i
     !end if
  end do
  !if(myid==0)then
  !   print *,'bending average',bendingSum/nTimeIntervalsSpectral
  !end if
  functionValue(costFuncBendingCoef)=bendingSum/nTimeIntervalsSpectral

  
  call computeAeroCoef(globalCFVals,sps)
  
 !  call farFieldDrag(value1)
!   call farFieldInducedDrag(value2)
!   if (myID==0)then
!      print *,'Total Drag:',value2-value1
!   end if
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
     !(cl0,cd0,cmz0,dcldalpha,dcddalpha,&
     !     dcmzdalpha,dcldalphadot,dcddalphadot,dcmzdalphadot,dcldq,&
     !     dcddq,dcmzdq,dcldqdot,dcddqdot,dcmzdqdot)

     functionValue( costFuncCl0  )       = coef0(1)!cl0
     functionValue( costFuncCd0 )        = coef0(2)!cd0
     functionValue( costFuncCFy0 )       = coef0(4)
     functionValue( costFuncCm0 )        = coef0(8)!cmz0

     functionValue( costFuncClAlpha)     = dcdalpha(1)!dcldalpha
     functionValue( costFuncCdAlpha)     = dcdalpha(2)!dcddalpha
     functionValue( costFuncCFyAlpha)    = dcdalpha(4)!dcddalpha
     functionValue( costFuncCmzAlpha)    = dcdalpha(8)!dcmzdalpha

     functionValue( costFuncClAlphaDot)     = dcdalphadot(1)!dcldalphadot
     functionValue( costFuncCdAlphaDot)     = dcdalphadot(2)!dcddalphadot
     functionValue( costFuncCFyAlphaDot)    = dcdalphadot(4)!dcddalphadot
     functionValue( costFuncCmzAlphaDot)    = dcdalphadot(8)! dcmzdalphadot
     
     functionValue( costFuncClq)         = dcdq(1)!dcldq
     functionValue( costFuncCdq)         = dcdq(2)!dcddq
     functionValue( costFuncCfyq)        = dcdq(4)!dcfydq
     functionValue( costFuncCmzq)        = dcdq(8)!dcmzdq

     functionValue( costFuncClqDot)         = dcdqdot(1)!dcldqDot
     functionValue( costFuncCdqDot)         = dcdqdot(2)!dcddqDot
     functionValue( costFuncCfyqDot)        = dcdqdot(4)!dcfydqDot
     functionValue( costFuncCmzqDot)        = dcdqdot(8)!dcmzdqDot
     !print *,'function value',functionValue,'dcdq',dcdq
  end if

 
end subroutine getSolution
