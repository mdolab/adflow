subroutine getCostFunction2(force, moment, sepSensor, Cavitation, alpha, beta, liftIndex)

  ! Compute the value of the actual objective function based on the
  ! (summed) forces and moments and any other "extra" design
  ! variables. The index of the objective is determined by 'iDV'. This
  ! function is intended to be AD'ed in reverse mode. 
  use inputTimeSpectral
  use costFunctions
  use inputPhysics
  use flowvarRefState
  use inputTSStabDeriv
  implicit none

  ! Input 
  integer(kind=intType), intent(in) :: liftIndex
  real(kind=realType), intent(in), dimension(3, nTimeIntervalsSpectral) :: force, moment
  real(kind=realType), intent(in) ::  sepSensor, Cavitation
  real(kind=realType), intent(in) :: alpha, beta

  ! Working
  real(kind=realType) :: fact, factMoment, scaleDim, ovrNTS
  real(kind=realType), dimension(3) :: cf, cm
  real(kind=realType) :: elasticMomentx, elasticMomenty, elasticMomentz
  real(kind=realType), dimension(nTimeIntervalsSpectral, 8) :: baseCoef
  real(kind=realType), dimension(8) :: coef0, dcdalpha, dcdalphadot, dcdq, dcdqdot
  real(kind=realType) :: bendingMoment
  integer(kind=intType) :: sps

  ! Generate constants
  scaleDim = pRef/pInf
  fact = two/(gammaInf*pInf*MachCoef**2 &
       *surfaceRef*LRef**2*scaleDim)
  factMoment = fact/(lengthRef*LRef)

  ovrNTS = one/nTimeIntervalsSpectral

  ! Pre-compute TS stability info if required:
  if (TSStability) then
     call computeTSDerivatives(force, moment, liftIndex, coef0, dcdalpha, &
          dcdalphadot, dcdq, dcdqdot)
  end if
  
  funcValues = zero
  ! Now we just compute each cost function:
  do sps=1,nTimeIntervalsSpectral
     funcValues(costFuncForceX) = funcValues(costFuncForceX) + ovrNTS*force(1, sps)
     funcValues(costFuncForceY) = funcValues(costFuncForceY) + ovrNTS*force(2, sps)
     funcValues(costFuncForceZ) = funcValues(costFuncForceZ) + ovrNTS*force(3, sps)

     funcValues(costFuncMomX) = funcValues(costFuncMomX) + ovrNTS*moment(1, sps)
     funcValues(costFuncMomY) = funcValues(costFuncMomY) + ovrNTS*moment(2, sps)
     funcValues(costFuncMomZ) = funcValues(costFuncMomZ) + ovrNTS*moment(3, sps)

     funcValues(costFuncSepSensor) = funcValues(costFuncSepSensor) + ovrNTS*sepSensor
     funcValues(costFuncCavitation) = funcValues(costFuncCavitation) + ovrNTS*Cavitation

     ! Bending moment calc
     cm = factMoment*moment(:, sps)
     cf = fact*force(:, sps)
     call computeRootBendingMoment(cf, cm, liftIndex, bendingMoment)
     funcValues(costFuncBendingCoef) = funcValues(costFuncBendingCoef) + ovrNTS*bendingMoment
  end do

  funcValues(costFuncForceXCoef) = funcValues(costFuncForceX) * fact
  funcValues(costFuncForceYCoef) = funcValues(costFuncForceY) * fact
  funcValues(costFuncForceZCoef) = funcValues(costFuncForceZ) * fact
  
  funcValues(costFuncMomXCoef) = funcValues(costFuncMomX) * factMoment
  funcValues(costFuncMomYCoef) = funcValues(costFuncMomY) * factMoment
  funcValues(costFuncMomZCoef) = funcValues(costFuncMomZ) * factMoment

  funcValues(costFuncLift) = &
       funcValues(costFuncForceX)*liftDirection(1) + &
       funcValues(costFuncForceY)*liftDirection(2) + &
       funcValues(costFuncForceZ)*liftDirection(3)

  funcValues(costFuncDrag) = &
       funcValues(costFuncForceX)*dragDirection(1) + &
       funcValues(costFuncForceY)*dragDirection(2) + &
       funcValues(costFuncForceZ)*dragDirection(3)

  funcValues(costFuncLiftCoef) = funcValues(costFuncLift)*fact
  funcValues(costFuncDragCoef) = funcValues(costFuncDrag)*fact

     ! -------------------- Time Spectral Objectives ------------------
  funcValues(costFuncCl0) = coef0(1)
  funcValues(costFuncCd0) = coef0(2)
  funcValues(costFuncCm0) = coef0(8)
  funcValues(costFuncClAlpha) = dcdalpha(1)
  funcValues(costFuncCdAlpha) = dcdalpha(2)
  funcValues(costFuncCmzAlpha) = dcdalpha(8)
  funcValues(costFuncClAlphaDot) = dcdalphadot(1)
  funcValues(costFuncCdAlphaDot) = dcdalphadot(2)
  funcValues(costFuncCmzAlphaDot) = dcdalphadot(8)
  funcValues(costFuncClq) = dcdq(1)
  funcValues(costFuncCdq) = dcdq(2)
  funcValues(costFuncCmzq) = dcdq(8)
  funcValues(costFuncClqDot) = dcdqdot(1)
  funcValues(costFuncCdqDot) = dcdqdot(2)
  funcValues(costFuncCmzqDot) = dcdqdot(8)

end subroutine getCostFunction2

