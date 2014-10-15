subroutine getCostFunction(costFunction, force, moment, sepSensor, Cavitation, &
  alpha, beta, liftIndex, objValue)

  ! Compute the value of the actual objective function based on the
  ! (summed) forces and moments and any other "extra" design
  ! variables. The index of the objective is determined by 'iDV'. This
  ! function is intended to be AD'ed in reverse mode. 
  use inputTimeSpectral
  use costFunctions
  use inputPhysics
  use flowvarRefState
  implicit none

  ! Input 
  integer(kind=intType), intent(in) :: costFunction
  integer(kind=intType), intent(in) :: liftIndex
  real(kind=realType), intent(in), dimension(3, nTimeIntervalsSpectral) :: force, moment
  real(kind=realType), intent(in), dimension(nTimeIntervalsSpectral) :: sepSensor, Cavitation
  real(kind=realType), intent(in) :: alpha, beta

  ! Output
  real(kind=realType) :: objValue

  ! Working
  real(kind=realType) :: fact, factMoment, scaleDim, ovrNTS
  real(kind=realType), dimension(3) :: cf, cm
  real(kind=realType) :: elasticMomentx, elasticMomenty, elasticMomentz
  real(kind=realType), dimension(nTimeIntervalsSpectral, 8) :: baseCoef
  real(kind=realType), dimension(8) :: coef0, dcdalpha, dcdalphadot, dcdq, dcdqdot
  real(kind=realType) :: bendingMoment
  integer(kind=intType) :: sps

  ! Compute the factor since we may need it below
  call adjustInflowAngle(alpha, beta, liftIndex)
  call referenceState()
  scaleDim = pRef/pInf

  fact = two/(gammaInf*pInf*MachCoef**2 &
       *surfaceRef*LRef**2*scaleDim)
  factMoment = fact/(lengthRef*LRef)

  objValue = zero
  ovrNTS = one/nTimeIntervalsSpectral

  ! Pre-compute TS stability info if required:
  select case(costFunction)
  case(costFuncCl0,costFuncCd0,costFuncCm0, &
       costFuncClAlpha,costFuncCdAlpha,costFuncCmzAlpha,&
       costFuncClAlphaDot,costFuncCdAlphaDot,costFuncCmzAlphaDot,&
       costFuncClq,costFuncCdq,costFuncCmzq,&
       costFuncClqDot,costFuncCdqDot,costFuncCmzqDot)

     call computeTSDerivatives(force, moment, liftIndex, coef0, dcdalpha, &
          dcdalphadot, dcdq, dcdqdot)
  end select

  ! Main cost function selection
  select case(costFunction)

     ! ------------------ Steady or Average Objectives ---------------
  case(costFuncForceX)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*force(1, sps)
     end do
  case(costFuncForceY)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*force(2, sps)
     end do
  case(costFuncForceZ)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*force(3, sps)
     end do
  case(costFuncMomX)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*moment(1, sps)
     end do
  case(costFuncMomY)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*moment(2, sps)
     end do
  case(costFuncMomZ)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*moment(3, sps)
     end do

  case(costFuncForceXCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*force(1, sps)*fact
     end do
  case(costFuncForceYCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*force(2, sps)*fact
     end do
  case(costFuncForceZCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*force(3, sps)*fact
     end do
  case(costFuncMomXCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*moment(1, sps)*factMoment
     end do
  case(costFuncMomYCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*moment(2, sps)*factMoment
     end do
  case(costFuncMomZCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*moment(3, sps)*factMoment
     end do

  case(costFuncLift)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*(&
             force(1, sps)*liftDirection(1) + &
             force(2, sps)*liftDirection(2) + &
             force(3, sps)*liftDIrection(3))
     end do
  case(costFuncDrag)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*(&
             force(1, sps)*dragDirection(1) + &
             force(2, sps)*dragDirection(2) + &
             force(3, sps)*dragDIrection(3))
     end do

  case(costFuncLiftCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*fact*(&
             force(1, sps)*liftDirection(1) + &
             force(2, sps)*liftDirection(2) + &
             force(3, sps)*liftDIrection(3))
     end do
  case(costFuncDragCoef)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*fact*(&
             force(1, sps)*dragDirection(1) + &
             force(2, sps)*dragDirection(2) + &
             force(3, sps)*dragDIrection(3))
     end do
     ! -------------------- Time Spectral Objectives ------------------
  case(costFuncCl0)
     objValue = coef0(1)
  case(costFuncCd0)
     objValue = coef0(2)
  case(costFuncCm0)
     objValue = coef0(8)

  case(costFuncClAlpha)
     objValue = dcdalpha(1)
  case(costFuncCdAlpha)
     objValue = dcdalpha(2)
  case(costFuncCmzAlpha)
     objValue = dcdalpha(8)
  case(costFuncClAlphaDot)
     objValue = dcdalphadot(1)
  case(costFuncCdAlphaDot)
     objValue = dcdalphadot(2)
  case(costFuncCmzAlphaDot)
     objValue = dcdalphadot(8)

  case(costFuncClq)
     objValue = dcdq(1)
  case(costFuncCdq)
     objValue = dcdq(2)
  case(costFuncCmzq)
     objValue = dcdq(8)
  case(costFuncClqDot)
     objValue = dcdqdot(1)
  case(costFuncCdqDot)
     objValue = dcdqdot(2)
  case(costFuncCmzqDot)
     objValue = dcdqdot(8)

     ! ----------------- Bending moment calc -------------------
  case(costFuncBendingCoef)
     do sps=1,nTimeIntervalsSpectral
        cm = factMoment*moment(:, sps)
        cf = fact*force(:, sps)
        call computeRootBendingMoment(cf, cm, liftIndex, bendingMoment)
        objValue = objValue + ovrNTS*bendingMoment
     end do

  case (costFuncSepSensor)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*sepSensor(sps)
     end do
        
  case (costFuncCavitation)
     do sps=1,nTimeIntervalsSpectral
        objValue = objValue + ovrNTS*Cavitation(sps)
     end do

  end select
end subroutine getCostFunction

