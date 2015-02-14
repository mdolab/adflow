
subroutine getCostFuncMat(alpha, beta, liftIndex)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the cost function vector needed to calculate           *
  !     * dIdw, dIdx and dIda.                                            *
  !     *                                                                *
  !     ******************************************************************
  !
  use constants
  use costFunctions
  use flowVarRefState 
  use inputPhysics
  use inputTimeSpectral
  implicit none
 
  !Subroutine Vars Input/Ouput
  real(kind=realType), intent(in) :: alpha, beta
  integer(kind=intType), intent(in) ::liftIndex

  ! Working vars
  integer(kind=intType) :: sps
  real(kind=realtype) :: fact, scaleDim

  call adjustInflowAngle(alpha, beta, liftIndex)
  scaleDim = pRef/pInf

  fact = two/(gammaInf*pInf*MachCoef**2 &
       *surfaceRef*LRef**2*scaleDim)

  costFuncMat = zero

  do sps=1,nTimeIntervalsSpectral

     ! Just set the ones we know
     costFuncMat(:, costFuncLift,sps) = (/liftDirection(1), liftDirection(2), &
          liftDirection(3), zero, zero, zero/)
     
     costFuncMat(:, costFuncDrag,sps) = (/dragDirection(1), dragDirection(2), &
          dragDirection(3), zero, zero, zero/)
     
     costFuncMat(:, costFuncLiftCoef,sps) = (/liftDirection(1)*fact, &
          liftDirection(2)*fact, liftDirection(3)*fact, zero, zero, zero/)
     
     costFuncMat(:, costFuncDragCoef,sps) = (/dragDirection(1)*fact, &
          dragDirection(2)*fact, dragDirection(3)*fact, zero, zero, zero/)
     
     costFuncMat(:, costFuncForceX,sps) = (/one, zero, zero, zero, zero, zero/)
     costFuncMat(:, costFuncForceY,sps) = (/zero, one, zero, zero, zero, zero/)
     costFuncMat(:, costFuncForceZ,sps) = (/zero, zero, one, zero, zero, zero/)
     
     costFuncMat(:, costFuncForceXCoef,sps) = (/fact, zero, zero, zero, zero, zero/)
     costFuncMat(:, costFuncForceYCoef,sps) = (/zero, fact, zero, zero, zero, zero/)
     costFuncMat(:, costFuncForceZCoef,sps) = (/zero, zero, fact, zero, zero, zero/)
     
     costFuncMat(:, costFuncMomX,sps) = (/zero, zero, zero, one, zero, zero/)
     costFuncMat(:, costFuncMomY,sps) = (/zero, zero, zero, zero, one, zero/)
     costFuncMat(:, costFuncMomZ,sps) = (/zero, zero, zero, zero, zero, one/)
     
     ! update fact to get the moment
     fact = fact/(lengthRef*LRef)
     
     costFuncMat(:, costFuncMomXCoef,sps) = (/zero, zero, zero, fact, zero, zero/)
     costFuncMat(:, costFuncMomYCoef,sps) = (/zero, zero, zero, zero, fact, zero/)
     costFuncMat(:, costFuncMomZCoef,sps) = (/zero, zero, zero, zero, zero, fact/)

  end do

  costFuncMat = costFuncMat / nTimeIntervalsSpectral

end subroutine getCostFuncMat
