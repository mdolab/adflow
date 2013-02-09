
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
  implicit none
 
  !Subroutine Vars Input/Ouput
  real(kind=realType), intent(in) :: alpha, beta
  integer(kind=intType), intent(in) ::liftIndex

  ! Working vars
  real(kind=realtype) :: fact

  call adjustInflowAngle(alpha, beta, liftIndex)

  fact = two/(gammaInf*pInf*MachCoef**2 &
       *surfaceRef*LRef**2)

  costFuncMat = zero

  ! Just set the ones we know
  costFuncMat(:, costFuncLift) = (/liftDirection(1), liftDirection(2), &
       liftDirection(3), zero, zero, zero/)
  
  costFuncMat(:, costFuncDrag) = (/dragDirection(1), dragDirection(2), &
       dragDirection(3), zero, zero, zero/)

  costFuncMat(:, costFuncLiftCoef) = (/liftDirection(1)*fact, &
       liftDirection(2)*fact, liftDirection(3)*fact, zero, zero, zero/)

  costFuncMat(:, costFuncDragCoef) = (/dragDirection(1)*fact, &
       dragDirection(2)*fact, dragDirection(3)*fact, zero, zero, zero/)

  costFuncMat(:, costFuncForceX) = (/one, zero, zero, zero, zero, zero/)
  costFuncMat(:, costFuncForceY) = (/zero, one, zero, zero, zero, zero/)
  costFuncMat(:, costFuncForceZ) = (/zero, zero, one, zero, zero, zero/)

  costFuncMat(:, costFuncForceXCoef) = (/fact, zero, zero, zero, zero, zero/)
  costFuncMat(:, costFuncForceYCoef) = (/zero, fact, zero, zero, zero, zero/)
  costFuncMat(:, costFuncForceZCoef) = (/zero, zero, fact, zero, zero, zero/)

  costFuncMat(:, costFuncMomX) = (/zero, zero, zero, one, zero, zero/)
  costFuncMat(:, costFuncMomY) = (/zero, zero, zero, zero, one, zero/)
  costFuncMat(:, costFuncMomZ) = (/zero, zero, zero, zero, zero, one/)

  costFuncMat(:, costFuncMomXCoef) = (/zero, zero, zero, fact, zero, zero/)
  costFuncMat(:, costFuncMomYCoef) = (/zero, zero, zero, zero, fact, zero/)
  costFuncMat(:, costFuncMomZCoef) = (/zero, zero, zero, zero, zero, fact/)

end subroutine getCostFuncMat
