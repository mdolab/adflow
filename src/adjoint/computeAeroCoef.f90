!
!     ******************************************************************
!     *                                                                *
!     * File:          computeAeroCoef.f90                             *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader,Gaetan Kenway  *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-13-2011                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine computeAeroCoef(globalCFVals,sps)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the aerodynamic coefficients from the force and moment *
  !     * produced by the pressure and shear stresses on the body walls: *
  !     *                                                                *
  !     ******************************************************************
  !

  use blockPointers 
  use communication 
  use inputPhysics   
  use iteration      
  use BCTypes
  use costFunctions
  use inputTimeSpectral
  use flowVarRefState
  use overset, only : oversetPresent
  use utils, only : EChk, setPointers
  implicit none

  ! Input/Ouput Variables
  integer(kind=intType), intent(in) :: sps
  real(kind=realType), intent(out), dimension(nCostFunction)::globalCFVals

  !      Local variables.
  integer(kind=intType) :: nn, ierr
  real(kind=realType) :: force(3), cforce(3), Lift, Drag, CL, CD
  real(kind=realType) :: Moment(3),cMoment(3), fact, scaleDim
  real(kind=realType), dimension(3) :: cFpLocal, cFvLocal, cMpLocal, cMvLocal, sepSensorAvgLocal
  real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yPlusMax, sepSensorLocal, cavitationLocal
  real(Kind=realType) :: sepSensor, sepSensorAvg(3), Cavitation
  real(kind=realType), dimension(nCostFunction)::localCFVals

  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  cFpLocal = zero
  cFvLocal = zero
  cMpLocal = zero
  cMvLocal = zero
  sepSensorLocal = zero
  sepSensorAvgLocal = zero
  cavitationLocal = zero
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     
     call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax, sepSensor, &
          sepSensorAvg, Cavitation)

     cFpLocal = cFpLocal + cFp
     cFvLocal = cFvLocal + cFv
     cMpLocal = cMpLocal + cMp
     cMvLocal = cMvLocal + cMv
     sepSensorLocal = sepSensorLocal + sepSensor
     sepSensorAvgLocal = sepSensorAvgLocal + sepSensorAvg
     cavitationLocal = cavitationLocal + cavitation
  end do domains

  if (oversetPresent) then 
     call forcesAndMomentsZipper(cfp, cfv, cmp, cmv, sps)
  
     ! Add the extra zipper component on the root proc
     if (myid == 0) then 
        cFpLocal = cFpLocal + cFp
        cFvLocal = cFvLocal + cFv
        cMpLocal = cMpLocal + cMp
        cMvLocal = cMvLocal + cMv
     end if
  end if

  ! Now compute the other variables
  cFp = cFpLocal
  cFv = cFvLocal
  cMp = cMpLocal
  cMv = cMvLocal
  sepSensor = sepSensorLocal
  cavitation = cavitationLocal
  sepSensorAvg = sepSensorAvgLocal
  

  scaleDim = pRef/pInf

  ! Sum pressure and viscous contributions
  cForce = cFp + cFv
  cMoment = cMp + cMv

  ! Get Lift coef and Drag coef
  CD =  cForce(1)*dragDirection(1) &
       + cForce(2)*dragDirection(2) &
       + cForce(3)*dragDirection(3)
  
  CL =  cForce(1)*liftDirection(1) &
       + cForce(2)*liftDirection(2) &
       + cForce(3)*liftDirection(3)
  
  ! Divide by fact to get the forces, Lift and Drag back
  fact = two/(gammaInf*pInf*MachCoef*MachCoef &
       *surfaceRef*LRef*LRef*scaleDim)
  Force = cForce / fact
  Lift  = CL / fact
  Drag  = CD / fact
  
  ! Moment factor has an extra lengthRef
  fact = fact/(lengthRef*LRef)
  Moment = cMoment / fact

  localCFVals(costFuncLift) = Lift
  localCFVals(costFuncDrag) = Drag
  localCFVals(costFuncLiftCoef) = Cl
  localCFVals(costFuncDragCoef) = Cd
  localCFVals(costFuncForceX) = Force(1)
  localCFVals(costFuncForceY) = Force(2)
  localCFVals(costFuncForceZ) = Force(3)
  localCFVals(costFuncForceXCoef) = cForce(1)
  localCFVals(costFuncForceYCoef) = cForce(2)
  localCFVals(costFuncForceZCoef) = cForce(3)
  localCFVals(costFuncMomX) = moment(1)
  localCFVals(costFuncMomY) = moment(2)
  localCFVals(costFuncMomZ) = moment(3)
  localCFVals(costFuncMomXCoef) = cmoment(1)
  localCFVals(costFuncMomYCoef) = cmoment(2)
  localCFVals(costFuncMomZCoef) = cmoment(3)
  localCFVals(costFuncSepSensor) = sepSensor
  localCFVals(costFuncSepSensorAvgX) = sepSensorAvg(1)
  localCFVals(costFuncSepSensorAvgY) = sepSensorAvg(2)
  localCFVals(costFuncSepSensorAvgZ) = sepSensorAvg(3)
  localCFVals(costFuncCavitation) = Cavitation

  ! Now we will mpi_allReduce them into globalCFVals
  call mpi_allreduce(localCFVals, globalCFVals, nCostFunction, sumb_real, &
       mpi_sum, SUmb_comm_world, ierr)

end subroutine computeAeroCoef
