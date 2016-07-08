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
  implicit none

  ! Input/Ouput Variables
  integer(kind=intType), intent(in) :: sps
  real(kind=realType), intent(out), dimension(nCostFunction)::globalCFVals

  !      Local variables.
  integer(kind=intType) :: nn, ierr
  real(kind=realType) :: force(3), cforce(3), Lift, Drag, CL, CD
  real(kind=realType) :: Moment(3),cMoment(3), fact, scaleDim
  real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yPlusMax
  real(kind=realType) :: cMpaxis, CMvaxis, cAxisMoment
  real(Kind=realType) :: sepSensor, sepSensorAvg(3), Cavitation, axisMoment
  real(kind=realType), dimension(nCostFunction)::localCFVals

  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  !Zero the summing variable
  localCFVals(:) = 0.0
  globalCFVals(:) = 0.0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     
     call forcesAndMoments(cFp, cFv, cMp, cMv, cMpaxis, cMvaxis, & 
          yplusMax, sepSensor, sepSensorAvg, Cavitation)
     scaleDim = pRef/pInf

     ! Sum pressure and viscous contributions
     cForce = cFp + cFv
     cMoment = cMp + cMv
     cAxisMoment = cMpaxis + cMvaxis

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
     axisMoment = cAxisMoment / fact

     localCFVals(costFuncLift) = localCFVals(costFuncLift) + Lift
     localCFVals(costFuncDrag) = localCFVals(costFuncDrag) + Drag
     localCFVals(costFuncLiftCoef) = localCFVals(costFuncLiftCoef) + Cl
     localCFVals(costFuncDragCoef) = localCFVals(costFuncDragCoef) + Cd
     localCFVals(costFuncForceX) = localCFVals(costFuncForceX) + Force(1)
     localCFVals(costFuncForceY) = localCFVals(costFuncForceY) + Force(2)
     localCFVals(costFuncForceZ) = localCFVals(costFuncForceZ) + Force(3)
     localCFVals(costFuncForceXCoef) = localCFVals(costFuncForceXCoef) + cForce(1)
     localCFVals(costFuncForceYCoef) = localCFVals(costFuncForceYCoef) + cForce(2)
     localCFVals(costFuncForceZCoef) = localCFVals(costFuncForceZCoef) + cForce(3)
     localCFVals(costFuncMomX) = localCFVals(costFuncMomX) + moment(1)
     localCFVals(costFuncMomY) = localCFVals(costFuncMomY) + moment(2)
     localCFVals(costFuncMomZ) = localCFVals(costFuncMomZ) + moment(3)
     localCFVals(costFuncMomXCoef) = localCFVals(costFuncMomXCoef) + cmoment(1)
     localCFVals(costFuncMomYCoef) = localCFVals(costFuncMomYCoef) + cmoment(2)
     localCFVals(costFuncMomZCoef) = localCFVals(costFuncMomZCoef) + cmoment(3)
     localCFVals(costFuncSepSensor) = localCFVals(costFuncSepSensor) + sepSensor
     localCFVals(costFuncSepSensorAvgX) = localCFVals(costFuncSepSensorAvgX) + sepSensorAvg(1)
     localCFVals(costFuncSepSensorAvgY) = localCFVals(costFuncSepSensorAvgY) + sepSensorAvg(2)
     localCFVals(costFuncSepSensorAvgZ) = localCFVals(costFuncSepSensorAvgZ) + sepSensorAvg(3)
     localCFVals(costFuncCavitation) = localCFVals(costFuncCavitation) + Cavitation
     localCFVals(costFuncAxisMoment) = localCFVals(costFuncAxisMoment) + AxisMoment
     
  end do domains

  ! Now we will mpi_allReduce them into globalCFVals
  call mpi_allreduce(localCFVals, globalCFVals, nCostFunction, sumb_real, &
       mpi_sum, SUmb_comm_world, ierr)

end subroutine computeAeroCoef
