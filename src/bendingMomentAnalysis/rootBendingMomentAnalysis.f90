!
!****************************************************
!
!   Filename: rootBendingMomentAnalysis.f90
!   Author: C.A.(Sandy) Mader
!   Date Started: July 10, 2011
!   Date Modified: July 10,2011
!
!****************************************************

subroutine computeRootBendingMoment(sol,bendingMoment)

  !*******************************************************
  !                                                      *
  ! Compute a normalized bending moment coefficient from *
  ! the force and moment coefficient. At the moment this *
  ! Routine only works for a half body. Additional logic *
  ! would be needed for a full body.                     *
  !                                                      *
  !*******************************************************

  use inputPhysics   ! liftDirection, dragDirection,pointref,pointrefec
  use costFunctions
  use communication
  use inputIteration
  implicit none

  !input/output variables
  real(kind=realType), dimension(nCostFunction):: Sol
  real(kind=realType)::bendingMoment
  !Subroutine Variables
  real(kind=realType):: cmx,cmy,cmz,cfx,cfy,cfz 
  real(kind=realType)::elasticMomentx,elasticMomenty,elasticMomentz
  integer(kind=intType)::liftIndex
  real(kind=realType)::alpha,beta

  !Begin execution
  !determine the liftIndex from the flow and liftdirection
  call getDirAngle(velDirFreestream,LiftDirection,&
       liftIndex,alpha,beta)

  !Set some basic variables from the solution, do the analysis based on
  !the coefficients so that it is well scaled
  cmx = sol(costFuncMomXCoef)
  cmy = sol(costFuncMomYCoef)
  cmz = sol(costFuncMomZCoef)
  cfx = sol(costFuncForceXCoef)
  cfy = sol(costFuncForceYCoef)
  cfz = sol(costFuncForceZCoef)
  !sum up the moments about the root elastic center to determine the effect of sweep on the moment
  
  if (liftIndex == 2) then
     !z out wing sum momentx,momentz
     elasticMomentx = cmx + cfy*(pointRefEC(3)-pointRef(3))/lengthref-cfz*(pointRefEC(2)-pointRef(2))/lengthref
     elasticMomentz = cmz - cfy*(pointRefEC(1)-pointref(1))/lengthref+cfx*(pointRefEC(2)-pointRef(2))/lengthref
    
     BendingMoment = sqrt(elasticMomentx**2+elasticMomentz**2)
     elasticMomenty = zero
  elseif (liftIndex == 3) then
     !y out wing sum momentx,momenty
     elasticMomentx = cmx + cfz*(pointrefEC(2)-pointRef(2))/lengthref+cfy*(pointrefEC(3)-pointref(3))/lengthref
     elasticMomenty = cmy + cfz*(pointRefEC(1)-pointRef(1))/lengthref+cfx*(pointrefEC(3)-pointRef(3))/lengthref
     
     BendingMoment = sqrt(elasticMomentx**2+elasticMomenty**2)
     elasticMomentz = zero
  end if
  if( myid==0 .and. printIterations) then
     print *,'Bending moment components',elasticMomentx,elasticMomenty,elasticMomentz
  end if

end subroutine computeRootBendingMoment

