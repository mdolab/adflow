!
!     ******************************************************************
!     *                                                                *
!     * File:          computeAeroCoef.f90                             *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-17-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine computeAeroCoef(sps)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the aerodynamic coefficients from the force and moment *
  !     * produced by the pressure and shear stresses on the body walls: *
  !     *                                                                *
  !     ******************************************************************
  !

  use blockPointers  ! nDom
  use communication  ! my_ID SUmb_comm_world
  use inputPhysics   ! liftDirection, dragDirection
  use iteration      ! groundLevel
  use BCTypes
  !use monitor        ! monLoc, monGlob, nMonSum
  use costFunctions
  implicit none
  !
  !     Subroutine arguments.
  integer(kind=intType) :: sps
  !      Local variables.
  !
  integer :: ierr, nn,mm
  integer :: iBeg,iEnd,jBeg,jEnd,ii,npts

  real(kind=realType) :: force(3),cforce(3),Lift,Drag,CL,CD
  real(kind=realType) :: Moment(3),cMoment(3)
  real(kind=realType) :: alpha,beta
  integer(kind=intType) :: liftIndex
  real(kind=realType), dimension(:),allocatable :: lVars
  real(kind=realType), dimension(:),allocatable :: gVars 
  real(kind=realType),dimension(:,:),allocatable :: pts

  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  call getForceSize(npts)
  allocate(pts(3,npts))
  call getForcePoints(pts,npts)

  ii = 0

  call getDirAngle(velDirFreestream,LiftDirection,&
       liftIndex,alpha,beta)
 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

           call computeForceAndMomentAdj(Force,cForce,Lift,Drag,Cl,Cd,&
                moment,cMoment,alpha,beta,liftIndex,MachCoef,&
                pointRef,pts,npts,w,rightHanded,bcfaceid(mm),&
                iBeg,iEnd,jBeg,jEnd,ii)
           ii = ii + (iEnd-iBeg+1)*(jEnd-jBeg+1)

           functionValue(costFuncLift) = functionValue(costFuncLift) + Lift
           functionValue(costFuncDrag) = functionValue(costFuncDrag) + Drag
           functionValue(costFuncLiftCoef) = functionValue(costFuncLiftCoef) + Cl
           functionValue(costFuncDragCoef) = functionValue(costFuncDragCoef) + Cd
           functionValue(costFuncForceX) = functionValue(costFuncForceX) + Force(1)
           functionValue(costFuncForceY) = functionValue(costFuncForceY) + Force(2)
           functionValue(costFuncForceZ) = functionValue(costFuncForceZ) + Force(3)
           functionValue(costFuncForceXCoef) = functionValue(costFuncForceXCoef) + cForce(1)
           functionValue(costFuncForceYCoef) = functionValue(costFuncForceYCoef) + cForce(2)
           functionValue(costFuncForceZCoef) = functionValue(costFuncForceZCoef) + cForce(3)
           functionValue(costFuncMomX) = functionValue(costFuncMomX) + moment(1)
           functionValue(costFuncMomY) = functionValue(costFuncMomY) + moment(2)
           functionValue(costFuncMomZ) = functionValue(costFuncMomZ) + moment(3)
           functionValue(costFuncMomXCoef) = functionValue(costFuncMomXCoef) + cmoment(1)
           functionValue(costFuncMomYCoef) = functionValue(costFuncMomYCoef) + cmoment(2)
           functionValue(costFuncMomZCoef) = functionValue(costFuncMomZCoef) + cmoment(3)
        end if
     end do bocos
  end do domains

  deallocate(pts)
end subroutine computeAeroCoef
