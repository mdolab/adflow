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

  use blockPointers  ! nDom
  use communication  ! my_ID SUmb_comm_world
  use inputPhysics   ! liftDirection, dragDirection
  use iteration      ! groundLevel
  use BCTypes
  !use monitor        ! monLoc, monGlob, nMonSum
  use costFunctions
  use inputTimeSpectral
  implicit none
  !
  !     Subroutine arguments.
  integer(kind=intType) :: sps
  real(kind=realType), dimension(nCostFunction)::globalCFVals
  !      Local variables.
  !
  integer :: ierr, nn,mm
  integer :: iBeg,iEnd,jBeg,jEnd,ii,npts,nTS

  real(kind=realType) :: force(3),cforce(3),Lift,Drag,CL,CD
  real(kind=realType) :: Moment(3),cMoment(3)
  real(kind=realType) :: alpha,beta
  integer(kind=intType) :: liftIndex
  real(kind=realType), dimension(nCostFunction)::localCFVals
  real(kind=realType),dimension(:,:,:),allocatable :: pts
  logical :: forcesTypeSave
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  call getForceSize(npts,nTS)
  allocate(pts(3,npts,nTimeIntervalsSpectral))
  call getForcePoints(pts,npts,nTS)
  ii = 0

  call getDirAngle(velDirFreestream,LiftDirection,&
       liftIndex,alpha,beta)
  
  forcesTypeSave = forcesAsTractions
  forcesAsTractions = .False.
  !Zero the summing variable
  localCFVals(:) = 0.0
  globalCFVals(:) = 0.0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

           call computeForceAndMomentAdj(Force,cForce,Lift,Drag,Cl,Cd,&
                moment,cMoment,alpha,beta,liftIndex,MachCoef,&
                pointRef,lengthRef,surfaceRef,pts(:,:,sps),npts,w,&
                rightHanded,bcfaceid(mm),iBeg,iEnd,jBeg,jEnd,ii,sps)
           ii = ii + (iEnd-iBeg+1)*(jEnd-jBeg+1)

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

        end if
     end do bocos
  end do domains

  ! Now we will mpi_allReduce them into globalCFVals
 
  call mpi_allreduce(localCFVals, globalCFVals, nCostFunction, sumb_real, &
       mpi_sum, SUmb_comm_world, ierr)

  forcesAsTractions = forcesTypeSave

  deallocate(pts)
end subroutine computeAeroCoef
