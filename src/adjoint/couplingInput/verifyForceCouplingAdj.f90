!
!     ******************************************************************
!     *                                                                *
!     * File:          verifyForcesAdj.f90                             *
!     * Author:        C.A.(Sandy) Mader, Andre C. Marta               *
!     *                Seongim Choi                                    *
!     * Starting date: 12-14-2007                                      *
!     * Last modified: 12-27-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifyForceCouplingAdj(level)
!
!     ******************************************************************
!     *                                                                *
!     * Computes the Force coefficients for the current configuration  *
!     * for the finest grid level over all time instances using the    *
!     * auxiliary routines modified for tapenade and compares them to  *
!     * the force coefficients computed with the original code.        *
!     *                                                                *
!     * This is only executed in debug mode.                           *
!     *                                                                *
!     ******************************************************************
!
      use blockPointers
      use communication       ! myID
      use inputPhysics        !
      use flowVarRefState     !
      use inputDiscretization ! spaceDiscr
      use iteration           ! currentLevel
      use monitor             ! monLoc, MonGlob, nMonSum
      use inputTimeSpectral   ! nTimeInstancesMax
      use section
      use bcTypes             ! EulerWall, NSWallAdiabatic, NSWallIsothermal
      use mddatalocal
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level
!
!     Local variables.
!
      integer(kind=intType) :: discr, nHalo
      integer(kind=intType) :: mm, nn, sps=1
      integer(kind=intType)::  famID=0
       integer(kind=intType) :: startInd, endInd

      logical :: fineGrid, correctForK, exchangeTurb

      integer(kind=intType)::liftIndex

      real(kind=realType), dimension(:,:,:,:), allocatable :: xAdj,wAdj
      real(kind=realType), dimension(:,:,:), allocatable :: pAdj
!      real(kind=realType), dimension(:,:,:,:), allocatable :: siAdj, sjAdj, skAdj

      integer(kind=intType) :: i2Beg, i2End, j2Beg, j2End
      integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd

      real(kind=realType), dimension(3) :: cFpAdj, cFvAdj
!!$      real(kind=realType), dimension(3) :: cMpAdj, cMvAdj
!!$
!!$      real(kind=realType), dimension(3) :: cFp, cFv
!!$      real(kind=realType), dimension(3) :: cMp, cMv
!!$
!!$      real(kind=realType), dimension(3) :: cFpFD, cFvFD
!!$      real(kind=realType), dimension(3) :: cMpFD, cMvFD
!!$
!!$      real(kind=realType) :: Cl,Cd,Cfx,Cfy,Cfz,Cmx,Cmy,Cmz
!!$      real(kind=realType), dimension(nTimeIntervalsSpectral) :: ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj, &
!!$                                                           CmxAdj,CmyAdj,CmzAdj
!!$!      real(kind=realType), dimension(:,:,:),allocatable:: normAdj
      real(kind=realType), dimension(3) :: refPoint
      real(kind=realType) :: yplusMax

      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,machCoefAdj
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj


      integer(kind=intType) :: nSurfNodesLoc, modFamID,ii
      real(kind=realType), dimension(:,:),allocatable :: forceLoc
      integer(kind=intType) :: nmonsum2 ,i
 
      real(kind=realType) :: fact!temporary

      logical :: contributeToForce, viscousSubface,secondhalo, righthanded

      integer :: ierr

!for debug
real(kind=realType), dimension(3) :: cfpadjout, cmpadjout
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      if(myID == 0) then
        write(*,*) "Running verifyForceCouplingAdj..."
      endif

      ! Set the grid level of the current MG cycle, the value of the
      ! discretization and the logical fineGrid.

      currentLevel = level
      discr        = spaceDiscr
      fineGrid     = .true.

      ! Determine whether or not the total energy must be corrected
      ! for the presence of the turbulent kinetic energy and whether
      ! or not turbulence variables should be exchanged.

      correctForK  = .false.
      exchangeTurb = .false.
      secondhalo = .true.
!
!     ******************************************************************
!     *                                                                *
!     * Exchange halo data to make sure it is up-to-date.              *
!     * (originally called inside "rungeKuttaSmoother" subroutine).    *
!     *                                                                *
!     ******************************************************************

!      print *,'exchanging halo data'
      ! Exchange the pressure if the pressure must be exchanged early.
      ! Only the first halo's are needed, thus whalo1 is called.
      ! Only on the fine grid.
      
      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
           .false., .false.)
      
      ! Apply all boundary conditions to all blocks on this level.
      
      call applyAllBC(secondHalo)
      
      ! Exchange the solution. Either whalo1 or whalo2
      ! must be called.
      
      if( secondHalo ) then
         call whalo2(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      else
         call whalo1(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      endif

      call mpi_barrier(SUmb_comm_world, ierr)      
!
!     ******************************************************************
!     *                                                                *
!     * Update the force coefficients using the usual flow solver      *
!     * routine.                                                       *
!     *                                                                *
!     ******************************************************************
!
      !print *,' Calling original routines',level
      call metric(level) 

      call mdCreateSurfForceListLocal(sps,famID,startInd,endInd)

      

!     ******************************************************************
!     *                                                                *
!     * Compute the forces using the tapenade routines.                *
!     *                                                                *
!     ******************************************************************
        
!     ******************************************************************
!     * From forcesAndMoments.f90                                      *
!     * Determine the reference point for the moment computation in    *
!     * meters.                                                        *
!     ******************************************************************
       ! Allocate the memory for the local forces and initialize it
       ! to zero.  ModFamID is introduced to take famID == 0
       ! into account.

      modFamID = max(famID, 1_intType)
      nSurfNodesLoc = mdNSurfNodesLocal(modFamID)
      
      allocate(forceLoc(3,nSurfNodesLoc), stat=ierr)
      if(ierr /= 0)                             &
           call terminate("verifyForceCoupling", &
           "Memory allocation failure for forceLoc")
       
      forceLoc = zero
 
      !print *,' starting Tapenade routines'
      refPoint(1) = LRef*pointRef(1)
      refPoint(2) = LRef*pointRef(2)
      refPoint(3) = LRef*pointRef(3)

      yplusMax = zero

      spectralLoopAdj: do sps=1,nTimeIntervalsSpectral

         ! Loop over the number of local blocks.
         ii=0.0
         domainForcesLoop: do nn=1,nDom
            !print *,'setting pointers'
            ! Set some pointers to make the code more readable.
            call setPointersAdj(nn,groundlevel,sps)
            !print *,'allocating'
            allocate(xAdj(0:ie,0:je,0:ke,3), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for xAdj.")

  
            allocate(wAdj(0:ib,0:jb,0:kb,nw), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for wAdj.")

            allocate(pAdj(0:ib,0:jb,0:kb), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("Memory allocation failure for pAdj.")

            righthanded = flowDoms(nn,level,sps)%righthanded

            !print *,'copying stencil'
            ! Copy the coordinates into xAdj and
            ! Compute the face normals on the subfaces
            call copyADjointForcesStencil(wAdj,xAdj,alphaAdj,betaAdj,&
           MachAdj,machCoefAdj,prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
           rhoinfAdj, pinfAdj,murefAdj, timerefAdj,pInfCorrAdj,nn,level,sps,&
           liftIndex)
            !call copyADjointForcesStencil(wAdj,xAdj,nn,level,sps)

            bocoLoop: do mm=1,nBocos
               ! Check if this BC needs force integration (Either EulerWall or NSWall)

!s               invForce: if(BCType(mm) == EulerWall        .or. &
!s                            BCType(mm) == NSWallAdiabatic .or.  &
!s                            BCType(mm) == NSWallIsothermal) then
               !print *,'setting indicies',ii
                  ! Determine the range of cell indices of the owned cells
                  ! Notice these are not the node indices
                  iiBeg = BCData(mm)%icBeg
                  iiEnd = BCData(mm)%icEnd
                  jjBeg = BCData(mm)%jcBeg
                  jjEnd = BCData(mm)%jcEnd

                  i2Beg= BCData(mm)%inBeg+1; i2End = BCData(mm)%inEnd
                  j2Beg= BCData(mm)%jnBeg+1; j2End = BCData(mm)%jnEnd

                  call computeForceCouplingAdj(xAdj,wAdj,pAdj, &
                       iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End, &
                       mm,yplusMax,refPoint,nSurfNodesLoc,forceLoc,  &
                       nn,level,sps,righthanded,secondhalo,&
                       alphaAdj,betaAdj,machAdj,machcoefAdj,prefAdj,&
                       rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                       rhoinfAdj, pinfAdj,murefAdj, timerefAdj,pInfCorrAdj,&
                       liftIndex,ii)


            end do bocoLoop
            !print *,'deallocating p'
            deallocate(pAdj)
            !print *,'deallocating w'
            deallocate(wAdj)
            !print *,'deallocating x'
            deallocate(xAdj)
!            print *,'deallocation finished'
!            deallocate(wAdj)
      

         end do domainForcesLoop


      end do spectralLoopAdj

      call mpi_barrier(SUmb_comm_world, ierr)
      ! Root processor outputs results.
!      print *,'printing results'
      do sps=1,nTimeIntervalsSpectral
         print *,'ForceLoc MDSurfForce'
         do i=1,nSurfNodesLoc
            print *,'forcex',forceLoc(1,i),mdsurfforcelocal(1,i),forceLoc(1,i)-mdsurfforcelocal(1,i)
            print *,'forcey',forceLoc(2,i),mdsurfforcelocal(2,i),forceLoc(2,i)-mdsurfforcelocal(2,i)
            print *,'forcez',forceLoc(3,i),mdsurfforcelocal(3,i),forceLoc(3,i)-mdsurfforcelocal(3,i)
         enddo
      end do
       
      !print *,'finished computing forces'
      ! Flush the output buffer and synchronize the processors.
       
!      call f77flush()

       
      ! Output format.
       
10    format(1x,4a14)
20    format(1x,a,4(1x,e13.6))
      
    end subroutine verifyForceCouplingAdj
     
