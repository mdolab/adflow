!________________________________________________________________________
!
!     File:          verifyResiduals.f90
!     Author:        Joaquim R. R. A. Martins, C.A.(Sandy) Mader
!     Starting date: 04-06-2006
!     Last modified: 06-11-2008
!________________________________________________________________________

      subroutine verifyResiduals(level)
!
!________________________________________________________________________
!
!     This subroutine computes all residuals cell by cell by calling
!     residualAdj and then checks them against the current residuals
!     that SUmb has calculated.
!________________________________________________________________________

      !&&use precision
      use blockPointers
      !&&use flowvarrefstate
      use communication
      !&&use constants
      ! use flowVarRefState
      !use inputIteration
      use inputPhysics
      !&&use iteration
      !use killSignals
      !use monitor
      !&&implicit none

!!$!from transfer to fine grid
      use BCTypes
!!$      !use blockPointers
!!$      !use flowVarRefState
!!$      !use inputIteration
      use inputTimeSpectral
!!$      !use iteration
!!$      ! implicit none


!*****************
!from execute mg cycle
      use flowVarRefState
      use iteration
      use inputIteration
      implicit none

!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level!, sps
!
!      Local variables.
!
       integer(kind=intType) :: nn,icell, jcell, kcell,ii, jj, kk
!************************
!
!     Local variables
!
      integer(kind=intType) :: i, j, k,q,sps=1
      real(kind=realType), dimension(nw) :: dwL2
     ! real(kind=realType), dimension(0:nx,0:ny,0:nz, nw) :: dwerr
      real(kind=realType), dimension(nx, ny, nz, nw) :: dwerr
!*****
      real(kind=realType), dimension(0:ib,0:jb,0:kb, 1:nw) :: wtemp,wtemp2,dwref
      real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ptemp,ptemp2
      real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdjtemp
!*******
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes,machref
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdj
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xAdj
      real(kind=realType), dimension(nw) :: dwAdjRef
      character fileName*32, dataName*32

      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,machCoefAdj
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj

      !********************************
      logical :: secondHalo =.true.! .false.!.true.
      ! Set the value of secondHalo, depending on the situation.
      !secondHalo = .true.

      logical :: corrections =  .false.!.true.
      logical :: fineGrid, correctForK, exchangeTurb

      integer(kind=intType) :: coarseLevel, nVarInt 
      integer :: ierr


!File Parameters
      integer :: unit = 8,unit1 = 30,unit3 = 12,unit4 = 13,ierror
      character(len = 9)::outfile
      
      outfile = "oldres.txt"
      
      open (UNIT=unit,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifyResiduals", &
           "Something wrong when &
           &calling open")
      
      outfile = "newres.txt"
      
      open (UNIT=unit1,File=outfile,status='replace',action='write',iostat=ierror)
       if(ierror /= 0)                        &
            call terminate("verifyResiduals", &
            "Something wrong when &
            &calling open2")

  outfile = "oldw.txt"
      
      open (UNIT=unit3,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifyResiduals", &
           "Something wrong when &
           &calling open")
      
!!$      outfile = "neww.txt"
!!$      
!!$      open (UNIT=unit4,File=outfile,status='replace',action='write',iostat=ierror)
!!$       if(ierror /= 0)                        &
!!$            call terminate("verifyResiduals", &
!!$            "Something wrong when &
!!$            &calling open2")










!----------------------
!initialize residuals
!----------------------
      ! Set the grid level of the current MG cycle, the value of the
      ! discretization and the logical fineGrid.

      currentLevel = level
!      discr        = spaceDiscr
      fineGrid     = .true.

      ! Determine whether or not the total energy must be corrected
      ! for the presence of the turbulent kinetic energy.

      if( kPresent ) then
         if((currentLevel <= groundLevel) .or. turbCoupled) then
            correctForK = .true.
         else
            correctForK = .false.
         endif
      else
         correctForK = .false.
      endif

      ! and whether or not turbulence variables should be exchanged
      exchangeTurb = .false.

      
      ! Set the value of secondHalo, depending on the situation.
      ! In the full MG (currentLevel < groundLevel) the second halo is
      ! always set; otherwise only on the finest mesh in the current mg
      ! cycle.

      if(currentLevel <= groundLevel) then
         secondHalo = .true.
      else
         secondHalo = .false.
      endif
      
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Exchange halo data to make sure it is up-to-date.              *
!!$!     * (originally called inside "rungeKuttaSmoother" subroutine).    *
!!$!     *                                                                *
!!$!     ******************************************************************
!!$!
!!$      ! Exchange the pressure if the pressure must be exchanged early.
!!$      ! Only the first halo's are needed, thus whalo1 is called.
!!$      ! Only on the fine grid.
!!$      
!!$      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
!!$           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
!!$           .false., .false.)
!!$      
!!$      ! Apply all boundary conditions to all blocks on this level.
!!$      
!!$      call applyAllBC(secondHalo)
!!$      
!!$      ! Exchange the solution. Either whalo1 or whalo2
!!$      ! must be called.
!!$      
!!$      if( secondHalo ) then
!!$         call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$              .true., .true.)
!!$      else
!!$         call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$              .true., .true.)
!!$      endif
!!$      
!!$      ! Reset the values of rkStage and currentLevel, such that
!!$      ! they correspond to a new iteration.
!!$
!!$      rkStage = 0
!!$      currentLevel = groundLevel
!!$
!!$      ! Compute the latest values of the skin friction velocity.
!!$      ! The currently stored values are of the previous iteration.
!!$      
!!$      call computeUtau
!!$      
!!$      ! Apply an iteration to the turbulent transport equations in
!!$      ! case these must be solved segregatedly.
!!$      
!!$      if( turbSegregated ) call turbSolveSegregated
!!$
!!$      ! Compute the time step.
!!$      
!!$      call timeStep(.false.)
!!$      
!!$      ! Compute the residual of the new solution on the ground level.
!!$      
!!$      if( turbCoupled ) then
!!$         call initres(nt1MG, nMGVar)
!!$         call turbResidual
!!$      endif
!!$
!!$      call initres(1_intType, nwf)
!!$      call residual
!!$
!!$    ! Get the final time for original routines.
!!$      call mpi_barrier(SUmb_comm_world, ierr)


!------------------------
! reference residual
!------------------------

       
      nVarInt = nw
      if( corrections ) nVarInt = nMGVar
!________________________________________________________________________
!
!     Begin execution
!________________________________________________________________________

write(*,*) 'in verify residuals'

!
!     Compute difference in residuals

!Baseline residual calculation

       !dw(:,:,:, :) =0.0

       wtemp(:,:,:,:) = w(:,:,:,:)
       ptemp(:,:,:) = p(:,:,:)
       call residual
       dwref(:,:,:,:) = dw(:,:,:,:)

      ! w(15,2,5,:) = w(15,2,5,:)+0.005
!!$       i = 33
!!$       j = 2
!!$       k = 12

       i = 12
       j = 4
       k = 3

       w(i,j,k,:) = w(i,j,k,:)+0.005
     !  machref = mach
      ! mach = machref+0.005
       wtemp2(:,:,:,:) = w(:,:,:,:)

       

      ! call referenceState
       
      ! call setFlowInfinityState

       call computePressureAdj(w(i-2:i+2, j-2:j+2, k-2:k+2,:), pAdjtemp)

       p(i, j, k) = pAdjtemp(0,0,0)
       !wtemp(:,:,:,:) = w(:,:,:,:)

       write(unit,*)'States0'
       write(unit1,*)'States0'
       do kk= 0,kb
          do jj= 0,jb
             do ii= 0,ib
                do nn = 1,5
                   write(unit,21) w(ii,jj,kk,nn),ii,jj,kk,nn
                   write(unit1,21) wtemp(ii,jj,kk,nn),ii,jj,kk,nn
!21                format(1x,'w ',f18.10,4I4)
                enddo
             enddo
          enddo
       enddo
       write(unit,*)'doneStates0'
       write(unit1,*)'doneStates0'

      ! Exchange the pressure if the pressure must be exchanged early.
      ! Only the first halo's are needed, thus whalo1 is called.
      ! Only on the fine grid.
      
      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
           .false., .false.)
       groundLevel = 1

 
!!$       write(unit3,*)'initres'
!!$       do kk= -2,2
!!$          do jj= -2,2
!!$             do ii= -2,2
!!$                do nn = 1,5
!!$                   icell = i+ii
!!$                   jcell = j+jj
!!$                   kcell = k+kk
!!$                   write(unit3,22)w(icell,jcell,kcell,nn),icell,jcell,kcell,nn
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo


       call applyAllBC(secondHalo)

!!$       write(unit3,*)'applybcs'
!!$       do kk= -2,2
!!$          do jj= -2,2
!!$             do ii= -2,2
!!$                do nn = 1,5
!!$                   icell = i+ii
!!$                   jcell = j+jj
!!$                   kcell = k+kk
!!$                   write(unit3,22)w(icell,jcell,kcell,nn),icell,jcell,kcell,nn
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo

       write(unit,*)'applybcs'
       write(unit1,*)'applybcs'
       do kk= 0,kb
          do jj= 0,jb
             do ii= 0,ib
                do nn = 1,5
                   write(unit,21) w(ii,jj,kk,nn),ii,jj,kk,nn
                   write(unit1,21) wtemp(ii,jj,kk,nn),ii,jj,kk,nn
!21                format(1x,'w ',f18.10,4I4)
                enddo
             enddo
          enddo
       enddo
       write(unit,*)'doneStates0'
       write(unit1,*)'doneStates0'

       ! Exchange the solution. Either whalo1 or whalo2
       ! must be called.

       if( secondHalo ) then
         write(*,*)'2ndHalo..........'
         call whalo2(currentLevel, 1_intType, nVarInt, .true., &
                     .true., .true.)
       else
          write(*,*)'1stHalo..........'
         call whalo1(currentLevel, 1_intType, nVarInt, .true., &
                     .true., .true.)
       endif

!!$       write(unit3,*)'whalo'
!!$       do kk= -2,2
!!$          do jj= -2,2
!!$             do ii= -2,2
!!$                do nn = 1,5
!!$                   icell = i+ii
!!$                   jcell = j+jj
!!$                   kcell = k+kk
!!$                   write(unit3,22)w(icell,jcell,kcell,nn),icell,jcell,kcell,nn
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo

       write(unit,*)'whalo'
       write(unit1,*)'whalo'
       do kk= 0,kb
          do jj= 0,jb
             do ii= 0,ib
                do nn = 1,5
                   write(unit,21) w(ii,jj,kk,nn),ii,jj,kk,nn
                   write(unit1,21) wtemp(ii,jj,kk,nn),ii,jj,kk,nn
!21                format(1x,'w ',f18.10,4I4)
                enddo
             enddo
          enddo
       enddo
       write(unit,*)'doneStates0'
       write(unit1,*)'doneStates0'

      print *, "Calling Residual ="

       call initres(1_intType, nwf)
       write(unit,*)'initres0'
       write(unit1,*)'initres0'
       do kk= 0,kb
          do jj= 0,jb
             do ii= 0,ib
                do nn = 1,5
                   write(unit,21) w(ii,jj,kk,nn),ii,jj,kk,nn
                   write(unit1,21) wtemp(ii,jj,kk,nn),ii,jj,kk,nn
!21                format(1x,'w ',f18.10,4I4)
                enddo
             enddo
          enddo
       enddo
       write(unit,*)'doneStates0'
       write(unit1,*)'doneStates0'

       call residual
       print *, "Called Residual =", w(i,j,k,1),w(4,3,3,1),w(i,j,k,1)-w(4,3,3,1)
       
       write(unit3,*)'residual'
       do kk= -2,2
          do jj= -2,2
             do ii= -2,2
                do nn = 1,5
                   icell = i+ii
                   jcell = j+jj
                   kcell = k+kk
                   write(unit3,22)w(icell,jcell,kcell,nn),icell,jcell,kcell,nn
                enddo
             enddo
          enddo
       enddo

       write(unit,*)'resStates'
       write(unit1,*)'resStates'
       do kk= 0,kb
          do jj= 0,jb
             do ii= 0,ib
                do nn = 1,5
                   write(unit,21) w(ii,jj,kk,nn),ii,jj,kk,nn
                   write(unit1,21) wtemp(ii,jj,kk,nn),ii,jj,kk,nn
21                format(1x,'w ',f18.10,4I4)
                enddo
             enddo
          enddo
       enddo
       write(unit,*)'doneStates'
       write(unit1,*)'doneStates'
       
       ptemp2(:,:,:) = p(:,:,:)
       p(:,:,:) = ptemp(:,:,:)

       w(:,:,:,:) = wtemp2(:,:,:,:)
       print *, "reset State =", w(i,j,k,1),w(5,1,5,1),w(i,j,k,1)-w(5,1,5,1)

!***********************************
                 
      dwL2(:) = 0.
      dwerr(:,:,:,:) = 0.
      call cpu_time(time(1))
      call setPointersAdj(1,level,sps)
      
      do k= 2, kl
         do j= 2, jl
            do i= 2, il
               !write(*,*)'adj loop'
               call copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,&
                    MachAdj,MachCoefAdj,i, j, k,prefAdj,&
                    rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                    rhoinfAdj, pinfAdj,&
                    murefAdj, timerefAdj,pInfCorrAdj)
               !(wAdj, xAdj, iCell, jCell, kCell)
               
               ! Compute the total residual.
               ! This includes inviscid and viscous fluxes, artificial
               ! dissipation, and boundary conditions.                   
               !print *,'Calling compute ADjoint'
               call computeRAdjoint(wAdj,xAdj,dwAdjRef,alphaAdj,&
                    betaAdj,MachAdj, MachCoefAdj,&
                    i, j,  k, &
                    nn,sps, correctForK,secondHalo,prefAdj,&
                    rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                    rhoinfAdj, pinfAdj,&
                    murefAdj, timerefAdj,pInfCorrAdj)
               
               !wAdj(-2:2,-2:2,-2:2,1:nw) = &
               !     w(i-2:i+2, j-2:j+2, k-2:k+2, 1:nw)               
               !call residualAdj(wAdj, dwAdjRef, i, j ,k)               
               dwerr(i-1, j-1, k-1, :) = (dwAdjRef(:) - dw(i, j, k, :))
!!$       if ((i==33).and.(j == 2).and.(k == 12))then
!!$          print *,'testing'
!!$
!!$       write(unit4,*)'residual'
!!$       do kk= -2,2
!!$          do jj= -2,2
!!$             do ii= -2,2
!!$                do nn = 1,5
!!$                   icell = i+ii
!!$                   jcell = j+jj
!!$                   kcell = k+kk
!!$                   write(unit4,22)wadj(ii,jj,kk,nn),icell,jcell,kcell,nn
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
!!$       endif


!!$               !dwerr(:, :,:, :) = 0.0
!!$               dwerr(i-1, j-1, k-1, :) = 0.0
!!$               !if ((k == 13).or.(j==5).or.(i==2))then
!!$               if ((i==14))then
!!$                  dwerr(i-1, j-1, k-1, :) = 25
!!$               end if
!----------------------------------------------
               
               do q = 1,nw
                  !write errors one by one to a file for xx diff
                  write(unit,20) dw(i,j,k,q),ptemp2(i,j,k),i,j,k,q
!                  write(unit1,20) dwref(i,j,k,q),ptemp2(i,j,k),i,j,k,q
                  write(unit1,20) dwAdjRef(q), p(i,j,k),i,j,k,q
20                format(1x,'dw ',f18.10,' p ',f18.10,4I4)
               end do
!----------------------------------------------
               dwL2(:) = dwL2(:) + dwerr(i-1, j-1, k-1, :)**2 
            end do
         end do
      end do

      w(:,:,:,:) = wtemp(:,:,:,:)

      !mach = machref

      call cpu_time(time(2))
      timeRes = (time(2)-time(1))/(nx*ny*nz)
      print *, "Average time for each cell residual calculation =", timeRes
      dwL2(:) = sqrt(dwL2(:))
      print *, "L-2 norm of differences ="
      print *, dwL2(:)

      write(*,*)'shape', shape(dwerr)
!
!     Write residuals to the CGNS file
!
      fileName = "bump_residual_errors.cgns"
      call writeCGNSMesh(fileName)    
      dataName = "ResErr"
      call writeCGNSData(fileName, dataName, dwerr(:,:,:,:))

22    format(1x,'w ',f18.10,4I4)

      end subroutine verifyResiduals
