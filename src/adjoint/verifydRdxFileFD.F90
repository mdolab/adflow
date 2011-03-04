!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdxfile.f90                              *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 06-12-2009                                      *
!     * Last modified: 06-12-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdxfileFD(level)
!
!     ******************************************************************
!     *                                                                *
!     *  This subroutine computes the values for dRdw and compares     *
!     *  the Tapenade result to the finite-difference results.         *
!     *                                                                *
!     ******************************************************************

#ifndef USE_NO_PETSC
      use blockPointers ! block (nDoms,flowDoms), globalCell
      use flowvarrefstate
      use communication
      use iteration     ! groundLevel
      use inputTimeSpectral ! spaceDiscr
      use inputIO

      !from old verify routine
      use ADjointPETSc, only: drdx,drdxfd,petscone,insert_values,petscierr,mat_final_assembly,petsc_viewer_draw_world,petsc_viewer_stdout_world,add_values
      !use FDPETSc, only: DRDWFD
      use precision
      !use blockPointers
      !use flowvarrefstate
      !use iteration
      use inputIteration
      !implicit none
     

      implicit none

!
!     Subroutine arguments
      integer(kind=intType), intent(in) :: level


!
!     Local variables 
!
      integer(kind=intType) :: i, j, k, n
      integer(kind=intType) :: iCell,jCell,kCell,liftIndex
      real(kind=realType), dimension(nw) :: dwL2
      real(kind=realType), dimension(nx, ny, nz, nw) :: dwerr
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes
!      real(kind=realType), dimension(4) :: time
      real(kind=realType) ::  timeOri
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdj
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjb
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wFD
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdj
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdjb
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xFD
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xAdj
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xAdjb
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xFD


      REAL(KIND=REALTYPE) :: machadj, machcoefadj, pinfcorradj,machgridadj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb, pinfcorradjb,machgridadjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: rhoinfadjb
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb

      character fileName*32, dataName*32
      real(kind=realType), dimension(nw) :: dwAdj,dwAdjb,dwAdjRef
      real(kind=realType), dimension(nw) :: dwAdjP, dwAdjM
      real(kind=realType) :: deltax
      real(kind=realType) :: timeAdj, timeFD, timeResAdj

      integer(kind=intType), dimension(nDom) :: maxglobalcell
      integer(kind=intType) :: idx, ii, jj, kk, idxnode, idxres, m, l,nnn
      integer(kind=intType) :: sps, nTime, max_nTime, nHalo, nn, discr
      !real(kind=realType), allocatable, dimension(:,:,:,:) :: dRdwAdj,dRdwFD1,dRdwFD2,dRdwErr
      real(kind=realType), allocatable, dimension(:,:,:,:) ::xtemp
      real(kind=realType), allocatable, dimension(:,:,:) ::ptemp

      REAL(KIND=REALTYPE), DIMENSION(3) :: rotcenteradj
      REAL(KIND=REALTYPE), DIMENSION(3):: rotrateadj
      REAL(KIND=REALTYPE) :: rotrateadjb(3)

      integer :: ierr, testnode
      logical :: fineGrid, correctForK, exchangeTurb,secondHalo
      REAL(KIND=REALTYPE) ::value,xref
 

      integer :: unitdRdx = 8,ierror
      character(len = 20)::outfile,testfile
      write(testfile,100) myid!12
100   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,101) trim(testfile)!testfile
101   format("FDdRdxfile",a,".out")
      unitdrdx = 8+myID

      
      open (UNIT=unitdRdx,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdxFile", &
           "Something wrong when &
           &calling open")

      !FULL FD Variables
      
!!$      !real(kind=realType), dimension(nx, ny, nz, nw,nx) :: dRdwErr, dRdwErrRel
!!$      real(kind=realType), dimension(nw) :: dRdwL2, dRdwL2Rel
!!$      real(kind=realtype), dimension(0:ib,0:jb,0:kb,1:nw)::dwp,dwm,dwtemp
!!$      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw) :: wtemp
!!$      real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ptemp
!!$      integer(kind=intType) :: istate, jstate, kstate,ires
!!$      real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdjtemp
!!$      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjtemp
!!$      real(kind=realType), dimension(nx,ny,nz,nw) :: wFD2
!!$      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdw
!!$      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdwFD

!     ******************************************************************
!     *                                                                *
!     *  Begin execution                                               *
!     *                                                                *
!     ******************************************************************

      if( myID==0 ) write(*,*) "Running verifydRdxFile..."

      currentLevel = level
      !discr        = spaceDiscr
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


!
!     ******************************************************************
!     *                                                                *
!     * Exchange halo data to make sure it is up-to-date.              *
!     * (originally called inside "rungeKuttaSmoother" subroutine).    *
!     *                                                                *
!     ******************************************************************
!
      ! Exchange the pressure if the pressure must be exchanged early.
      ! Only the first halo's are needed, thus whalo1 is called.
      ! Only on the fine grid.
      
      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
           .false., .false.)
      
      ! Apply all boundary conditions to all blocks on this level.
      print *,'second halo',secondhalo,exchangePressureEarly
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

     ! Reset the values of rkStage and currentLevel, such that
      ! they correspond to a new iteration.

      rkStage = 0
      currentLevel = groundLevel

      ! Compute the latest values of the skin friction velocity.
      ! The currently stored values are of the previous iteration.
      
      call computeUtau
      
      ! Apply an iteration to the turbulent transport equations in
      ! case these must be solved segregatedly.
      
      if( turbSegregated ) call turbSolveSegregated

      ! Compute the time step.
      
      call timeStep(.false.)
      
      ! Compute the residual of the new solution on the ground level.
      
      if( turbCoupled ) then
         call initres(nt1MG, nMGVar)
         call turbResidual
      endif

      call initres(1_intType, nwf)
      call residual
      

      ! Get the final time for original routines.
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
        call cpu_time(time(2))
        timeOri = time(2)-time(1)
      endif


!     ******************************************************************
!     *                                                                *
!     *  Compute d(dw)/d(w) using Tapenade                             *
!     *                                                                *
!     ******************************************************************

      ! Get the initial AD time.

      call mpi_barrier(SUmb_comm_world, ierr)
      if( myID==0 ) call cpu_time(time(1))
      
      !allocate memory for FD
      allocatedomains: do nn = 1,ndom
         print *,'domain',nn
         groundLevel = 1
         sps = 1
         call setPointersAdj(nn,1,sps)
         allocate(flowDoms(nn,level,sps)%dwp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
         allocate(flowDoms(nn,level,sps)%dwm(0:ib,0:jb,0:kb,1:nw),stat=ierr)
         allocate(flowDoms(nn,level,sps)%dwtemp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
         dwtemp = zero
         dwp = zero
         dwm = zero
      end do allocatedomains


      storedomains: do nn = 1,ndom
         groundLevel = 1
         sps = 1
         call setPointersAdj(nn,1,sps)
         dwtemp = dw
      end do storedomains
      !zero the matrix for dRdW Insert call
      call MatZeroEntries(dRdxFD,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("verifydRdxFileFD", "Error in MatZeroEntries drdxfd")

      deltax = 1e-5!1.d-5
      !print *, "deltaw=", deltaw
      !wFD2(:,:,:,:) = 0. 
      call cpu_time(time(3))
      domains: do nn = 1,ndom
         print *,'domain',nn
         groundLevel = 1
         sps = 1
         call setPointersAdj(nn,1,sps)
         allocate(xtemp(0:ie,0:je,0:ke,1:3),ptemp(0:ib,0:jb,0:kb))
         xtemp(:,:,:,:) = x(:,:,:,:)
         ptemp = p
         do kCell = 1,kl!0, ke
            do jCell = 1,jl!0, je
               do iCell = 1,il!0, ie
                  !print *,'ie',icell,ie,jcell,je,kcell,ke
                  do m = 1, 3
                     do nnn = 1,ndom
                        call setPointersAdj(nnn,1,sps)
                        dw = dwtemp
                     enddo
                     call setPointersAdj(nn,1,sps)
                     
                     xref =x(icell,jcell,kcell,m) 
                     x(icell,jcell,kcell,m)  = xref+ deltax
                     p = ptemp
                     call xhalo(groundlevel)
                     call metric(groundlevel)
                     !call setPointers(nn,level,sps)
                     call checkSymmetry(groundlevel)
                     !call referenceState
                      
                     !call setFlowInfinityState
                     
                     call initres(1_intType, nwf)
                     
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
                     !if( secondHalo ) then
                     !   !  write(*,*)'2ndHalo..........'
                     !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
                     !        .true., .true.)
                     !else
                     !   ! write(*,*)'1stHalo..........'
                     !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
                     !        .true., .true.)
                     !endif
                     
                     call residual
                     !save residual
                     do nnn = 1,ndom
                        groundLevel = 1
                        sps = 1
                        call setPointersAdj(nnn,1,sps)
                        dwp = dw
                     end do
                     !print *, "Called Residual =",nn,domain
                     do nnn = 1,ndom
                        call setPointersAdj(nnn,1,sps)
                        dw = dwtemp
                        
                     enddo
                     !nn=domain
                     groundlevel = 1
                     sps = 1
                     call setPointersAdj(nn,1,sps)
                     p(:,:,:) = ptemp(:,:,:)
                     x(icell,jcell,kcell,m)  = xref- deltax
                     
                     call xhalo(groundlevel)
                     call metric(groundlevel)
                     !call setPointers(nn,level,sps)
                     call checkSymmetry(groundlevel)
                     !call referenceState
                      
                     !call setFlowInfinityState
                     
                     call initres(1_intType, nwf)
                     
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
                     !if( secondHalo ) then
                     !   !  write(*,*)'2ndHalo..........'
                     !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
                     !        .true., .true.)
                     !else
                     !   ! write(*,*)'1stHalo..........'
                     !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
                     !        .true., .true.)
                     !endif
                     
                     call residual
                     !save residual
                     do nnn = 1,ndom
                        groundLevel = 1
                        sps = 1
                        call setPointersAdj(nnn,1,sps)
                        dwm = dw
                     end do
                     !print *, "Called Residual =",nn,domain
                     
                     !nn=domain
                     groundlevel = 1
                     sps = 1
                     call setPointersAdj(nn,1,sps)

                     idxnode   = globalnode(iCell,jCell,kCell)*3+m 
                     !testnode =  globalnode(iCell,jCell,kCell)
                     do nnn = 1,ndom
                        call setPointersAdj(nnn,1,sps)
                        DO I=2,Il
                           DO J=2,Jl
                              DO K=2,Kl
                                 do n = 1,nw
                                    idxres = globalCell(i,j,k)*nw+n
                                    
                                    ! Loop over output cell residuals (R)
                                    ! Loop over location of output (R) cell of residual
                                    value = (dwp(i,j,k,n)-dwm(i,j,k,n))/(2*deltax)
!                                    value = (dwp(i,j,k,n)-dwtemp(i,j,k,n))/(deltax)
                                    if ((idxres-1)>=0 .and. (idxnode-1)>=0)then
                                       if(value>1e-10)then
                                          print *,'dx',value,dwp(i,j,k,n),dwm(i,j,k,n),dwtemp(i,j,k,n),(2*deltax)
                                          write(unitdrdx,13) idxnode,idxres,m,icell,jcell,kcell,nn,n,k,j,i,nnn,value
                                          !write(unitWarp,13) xderiv,i,j,k,n,nnn,nn,mm,ll
13                                        format(1x,'drdx',12I8,f18.10)
                                       endif
                                    end if
                                 end do
                              end do
                           end do
                        end do
                        !call setPointersAdj(nn,1,sps)
                     enddo
                     call setPointersAdj(nn,1,sps)
                     x(icell,jcell,kcell,m) =xref
                     p = ptemp
                  end do
               enddo
            end do
         enddo
         call setPointersAdj(nn,1,sps)
         x(:,:,:,:) = xtemp(:,:,:,:)
         deallocate(xtemp,ptemp)
      enddo domains


      !deallocate memory for FD
      deallocatedomains: do nn = 1,ndom
         print *,'domain',nn
         groundLevel = 1
         sps = 1
         call setPointersAdj(nn,1,sps)
         deallocate(flowDoms(nn,level,sps)%dwp)
         deallocate(flowDoms(nn,level,sps)%dwm)
         deallocate(flowDoms(nn,level,sps)%dwtemp)
      end do deallocatedomains

      ! Get new time and compute the elapsed AD time.
      !      stop
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
         call cpu_time(time(2))
         timeAdj = time(2)-time(1)
      endif
      
      
      !Print *,'barriercall',myID
      call mpi_barrier(SUmb_comm_world, ierr)
      
      close(unitdrdx)
      
      call mpi_barrier(SUmb_comm_world, ierr)
      
      
111   format(4I4, 3ES22.14)
      
#endif      
      
    end subroutine verifydRdxfileFD
