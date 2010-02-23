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
      use inputTimeSpectral!nTimeIntervalsSpectral
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
      integer(kind=intType) :: i, j, k,n,sps=1,liftIndex,nnn,ssps,sps2
      real(kind=realType), dimension(nw) :: dwL2
     ! real(kind=realType), dimension(0:nx,0:ny,0:nz, nw) :: dwerr
      real(kind=realType), dimension(nx, ny, nz, nw) :: dwerr
!*****
      real(kind=realType), dimension(0:ib,0:jb,0:kb, 1:nw) :: wtemp,wtemp2,dwref
      real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ptemp,ptemp2
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) :: pAdjtemp
      real(kind=realType), dimension(0:ie,0:je,0:ke, 1:3) :: xtemp,xtemp2
!*******
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes,machref
      real(kind=realType),dimension(3)::xref,deltax
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral) :: wAdj
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral)  :: xAdj
      real(kind=realType), dimension(nw,nTimeIntervalsSpectral) :: dwAdj
      character fileName*32, dataName*32

      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,machCoefAdj,machGridAdj
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
      real(kind=realType), dimension(3) ::rotRateAdj,rotCenterAdj,rotrateadjb
      real(kind=realType), dimension(2,2,2,3,nTimeIntervalsSpectral)::xBlockCornerAdj

      !********************************
      logical :: secondHalo =.true.! .false.!.true.
      ! Set the value of secondHalo, depending on the situation.
      !secondHalo = .true.

      logical :: corrections =  .false.!.true.
      logical :: fineGrid, correctForK, exchangeTurb

      integer(kind=intType) :: coarseLevel, nVarInt 
      integer :: ierr


!File Parameters
      integer :: unitRes = 8,unitresAD = 30,unitx = 12,unitxAD = 13,ierror
      integer ::iii,iiii,jjj,jjjj,kkk,kkkk,nnnn,istart,jstart,kstart,iend,jend,kend
      character(len = 16)::outfile
      
      outfile = "originalres.txt"
      
      open (UNIT=unitRes,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifyResiduals", &
           "Something wrong when &
           &calling open")
      
      outfile = "ADres.txt"
      
      open (UNIT=unitresAD,File=outfile,status='replace',action='write',iostat=ierror)
       if(ierror /= 0)                        &
            call terminate("verifyResiduals", &
            "Something wrong when &
            &calling open2")

  outfile = "xoriginal.txt"
      
      open (UNIT=unitx,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifyResiduals", &
           "Something wrong when &
           &calling open")
      
      outfile = "xAD.txt"
      
      open (UNIT=unitxAD,File=outfile,status='replace',action='write',iostat=ierror)
       if(ierror /= 0)                        &
            call terminate("verifyResiduals", &
            "Something wrong when &
            &calling open2")




!close(unitxAD)
!close(unitx)






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
      
      call applyAllBC(secondHalo)
      
      ! If case this routine is called in full mg mode call the mean
      ! flow boundary conditions again such that the normal momentum
      ! boundary condition is treated correctly.
      
      if(.not. corrections) call applyAllBC(secondHalo)
      
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


!------------------------
! reference residual
!------------------------
      print *,'storing reference states'
      groundLevel = 1
      do sps = 1,nTimeIntervalsSpectral!sps = 1
         !Store current state values and residuals
         do nn=1,ndom
            !print *,'domains 1',nn
            call setPointersAdj(nn,level,sps)
            wtmp = w
            ptmp = p
            dwtmp = dw
         end do
      end do
      nVarInt = nw
      if( corrections ) nVarInt = nMGVar
!________________________________________________________________________
!
!     Begin execution
!________________________________________________________________________

      write(*,*) 'in verify residuals'

!
!     Compute difference in residuals
      deltax(:) = 0.!0.05
      deltax(3) = 0.0
       ii = 6!11!1!5
       jj = 2!6!1!4
       kk = 2!6!1!3
       ssps=1
       nn = 1
       call setPointers(nn,level,ssps)
       xref = x(ii,jj,kk,:)
       x(ii,jj,kk,:) = xref+deltax
!!$       x(ii,jj,kk,1) = xref(1)+deltax
!!$       x(ii,jj,kk,2) = xref(2)+deltax
!!$       x(ii,jj,kk,3) = xref(3)+deltax
       w(ii,jj,kk,:) = w(ii,jj,kk,:)!+0.005
       machref = mach
       mach = machref!+0.005
       
       call referenceState
       
       call setFlowInfinityState
!!$     !print out x
!!$       do sps = 1,nTimeIntervalsSpectral
!!$          do nnnn=1,ndom
!!$          call setPointersAdj(nnnn,1,sps)
!!$          do iii = 2,il
!!$             do jjj = 2,jl
!!$                do kkk = 2,kl
!!$                   istart = -3
!!$                   jstart = -3
!!$                   kstart = -3
!!$                   iend = 2
!!$                   jend = 2
!!$                   kend = 2
!!$                   if(iii==2) istart=-2
!!$                   if(jjj==2) jstart=-2
!!$                   if(kkk==2) kstart=-2
!!$                   if(iii==il) iend=1
!!$                   if(jjj==jl) jend=1
!!$                   if(kkk==kl) kend=1
!!$                   do iiii = istart,iend
!!$                      do jjjj = jstart,jend
!!$                         do kkkk = kstart,kend
!!$                            do n = 1,3
!!$                               do sps2 = 1,nTimeIntervalsSpectral
!!$                                  call setPointersAdj(nnnn,1,sps2)
!!$                                  i = iii+iiii
!!$                                  j = jjj+jjjj
!!$                                  k = kkk+kkkk
!!$                                  write(unitx,11) i,j,k,n,nnnn,sps,sps2,x(i,j,k,n)
!!$11                                format(1x,'xadj',7I8,f20.14)
!!$                               enddo
!!$                            enddo
!!$                         enddo
!!$                      enddo
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          enddo 
!!$          call setPointersAdj(nn,1,sps)
!!$       enddo
!!$    end do
       call xhalo(level)
!!$       do sps = 1,nTimeIntervalsSpectral
!!$       do nnnn=1,ndom
!!$          call setPointersAdj(nnnn,1,sps)
!!$          do iii = 2,il
!!$             do jjj = 2,jl
!!$                do kkk = 2,kl
!!$                   istart = -3
!!$                   jstart = -3
!!$                   kstart = -3
!!$                   iend = 2
!!$                   jend = 2
!!$                   kend = 2
!!$                   if(iii==2) istart=-2
!!$                   if(jjj==2) jstart=-2
!!$                   if(kkk==2) kstart=-2
!!$                   if(iii==il) iend=1
!!$                   if(jjj==jl) jend=1
!!$                   if(kkk==kl) kend=1
!!$                   do iiii = istart,iend
!!$                      do jjjj = jstart,jend
!!$                         do kkkk = kstart,kend
!!$                            do n = 1,3
!!$                               do sps2 = 1,nTimeIntervalsSpectral
!!$                                  call setPointersAdj(nnnn,1,sps2)
!!$                                  i = iii+iiii
!!$                                  j = jjj+jjjj
!!$                                  k = kkk+kkkk
!!$                                  write(unitx,11) i,j,k,n,nnnn,sps,sps2,x(i,j,k,n)
!!$11                                format(1x,'xadj',7I8,f20.14)
!!$                               end do
!!$                            enddo
!!$                         enddo
!!$                      enddo
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          enddo
!!$          call setPointersAdj(nn,1,sps)
!!$       enddo
!!$    end do
       call metric(level)
       call checkSymmetry(level)
       !print out x
       do sps = 1,nTimeIntervalsSpectral
          do nnnn=1,1!ndom
          call setPointersAdj(nnnn,1,sps)
          do iii = 2,il
             do jjj = 2,jl
                do kkk = 2,kl
                   istart = -3
                   jstart = -3
                   kstart = -3
!!$                   istart = -2!3
!!$                   jstart = -2!3
!!$                   kstart = -2!3
                   iend = 2
                   jend = 2
                   kend = 2
                   if(iii==2) istart=-2
                   if(jjj==2) jstart=-2
                   if(kkk==2) kstart=-2
                   if(iii==il) iend=1
                   if(jjj==jl) jend=1
                   if(kkk==kl) kend=1
                   do iiii = istart,iend
                      do jjjj = jstart,jend
                         do kkkk = kstart,kend
                            do n = 1,3
                            !do n = 1,nw
                               do sps2 = 1,nTimeIntervalsSpectral
                                  call setPointersAdj(nnnn,1,sps2)
                                  i = iii+iiii
                                  j = jjj+jjjj
                                  k = kkk+kkkk
                                  write(unitx,11) i,j,k,n,nnnn,sps,sps2,x(i,j,k,n)
                                  !write(unitx,11) i,j,k,n,nnnn,sps,sps2,w(i,j,k,n)
11                                format(1x,'res',7I8,f20.14)
                                  call setPointersAdj(nnnn,1,sps)
                               end do
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          call setPointersAdj(nn,1,sps)
       enddo
    end do
        do sps = 1,nTimeIntervalsSpectral
           do nnn=1,ndom
              call setPointersAdj(nnn,1,sps)
              call computeForcesPressureAdj(w,p)
           end do
        end do
       call setPointersAdj(nn,1,ssps)
  
!!$     !print out w
!!$       do sps = 1,nTimeIntervalsSpectral
!!$          do nnnn=1,ndom
!!$          call setPointersAdj(nnnn,1,sps)
!!$          do iii = 2,il
!!$             do jjj = 2,jl
!!$                do kkk = 2,kl
!!$                   istart = -2!3
!!$                   jstart = -2!3
!!$                   kstart = -2!3
!!$                   iend = 2
!!$                   jend = 2
!!$                   kend = 2
    !              if(iii==2) istart=-1
    !              if(jjj==2) jstart=-1
    !              if(kkk==2) kstart=-1
    !              if(iii==il) iend=1
    !              if(jjj==jl) jend=1
    !              if(kkk==kl) kend=1
!!$                   !do sps2 = 1,nTimeIntervalsSpectral
!!$                   do iiii = istart,iend
!!$                      do jjjj = jstart,jend
!!$                         do kkkk = kstart,kend
!!$                            !do n = 1,3
!!$                               !do n = 1,1!nw
!!$                               do n = 1,nw
!!$                               do sps2 = 1,nTimeIntervalsSpectral
!!$                                  call setPointersAdj(nnnn,1,sps2)
!!$                                  i = iii+iiii
!!$                                  j = jjj+jjjj
!!$                                  k = kkk+kkkk
!!$                                 
!!$                                  write(unitx,11) i,j,k,n,nnnn,sps,sps2,w(i,j,k,n)
!!$                                  !write(unitx,11) i,j,k,n,nnnn,sps,sps2,p(i,j,k)
!!$                                  !write(unitx,11) i,j,k,n,nnnn,sps,sps2,sFaceI(i,j,k)
!!$                                  !write(unitx,11) i,j,k,n,nnnn,sps,sps2,s(i,j,k,n)
!!$11                                format(1x,'wadj',7I8,f20.14)
!!$                                  call setPointersAdj(nnnn,1,sps)
!!$                               end do
!!$                            enddo
!!$                         enddo
!!$                      enddo
!!$                   enddo
!!$                !end do
!!$                enddo
!!$             enddo
!!$          enddo
!!$          call setPointersAdj(nn,1,sps)
!!$       enddo
!!$    end do
       ! Exchange the pressure if the pressure must be exchanged early.
       ! Only the first halo's are needed, thus whalo1 is called.
       ! Only on the fine grid.                
       
       if(exchangePressureEarly .and. currentLevel <= groundLevel) &
            call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
            .false., .false.)
       
       
       call applyAllBC(secondHalo)
!!$      !print out w
!!$       do sps = 1,nTimeIntervalsSpectral
!!$          do nnnn=1,ndom
!!$          call setPointersAdj(nnnn,1,sps)
!!$          do iii = 2,il
!!$             do jjj = 2,jl
!!$                do kkk = 2,kl
!!$                   istart = -2!3
!!$                   jstart = -2!3
!!$                   kstart = -2!3
!!$                   iend = 2
!!$                   jend = 2
!!$                   kend = 2
!!$                   do iiii = istart,iend
!!$                      do jjjj = jstart,jend
!!$                         do kkkk = kstart,kend
!!$                            do n = 1,nw
!!$                               do sps2 = 1,nTimeIntervalsSpectral
!!$                                  call setPointersAdj(nnnn,1,sps2)
!!$                                  i = iii+iiii
!!$                                  j = jjj+jjjj
!!$                                  k = kkk+kkkk
!!$                                 
!!$                                  write(unitx,11) i,j,k,n,nnnn,sps,sps2,iii,jjj,kkk,w(i,j,k,n)
!!$11                                format(1x,'wadj',10I8,f20.12)
!!$                                  call setPointersAdj(nnnn,1,sps)
!!$                               end do
!!$                            enddo
!!$                         enddo
!!$                      enddo
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          enddo
!!$          call setPointersAdj(nn,1,sps)
!!$       enddo
!!$    end do      
       ! Exchange the solution. Either whalo1 or whalo2
       ! must be called.
       if( secondHalo ) then
          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
               .true., .true.)
       else
          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
               .true., .true.)
       endif

       print *, "Calling Residual ="
       call timeStep(.false.)

       !print out x
!!$       do sps = 1,nTimeIntervalsSpectral
!!$          do nnnn=1,ndom
!!$             call setPointersAdj(nnnn,1,sps)
!!$             do iii = 2,il
!!$                do jjj = 2,jl
!!$                   do kkk = 2,kl
!!$                      istart = -2
!!$                      jstart = -2
!!$                      kstart = -2
!!$                      iend = 2
!!$                      jend = 2
!!$                      kend = 2
!!$                      do iiii = istart,iend
!!$                         do jjjj = jstart,jend
!!$                            do kkkk = kstart,kend
!!$                               do n = 1,nw
!!$                                  do sps2 = 1,nTimeIntervalsSpectral
!!$                                     call setPointersAdj(nnnn,1,sps2)
!!$                                     i = iii+iiii
!!$                                     j = jjj+jjjj
!!$                                     k = kkk+kkkk
!!$                                     write(unitx,11) i,j,k,n,nnnn,sps,sps2,w(i,j,k,n)
!!$11                                   format(1x,'res',7I8,f20.14)
!!$                                     call setPointersAdj(nnnn,1,sps)
!!$                                  end do
!!$                               enddo
!!$                            enddo
!!$                         enddo
!!$                      enddo
!!$                   enddo
!!$                enddo
!!$             enddo
!!$             call setPointersAdj(nn,1,sps)
!!$          enddo
!!$       end do    
  !print out w
!!$       do sps = 1,nTimeIntervalsSpectral
!!$          do nnnn=1,1!ndom
!!$          call setPointersAdj(nnnn,1,sps)
!!$          do iii = 2,il
!!$             do jjj = 2,jl
!!$                do kkk = 2,kl
!!$                   istart = -2!3
!!$                   jstart = -2!3
!!$                   kstart = -2!3
!!$                   iend = 2
!!$                   jend = 2
!!$                   kend = 2
!!$                   do iiii = istart,iend
!!$                      do jjjj = jstart,jend
!!$                         do kkkk = kstart,kend
!!$                            do n = 1,nw
!!$                               do sps2 = 1,nTimeIntervalsSpectral
!!$                                  call setPointersAdj(nnnn,1,sps2)
!!$                                  i = iii+iiii
!!$                                  j = jjj+jjjj
!!$                                  k = kkk+kkkk
!!$                                 
!!$                                  write(unitx,11) i,j,k,n,nnnn,sps,sps2,iii,jjj,kkk,w(i,j,k,n)
!!$11                                format(1x,'wadj',10I8,f20.12)
!!$                                  call setPointersAdj(nnnn,1,sps)
!!$                               end do
!!$                            enddo
!!$                         enddo
!!$                      enddo
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          enddo
!!$          call setPointersAdj(nn,1,sps)
!!$       enddo
!!$    end do      
       if( turbCoupled ) then
          call initres(nt1MG, nMGVar)
          call turbResidual
       endif
       
       call initres(1_intType, nwf) 
!!$     do sps = 1,nTimeIntervalsSpectral
!!$        do nnn =1,ndom
!!$           call setPointers(nnn,level,sps)
!!$           dw(:,:,:,:) = 0.0
!!$        enddo
!!$     end do


!!$       !print dw 
!!$       do sps = 1,nTimeIntervalsSpectral
!!$          do nnn =1,ndom
!!$             call setPointers(nnn,level,sps)
!!$             do i = 2,il
!!$                do j = 2,jl
!!$                   do k = 2,kl
!!$                      do n = 1,nw
!!$                         write(unitRes,10) i,j,k,n,nnn,sps,dw(i,j,k,n)
!!$10                       format(1x,'res',6I8,f20.14)
!!$                      enddo
!!$                   enddo
!!$                end do
!!$             end do
!!$          end do
       call residual
       

       !print dw 
       do sps = 1,nTimeIntervalsSpectral
          do nnn =1,ndom
             call setPointers(nnn,level,sps)
             do i = 2,il
                do j = 2,jl
                   do k = 2,kl
                      do n = 1,nw
                         write(unitRes,10) i,j,k,n,nnn,sps,dw(i,j,k,n)
10                       format(1x,'res',6I8,f20.14)
                      enddo
                   enddo
                end do
             end do
          end do
          
          do nnn=1,ndom
          !print *,'domains reset',nnn
             call setPointersAdj(nnn,1,sps)
             w = wtmp
             p = ptmp
             dw = dwtmp
          end do
       end do
 
       print *,'nn',nn! = 1
       call setPointers(nn,level,ssps)
       
       x(ii,jj,kk,:) = xref
       !reset initial conditions
       call xhalo(level)
       call metric(level)
       call checkSymmetry(level)
    
       close(unitx)

       call cpu_time(time(1))
    do sps = 1,nTimeintervalsSpectral
       print *,'nn',nn! = 1
       call setPointers(nn,level,ssps)
       x(ii,jj,kk,:) = xref+deltax

!***********************************
!Now compute residuals using ADjoint routines
!***********************************

      do nnn = 1,1!ndom
         call setPointersAdj(nnn ,level,sps)
         !print *,'in AD loop',nnn 
         !write(unitxAD,*) 'block Num',nnn
         do icell= 2, il
            do jcell= 2, jl
               do kcell= 2, kl
                  
                  !print *,'index',i,j,k
                  call  copyADjointStencil(wAdj, xAdj,xBlockCornerAdj,alphaAdj,&
                          betaAdj,MachAdj,machCoefAdj,machGridAdj,iCell, jCell, kCell,&
                          nnn,level,sps,&
                          prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)


!!$                  !print out xAdj
!!$                  istart = -3
!!$                  jstart = -3
!!$                  kstart = -3
!!$                  iend = 2
!!$                  jend = 2
!!$                  kend = 2
!!$                  if(icell==2) istart=-2
!!$                  if(jcell==2) jstart=-2
!!$                  if(kcell==2) kstart=-2
!!$                  if(icell==il) iend=1
!!$                  if(jcell==jl) jend=1
!!$                  if(kcell==kl) kend=1
!!$                  do iiii = istart,iend
!!$                     do jjjj = jstart,jend
!!$                        do kkkk = kstart,kend
!!$                           do n = 1,3 
!!$                              do sps2 = 1,nTimeIntervalsSpectral
!!$                                 i = icell+iiii
!!$                                 j = jcell+jjjj
!!$                                 k = kcell+kkkk
!!$                                 write(unitxAD,11) i,j,k,n,nnn,sps,sps2,xAdj(iiii,jjjj,kkkk,n,sps2)
!!$                              end do
!!$                           enddo
!!$                        enddo
!!$                     enddo
!!$                  enddo
        
                  
                  ! Compute the total residual.
                  ! This includes inviscid and viscous fluxes, artificial
                  ! dissipation, and boundary conditions.                   
                  !print *,'Calling compute ADjoint' 

                  call computeRAdjoint(wAdj,xAdj,xBlockCornerAdj,dwAdj,alphaAdj,&
                          betaAdj,MachAdj, &
                          MachCoefAdj,machGridAdj,iCell, jCell,  kCell, &
                          nnn,level,sps, correctForK,secondHalo,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                  !print out xAdj
                  istart = -3
                  jstart = -3
                  kstart = -3
                  iend = 2
                  jend = 2
                  kend = 2
                  if(icell==2) istart=-2
                  if(jcell==2) jstart=-2
                  if(kcell==2) kstart=-2
                  if(icell==il) iend=1
                  if(jcell==jl) jend=1
                  if(kcell==kl) kend=1
                  do iiii = istart,iend
                     do jjjj = jstart,jend
                        do kkkk = kstart,kend
                           do n = 1,3
                              do sps2 = 1,nTimeIntervalsSpectral
                                 i = icell+iiii
                                 j = jcell+jjjj
                                 k = kcell+kkkk
                                 write(unitxAD,11) i,j,k,n,nnn,sps,sps2,xAdj(iiii,jjjj,kkkk,n,sps2)
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
!!$                 !print out wAdj
!!$                  istart = -2
!!$                  jstart = -2
!!$                  kstart = -2
!!$                  iend = 2
!!$                  jend = 2
!!$                  kend = 2
!!$                  do iiii = istart,iend
!!$                     do jjjj = jstart,jend
!!$                        do kkkk = kstart,kend
!!$                           !do n = 1,1!nw
!!$                              do n = 1,nw 
!!$                              do sps2 = 1,nTimeIntervalsSpectral
!!$                                 i = icell+iiii
!!$                                 j = jcell+jjjj
!!$                                 k = kcell+kkkk
!!$                                 write(unitxAD,11) i,j,k,n,nnn,sps,sps2,icell,jcell,kcell,wAdj(iiii,jjjj,kkkk,n,sps2) 
!!$                                 !write(unitxAD,11) i,j,k,n,nnn,sps,sps2,pAdj(iiii,jjjj,kkkk,sps2)
!!$                              enddo
!!$                           enddo
!!$                        enddo
!!$                     enddo
!!$                  enddo
        
                  !----------------------------------------------
                  
                  do n = 1,nw
                     !print *,'dwadj2',dwadj,'single',dwadj(n,sps),n,sps,icell,jcell,kcell
                     write(unitResAD,10) icell,jcell,kcell,n,nnn,sps,dwAdj(n,sps)
!10                    format(1x,'drdx',5I8,f18.10) 
                  end do
                  !----------------------------------------------
                  
               end do
            end do
         end do
      end do
      call setPointersAdj(nn ,level,ssps)
      x(ii,jj,kk,:) = xref
   end do
      call cpu_time(time(2))
 !     timeRes = (time(2)-time(1))/(nx*ny*nz)
 !     print *, "Average time for each cell residual calculation =", timeRes
   
      close(unitxAD)
      close(unitResAD)
      close(unitRes)
 


    end subroutine verifyResiduals
