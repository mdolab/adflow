!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdw.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-05-2008                                      *
!     * Last modified: 05-06-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdx(level)
!
!     ******************************************************************
!     *                                                                *
!     *  This subroutine computes the values for dRdw and compares     *
!     *  the Tapenade result to the finite-difference results.         *
!     *                                                                *
!     ******************************************************************


      use blockPointers ! block (nDoms,flowDoms), globalCell
      use flowvarrefstate
      use communication
      use iteration     ! groundLevel
      use inputTimeSpectral ! spaceDiscr
      use inputIO
     

      implicit none

!
!     Subroutine arguments
      integer(kind=intType), intent(in) :: level


!
!     Local variables 
!
      integer(kind=intType) :: i, j, k, n
      integer(kind=intType) :: iCell,jCell,kCell
      real(kind=realType), dimension(nw) :: dwL2
      real(kind=realType), dimension(nx, ny, nz, nw) :: dwerr
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes
!      real(kind=realType), dimension(4) :: time
      real(kind=realType) ::  timeOri
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdj
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjb
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wFD
      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdj
      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdjb
      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xFD

      character fileName*32, dataName*32
      real(kind=realType), dimension(nw) :: dwAdj,dwAdjb,dwAdjRef
      real(kind=realType), dimension(nw) :: dwAdjP, dwAdjM
      real(kind=realType) :: deltax, xAdjRef
      real(kind=realType) :: timeAdj, timeFD, timeResAdj

      integer(kind=intType), dimension(nDom) :: maxglobalcell
      integer(kind=intType) :: idx, ii, jj, kk, idxstate, idxres, m, l
      integer(kind=intType) :: sps, nTime, max_nTime, nHalo, nn, discr
      real(kind=realType), allocatable, dimension(:,:,:,:) :: dRdxAdj,dRdxFD,dRdxErr

      integer :: ierr
      logical :: fineGrid, correctForK, exchangeTurb,secondHalo


!     ******************************************************************
!     *                                                                *
!     *  Begin execution                                               *
!     *                                                                *
!     ******************************************************************

      if( myID==0 ) write(*,*) "Running verifydRdW..."

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

      !allocate memory for error arrays
      max_nTime=-100
      do i = 1,nDom
         print *,'allocating i'
         maxglobalcell(i) = maxval(flowDoms(i,currentLevel,1)%globalCell(:,:,:))
         print *,'maxglobalcell',maxglobalcell(i) 
         nTime     = nTimeIntervalsSpectral!sections(sectionID)%nTimeInstances
         if(nTime>=max_nTime) max_nTime = nTime
      enddo
      idx = maxval(maxglobalcell(:))
      print *,'allocating',idx,nw*(idx+1),3*(idx+1),ndom,max_nTime
      allocate(dRdxErr(nw*(idx+1),3*(idx+1),ndom,max_nTime), &
               dRdxAdj(nw*(idx+1),3*(idx+1),ndom,max_nTime), &
               dRdxFD(nw*(idx+1),3*(idx+1),ndom,max_nTime))
      print *,' allocated'

!           if(ierr /= 0)                       &
!                call terminate("memory?") 
      
      ! Initialize the temporary arrays.
      dRdxErr = 0
      dRdxAdj = 0
      dRdxFD  = 0


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
!     *  Compute d(dw)/d(x) using Tapenade                             *
!     *                                                                *
!     ******************************************************************

      ! Get the initial AD time.

      call mpi_barrier(SUmb_comm_world, ierr)
      if( myID==0 ) call cpu_time(time(1))
      
      print *,'Entering Domain loop'
      domainLoopAD: do nn=1,nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop: do sps=1,nTimeIntervalsSpectral
  !          print *,'Setting Pointers'
            call setPointersAdj(nn,level,sps)

            ! Loop over location of output (R) cell of residual
            do kCell = 2, kl
               do jCell = 2, jl
                  do iCell = 2, il
 !                    print *,'indices',icell,jcell,kcell
                     ! Copy the state w to the wAdj array in the stencil
                     call copyADjointStencil(wAdj, xAdj, iCell, jCell, kCell)                  
!                     print *,'Stencil Copied'

                     mLoop: do m = 1, nw           ! Loop over output cell residuals (R)
!                        print *,'initializing variables'
                        ! Initialize the seed for the reverse mode
                        dwAdjb(:) = 0.; dwAdjb(m) = 1.
                        dwAdj(:)  = 0.
                        wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                        xAdjb(:,:,:,:)  = 0.  !dR(m)/dx
                        

  !                      print *,'calling reverse mode'
!                        print *,'secondhalo',secondhalo
                        
                        ! Call reverse mode of residual computation
                        call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb,&
                             dwadj, dwadjb, icell, jcell, kcell, nn, sps,&
                             correctfork, secondhalo)
                        
                        ! Store the block Jacobians (by rows).
 !                       print *,'entering storage loop'
                        do ii=-2,2!1,il-1
                           do jj = -2,2!1,jl-1
                              do kk = -2,2!1,kl-1
                                 do l = 1,3
                                    i = iCell + ii
                                    j = jCell + jj
                                    k = kCell + kk
                                    !print *,'secondaryindicies',i,j,k,ii,jj,kk
                                    if(i>zero .and. j>zero .and. k>zero .and. i<il .and. j<jl .and. k<kl)then
                                       idxstate = globalCell(i,j,k)*3+l
                                       idxres   = globalCell(iCell,jCell,kCell)*nw+m
                                       !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                                       dRdxAdj(idxres,idxstate,nn,sps) = xAdjb(ii,jj,kk,l)
                                      ! if ( wAdjb(ii,jj,kk,l).ne. 0.0) then
                                      !    print *,'wAdjb', wAdjb(ii,jj,kk,l),i,j,k,l
                                      ! endif
                                    end if
                                 enddo !l
                              enddo !kk
                           enddo !jj
                        enddo !ii
                        
                     end do mLoop
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop
      end do domainLoopAD
      print *,'AD Completed'
      ! Get new time and compute the elapsed AD time.
!      stop
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
         call cpu_time(time(2))
         timeAdj = time(2)-time(1)
      endif


!     ******************************************************************
!     *                                                                *
!     *  Compute d(dw)/d(x) using central finite differences           *
!     *                                                                *
!     ******************************************************************

      ! Get the initial FD time.
      print *,'starting FD calculation'
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) call cpu_time(time(3))


      deltax = 1.d-6
      xFD(:,:,:,:) = 0. 
      call cpu_time(time(3))

      print *,'entering FD loop'
      domainResidualLoopFDorig: do nn=1,nDom         

         ! Determine the number of time instances for this block
        
         nTime     = nTimeIntervalsSpectral

         ! Loop over the number of time instances for this block.

         spectralLoop2: do sps=1,nTime
!            print *,'setting Pointers'
            call setPointersAdj(nn,level,sps)


               ! Loop over location of output (R) cell of residual
               do kCell = 2, kl
                  !print *, "iRes, kCell =", iRes, kCell
                  do jCell = 2, jl
                     do iCell = 2, il
                        
                        do m = 1, nw           ! Loop over output cell residuals (R)
!                           print *,'indices',icell,jcell,kcell,m
                        ! Copy the state in the stencil
                        ! actually no need to call again, but ...
                           call copyADjointStencil(wAdj,xAdj, iCell, jCell, kCell)
!                           print *,'stencil copied'
                           ! Loop over the inputs (w) 
                           !print *,'secondhalo',secondhalo
                        
                           do ii = -2, 2
                              do jj = -2, 2
                                 do kk = -2, 2
                                    do n = 1, 3
                                       i = iCell+ii
                                       j = jCell+jj
                                       k = kCell+kk
 !                                      print *,'secondary FD indices',i,j,k,n,ii,jj,kk
                                       if (ii==0 .or. jj==0 .or. kk==0) then
                                          xAdjRef = xAdj(ii,jj,kk,n)
                                          xAdj(ii,jj,kk,n) = xAdjRef + deltax
                                          call computeRAdjoint(wAdj,xAdj,dwAdjP,   &
                                               iCell, jCell,  kCell, &
                                               nn,sps, correctForK,secondHalo)
                                          
                                          xAdj(ii,jj,kk,n) = xAdjRef - deltax
                                          call computeRAdjoint(wAdj,xAdj,dwAdjM,   &
                                               iCell, jCell,  kCell, &
                                               nn,sps, correctForK,secondHalo)

                                          
                                          xFD(ii,jj,kk,n) = (dwAdjP(m)-dwAdjM(m))/(2.0*deltax)
                                          xAdj(ii,jj,kk,n) = xAdjRef

                                          if(i>zero .and. j>zero .and. k>zero .and. i<il .and. j<jl .and. k<kl)then
!                                             print *,'global FD indicies',idxstate,idxres
                                             idxstate = globalCell(i,j,k)*3+n
                                             idxres   = globalCell(iCell,jCell,kCell)*nw+m
                                             
                                             dRdxFD(idxres,idxstate,nn,sps) = xFD(ii,jj,kk,n)
                                          end if
                                       end if
                                    end do ! n
                                 end do ! ii
                              end do !jj
                           end do !kk
                        
                        !print 111, iRes, iCell-1, jCell-1, kCell-1, dRdw(iCell-1, jCell-1, kCell-1, iRes), dRdwFD(iCell-1, jCell-1, kCell-1, iRes)
                     end do !iRes
                  end do !iCell
               end do !jCell
            end do !kCell
         end do spectralLoop2
      end do domainResidualLoopFDorig
      
      call cpu_time(time(4))
      timeFD = time(4)-time(3)

!_______________________________________________________
!
!     Compute the errors in dR/dx
!_______________________________________________________

      print *, "Computing the the error..."

    do m=1,nw
      do iCell=2,il
        do jCell=2,jl
          do kCell=2,kl
            do ii=-2,2
              do jj=-2,2
                do kk=-2,2
                 do n=1,3
                  i=iCell+ii; j=jCell+jj; k=kCell+kk
                     if(i>zero .and. j>zero .and. k>zero .and. i<il .and. j<jl .and. k<kl)then
                        idxstate = globalCell(i,j,k)*3+l
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        dRdxErr(idxres, idxstate, 1, 1) = dRdxAdj(idxres, idxstate, 1, 1) - dRdxFD(idxres, idxstate, 1, 1)
                         if(dRdxFD(idxres,idxstate,1,1).ne.0) then
                         write(*,*) iCell,jCell,kCell,ii,jj,kk, &
                         dRdxAdj(idxres,idxstate,1,1), dRdxFD(idxres,idxstate,1,1), dRdxErr(idxres,idxstate,1,1)
                         end if
                      end if
                end do
              end do
            end do
          end do
        end do
      end do
     end do
    end do
stop

!      dRdwErrRel_1(:, :, :, :) = abs(dRdwErr_q(:, :, :, :)/dRdw_q(:,:,:,:))
!      print *
      print *, "max, min of dRdwErr =", maxval(dRdxErr(:,:,1,1)), minval(dRdxErr(:,:,1,1))
!      print *, "max, min of dRdwErrRel =", maxval(dRdwErrRel_q(:,:,:,1)), minval(dRdwErrRel_q(:,:,:,1))
!      print *


 111  format(4I4, 3ES22.14)



    end subroutine verifydRdx
