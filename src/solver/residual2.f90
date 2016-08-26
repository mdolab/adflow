       subroutine residual2
!
!       residual computes the residual of the mean flow equations on   
!       the current MG level.                                          
!
       use blockPointers
       use cgnsGrid
       use flowVarRefState
       use inputIteration
       use inputDiscretization
       use inputTimeSpectral
       use iteration
       use communication
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, discr
       integer(kind=intType) :: i, j, k, l
       logical :: fineGrid

       ! Add the source terms from the level 0 cooling model.

       call level0CoolingModel

       ! Set the value of rFil, which controls the fraction of the old
       ! dissipation residual to be used. This is only for the runge-kutta
       ! schemes; for other smoothers rFil is simply set to 1.0.
       ! Note the index rkStage+1 for cdisRK. The reason is that the
       ! residual computation is performed before rkStage is incremented.

       if(smoother == RungeKutta) then
         rFil = cdisRK(rkStage+1)
       else
         rFil = one
       endif

       ! Initialize the local arrays to monitor the massflows to zero.

       massFlowFamilyInv  = zero
       massFlowFamilyDiss = zero

       ! Set the value of the discretization, depending on the grid level,
       ! and the logical fineGrid, which indicates whether or not this
       ! is the finest grid level of the current mg cycle.

       discr = spaceDiscrCoarse
       if(currentLevel == 1) discr = spaceDiscr

       fineGrid = .false.
       if(currentLevel == groundLevel) fineGrid = .true.

       ! Loop over the number of spectral solutions and local
       ! number of blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domainLoop: do nn=1,nDom

           ! Set the pointers to this block and compute the central
           ! inviscid flux.

           call setPointers(nn, currentLevel, sps)
           
           call inviscidCentralFlux
!            do l=1,nw
!               do k=0,kb
!                  do j=0,jb
!                     do i=0,ib
!                        if (myisnan(w(i,j,k,l))) then
!                           print *,'Error in w at:',myid, nbkGlobal,i,j,k,l
!                        end if
!                        if (myisnan(dw(i,j,k,l))) then
!                           print *,'Error in dw at:',myid, nbkGlobal,i,j,k,l
!                        end if
!                        if (myisnan(p(i,j,k))) then
!                           print *,'Error in p at:',myid, nbkGlobal,i,j,k
!                        end if
!                     end do
!                  end do
!               end do
!            end do

           ! Compute the artificial dissipation fluxes.
           ! This depends on the parameter discr.

           select case (discr)

             case (dissScalar) ! Standard scalar dissipation scheme.

               if( fineGrid ) then
                 call inviscidDissFluxScalar
               else
                  call inviscidDissFluxScalarCoarse
               endif

             !===========================================================

             case (dissMatrix) ! Matrix dissipation scheme.

               if( fineGrid ) then
                 call inviscidDissFluxMatrix
               else
                 call inviscidDissFluxMatrixCoarse
               endif

             !===========================================================

             case (dissCusp) ! Cusp dissipation scheme.

               if( fineGrid ) then
                 call inviscidDissFluxCusp
               else
                 call inviscidDissFluxCuspCoarse
               endif

             !===========================================================

             case (upwind) ! Dissipation via an upwind scheme.
!                 print *,'Before upwind check'
!                 do l=1,nw
!                    do k=0,kb
!                       do j=0,jb
!                          do i=0,ib
!                             if (myisnan(w(i,j,k,l))) then
!                                print *,'Error(2) in w at:',myid, nbkGlobal,i,j,k,l
!                             end if
!                             if (myisnan(p(i,j,k))) then
!                                print *,'Error(2) in p at:',myid, nbkGlobal,i,j,k
!                             end if

!                             if (myisnan(dw(i,j,k,l))) then
!                                print *,'Error(2) in dw at:',myid, nbkGlobal,i,j,k
!                             end if


!                          end do
!                       end do
!                    end do
!                 end do


               call inviscidUpwindFlux(fineGrid)
!                print *,'after upwind check'
!                do l=1,nw
!                   do k=0,kb
!                      do j=0,jb
!                         do i=0,ib
!                            if (myisnan(fw(i,j,k,l))) then
!                               print *,'Error(2) in fw at:',myid, nbkGlobal,i,j,k,l
!                            end if
!                            if (myisnan(dw(i,j,k,l))) then
!                               print *,'Error(2) in dw at:',myid, nbkGlobal,i,j,k,l
!                            end if


!                         end do
!                      end do
!                   end do
!                end do
               



           end select

           ! Compute the viscous flux in case of a viscous computation.

           if( viscous ) call viscousFlux

           ! Add the dissipative and possibly viscous fluxes to the
           ! Euler fluxes. Loop over the owned cells and add fw to dw.
           ! Also multiply by iblank so that no updates occur in holes

           do l=1,nwf
             do k=2,kl
               do j=2,jl
                 do i=2,il
                   dw(i,j,k,l) = (dw(i,j,k,l) + fw(i,j,k,l)) &
                               * max(real(iblank(i,j,k), realType), zero)
                 enddo
               enddo
             enddo
           enddo

         enddo domainLoop
       enddo spectralLoop

       end subroutine residual2
