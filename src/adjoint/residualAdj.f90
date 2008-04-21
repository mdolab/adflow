!
!      ******************************************************************
!      *                                                                *
!      * File:          residual.f90                                    *
!      * Author:        Edwin van der Weide, Steve Repsher (blanking)   *
!      * Starting date: 04-21-2008                                      *
!      * Last modified: 04-21-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine residualAdj(wAdj,siAdj,sjAdj,skAdj,volAdj,normAdj,&
                              dwAdj, iCell, jCell, kCell,  &  
                              correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * residual computes the residual of the mean flow equations on   *
!      * the current MG level.                                          *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use flowVarRefState
       use inputIteration
       use inputDiscretization
       use inputTimeSpectral
       use iteration
       implicit none
       
!       Subroutine Variables
       integer(kind=intType) :: iCell, jCell, kCell
       real(kind=realType), dimension(-2:2,-2:2,-2:2,3), intent(inout) :: siAdj, sjAdj, skAdj
       real(kind=realType), dimension(0:0,0:0,0:0), intent(inout) :: volAdj
       real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(inout) :: normAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
                                                 intent(inout) :: wAdj
       real(kind=realType), dimension(nw), intent(inout) :: dwAdj
       !integer(kind=intType), intent(in) :: discr
       logical, intent(in) :: correctForK

       real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdj
        real(kind=realType), dimension(nw) :: fwAdj
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, discr
       integer(kind=intType) :: i, j, k, l

       logical :: fineGrid
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

!   Come back to this later....
!!$       ! Add the source terms from the level 0 cooling model.
!!$
!!$       call level0CoolingModel

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

!moved outside...
!!$
!!$       ! Loop over the number of spectral solutions and local
!!$       ! number of blocks.
!!$
!!$       spectralLoop: do sps=1,nTimeIntervalsSpectral
!!$         domainLoop: do nn=1,nDom
!!$
!!$           ! Set the pointers to this block and compute the central
!!$           ! inviscid flux.
!!$
!!$           call setPointers(nn, currentLevel, sps)

           call inviscidCentralFluxAdj(wAdj,  pAdj,  dwAdj,         &
                                         siAdj, sjAdj, skAdj, volAdj, &
                                         iCell, jCell, kCell)

           ! Compute the artificial dissipation fluxes.
           ! This depends on the parameter discr.

           select case (discr)

!!$             case (dissScalar) ! Standard scalar dissipation scheme.

!!$               if( fineGrid ) then
!!$                 call inviscidDissFluxScalarAdj()
!!$               else
!!$                  call terminate("residualAdj", &
!!$                        "ADjoint does not function on coarse grid level")
!!$                  !call inviscidDissFluxScalarCoarse
!!$               endif

             !===========================================================

!!$             case (dissMatrix) ! Matrix dissipation scheme.
!!$
!!$               if( fineGrid ) then
!!$                 call inviscidDissFluxMatrixAdj()
!!$               else
!!$                 call terminate("residualAdj", &
!!$                        "ADjoint does not function on coarse grid level")
!!$                 !call inviscidDissFluxMatrixCoarse
!!$               endif

             !===========================================================

!!$             case (dissCusp) ! Cusp dissipation scheme.
!!$
!!$               if( fineGrid ) then
!!$                 call inviscidDissFluxCuspAdj()
!!$               else
!!$                 call terminate("residualAdj", &
!!$                        "ADjoint does not function on coarse grid level")
!!$                 !call inviscidDissFluxCuspCoarse
!!$               endif

             !===========================================================

             case (upwind) ! Dissipation via an upwind scheme.

               call inviscidUpwindFluxAdj(wAdj,  pAdj,  dwAdj, &
                                        siAdj, sjAdj, skAdj, &
                                        iCell, jCell, kCell)

           end select

!come back to later
!!$           ! Compute the viscous flux in case of a viscous computation.
!!$
!!$           if( viscous ) call viscousFlux

           ! Add the dissipative and possibly viscous fluxes to the
           ! Euler fluxes. Loop over the owned cells and add fw to dw.
           ! Also multiply by iblank so that no updates occur in holes
           ! or on the overset boundary.

!!$           do l=1,nwf
!!$             do k=2,kl
!!$               do j=2,jl
!!$                 do i=2,il
!!$                   dw(i,j,k,l) = (dw(i,j,k,l) + fw(i,j,k,l)) &
!!$                               * real(iblank(i,j,k), realType)
!!$                 enddo
!!$               enddo
!!$             enddo
!!$           enddo


           do l=1,nwf
              dw(i,j,k,l) = (dwAdj(l) + fwAdj(l)) &
                               * real(iblank(iCell,jCell,kCell), realType)
           enddo


         end subroutine residualAdj
