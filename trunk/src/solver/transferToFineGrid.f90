!
!      ******************************************************************
!      *                                                                *
!      * File:          transferToFineGrid.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-13-2003                                      *
!      * Last modified: 08-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine transferToFineGrid(corrections)
!
!      ******************************************************************
!      *                                                                *
!      * transferToFineGrid interpolates either the corrections or the  *
!      * solution to the next finer grid level. A standard trilinear    *
!      * interpolation is used.                                         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       use inputIteration
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: corrections
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, i, j, k, l
       integer(kind=intType) :: ii, jj, kk, ii1, jj1, kk1
       integer(kind=intType) :: coarseLevel, nVarInt

       real(kind=realType) :: fact
       real(kind=realType), dimension(:,:,:,:), pointer :: ww, ww1, res
       real(kind=realType), dimension(:,:,:),   pointer :: pp, pp1

       logical :: secondHalo, correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the coarse grid level in coarseLevel.

       coarseLevel = currentLevel +1

       ! Set the number of variables for which either the corrections
       ! or the solution must be interpolated. If the corrections must
       ! be interpolated, this value is set to the number of variables
       ! to which multigrid must be applied, otherwise all conservative
       ! variables are interpolated.

       nVarInt = nw
       if( corrections ) nVarInt = nMGVar

       ! Set the value of secondHalo, depending on the situation.
       ! In the full MG (currentLevel < groundLevel) the second halo is
       ! always set; otherwise only on the finest mesh in the current mg
       ! cycle.

       if(currentLevel <= groundLevel) then
         secondHalo = .true.
       else
         secondHalo = .false.
       endif

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

       ! Set fact to either 0.0 or 1.0, depending whether neumann or
       ! dirichlet boundary conditions should be used for the boundary
       ! halo's when interpolating.
 
       fact = one
       if(corrections .and. mgBoundCorr == bcDirichlet0) fact = zero

       ! Loop over the number of spectral solutions and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers to the fine level, i.e. currentLevel.
           ! Also set the pointers for ww, pp, ww1 and pp1 to the
           ! coarse grid values.

           call setPointers(nn, currentLevel,sps)
           ww  => flowDoms(nn,coarseLevel,sps)%w
           pp  => flowDoms(nn,coarseLevel,sps)%p
           ww1 => flowDoms(nn,coarseLevel,sps)%w1
           pp1 => flowDoms(nn,coarseLevel,sps)%p1

           ! Store the correction in ww if the corrections must be
           ! interpolated. The 1st level halo's are included, because
           ! these values are needed for the interpolation. Note that
           ! flowDoms(nn,coarseLevel)%ie, etc. must be used, because
           ! ie is equal to the fine grid value.

           testCorrections1: if( corrections ) then

             ! Flow field variables. Have res point to dw and compute the
             ! corrections. Note that the pressure correction must be
             ! stored instead of the total energy.

             res => dw

             do k=1,flowDoms(nn,coarseLevel,sps)%ke
               do j=1,flowDoms(nn,coarseLevel,sps)%je
                 do i=1,flowDoms(nn,coarseLevel,sps)%ie
                   ww(i,j,k,irho)  = ww(i,j,k,irho) - ww1(i,j,k,irho)
                   ww(i,j,k,ivx)   = ww(i,j,k,ivx)  - ww1(i,j,k,ivx)
                   ww(i,j,k,ivy)   = ww(i,j,k,ivy)  - ww1(i,j,k,ivy)
                   ww(i,j,k,ivz)   = ww(i,j,k,ivz)  - ww1(i,j,k,ivz)
                   ww(i,j,k,irhoE) = pp(i,j,k)      - pp1(i,j,k)
                 enddo
               enddo
             enddo

             ! The possible turbulent variables.

             do l=nt1,nVarInt
               do k=1,flowDoms(nn,coarseLevel,sps)%ke
                 do j=1,flowDoms(nn,coarseLevel,sps)%je
                   do i=1,flowDoms(nn,coarseLevel,sps)%ie
                     ww(i,j,k,l) = ww(i,j,k,l) - ww1(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo

           else testCorrections1

             ! The solution must be interpolated. Have res point to w and
             ! store the pressure instead of the total energy.

             res => w

             do k=1,flowDoms(nn,coarseLevel,sps)%ke
               do j=1,flowDoms(nn,coarseLevel,sps)%je
                 do i=1,flowDoms(nn,coarseLevel,sps)%ie
                   ww(i,j,k,irhoE) = pp(i,j,k)
                 enddo
               enddo
             enddo

           endif testCorrections1

           ! Set the values of the coarse grid boundary halo cells.

           call setCorrectionsCoarseHalos(sps, nn, coarseLevel, &
                                          fact, nVarInt)

           ! Loop over the owned fine grid cells and determine res
           ! by trilinear interpolation. Note that the coarse grid cell
           ! ii (and jj and kk) are such that they are closest to the
           ! fine grid cell center and ii1 is further away. This means
           ! that in 1d ii get the weight 3/4 and ii1 1/4.

           do k=2,kl

             ! Determine the coarse grid cells kk and kk1.

             kk  = mgKCoarse(k,1)
             kk1 = mgKCoarse(k,2)

             do j=2,jl

               ! Determine the coarse grid cells jj and jj1.

               jj  = mgJCoarse(j,1)
               jj1 = mgJCoarse(j,2)

               do i=2,il

                 ! Determine the coarse grid cells ii and ii1.

                 ii  = mgICoarse(i,1)
                 ii1 = mgICoarse(i,2)

                 ! Loop over the number of variables and interpolate them.
                 ! The weights involved are 27/64, 9/64, 3/64 and 1/64.
                 ! For computational efficiency their (exact) decimal
                 ! counterparts are used in the loop below, which are
                 ! 0.421875, 0.140625, 0.046875 and 0.015625 respectively.

                 do l=1,nVarInt
                   res(i,j,k,l) = 0.421875_realType* ww(ii,jj,kk,l)    &
                                + 0.140625_realType*(ww(ii1,jj,kk,l)   &
                                +                    ww(ii,jj1,kk,l)   &
                                +                    ww(ii,jj,kk1,l))  &
                                + 0.046875_realType*(ww(ii1,jj1,kk,l)  &
                                +                    ww(ii1,jj,kk1,l)  &
                                +                    ww(ii,jj1,kk1,l)) &
                                + 0.015625_realType* ww(ii1,jj1,kk1,l)
                 enddo

               enddo
             enddo
           enddo

           ! Possibility to do smoothing on the corrections, if desired.


           ! Compute the new state vector on the fine mesh in case the
           ! corrections have just been interpolated.

           testCorrections2: if( corrections ) then

             ! Flow field variables. Again the pressure is updated
             ! and not the total energy. Make sure that the pressure and
             ! density do not become negative.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   w(i,j,k,irho) = w(i,j,k,irho) + dw(i,j,k,irho)
                   w(i,j,k,ivx)  = w(i,j,k,ivx)  + dw(i,j,k,ivx)
                   w(i,j,k,ivy)  = w(i,j,k,ivy)  + dw(i,j,k,ivy)
                   w(i,j,k,ivz)  = w(i,j,k,ivz)  + dw(i,j,k,ivz)
                   p(i,j,k)      = p(i,j,k)      + dw(i,j,k,irhoE)

                   w(i,j,k,irho) = max(w(i,j,k,irho), &
                                       1.e-4_realType*rhoInf)
                   p(i,j,k)      = max(p(i,j,k), &
                                       1.e-4_realType*pInfCorr)
                 enddo
               enddo
             enddo

             ! The possible turbulent variables.

             do l=nt1,nVarInt
               do k=2,kl
                 do j=2,jl
                   do i=2,il
                     w(i,j,k,l) = w(i,j,k,l) + dw(i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo

           else testCorrections2

             ! The solution must be interpolated. At the moment the
             ! pressure is stored at the place of the total energy.
             ! Copy this value in the pressure array.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   p(i,j,k) = w(i,j,k,irhoE)
                 enddo
               enddo
             enddo

           endif testCorrections2

           ! Compute the total energy for the owned cells of this block.
           ! If the solution must be interpolated, extrapolate the
           ! solution in the halo's.

           call computeEtot(2_intType,il, 2_intType,jl, &
                            2_intType,kl, correctForK)
           if(.not. corrections) call extrapolateSolution

           ! Compute the laminar viscosity and eddy viscosity for the
           ! owned cells of this block. If the solution must be
           ! interpolated, extrapolate the viscosities in the halo's.

           call computeLamViscosity
           call computeEddyViscosity
           if(.not. corrections) call extrapolateViscosities

         enddo domains
       enddo spectralLoop

       ! Exchange the pressure if the pressure must be exchanged early.
       ! Only the first halo's are needed, thus whalo1 is called.
       ! On the finest mesh only.

       if(exchangePressureEarly .and. currentLevel <= groundLevel) &
         call whalo1(currentLevel, 1_intType, 0_intType, &
                     .true., .false., .false.)

       ! Apply all boundary conditions to all blocks on this level.
       ! In case of a full mg mode, and a segegated turbulent solver,
       ! first call the turbulent boundary conditions, such that the
       ! turbulent kinetic energy is properly initialized in the halo's.

       if(turbSegregated .and. (.not. corrections)) &
         call applyAllTurbBC(secondHalo)

       ! Apply all boundary conditions of the mean flow.

       call applyAllBC(secondHalo)

       ! If case this routine is called in full mg mode call the mean
       ! flow boundary conditions again such that the normal momentum
       ! boundary condition is treated correctly.

       if(.not. corrections) call applyAllBC(secondHalo)

       ! Exchange the solution. Either whalo1 or whalo2
       ! must be called.

       if( secondHalo ) then
         call whalo2(currentLevel, 1_intType, nVarInt, .true., &
                     .true., .true.)
       else
         call whalo1(currentLevel, 1_intType, nVarInt, .true., &
                     .true., .true.)
       endif

       ! For full multigrid mode the bleeds must be determined, the
       ! boundary conditions must be applied one more time and the
       ! solution must be exchanged again.

       if(.not. corrections) then
         call BCDataMassBleedOutflow(.true., .true.)
         call applyAllBC(secondHalo)

         if( secondHalo ) then
           call whalo2(currentLevel, 1_intType, nVarInt, .true., &
                       .true., .true.)
         else
           call whalo1(currentLevel, 1_intType, nVarInt, .true., &
                       .true., .true.)
         endif
       endif

       end subroutine transferToFineGrid

!      ==================================================================

       subroutine extrapolateSolution
!
!      ******************************************************************
!      *                                                                *
!      * extrapolateSolution sets the solution of the cell halos by a   *
!      * constant extrapolation. This routine is called after the       *
!      * solution has been interpolated to the next finer grid and this *
!      * routine makes sure that the halo's are initialized.            *
!      * Only the block to which the pointers in blockPointers          *
!      * currently point is treated.                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Constant extrapolation in i-direction.

       do k=2,kl
         do j=2,jl

           do l=1,nw
             w(0,j,k,l)  = w(2,j,k,l)
             w(1,j,k,l)  = w(2,j,k,l)
             w(ie,j,k,l) = w(il,j,k,l)
             w(ib,j,k,l) = w(il,j,k,l)
           enddo

           p(0,j,k)  = p(2,j,k)
           p(1,j,k)  = p(2,j,k)
           p(ie,j,k) = p(il,j,k)
           p(ib,j,k) = p(il,j,k)

         enddo
       enddo

       ! Constant extrapolation in the j-direction. Take the just
       ! interpolated values in i-direction into account.

       do k=2,kl
         do i=0,ib

           do l=1,nw
             w(i,0,k,l)  = w(i,2,k,l)
             w(i,1,k,l)  = w(i,2,k,l)
             w(i,je,k,l) = w(i,jl,k,l)
             w(i,jb,k,l) = w(i,jl,k,l)
           enddo

           p(i,0,k)  = p(i,2,k)
           p(i,1,k)  = p(i,2,k)
           p(i,je,k) = p(i,jl,k)
           p(i,jb,k) = p(i,jl,k)

         enddo
       enddo

       ! Constant extrapolation in the k-direction. Take the just
       ! interpolated values in i- and j-direction into account.

       do j=0,jb
         do i=0,ib

           do l=1,nw
             w(i,j,0,l)  = w(i,j,2,l)
             w(i,j,1,l)  = w(i,j,2,l)
             w(i,j,ke,l) = w(i,j,kl,l)
             w(i,j,kb,l) = w(i,j,kl,l)
           enddo

           p(i,j,0)  = p(i,j,2)
           p(i,j,1)  = p(i,j,2)
           p(i,j,ke) = p(i,j,kl)
           p(i,j,kb) = p(i,j,kl)

         enddo
       enddo

       end subroutine extrapolateSolution

!      ==================================================================

       subroutine extrapolateViscosities
!
!      ******************************************************************
!      *                                                                *
!      * extrapolateViscosities sets the laminar and eddy viscosities   *
!      * of the cell halos by a constant extrapolation. This routine is *
!      * called after the solution has been interpolated to the next    *
!      * finer grid and this routine makes sure that the halo's are     *
!      * initialized.                                                   *
!      * Only the block to which the pointers in blockPointers          *
!      * currently point is treated and only for a viscous problem.     *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is not a viscous problem.

       if(.not. viscous) return

       ! Constant extrapolation in i-direction.

       do k=2,kl
         do j=2,jl

           rlv(0,j,k)  = rlv(2,j,k)
           rlv(1,j,k)  = rlv(2,j,k)
           rlv(ie,j,k) = rlv(il,j,k)
           rlv(ib,j,k) = rlv(il,j,k)

           if( eddyModel ) then
             rev(0,j,k)  = rev(2,j,k)
             rev(1,j,k)  = rev(2,j,k)
             rev(ie,j,k) = rev(il,j,k)
             rev(ib,j,k) = rev(il,j,k)
           endif

         enddo
       enddo

       ! Constant extrapolation in the j-direction. Take the just
       ! interpolated values in i-direction into account.

       do k=2,kl
         do i=0,ib

           rlv(i,0,k)  = rlv(i,2,k)
           rlv(i,1,k)  = rlv(i,2,k)
           rlv(i,je,k) = rlv(i,jl,k)
           rlv(i,jb,k) = rlv(i,jl,k)

           if( eddyModel ) then
             rev(i,0,k)  = rev(i,2,k)
             rev(i,1,k)  = rev(i,2,k)
             rev(i,je,k) = rev(i,jl,k)
             rev(i,jb,k) = rev(i,jl,k)
           endif

         enddo
       enddo

       ! Constant extrapolation in the k-direction. Take the just
       ! interpolated values in i- and j-direction into account.

       do j=0,jb
         do i=0,ib

           rlv(i,j,0)  = rlv(i,j,2)
           rlv(i,j,1)  = rlv(i,j,2)
           rlv(i,j,ke) = rlv(i,j,kl)
           rlv(i,j,kb) = rlv(i,j,kl)

           if( eddyModel ) then
             rev(i,j,0)  = rev(i,j,2)
             rev(i,j,1)  = rev(i,j,2)
             rev(i,j,ke) = rev(i,j,kl)
             rev(i,j,kb) = rev(i,j,kl)
           endif

         enddo
       enddo

       end subroutine extrapolateViscosities
