       subroutine transferToCoarseGrid
!
!       transferToCoarseGrid restricts both the solution and the       
!       residual to the next coarser grid level and computes the       
!       residual forcing term on this level.                           
!
       use blockPointers
       use flowVarRefState
       use inputIteration
       use inputTimeSpectral
       use iteration
       use utils, only : setPointers
       use haloExchange, only : whalo1
       use flowUtils, only : computeEtotBlock

       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn, sps, i, j, k, l
       integer(kind=intType) :: ii, jj, kk, ii1, jj1, kk1
       integer(kind=intType) :: fineLevel

       integer(kind=intType), dimension(:,:,:), pointer :: iiblank

       real(kind=realType) :: vola, tmp, weigth, blankFact

       real(kind=realType), dimension(:,:,:,:), pointer :: ww
       real(kind=realType), dimension(:,:,:), pointer :: pp, vvol, rrev

       logical :: correctForK

       ! Compute the residual on the fine grid. It is assumed that the
       ! halo's contain the correct values. The time step is computed,
       ! because this routine also computes the spectral radii for the
       ! artificial dissipation terms.

       rkStage = 0
       call timeStep(.true.)

       call initres(1_intType, nwf)
       call residual

       ! Store the fine grid level and update the current level such
       ! that it corresponds to the coarse grid

       fineLevel    = currentLevel
       currentLevel = currentLevel +1

       ! Set the logical correctForK. 
       correctForK = .false.

       ! Set the value of the blanking factor for the restricted
       ! residual to 1. This will be overwritten below if needed.

       blankFact = one

       ! Loop over the number of spectral solutions and domains.

       spectralLoop1: do sps=1,nTimeIntervalsSpectral
         domains1: do nn=1,nDom

           ! Set the pointers to the coarse block and to the fine grid
           ! solution and volumes. Note that this is not needed for the,
           ! residual, because only the fine grid residual is allocated.

           call setPointers(nn, currentLevel, sps)

           ww      => flowDoms(nn,fineLevel,sps)%w
           pp      => flowDoms(nn,fineLevel,sps)%p
           vvol    => flowDoms(nn,fineLevel,sps)%vol
           rrev    => flowDoms(nn,fineLevel,sps)%rev
           iiblank => flowDoms(nn,fineLevel,sps)%iblank

           ! Restrict the solution and the residual to the coarser
           ! meshes. The solution is done in a volume averaged way. Note
           ! that the sum of the fine grid volumes is used to divide and
           ! not the coarse grid volume; these are not necessarily the
           ! same, especially for the flexible mg used.
 
           do k=2,kl
             kk  = mgKFine(k,1)
             kk1 = mgKFine(k,2)
             do j=2,jl
               jj  = mgJFine(j,1)
               jj1 = mgJFine(j,2)
               do i=2,il
                 ii  = mgIFine(i,1)
                 ii1 = mgIFine(i,2)

                 ! Determine the weight for the restricted residual. This
                 ! weight is less than 1.0 if in at least 1 direction an
                 ! irregular coarsening is used. This weight is not
                 ! applied to the solution, because volume weighting is
                 ! used there.

                 weigth = mgKWeight(k)*mgJWeight(j)*mgIWeight(i)

                 ! Compute the sum of the fine grid volumes.

                 vola = vvol(ii,jj,kk)   + vvol(ii1,jj,kk)   &
                      + vvol(ii,jj1,kk)  + vvol(ii1,jj1,kk)  &
                      + vvol(ii,jj,kk1)  + vvol(ii1,jj,kk1)  &
                      + vvol(ii,jj1,kk1) + vvol(ii1,jj1,kk1)

                 ! Invert the sum of the fine grid volumes.

                 vola = one/vola

                 ! Store the restricted residual in wr, the residual
                 ! forcing term, and the solution in ww.

                 do l=1,nwf
                   wr(i,j,k,l) = (dw(ii, jj,kk, l) + dw(ii, jj1,kk, l)  &
                               +  dw(ii1,jj,kk, l) + dw(ii1,jj1,kk, l)  &
                               +  dw(ii, jj,kk1,l) + dw(ii, jj1,kk1,l)  &
                               +  dw(ii1,jj,kk1,l) + dw(ii1,jj1,kk1,l)) &
                               *  weigth*blankFact
                 enddo

                 ! Restrict the solution.

                 ! Density.

                 w(i,j,k,irho) = (vvol(ii, jj, kk) *ww(ii, jj, kk, irho)  &
                               +  vvol(ii, jj1,kk) *ww(ii, jj1,kk, irho)  &
                               +  vvol(ii1,jj, kk) *ww(ii1,jj, kk, irho)  &
                               +  vvol(ii1,jj1,kk) *ww(ii1,jj1,kk, irho)  &
                               +  vvol(ii, jj, kk1)*ww(ii, jj, kk1,irho)  &
                               +  vvol(ii, jj1,kk1)*ww(ii, jj1,kk1,irho)  &
                               +  vvol(ii1,jj, kk1)*ww(ii1,jj, kk1,irho)  &
                               +  vvol(ii1,jj1,kk1)*ww(ii1,jj1,kk1,irho)) &
                               * vola

                 ! X-velocity.

                 w(i,j,k,ivx) = (vvol(ii, jj, kk) *ww(ii, jj, kk, ivx)  &
                              +  vvol(ii, jj1,kk) *ww(ii, jj1,kk, ivx)  &
                              +  vvol(ii1,jj, kk) *ww(ii1,jj, kk, ivx)  &
                              +  vvol(ii1,jj1,kk) *ww(ii1,jj1,kk, ivx)  &
                              +  vvol(ii, jj, kk1)*ww(ii, jj, kk1,ivx)  &
                              +  vvol(ii, jj1,kk1)*ww(ii, jj1,kk1,ivx)  &
                              +  vvol(ii1,jj, kk1)*ww(ii1,jj, kk1,ivx)  &
                              +  vvol(ii1,jj1,kk1)*ww(ii1,jj1,kk1,ivx)) &
                              * vola

                 ! Y-velocity.

                 w(i,j,k,ivy) = (vvol(ii, jj, kk) *ww(ii, jj, kk, ivy)  &
                              +  vvol(ii, jj1,kk) *ww(ii, jj1,kk, ivy)  &
                              +  vvol(ii1,jj, kk) *ww(ii1,jj, kk, ivy)  &
                              +  vvol(ii1,jj1,kk) *ww(ii1,jj1,kk, ivy)  &
                              +  vvol(ii, jj, kk1)*ww(ii, jj, kk1,ivy)  &
                              +  vvol(ii, jj1,kk1)*ww(ii, jj1,kk1,ivy)  &
                              +  vvol(ii1,jj, kk1)*ww(ii1,jj, kk1,ivy)  &
                              +  vvol(ii1,jj1,kk1)*ww(ii1,jj1,kk1,ivy)) &
                              * vola

                 ! Z-velocity.

                 w(i,j,k,ivz) = (vvol(ii, jj, kk) *ww(ii, jj, kk, ivz)  &
                              +  vvol(ii, jj1,kk) *ww(ii, jj1,kk, ivz)  &
                              +  vvol(ii1,jj, kk) *ww(ii1,jj, kk, ivz)  &
                              +  vvol(ii1,jj1,kk) *ww(ii1,jj1,kk, ivz)  &
                              +  vvol(ii, jj, kk1)*ww(ii, jj, kk1,ivz)  &
                              +  vvol(ii, jj1,kk1)*ww(ii, jj1,kk1,ivz)  &
                              +  vvol(ii1,jj, kk1)*ww(ii1,jj, kk1,ivz)  &
                              +  vvol(ii1,jj1,kk1)*ww(ii1,jj1,kk1,ivz)) &
                              * vola

                 ! Pressure.

                 p(i,j,k) = (vvol(ii, jj, kk) *pp(ii, jj, kk)  &
                          +  vvol(ii, jj1,kk) *pp(ii, jj1,kk)  &
                          +  vvol(ii1,jj, kk) *pp(ii1,jj, kk)  &
                          +  vvol(ii1,jj1,kk) *pp(ii1,jj1,kk)  &
                          +  vvol(ii, jj, kk1)*pp(ii, jj, kk1) &
                          +  vvol(ii, jj1,kk1)*pp(ii, jj1,kk1) &
                          +  vvol(ii1,jj, kk1)*pp(ii1,jj, kk1) &
                          +  vvol(ii1,jj1,kk1)*pp(ii1,jj1,kk1))*vola

                 ! Restrict the eddy viscosity if needed.
                 rev(i,j,k) = (vvol(ii, jj, kk) *rrev(ii, jj, kk)  &
                      +  vvol(ii, jj1,kk) *rrev(ii, jj1,kk)  &
                      +  vvol(ii1,jj, kk) *rrev(ii1,jj, kk)  &
                      +  vvol(ii1,jj1,kk) *rrev(ii1,jj1,kk)  &
                      +  vvol(ii, jj, kk1)*rrev(ii, jj, kk1) &
                      +  vvol(ii, jj1,kk1)*rrev(ii, jj1,kk1) &
                      +  vvol(ii1,jj, kk1)*rrev(ii1,jj, kk1) &
                      +  vvol(ii1,jj1,kk1)*rrev(ii1,jj1,kk1))*vola
               enddo
             enddo
           enddo

           ! Compute the total energy, laminar viscosity and eddy viscosity
           ! for the owned cells of this block.

           call computeEtotBlock(2_intType,il, 2_intType,jl, &
                2_intType,kl, correctForK)
           call computeLamViscosity
           call computeEddyViscosity

           ! Set the values of the 1st layer of corner row halo's to avoid
           ! divisions by zero and uninitialized variables.

           call setCornerRowHalos(nwf)

         enddo domains1
       enddo spectralLoop1

       ! Apply all boundary conditions to all blocks on this level.
       ! No need to exchange the pressure before, because on the coarser
       ! grids a constant pressure boundary condition is used for the
       ! inviscid walls.

       call applyAllBC(.false.)

       ! Exchange the solution. As on the coarse grid only the first
       ! layer of halo's is needed, whalo1 is called.

       call whalo1(currentLevel, 1_intType, nwf, .true., &
                  .true., .true.)

       ! The second part of the residual forcing term is the residual
       ! of the just restricted solution.
       ! First compute the time step.

       rkStage = 0
       call timeStep(.false.)

       ! The second part of the residual forcing term for the mean
       ! flow equations. Furthermore the solution, primitive variables,
       ! is stored in w1 and p1, also of the 1st level halo's. These
       ! may be needed to determine the fine grid corrections.

       spectralLoop3: do sps=1,nTimeIntervalsSpectral
         domains3: do nn=1,nDom

           ! Have the pointers point to this block.

           call setPointers(nn, currentLevel, sps)

           ! Initialize the mean flow residual to zero.

           do k=1,ke
             do j=1,je
               do i=1,ie
                 dw(i,j,k,irho)  = zero
                 dw(i,j,k,imx)   = zero
                 dw(i,j,k,imy)   = zero
                 dw(i,j,k,imz)   = zero
                 dw(i,j,k,irhoE) = zero
               enddo
             enddo
           enddo

           ! Store the restricted solution in w1 and p1.

           ! Flow field variables.

           do k=1,ke
             do j=1,je
               do i=1,ie
                 w1(i,j,k,irho)  = w(i,j,k,irho)
                 w1(i,j,k,ivx)   = w(i,j,k,ivx)
                 w1(i,j,k,ivy)   = w(i,j,k,ivy)
                 w1(i,j,k,ivz)   = w(i,j,k,ivz)
                 w1(i,j,k,irhoE) = w(i,j,k,irhoE)

                 p1(i,j,k)       = p(i,j,k)
               enddo
             enddo
           enddo

         enddo domains3
       enddo spectralLoop3

       ! Compute the mean flow residual.

       call residual

       ! Substract the residual from the restricted residual and form
       ! the residual forcing term, where the relaxation factor fcoll is
       ! taken into account. Store the restricted residual, currently
       ! stored in wr, in residual afterwards. This is the net result
       ! of adding the residual to the residual forcing term.

       spectralLoop4: do sps=1,nTimeIntervalsSpectral
         domains4: do nn=1,nDom

           ! Have the pointers point to this block.

           call setPointers(nn, currentLevel, sps)

           ! Loop over the owned cells. No need to do anything on the
           ! halo's.

           do l=1,nwf
             do k=2,kl
               do j=2,jl
                 do i=2,il
                   tmp = fcoll*wr(i,j,k,l)
                   wr(i,j,k,l) = tmp - dw(i,j,k,l)
                   dw(i,j,k,l) = tmp
                 enddo
               enddo
             enddo
           enddo

         enddo domains4
       enddo spectralLoop4

       end subroutine transferToCoarseGrid
