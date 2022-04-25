module multigrid

contains

  subroutine transferToCoarseGrid
    !
    !       transferToCoarseGrid restricts both the solution and the
    !       residual to the next coarser grid level and computes the
    !       residual forcing term on this level.
    !
    use constants
    use blockPointers, only : flowDoms, dw, il, jl, kl, ie, je, ke, w, &
         p1, p, rev, w1, mgIFine, mgJFine, mgKFine, nDom, wr, mgIWeight, &
         mgJWeight, mgKWeight, iblank
    use flowVarRefState, only : nwf, kPresent
    use inputIteration, only: fcoll
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : currentLevel, rkStage
    use utils, only : setPointers
    use haloExchange, only : whalo1
    use flowUtils, only : computeEtotBlock, computeLamViscosity
    use turbutils, only : computeeddyviscosity
    use solverUtils, only : timeStep
    use residuals, only : residual, initRes, sourceTerms
    use BCRoutines, only : applyAllBC
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
    call sourceTerms()
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
          call computeLamViscosity(.False.)
          call computeEddyViscosity(.False.)

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

  subroutine transferToFineGrid(corrections)
    !
    !       transferToFineGrid interpolates either the corrections or the
    !       solution to the next finer grid level. A standard trilinear
    !       interpolation is used.
    !
    use constants
    use blockPointers, only : flowDoms, dw, il, jl, kl, ie, je, ke, w, &
         p1, p, rev, w1, mgICoarse, mgJCoarse, mgKCoarse, nDom, wr, mgIWeight, &
         mgJWeight, mgKWeight, iblank
    use flowVarRefState, only : nwf, kPresent, pInfCorr, nw, rhoInf, nt1
    use inputPhysics, only : equations
    use inputIteration, only: fcoll, mgBoundCorr
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : currentLevel, rkStage, groundLevel, exchangePressureEarly
    use utils, only : setPointers, getCorrectForK
    use haloExchange, only : whalo1, whalo2
    use flowUtils, only : computeEtotBlock, computeLamViscosity
    use turbutils, only : computeeddyviscosity
    use turbBCRoutines, only : applyAllTurbBC
    use BCRoutines, only : applyAllBC
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

    ! Store the coarse grid level in coarseLevel.

    coarseLevel = currentLevel +1

    ! Set the number of variables for which either the corrections
    ! or the solution must be interpolated. If the corrections must
    ! be interpolated, this value is set to the number of variables
    ! to which multigrid must be applied, otherwise all conservative
    ! variables are interpolated.

    nVarInt = nw
    if( corrections ) nVarInt = nwf

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
    correctForK = getCorrectForK()

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

          call computeEtotBlock(2_intType,il, 2_intType,jl, &
               2_intType,kl, correctForK)
          if(.not. corrections) call extrapolateSolution

          ! Compute the laminar viscosity and eddy viscosity for the
          ! owned cells of this block. If the solution must be
          ! interpolated, extrapolate the viscosities in the halo's.

          call computeLamViscosity(.False.)
          call computeEddyViscosity(.False.)
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

    if(.not. corrections .and. equations==RANSEquations) then
       call applyAllTurbBC(secondHalo)
    end if

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
    !       extrapolateSolution sets the solution of the cell halos by a
    !       constant extrapolation. This routine is called after the
    !       solution has been interpolated to the next finer grid and this
    !       routine makes sure that the halo's are initialized.
    !       Only the block to which the pointers in blockPointers
    !       currently point is treated.
    !
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, l

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
    !       extrapolateViscosities sets the laminar and eddy viscosities
    !       of the cell halos by a constant extrapolation. This routine is
    !       called after the solution has been interpolated to the next
    !       finer grid and this routine makes sure that the halo's are
    !       initialized.
    !       Only the block to which the pointers in blockPointers
    !       currently point is treated and only for a viscous problem.
    !
    use blockPointers
    use flowVarRefState
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k

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

  subroutine executeMGCycle
    !
    !       executeMGCycle performs a multigrid cycle defined by
    !       cycling, see the module iteration.
    !
    use flowVarRefState
    use iteration
    use inputIteration
    use inputPhysics
    use utils, only : terminate
    use turbAPI, only : turbSolveDDADI
    use solverUtils, only : timeStep, computeUtau
    use smoothers, only : rungeKuttaSmoother, DADISmoother
    use residuals, only : residual, initRes, sourceTerms
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    ! Initialize currentLevel to groundLevel, the ground level
    ! of this multigrid cycle.

    currentLevel = groundLevel

    ! Loop over the number of steps in cycling.

    mgLoop: do nn=1,nstepsCycling

       ! Determine what must be done.

       select case (cycling(nn))

       case (-1_intType)

          ! Set the new currentLevel and and prolongate the
          ! the corrections to this grid level.

          currentLevel = currentLevel -1
          call transferToFineGrid(.true.)

       case ( 0_intType)

          ! Perform a smoothing iteration on this grid level.
          ! First determine the current situation. If this is the
          ! first entry in cycling the residual already contains the
          ! correct values and therefore it needs not to be computed.

          if(nn > 1) then

             ! Compute the residual if the previous action was not a
             ! restriction. In that case the residual already contains
             ! the correct values.

             if(cycling(nn-1) /= 1_intType) then

                ! Initialize and compute the residual.

                rkStage = 0
                call timeStep(.false.)

                call initres(1_intType, nwf)
                call sourceTerms()
                call residual

             endif
          endif

          ! Perform a smoothing step. Determine the smoother to
          ! be used and call the appropriate routine.

          select case (smoother)
          case (RungeKutta)
             call RungeKuttaSmoother
             iterType = "      RK"
          case (DADI)
             call DADISmoother
             iterType = "    DADI"
          case (nlLusgs)
             call terminate("executeMGCycle", &
                  "nlLusgs smoother not implemented yet")
          case (nlLusgsLine)
             call terminate("executeMGCycle", &
                  "nlLusgsLine smoother not implemented &
                  &yet")
          end select

       case ( 1_intType)

          ! Restrict the solution and residual to the next coarser
          ! grid level. Inside transferToCoarseGrid currentLevel
          ! is updated.

          call transferToCoarseGrid

       end select

    enddo mgLoop

    ! Reset the values of rkStage and currentLevel, such that
    ! they correspond to a new iteration.

    rkStage = 0
    currentLevel = groundLevel

    ! Compute the latest values of the skin friction velocity.
    ! The currently stored values are of the previous iteration.

    call computeUtau

    ! Apply an iteration to the turbulent transport equations in
    ! case these must be solved separately.

    if (equations == RANSEquations) then
       call turbSolveDDADI
    end if

    ! Compute the time step.

    call timeStep(.false.)

    ! Compute the residual of the new solution on the ground level.

    call initres(1_intType, nwf)
    call sourceTerms()
    call residual

    ! Set some information for monitoring purposes
    approxTotalIts = approxTotalIts + 1

  end subroutine executeMGCycle

  subroutine setCycleStrategy
    !
    !       setCycleStrategy sets the multigrid cycling strategy for the
    !       multigrid level groundLevel. It is cycle strategy for the
    !       fine grid cut off at the current grid level. If the grid level
    !       is not in the range of the fine grid cycle strategy, cycling
    !       will be set to a single grid strategy.
    !
    use constants
    use inputIteration, only : nMGSteps, cycleStrategy
    use iteration, only : cycling, groundLevel, nStepsCycling
    use utils, only : returnFail
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i
    integer(kind=intType) :: thisLevel, maxLevel

    ! Initialize thisLevel and maxLevel to 1, i.e. the finest grid.

    thisLevel = 1
    maxLevel  = 1

    ! Determine the cycling strategy for groundLevel by looping over
    ! the fine grid cycling strategy and picking the correct entries.

    nStepsCycling = 0
    do i=1,nMGSteps
       thisLevel = thisLevel + cycleStrategy(i)
       maxLevel  = max(maxLevel, thisLevel)

       ! Store this entry in cycling if a) we are on a coarser grid
       ! than groundLevel or b) if we are on groundLevel and
       ! cycleStrategy(i) does not correspond to a restriction,
       ! i.e. 1.

       if((thisLevel == groundLevel .and. cycleStrategy(i) /= 1) .or. &
            thisLevel > groundLevel) then
          nStepsCycling = nStepsCycling + 1
          cycling(nstepsCycling) = cycleStrategy(i)
       endif

       ! Break the loop if a cycle on the current grid level has
       ! been completed.

       if(thisLevel == groundLevel .and. cycleStrategy(i) == -1) exit

    enddo

    ! Take care of the case that groundLevel >= maxLevel.
    ! In this case a single grid strategy is used.

    if(groundLevel >= maxLevel) then
       nStepsCycling = 1
       cycling(1)    = 0
    endif

    ! Check in debug mode if the multigrid strategy created is
    ! a valid one.

    if( debug ) then

       thisLevel = 0
       do i=1,nstepsCycling
          thisLevel = thisLevel + cycling(i)
       enddo

       if(thisLevel /= 0) &
            call returnFail("setCyleStrategy", "Invalid strategy created")

    endif

  end subroutine setCycleStrategy

  subroutine setCornerRowHalos(nVar)
    !
    !       setCornerRowHalos initializes the halo's next to corner row
    !       halo's, such that it contains some values. Otherwise it may
    !       be uninitialized or cause a floating point exception, as this
    !       memory is also used to compute the mg corrections.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use blockPointers, only : w, p, rlv, rev, nx, ny, nz,  &
         il, ie, jl, je, kl, ke
    use flowVarRefState, only : eddyModel, viscous
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in) :: nVar
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, l, mm, ll

    !       Halo's on the i=iMin and i=iMax plane.
    !
    ! K-rows.

    mm = min(3_intType,jl)
    ll = max(2_intType,ny)

    do k=2,kl
       do l=1,nVar
          w(1, 2, k,l) = w(2, 2, k,l)
          w(1, mm,k,l) = w(2, mm,k,l)
          w(1, jl,k,l) = w(2, jl,k,l)
          w(1, ll,k,l) = w(2, ll,k,l)
          w(ie,2, k,l) = w(il,2, k,l)
          w(ie,mm,k,l) = w(il,mm,k,l)
          w(ie,jl,k,l) = w(il,jl,k,l)
          w(ie,ll,k,l) = w(il,ll,k,l)
       enddo

       p(1, 2, k) = p(2, 2, k)
       p(1, mm,k) = p(2, mm,k)
       p(1, jl,k) = p(2, jl,k)
       p(1, ll,k) = p(2, ll,k)
       p(ie,2, k) = p(il,2, k)
       p(ie,mm,k) = p(il,mm,k)
       p(ie,jl,k) = p(il,jl,k)
       p(ie,ll,k) = p(il,ll,k)

       if( viscous) then
          rlv(1, 2, k) = rlv(2, 2, k)
          rlv(1, mm,k) = rlv(2, mm,k)
          rlv(1, jl,k) = rlv(2, jl,k)
          rlv(1, ll,k) = rlv(2, ll,k)
          rlv(ie,2, k) = rlv(il,2, k)
          rlv(ie,mm,k) = rlv(il,mm,k)
          rlv(ie,jl,k) = rlv(il,jl,k)
          rlv(ie,ll,k) = rlv(il,ll,k)
       endif

       if( eddyModel ) then
          rev(1, 2, k) = rev(2, 2, k)
          rev(1, mm,k) = rev(2, mm,k)
          rev(1, jl,k) = rev(2, jl,k)
          rev(1, ll,k) = rev(2, ll,k)
          rev(ie,2, k) = rev(il,2, k)
          rev(ie,mm,k) = rev(il,mm,k)
          rev(ie,jl,k) = rev(il,jl,k)
          rev(ie,ll,k) = rev(il,ll,k)
       endif
    enddo

    ! J-rows; no need to include the corners. These have been set in
    ! the previous k-loop.

    mm = min(3_intType,kl)
    ll = max(2_intType,nz)

    do j=3,ny
       do l=1,nVar
          w(1, j,2, l) = w(2, j,2, l)
          w(1, j,mm,l) = w(2, j,mm,l)
          w(1, j,kl,l) = w(2, j,kl,l)
          w(1, j,ll,l) = w(2, j,ll,l)
          w(ie,j,2, l) = w(il,j,2, l)
          w(ie,j,mm,l) = w(il,j,mm,l)
          w(ie,j,kl,l) = w(il,j,kl,l)
          w(ie,j,ll,l) = w(il,j,ll,l)
       enddo

       p(1, j, 2) = p(2, j, 2)
       p(1, j,mm) = p(2, j,mm)
       p(1, j,kl) = p(2, j,kl)
       p(1, j,ll) = p(2, j,ll)
       p(ie,j, 2) = p(il,j, 2)
       p(ie,j,mm) = p(il,j,mm)
       p(ie,j,kl) = p(il,j,kl)
       p(ie,j,ll) = p(il,j,ll)

       if( viscous) then
          rlv(1, j, 2) = rlv(2, j, 2)
          rlv(1, j,mm) = rlv(2, j,mm)
          rlv(1, j,kl) = rlv(2, j,kl)
          rlv(1, j,ll) = rlv(2, j,ll)
          rlv(ie,j, 2) = rlv(il,j, 2)
          rlv(ie,j,mm) = rlv(il,j,mm)
          rlv(ie,j,kl) = rlv(il,j,kl)
          rlv(ie,j,ll) = rlv(il,j,ll)
       endif

       if( eddyModel ) then
          rev(1, j, 2) = rev(2, j, 2)
          rev(1, j,mm) = rev(2, j,mm)
          rev(1, j,kl) = rev(2, j,kl)
          rev(1, j,ll) = rev(2, j,ll)
          rev(ie,j, 2) = rev(il,j, 2)
          rev(ie,j,mm) = rev(il,j,mm)
          rev(ie,j,kl) = rev(il,j,kl)
          rev(ie,j,ll) = rev(il,j,ll)
       endif
    enddo
    !
    !       Halo's on the j=jMin and j=jMax plane.
    !
    ! K-rows; no need to include the corners; this is done in the
    ! next i-loop.

    mm = min(3_intType,il)
    ll = max(2_intType,nx)

    do k=3,nz
       do l=1,nVar
          w(2, 1, k,l) = w(2, 2, k,l)
          w(mm,1, k,l) = w(mm,2, k,l)
          w(il,1, k,l) = w(il,2, k,l)
          w(ll,1, k,l) = w(ll,2, k,l)
          w(2, je,k,l) = w(2, jl,k,l)
          w(mm,je,k,l) = w(mm,jl,k,l)
          w(il,je,k,l) = w(il,jl,k,l)
          w(ll,je,k,l) = w(ll,jl,k,l)
       enddo

       p(2, 1, k) = p(2, 2, k)
       p(mm,1, k) = p(mm,2, k)
       p(il,1, k) = p(il,2, k)
       p(ll,1, k) = p(ll,2, k)
       p(2, je,k) = p(2, jl,k)
       p(mm,je,k) = p(mm,jl,k)
       p(il,je,k) = p(il,jl,k)
       p(ll,je,k) = p(ll,jl,k)

       if( viscous) then
          rlv(2, 1, k) = rlv(2, 2, k)
          rlv(mm,1, k) = rlv(mm,2, k)
          rlv(il,1, k) = rlv(il,2, k)
          rlv(ll,1, k) = rlv(ll,2, k)
          rlv(2, je,k) = rlv(2, jl,k)
          rlv(mm,je,k) = rlv(mm,jl,k)
          rlv(il,je,k) = rlv(il,jl,k)
          rlv(ll,je,k) = rlv(ll,jl,k)
       endif

       if( eddyModel ) then
          rev(2, 1, k) = rev(2, 2, k)
          rev(mm,1, k) = rev(mm,2, k)
          rev(il,1, k) = rev(il,2, k)
          rev(ll,1, k) = rev(ll,2, k)
          rev(2, je,k) = rev(2, jl,k)
          rev(mm,je,k) = rev(mm,jl,k)
          rev(il,je,k) = rev(il,jl,k)
          rev(ll,je,k) = rev(ll,jl,k)
       endif
    enddo

    ! I-rows, including halo's set on the iMin and iMax plane.

    mm = min(3_intType,kl)
    ll = max(2_intType,nz)

    do i=1,ie
       do l=1,nVar
          w(i,1, 2, l) = w(i,2, 2, l)
          w(i,1, mm,l) = w(i,2, mm,l)
          w(i,1, kl,l) = w(i,2, kl,l)
          w(i,1, ll,l) = w(i,2, ll,l)
          w(i,je,2, l) = w(i,jl,2, l)
          w(i,je,mm,l) = w(i,jl,mm,l)
          w(i,je,kl,l) = w(i,jl,kl,l)
          w(i,je,ll,l) = w(i,jl,ll,l)
       enddo

       p(i, 1, 2) = p(i, 2, 2)
       p(i, 1,mm) = p(i, 2,mm)
       p(i, 1,kl) = p(i, 2,kl)
       p(i, 1,ll) = p(i, 2,ll)
       p(i,je, 2) = p(i,jl, 2)
       p(i,je,mm) = p(i,jl,mm)
       p(i,je,kl) = p(i,jl,kl)
       p(i,je,ll) = p(i,jl,ll)

       if( viscous) then
          rlv(i, 1, 2) = rlv(i, 2, 2)
          rlv(i, 1,mm) = rlv(i, 2,mm)
          rlv(i, 1,kl) = rlv(i, 2,kl)
          rlv(i, 1,ll) = rlv(i, 2,ll)
          rlv(i,je, 2) = rlv(i,jl, 2)
          rlv(i,je,mm) = rlv(i,jl,mm)
          rlv(i,je,kl) = rlv(i,jl,kl)
          rlv(i,je,ll) = rlv(i,jl,ll)
       endif

       if( eddyModel ) then
          rev(i, 1, 2) = rev(i, 2, 2)
          rev(i, 1,mm) = rev(i, 2,mm)
          rev(i, 1,kl) = rev(i, 2,kl)
          rev(i, 1,ll) = rev(i, 2,ll)
          rev(i,je, 2) = rev(i,jl, 2)
          rev(i,je,mm) = rev(i,jl,mm)
          rev(i,je,kl) = rev(i,jl,kl)
          rev(i,je,ll) = rev(i,jl,ll)
       endif
    enddo
    !
    !       Halo's on the k=kMin and k=kMax plane.
    !
    ! J-rows, including halo's set on the jMin and jMax plane.

    mm = min(3_intType,il)
    ll = max(2_intType,nx)

    do j=1,je
       do l=1,nVar
          w(2, j,1, l) = w(2, j,2, l)
          w(mm,j,1, l) = w(mm,j,2, l)
          w(il,j,1, l) = w(il,j,2, l)
          w(ll,j,1, l) = w(ll,j,2, l)
          w(2, j,ke,l) = w(2, j,kl,l)
          w(mm,j,ke,l) = w(mm,j,kl,l)
          w(il,j,ke,l) = w(il,j,kl,l)
          w(ll,j,ke,l) = w(ll,j,kl,l)
       enddo

       p( 2,j, 1) = p( 2,j, 2)
       p(mm,j, 1) = p(mm,j, 2)
       p(il,j, 1) = p(il,j, 2)
       p(ll,j, 1) = p(ll,j, 2)
       p( 2,j,ke) = p( 2,j,kl)
       p(mm,j,ke) = p(mm,j,kl)
       p(il,j,ke) = p(il,j,kl)
       p(ll,j,ke) = p(ll,j,kl)

       if( viscous) then
          rlv( 2,j, 1) = rlv( 2,j, 2)
          rlv(mm,j, 1) = rlv(mm,j, 2)
          rlv(il,j, 1) = rlv(il,j, 2)
          rlv(ll,j, 1) = rlv(ll,j, 2)
          rlv( 2,j,ke) = rlv( 2,j,kl)
          rlv(mm,j,ke) = rlv(mm,j,kl)
          rlv(il,j,ke) = rlv(il,j,kl)
          rlv(ll,j,ke) = rlv(ll,j,kl)
       endif

       if( eddyModel ) then
          rev( 2,j, 1) = rev( 2,j, 2)
          rev(mm,j, 1) = rev(mm,j, 2)
          rev(il,j, 1) = rev(il,j, 2)
          rev(ll,j, 1) = rev(ll,j, 2)
          rev( 2,j,ke) = rev( 2,j,kl)
          rev(mm,j,ke) = rev(mm,j,kl)
          rev(il,j,ke) = rev(il,j,kl)
          rev(ll,j,ke) = rev(ll,j,kl)
       endif
    enddo

    ! I-rows, including halo's set on the iMin and iMax plane.

    mm = min(3_intType,jl)
    ll = max(2_intType,ny)

    do i=1,ie
       do l=1,nVar
          w(i, 2, 1,l) = w(i, 2, 2,l)
          w(i,mm, 1,l) = w(i,mm, 2,l)
          w(i,jl, 1,l) = w(i,jl, 2,l)
          w(i,ll, 1,l) = w(i,ll, 2,l)
          w(i, 2,ke,l) = w(i, 2,kl,l)
          w(i,mm,ke,l) = w(i,mm,kl,l)
          w(i,jl,ke,l) = w(i,jl,kl,l)
          w(i,ll,ke,l) = w(i,ll,kl,l)
       enddo

       p(i, 2, 1) = p(i, 2, 2)
       p(i,mm, 1) = p(i,mm, 2)
       p(i,jl, 1) = p(i,jl, 2)
       p(i,ll, 1) = p(i,ll, 2)
       p(i, 2,ke) = p(i, 2,kl)
       p(i,mm,ke) = p(i,mm,kl)
       p(i,jl,ke) = p(i,jl,kl)
       p(i,ll,ke) = p(i,ll,kl)

       if( viscous) then
          rlv(i, 2, 1) = rlv(i, 2, 2)
          rlv(i,mm, 1) = rlv(i,mm, 2)
          rlv(i,jl, 1) = rlv(i,jl, 2)
          rlv(i,ll, 1) = rlv(i,ll, 2)
          rlv(i, 2,ke) = rlv(i, 2,kl)
          rlv(i,mm,ke) = rlv(i,mm,kl)
          rlv(i,jl,ke) = rlv(i,jl,kl)
          rlv(i,ll,ke) = rlv(i,ll,kl)
       endif

       if( eddyModel ) then
          rev(i, 2, 1) = rev(i, 2, 2)
          rev(i,mm, 1) = rev(i,mm, 2)
          rev(i,jl, 1) = rev(i,jl, 2)
          rev(i,ll, 1) = rev(i,ll, 2)
          rev(i, 2,ke) = rev(i, 2,kl)
          rev(i,mm,ke) = rev(i,mm,kl)
          rev(i,jl,ke) = rev(i,jl,kl)
          rev(i,ll,ke) = rev(i,ll,kl)
       endif
    enddo

  end subroutine setCornerRowHalos

  subroutine setCorrectionsCoarseHalos(sps, nn, coarseLevel, &
       fact, nVarInt)
    !
    !       setCorrectionsCoarseHalos sets the values of the coarse
    !       grid boundary halo corrections. For all boundaries, either a
    !       homogeneous Dirichlet condition (fact = 0.0) or a Neumann
    !       condition (fact = 1.0) is used. Exception are symmetry planes,
    !       where a mirroring takes place.
    !
    use constants
    use block, only : BCDataType, flowDoms
    use flowVarRefState, only : nt1
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps, nn, coarseLevel
    integer(kind=intType), intent(in) :: nVarInt
    real(kind=realType), intent(in)   :: fact
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, l, mm
    integer(kind=intType) :: il, jl, kl, ie, je, ke

    real(kind=realType) :: nnx, nny, nnz, vn

    real(kind=realType), dimension(:,:,:,:), pointer :: ww
    real(kind=realType), dimension(:,:,:),   pointer :: ww1, ww2

    type(BCDataType), dimension(:), pointer :: BCData

    ! Set the pointer ww to the coarse grid variables. At the moment
    ! when this routine is called, these contain the corrections in a
    ! normal grid cycle and the true variables when used to
    ! interpolate the solution to the next finer mesh level in full
    ! multigrid mode. Also set the pointer for BCData, such that the
    ! unit normals are accessed easier.

    ww     => flowDoms(nn,coarseLevel,sps)%w
    BCData => flowDoms(nn,coarseLevel,sps)%BCData

    ! Easier storage of the upper coarse block range.

    il = flowDoms(nn,coarseLevel,sps)%il
    jl = flowDoms(nn,coarseLevel,sps)%jl
    kl = flowDoms(nn,coarseLevel,sps)%kl

    ie = flowDoms(nn,coarseLevel,sps)%ie
    je = flowDoms(nn,coarseLevel,sps)%je
    ke = flowDoms(nn,coarseLevel,sps)%ke

    ! Loop over the number of boundary subfaces to correct the
    ! boundary halo values.

    subfacesCoarse: do mm=1,flowDoms(nn,coarseLevel,sps)%nBocos

       ! Set the pointers for ww1 and ww2, depending on the block face
       ! on which this subface is located. Note that bcFaceID is the
       ! same for all spectral modes and only the 1st is allocated.
       ! Therefore the value of the 1st spectral mode is used here.

       select case (flowDoms(nn,coarseLevel,1)%BCFaceID(mm))

       case (iMin)
          ww1 => ww(1, 1:,1:,:); ww2 => ww(2, 1:,1:,:)
       case (iMax)
          ww1 => ww(ie,1:,1:,:); ww2 => ww(il,1:,1:,:)
       case (jMin)
          ww1 => ww(1:,1 ,1:,:); ww2 => ww(1:,2 ,1:,:)
       case (jMax)
          ww1 => ww(1:,je,1:,:); ww2 => ww(1:,jl,1:,:)
       case (kMin)
          ww1 => ww(1:,1:,1 ,:); ww2 => ww(1:,1:,2 ,:)
       case (kMax)
          ww1 => ww(1:,1:,ke,:); ww2 => ww(1:,1:,kl,:)

       end select

       ! Choose what to do based on the BC type. Note that BCType is
       ! the same for all spectral modes and only the 1st is allocated.
       ! Therefore the value of the 1st spectral mode is used here.

       select case (flowDoms(nn,coarseLevel,1)%BCType(mm))

       case (SlidingInterface,   OversetOuterBound,     &
            DomainInterfaceAll, DomainInterfaceRhoUVW, &
            DomainInterfaceP,   DomainInterfaceRho,    &
            DomainInterfaceTotal)

          ! None of these are physical boundary conditions and thus
          ! nothing needs to be done. The halos already contain the
          ! corrections or solution.

       case (Symm)

          ! This is a symmetry plane. Loop over the faces of this
          ! subface.

          do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
             do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                ! Compute twice the normal velocity component of the
                ! internal cell.

                nnx = BCData(mm)%norm(i,j,1)
                nny = BCData(mm)%norm(i,j,2)
                nnz = BCData(mm)%norm(i,j,3)

                vn = two*(ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny &
                     +      ww2(i,j,ivz)*nnz)

                ! Compute the halo state. Make sure that the average
                ! normal velocity component is zero.

                ww1(i,j,irho)  = ww2(i,j,irho)
                ww1(i,j,ivx)   = ww2(i,j,ivx) - vn*nnx
                ww1(i,j,ivy)   = ww2(i,j,ivy) - vn*nny
                ww1(i,j,ivz)   = ww2(i,j,ivz) - vn*nnz
                ww1(i,j,irhoE) = ww2(i,j,irhoE)

                do l=nt1,nVarInt
                   ww1(i,j,l) = ww2(i,j,l)
                enddo

             enddo
          enddo

       case default

          ! Other type of boundary subface. Set the boundary halo's
          ! as fact times the internal cell value.

          do l=1,nVarInt
             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                   ww1(i,j,l) = fact*ww2(i,j,l)
                enddo
             enddo
          enddo

       end select
    enddo subfacesCoarse

  end subroutine setCorrectionsCoarseHalos

end module multigrid
