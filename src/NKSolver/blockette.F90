module blockette

  use constants
  ! This temporary module contains all cache-blocked code. It also
  ! contains the statically allocated variables on which the blocked
  ! code operates.

  ! Dummy Block dimensions
  integer(kind=intType), parameter :: BS=8
  integer(kind=intType), parameter :: bbil=BS+1, bbjl=BS+1, bbkl=BS+1
  integer(kind=intType), parameter :: bbie=BS+2, bbje=BS+2, bbke=BS+2
  integer(kind=intType), parameter :: bbib=BS+3, bbjb=BS+3, bbkb=BS+3

  ! Actual dimensions to execute
  integer(kind=intType) :: nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb

  ! Variables to track transferring variables between blockettes
  integer(kind=intType) :: singleHaloStart, doubleHaloStart, nodeStart

  ! Current indices into the original block
  integer(kind=intType) :: ii, jj, kk

  ! Double halos
  real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb, 1:6) :: w
  real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb) :: P, gamma
  real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb) :: ss ! Entropy

  ! Single halos
  real(kind=realType), dimension(0:bbie, 0:bbje, 0:bbke, 3) :: x
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke):: rlv, rev, vol, aa
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke) :: radI, radJ, radK, dtl
  real(kind=realType),dimension(1:bbie, 1:bbje, 1:bbke, 3) :: dss ! Shock sensor

  ! No halos
  real(kind=realType), dimension(2:bbil, 2:bbjl, 2:bbkl) :: volRef, d2wall
  integer(kind=intType), dimension(2:bbil, 2:bbjl, 2:bbkl) :: iblank

  ! Face Porosities
  integer(kind=porType), dimension(1:bbil, 2:bbjl, 2:bbkl) :: porI
  integer(kind=porType), dimension(2:bbil, 1:bbjl, 2:bbkl) :: porJ
  integer(kind=porType), dimension(2:bbil, 2:bbjl, 1:bbkl) :: porK

  ! Single halos (only owned cells significant)
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:5) :: fw
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:6) :: dw

  ! Face projected areas
  real(kind=realType), dimension(0:bbie, 1:bbje, 1:bbke, 3) :: sI
  real(kind=realType), dimension(1:bbie, 0:bbje, 1:bbke, 3) :: sJ
  real(kind=realType), dimension(1:bbie, 1:bbje, 0:bbke, 3) :: sK

  ! Face velocities
  real(kind=realType), dimension(0:bbie, 1:bbje, 1:bbke) :: sFaceI
  real(kind=realType), dimension(1:bbie, 0:bbje, 1:bbke) :: sFaceJ
  real(kind=realType), dimension(1:bbie, 1:bbje, 0:bbke) :: sFaceK

  ! Nodal gradients
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: ux, uy, uz
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: vx, vy, vz
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: wx, wy, wz
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: qx, qy, qz

  ! Make *all* of these variables tread-private
  !$OMP THREADPRIVATE(nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb)
  !$OMP THREADPRIVATE(w, p, gamma, ss, x, rlv, rev, vol, aa, radI, radJ, radK)
  !$OMP THREADPRIVATE(dss, volRef, d2wall, iblank, porI, porJ, porK, fw, dw)
  !$OMP THREADPRIVATE(sI, sJ, sK, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz)
contains

  subroutine blocketteRes(useDissApprox, useViscApprox, useUpdateIntermed, useFlowRes, useTurbRes, useSpatial, &
       useStoreWall, famLists, funcValues, forces, bcDataNames, bcDataValues, bcDataFamLists)

    ! Copy the values from blockPointers (assumed set) into the
    ! blockette

    use constants
    use block, only : nDom
    use BCRoutines, only : applyallBC_block
    use bcdata, only : setBCData, setBCDataFineGrid
    use turbbcRoutines, only : applyallTurbBCthisblock, bcTurbTreatment
    use inputPhysics , only : turbProd, equationMode, equations, turbModel
    use inputDiscretization, only : lowSpeedPreconditioner, useApproxWallDistance, useBlockettes
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowUtils, only : computeLamViscosity, computePressureSimple, adjustInflowAngle
    use flowVarRefState, only : nwf, nw, nt1, nt2
    use initializeFlow, only : referenceState
    use section, only: sections, nSections
    use iteration, only : rFil, currentLevel
    use haloExchange, only : exchangeCoor, whalo2
    use wallDistance, only : updateWallDistancesQuickly
    use utils, only : setPointers, EChk
    use turbUtils, only : computeEddyViscosity
    use residuals, only : sourceTerms_block
    use surfaceIntegrations, only : getSolution
    use adjointExtra, only : volume_block, metric_block, boundaryNormals, xhalo_block
    use oversetData, only : oversetPresent
    use inputOverset, only : oversetUpdateMode
    use oversetCommUtilities, only : updateOversetConnectivity
    use actuatorRegionData, only : nActuatorRegions
    implicit none

    ! Input/Output
    logical, intent(in), optional :: useDissApprox, useViscApprox, useUpdateIntermed, useFlowRes
    logical, intent(in), optional :: useTurbRes, useSpatial, useStoreWall
    integer(kind=intType), optional, dimension(:, :), intent(in) :: famLists
    real(kind=realType), optional, dimension(:, :), intent(out) :: funcValues
    character, optional, dimension(:, :), intent(in) :: bcDataNames
    real(kind=realType), optional, dimension(:), intent(in) :: bcDataValues
    integer(kind=intType), optional, dimension(:, :) :: bcDataFamLists
    real(kind=realType), intent(out), optional, dimension(:, :, :) :: forces

    ! Misc
    logical :: dissApprox, viscApprox, updateIntermed, flowRes, turbRes, spatial, storeWall
    integer(kind=intType) :: nn, sps, fSize, lstart, lend, iRegion
    real(kind=realType) ::  pLocal

    ! Set the defaults. The default is to compute the full, exact,
    ! RANS residual without updating the spatial values or the local
    ! timeStep.
    dissApprox = .False.
    viscApprox = .False.
    ! Update intermediate flag is to copy out intermediate variables
    ! that are computed during the blockette residual computation from
    ! blockette memory back to the main memory. These are the time
    ! step, spectral radii for all cases, and nodal gradients and
    ! speed of sound squared for viscous simulations. The regular
    ! "block" residuals do not need to copy out these since they
    ! are already computed in place. For the block residual, this
    ! flag only determines if we update the time step along with
    ! the spectral radii.
    updateIntermed = .False.
    flowRes = .True.
    turbRes = .True.
    spatial = .False.
    storeWall = .True.

    ! Parse the input variables
    if (present(useDissApprox)) then
       dissApprox = useDissApprox
    end if

    if (present(useViscApprox)) then
       viscApprox = useViscApprox
    end if

    if (present(useUpdateIntermed)) then
      updateIntermed = useUpdateIntermed
    end if

    if (present(useFlowRes)) then
       flowRes = useFlowRes
    end if

    if (present(useTurbRes)) then
       turbRes = useTurbRes
    end if

    if (present(useSpatial)) then
       spatial = useSpatial
    end if

    if (present(useStoreWall)) then
       storeWall = useStoreWall
    end if

    ! Spatial-only updates first
    if (spatial) then
       call adjustInflowAngle()

       ! Update all the BCData
       call referenceState
       if (present(bcDataNames)) then
          do sps=1,nTimeIntervalsSpectral
             call setBCData(bcDataNames, bcDataValues, bcDataFamLists, sps, &
                  size(bcDataValues), size(bcDataFamLIsts, 2))
          end do
          call setBCDataFineGrid(.true.)
       end if

       do sps=1, nTimeIntervalsSpectral
          do nn=1, nDom
             call setPointers(nn, currentLevel, sps)
             call xhalo_block()
          end do
       end do

       ! Now exchange the coordinates (fine level only)
       call exchangecoor(1)

       do sps=1, nTimeIntervalsSpectral
          ! Update overset connectivity if necessary
          if (oversetPresent .and. &
               (oversetUpdateMode == updateFast .or. &
               oversetUpdateMode == updateFull)) then
             call updateOversetConnectivity(1_intType, sps)
          end if
       end do
    end if

    ! Compute the required derived values and apply the BCs
    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn, currentLevel, sps)

          if (spatial) then
             call volume_block
             call metric_block
             call boundaryNormals

             if (equations == RANSEquations .and. useApproxWallDistance) then
                call updateWallDistancesQuickly(nn, 1, sps)
             end if
          end if

          ! Compute the pressures/viscositites
          call computePressureSimple(.False.)

          ! Compute Laminar/eddy viscosity if required
          call computeLamViscosity(.False.)
          call computeEddyViscosity(.False.)

          ! Make sure to call the turb BC's first incase we need to
          ! correct for K
          if( equations == RANSEquations .and. turbRes) then
             call BCTurbTreatment
             call applyAllTurbBCthisblock(.True.)
          end if
          call applyAllBC_block(.True.)
       end do
    end do

    ! Compute the ranges of the residuals we are dealing with:
    if (flowRes .and. turbRes) then
       lStart = 1
       lEnd =  nw

    else if (flowRes .and. (.not. turbRes)) then
       lStart = 1
       lEnd = nwf

    else if ((.not. flowRes) .and. turbres) then
       lStart = nt1
       lEnd   = nt2
    end if

    ! Exchange values
    call whalo2(1_intType, lStart, lEnd, .True., .True., .True.)

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       do sps=1,nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers(nn, currentLevel, sps)
             if( equations == RANSEquations .and. turbRes) then
                call BCTurbTreatment
                call applyAllTurbBCthisblock(.True.)
             end if
             call applyAllBC_block(.True.)
          end do
       end do
    end if

    ! Main loop for the residual...This is where the blockette magic happens.
    spsLoop: do sps=1, nTimeIntervalsSpectral
       blockLoop: do  nn=1, nDom
          call setPointers(nn, currentLevel, sps)

          rFil = one
          blockettes: if (useBlockettes) then
             call blocketteResCore(dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall)
          else
             call blockResCore(dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall, nn, sps)
          end if blockettes

          if (currentLevel == 1) then
             do iRegion=1, nActuatorRegions
                call sourceTerms_block(nn, .True., iRegion, pLocal)
             end do
          end if
       end do blockLoop
    end do spsLoop

    ! Compute the final solution values
    if (present(famLists)) then
       call getSolution(famLists, funcValues)
    end if

    if (present(forces)) then
       do sps=1, nTimeIntervalsSpectral
          ! Now we can retrieve the forces/tractions for this spectral instance
          fSize = size(forces, 2)
          call getForces(forces(:, :, sps), fSize, sps)
       end do
    end if
  end subroutine blocketteRes

  subroutine blocketteResCore(dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall)

    ! Main subroutine for computing the reisdual for the given block using blockettes
    use constants

    use constants
    use blockPointers, only : &
         bnx=>nx, bny=>ny, bnz=>nz, &
         bil=>il, bjl=>jl, bkl=>kl, &
         bie=>ie, bje=>je, bke=>ke, &
         bib=>ib, bjb=>jb, bkb=>kb, &
         bw=>w, bp=>p, bgamma=>gamma, &
         bradi=>radi, bradj=>radj, bradk=>radk, &
         bux=>ux, buy=>uy, buz=>uz, &
         bvx=>vx, bvy=>vy, bvz=>vz, &
         bwx=>wx, bwy=>wy, bwz=>wz, &
         bqx=>qx, bqy=>qy, bqz=>qz, &
         bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
         biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw, &
         bShockSensor=>shockSensor, &
         bsi=>si, bsj=>sj, bsk=>sk, &
         bsFaceI=>sFaceI, bsFaceJ=>sFaceJ, bsFaceK=>sFaceK , &
         bdtl=>dtl, baa=>aa, &
         addGridVelocities
    use flowVarRefState, only : nwf, nw, viscous, nt1, nt2
    use iteration, only : currentLevel
    use inputPhysics , only : equationMode, equations, turbModel
    use inputDiscretization, only : spaceDiscr
    use utils, only : setPointers, EChk
    use turbUtils, only : computeEddyViscosity
    use oversetData, only : oversetPresent

    implicit none

    ! Input
    logical, intent(in) :: dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall

    ! Working:
    integer(kind=intType) :: i, j, k, l, lStart, lEnd

    ! Compute the ranges of the residuals we are dealing with:
    if (flowRes .and. turbRes) then
       lStart = 1
       lEnd =  nw

    else if (flowRes .and. (.not. turbRes)) then
       lStart = 1
       lEnd = nwf

    else if ((.not. flowRes) .and. turbres) then
       lStart = nt1
       lEnd   = nt2
    end if

    ! Block loop over the owned cells
    !$OMP parallel do private(i,j,k,l) collapse(2)
    do kk=2, bkl, BS
       do jj=2, bjl, BS
          do ii=2, bil, BS

             ! Determine the actual size this block will be and set
             ! the sizes in the blockette module for each of the
             ! subroutines.

             nx = min(ii+BS-1, bil) - ii + 1
             ny = min(jj+BS-1, bjl) - jj + 1
             nz = min(kk+BS-1, bkl) - kk + 1

             il = nx + 1; jl = ny + 1; kl = nz + 1
             ie = nx + 2; je = ny + 2; ke = nz + 2
             ib = nx + 3; jb = ny + 3; kb = nz + 3

             firstBlockette: if (ii==2) then

                ! First loop. Need to compute the extra stuff. Set
                ! the generic starts and copy the extra
                ! variables in to the starting slots
                singleHaloStart = 1
                doubleHaloStart = 0
                nodeStart = 1

                ! Double halos
                do k=0, kb
                   do j=0, jb
                      do i=0, 3
                         w(i,j,k,1:nw) = bw(i+ii-2, j+jj-2, k+kk-2, 1:nw)
                         p(i,j,k) = bP(i+ii-2, j+jj-2, k+kk-2)
                         gamma(i,j,k) = bgamma(i+ii-2, j+jj-2, k+kk-2)
                         if (currentLevel == 1) then
                            ss(i,j,k) = bShockSensor(i+ii-2, j+jj-2,k+kk-2)
                         end if
                      end do
                   end do
                end do

                ! Single halos
                do k=1, ke
                   do j=1, je
                      do i=1, 2
                         rlv(i,j,k) = brlv(i+ii-2, j+jj-2, k+kk-2)
                         rev(i,j,k) = brev(i+ii-2, j+jj-2, k+kk-2)
                         vol(i,j,k) = bvol(i+ii-2, j+jj-2, k+kk-2)
                      end do
                   end do
                end do

                ! X
                do k=0, ke
                   do j=0, je
                      do i=0, 1
                         x(i,j,k,:) = bx(i+ii-2, j+jj-2, k+kk-2, :)
                      end do
                   end do
                end do
             else

                ! Subsequent loop. We can save a bunch of work by
                ! copying some of the pre-computed values from the
                ! previous blockette to this blockette. Basically the
                ! values that are at the "I end" get shuffled back to
                ! the I-start. We *also* do this for some of the
                ! intermediate variables that are costly to compute
                ! like the nodal gradients, and spectral radius which
                ! helps cut back on the amount of data duplication.

                ! Important Note: This cell is not the first cell. If
                ! this code is being executed, the previous blockette
                ! was copied fully in the i direction.
                ! Therefore, we can just copy the values from
                ! the end of the blockette as it is allocated.
                ! To do this, we ignore the dimensions of the "current"
                ! blockette, and just take the baseline BS dimensions
                ! as the current blockette might be partially filled
                ! in the i direction.

                singleHaloStart = 3
                doubleHaloStart = 4
                nodeStart = 2

                ! Double halos
                do k=0, kb
                   do j=0, jb
                      do i=0, 3
                         w(i,j,k,1:nw) = w(BS+i, j, k, 1:nw)
                         p(i,j,k) = p(BS+i, j, k)
                         gamma(i,j,k) = gamma(BS+i, j, k)
                         ss(i,j,k) = ss(BS+i, j, k)
                      end do
                   end do
                end do

                ! Single halos
                do k=1, ke
                   do j=1, je
                      do i=1, 2
                         rlv(i,j,k) = rlv(BS+i, j, k)
                         rev(i,j,k) = rev(BS+i, j, k)
                         vol(i,j,k) = vol(BS+i, j, k)

                         ! Computed variables

                         ! DONT Copy the spectral-radii. The loop that calculates
                         ! spectral radii also calculates portion of the time step,
                         ! so we don't want to mess with its boundaries to keep
                         ! it simple.
                         aa(i,j,k) = aa(BS+i, j, k)
                         dss(i,j,k,:) = dss(BS+i, j, k, :)
                      end do
                   end do
                end do

                ! X
                do k=0, ke
                   do j=0, je
                      do i=0, 1
                         x(i,j,k,:) = x(BS+i, j, k, :)
                      end do
                   end do
                end do

                ! Nodal gradients
                do k=1, kl
                   do j=1, jl
                      ux(1, j, k) = ux(BS+1, j, k)
                      uy(1, j, k) = uy(BS+1, j, k)
                      uz(1, j, k) = uz(BS+1, j, k)

                      vx(1, j, k) = vx(BS+1, j, k)
                      vy(1, j, k) = vy(BS+1, j, k)
                      vz(1, j, k) = vz(BS+1, j, k)

                      wx(1, j, k) = wx(BS+1, j, k)
                      wy(1, j, k) = wy(BS+1, j, k)
                      wz(1, j, k) = wz(BS+1, j, k)

                      qx(1, j, k) = qx(BS+1, j, k)
                      qy(1, j, k) = qy(BS+1, j, k)
                      qz(1, j, k) = qz(BS+1, j, k)
                   end do
                end do
             end if firstBlockette

             ! -------------------------------------
             !      Fill in the remaining values
             ! -------------------------------------

             ! Double halos
             do k=0, kb
                do j=0, jb
                   do i=4, ib
                      w(i,j,k,1:nw) = bw(i+ii-2, j+jj-2, k+kk-2, 1:nw)
                      p(i,j,k) = bP(i+ii-2, j+jj-2, k+kk-2)
                      gamma(i,j,k) = bgamma(i+ii-2, j+jj-2, k+kk-2)
                      if (currentLevel == 1) then
                         ss(i,j,k) = bShockSensor(i+ii-2, j+jj-2,k+kk-2)
                      end if
                   end do
                end do
             end do

             ! Single halos
             do k=1, ke
                do j=1, je
                   do i=3, ie
                      rlv(i,j,k) = brlv(i+ii-2, j+jj-2, k+kk-2)
                      rev(i,j,k) = brev(i+ii-2, j+jj-2, k+kk-2)
                      vol(i,j,k) = bvol(i+ii-2, j+jj-2, k+kk-2)
                   end do
                end do
             end do

             ! X
             do k=0, ke
                do j=0, je
                   do i=2, ie
                      x(i,j,k,:) = bx(i+ii-2, j+jj-2, k+kk-2, :)
                   end do
                end do
             end do

             ! No Halos (no change)
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      iblank(i,j,k) = biblank(i+ii-2,j+jj-2,k+kk-2)
                      if (equations .eq. ransequations) &
                      d2wall(i,j,k) = bd2wall(i+ii-2,j+jj-2,k+kk-2)
                      volRef(i,j,k) = bvolRef(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             ! Porosities (no change)
             do k=2, kl
                do j=2, jl
                   do i=1, il
                      porI(i,j,k) = bporI(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             do k=2, kl
                do j=1, jl
                   do i=2, il
                      PorJ(i,j,k) = bporJ(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             do k=1, kl
                do j=2, jl
                   do i=2, il
                      PorK(i,j,k) = bporK(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             ! Face velocities if necessary
             if (addGridVelocities) then
                do k=1, ke
                   do j=1, je
                      do i=0, ie
                         sFaceI(i, j, k) = bsFaceI(ii+ii-2, j+jj-2, k+kk-2)
                      end do
                   end do
                end do

                do k=1, ke
                   do j=0, je
                      do i=1, ie
                         sFaceJ(i, j, k) = bsFaceJ(ii+ii-2, j+jj-2, k+kk-2)
                      end do
                   end do
                end do

                do k=0, ke
                   do j=1, je
                      do i=1, ie
                         sFaceK(i, j, k) = bsFaceK(ii+ii-2, j+jj-2, k+kk-2)
                      end do
                   end do
                end do
             else
                sFaceI = zero
                sFaceJ = zero
                sFaceK = zero
             end if

             ! Clear the viscous flux before we start.
             fw = zero

             ! Call the routines in order:
             call metrics
             call initRes(lStart, lEnd)

             ! Compute turbulence residual for RANS equations
             if( equations == RANSEquations .and. turbRes) then

                ! Initialize only the Turblent Variables
                !call unsteadyTurbSpectral_block(itu1, itu1, nn, sps)

                select case (turbModel)

                case (spalartAllmaras)
                   call saSource
                   call saAdvection
                   !call unsteadyTurbTerm(1_intType, 1_intType, itu1-1, qq)
                   call saViscous
                   call saResScale
                end select
             endif

             call timeStep(updateIntermed)

             if (flowRes) then
                call inviscidCentralFlux

                if (dissApprox) then
                   select case (spaceDiscr)
                   case (dissScalar)
                      call inviscidDissFluxScalarApprox
                   case (dissMatrix)
                      call inviscidDissFluxMatrixApprox
                   case (upwind)
                      call inviscidUpwindFlux(.False.)
                   end select
                else
                   select case (spaceDiscr)
                   case (dissScalar)
                      call inviscidDissFluxScalar
                   case (dissMatrix)
                      call inviscidDissFluxMatrix
                   case (upwind)
                      call inviscidUpwindFlux(.True.)
                   end select
                end if

                if (viscous) then
                   call computeSpeedOfSoundSquared
                   if (viscApprox) then
                      call viscousFluxApprox
                   else
                      call allNodalGradients
                      call viscousFlux(storeWall)
                   end if
                end if

                call sumDwAndFw
             end if

             ! Now we can just set the part of dw we computed
             ! (owned cells only) and we're done!
             do l=lStart, lEnd
                do k=2, kl
                   do j=2, jl
                      do i=2, il
                         bdw(i+ii-2,j+jj-2,k+kk-2,l) = dw(i,j,k,l)
                      end do
                   end do
                end do
             end do

             ! Also copy out the intermediate variables if asked for them
             ! we need these to be updated in main memory because
             ! the reverse mode AD routines do use these variables.
             ! after every ANK and NK step, blocketteRes is called
             ! with updateIntermed = True, and it will update these
             ! arrays in main memory. The time step is required
             ! for the ANK and MG solver steps.
             intermed: if (updateIntermed) then
               ! time step
               do k=2, kl
                  do j=2, jl
                     do i=2, il
                        bdtl(i+ii-2, j+jj-2, k+kk-2) = dtl(i, j, k)
                     end do
                  end do
               end do

               ! Spectral radii
               do k=1, ke
                  do j=1, je
                     do i=1, ie
                        bradi(i+ii-2, j+jj-2, k+kk-2) = radi(i, j, k)
                        bradj(i+ii-2, j+jj-2, k+kk-2) = radj(i, j, k)
                        bradk(i+ii-2, j+jj-2, k+kk-2) = radk(i, j, k)
                     end do
                  end do
               end do

               ! need aa and nodal gradients if we have viscous fluxes
               visc: if (viscous .and. flowRes) then

                  ! speed of sound squared
                  do k=1, ke
                     do j=1, je
                        do i=1, ie
                           baa(i+ii-2, j+jj-2, k+kk-2) = aa(i, j, k)
                        end do
                     end do
                  end do

                  ! nodal gradients
                  do k=1, kl
                     do j=1, jl
                        do i=1, il

                           bux(i+ii-2, j+jj-2, k+kk-2) = ux(i, j, k)
                           buy(i+ii-2, j+jj-2, k+kk-2) = uy(i, j, k)
                           buz(i+ii-2, j+jj-2, k+kk-2) = uz(i, j, k)

                           bvx(i+ii-2, j+jj-2, k+kk-2) = vx(i, j, k)
                           bvy(i+ii-2, j+jj-2, k+kk-2) = vy(i, j, k)
                           bvz(i+ii-2, j+jj-2, k+kk-2) = vz(i, j, k)

                           bwx(i+ii-2, j+jj-2, k+kk-2) = wx(i, j, k)
                           bwy(i+ii-2, j+jj-2, k+kk-2) = wy(i, j, k)
                           bwz(i+ii-2, j+jj-2, k+kk-2) = wz(i, j, k)

                           bqx(i+ii-2, j+jj-2, k+kk-2) = qx(i, j, k)
                           bqy(i+ii-2, j+jj-2, k+kk-2) = qy(i, j, k)
                           bqz(i+ii-2, j+jj-2, k+kk-2) = qz(i, j, k)

                        end do
                     end do
                  end do
               end if visc

             end if intermed

          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine blocketteResCore

  subroutine blockResCore(dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall, nn, sps)

    use constants
    use fluxes, only : inviscidCentralFlux_block=>inviscidCentralFlux, &
         inviscidDissFluxScalar_block=>inviscidDissFluxScalar, &
         inviscidDissFluxMatrix_block=>inviscidDissFluxMatrix, &
         inviscidUpwindFlux_block=>inviscidUpwindFlux, &
         inviscidDissFluxScalarApprox_block=>inviscidDissFluxScalarApprox, &
         inviscidDissFluxMatrixApprox_block=>inviscidDissFluxMatrixApprox, &
         viscousFlux_block=>viscousFlux, &
         viscousFluxApprox_block=>viscousFluxApprox
    use solverUtils, only : timeStep_block
    use flowVarRefState, only : nwf, nw, viscous, nt1, nt2
    use inputPhysics , only : equationMode, equations, turbModel
    use residuals, only : initres_block
    use sa, only : sa_block
    use adjointExtra, only : sumDwAndFw_block=>sumDwAndFw
    use inputDiscretization, only : spaceDiscr
    use flowUtils, only : allNodalGradients_block=>allNodalGradients, &
         computeSpeedOfSoundSquared_block=>computeSpeedOfSoundSquared

    implicit none
    ! Input
    logical, intent(in) :: dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall
    integer(kind=intType), intent(in) :: nn, sps

    ! Working:
    integer(kind=intType) :: i, j, k, lStart, lEnd

    ! Compute the ranges of the residuals we are dealing with:
    if (flowRes .and. turbRes) then
       lStart = 1
       lEnd =  nw

    else if (flowRes .and. (.not. turbRes)) then
       lStart = 1
       lEnd = nwf

    else if ((.not. flowRes) .and. turbres) then
       lStart = nt1
       lEnd   = nt2
    end if

    ! Compute time step
    call timestep_block(.not. updateIntermed)

    call initres_block(lStart, lEnd, nn, sps) ! Initialize only the Turblent Variables

    fw = zero

    ! Possible Turblent Equations
    if(equations == RANSEquations .and. turbRes) then
       ! Compute the skin-friction velocity (wall functions only)
       !call computeUtau_block

       ! Now call the selected turbulence model
       select case (turbModel)
       case (spalartAllmaras)
          call sa_block(.true.)
       end select
    endif

    if (flowRes) then

       call inviscidCentralFlux_block
       if (dissApprox) then
          select case (spaceDiscr)
          case (dissScalar)
             call inviscidDissFluxScalarApprox_block
          case (dissMatrix)
             call inviscidDissFluxMatrixApprox_block
          case (upwind)
             call inviscidUpwindFlux_block(.True.)
          end select
       else
          select case (spaceDiscr)
          case (dissScalar)
             call inviscidDissFluxScalar_block
          case (dissMatrix)
             call inviscidDissFluxMatrix_block
          case (upwind)
             call inviscidUpwindFlux_block(.True.)
          end select
       end if

       if (viscous) then
          call computeSpeedOfSoundSquared_block
          if (viscApprox) then
             call viscousFluxApprox_block
          else
             call allNodalGradients_block
             call viscousFlux_block
          end if
       end if

       call sumDwAndFw_block
    end if
  end subroutine blockResCore

  subroutine metrics
    ! ---------------------------------------------
    !              Metric computation
    ! ---------------------------------------------

    use constants
    use blockPointers, only : rightHanded
    implicit none

    integer(kind=intType) :: i, j, k, l, m, n
    real(kind=realType), dimension(3) :: v1, v2
    real(kind=realType) :: fact

    ! Projected areas of cell faces in the i direction.
    if (rightHanded) then
       fact = half
    else
       fact = -half
    end if
    do k=1,ke
       n = k -1
       do j=1,je
          m = j -1
          do i=0,ie

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,n,1) - x(i,m,k,1)
             v1(2) = x(i,j,n,2) - x(i,m,k,2)
             v1(3) = x(i,j,n,3) - x(i,m,k,3)

             v2(1) = x(i,j,k,1) - x(i,m,n,1)
             v2(2) = x(i,j,k,2) - x(i,m,n,2)
             v2(3) = x(i,j,k,3) - x(i,m,n,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo

    ! Projected areas of cell faces in the j direction.

    do k=1,ke
       n = k -1
       do j=0,je
          do i=1,ie
             l = i -1

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,n,1) - x(l,j,k,1)
             v1(2) = x(i,j,n,2) - x(l,j,k,2)
             v1(3) = x(i,j,n,3) - x(l,j,k,3)

             v2(1) = x(l,j,n,1) - x(i,j,k,1)
             v2(2) = x(l,j,n,2) - x(i,j,k,2)
             v2(3) = x(l,j,n,3) - x(i,j,k,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo

    ! Projected areas of cell faces in the k direction.

    do k=0,ke
       do j=1,je
          m = j -1
          do i=1,ie
             l = i -1

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,k,1) - x(l,m,k,1)
             v1(2) = x(i,j,k,2) - x(l,m,k,2)
             v1(3) = x(i,j,k,3) - x(l,m,k,3)

             v2(1) = x(l,j,k,1) - x(i,m,k,1)
             v2(2) = x(l,j,k,2) - x(i,m,k,2)
             v2(3) = x(l,j,k,3) - x(i,m,k,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo
  end subroutine metrics

  subroutine initRes(varStart, varEnd)
    ! ---------------------------------------------
    !                     Init Res
    ! ---------------------------------------------

    use constants
    implicit none

    integer(kind=intType) :: varStart, varEnd
    ! Obviously this needs to be more complex for the actual code.
    dw(:, :, :, varStart:varEnd) = zero

  end subroutine initRes

  subroutine saSource
    ! ---------------------------------------------
    !                    SA Source Term
    ! ---------------------------------------------

    use constants
    use paramTurb
    use blockPointers, only : sectionID
    use inputPhysics, only :useft2SA, useRotationSA, turbProd, equations
    use inputDiscretization, only : approxSA
    use section, only : sections
    use sa, only : cv13, kar2Inv, cw36, cb3Inv
    use flowvarRefState, only : timeRef

    implicit none

    ! Variables for sa Souce
    real(kind=realType) :: fv1, fv2, ft2
    real(kind=realType) :: sst, nu, dist2Inv, chi, chi2, chi3
    real(kind=realType) :: rr, gg, gg6, termFw, fwSa, term1, term2
    real(kind=realType) :: dfv1, dfv2, dft2, drr, dgg, dfw, sqrtProd
    real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
    real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
    real(kind=realType) :: vortx, vorty, vortz
    real(kind=realType) :: omegax, omegay, omegaz
    real(kind=realType) :: strainMag2, prod
    real(kind=realType), parameter :: xminn = 1.e-10_realType
    real(kind=realType), parameter :: f23 = two*third
    integer(kind=intType) :: i, j, k
    real(kind=realType) :: term1Fact

    ! Set model constants
    cv13    = rsaCv1**3
    kar2Inv = one/(rsaK**2)
    cw36    = rsaCw3**6
    cb3Inv  = one/rsaCb3

    ! set the approximate multiplier here
    term1Fact = one
    if (approxSA) term1Fact = zero

    ! Determine the non-dimensional wheel speed of this block.

    omegax = timeRef*sections(sectionID)%rotRate(1)
    omegay = timeRef*sections(sectionID)%rotRate(2)
    omegaz = timeRef*sections(sectionID)%rotRate(3)
    do k=2, kl
       do j=2, jl
          do i=2, il

             ! Compute the gradient of u in the cell center. Use is made
             ! of the fact that the surrounding normals sum up to zero,
             ! such that the cell i,j,k does not give a contribution.
             ! The gradient is scaled by the factor 2*vol.

             uux = w(i+1,j,k,ivx)*si(i,j,k,1) - w(i-1,j,k,ivx)*si(i-1,j,k,1) &
                  + w(i,j+1,k,ivx)*sj(i,j,k,1) - w(i,j-1,k,ivx)*sj(i,j-1,k,1) &
                  + w(i,j,k+1,ivx)*sk(i,j,k,1) - w(i,j,k-1,ivx)*sk(i,j,k-1,1)
             uuy = w(i+1,j,k,ivx)*si(i,j,k,2) - w(i-1,j,k,ivx)*si(i-1,j,k,2) &
                  + w(i,j+1,k,ivx)*sj(i,j,k,2) - w(i,j-1,k,ivx)*sj(i,j-1,k,2) &
                  + w(i,j,k+1,ivx)*sk(i,j,k,2) - w(i,j,k-1,ivx)*sk(i,j,k-1,2)
             uuz = w(i+1,j,k,ivx)*si(i,j,k,3) - w(i-1,j,k,ivx)*si(i-1,j,k,3) &
                  + w(i,j+1,k,ivx)*sj(i,j,k,3) - w(i,j-1,k,ivx)*sj(i,j-1,k,3) &
                  + w(i,j,k+1,ivx)*sk(i,j,k,3) - w(i,j,k-1,ivx)*sk(i,j,k-1,3)

             ! Idem for the gradient of v.

             vvx = w(i+1,j,k,ivy)*si(i,j,k,1) - w(i-1,j,k,ivy)*si(i-1,j,k,1) &
                  + w(i,j+1,k,ivy)*sj(i,j,k,1) - w(i,j-1,k,ivy)*sj(i,j-1,k,1) &
                  + w(i,j,k+1,ivy)*sk(i,j,k,1) - w(i,j,k-1,ivy)*sk(i,j,k-1,1)
             vvy = w(i+1,j,k,ivy)*si(i,j,k,2) - w(i-1,j,k,ivy)*si(i-1,j,k,2) &
                  + w(i,j+1,k,ivy)*sj(i,j,k,2) - w(i,j-1,k,ivy)*sj(i,j-1,k,2) &
                  + w(i,j,k+1,ivy)*sk(i,j,k,2) - w(i,j,k-1,ivy)*sk(i,j,k-1,2)
             vvz = w(i+1,j,k,ivy)*si(i,j,k,3) - w(i-1,j,k,ivy)*si(i-1,j,k,3) &
                  + w(i,j+1,k,ivy)*sj(i,j,k,3) - w(i,j-1,k,ivy)*sj(i,j-1,k,3) &
                  + w(i,j,k+1,ivy)*sk(i,j,k,3) - w(i,j,k-1,ivy)*sk(i,j,k-1,3)

             ! And for the gradient of w.

             wwx = w(i+1,j,k,ivz)*si(i,j,k,1) - w(i-1,j,k,ivz)*si(i-1,j,k,1) &
                  + w(i,j+1,k,ivz)*sj(i,j,k,1) - w(i,j-1,k,ivz)*sj(i,j-1,k,1) &
                  + w(i,j,k+1,ivz)*sk(i,j,k,1) - w(i,j,k-1,ivz)*sk(i,j,k-1,1)
             wwy = w(i+1,j,k,ivz)*si(i,j,k,2) - w(i-1,j,k,ivz)*si(i-1,j,k,2) &
                  + w(i,j+1,k,ivz)*sj(i,j,k,2) - w(i,j-1,k,ivz)*sj(i,j-1,k,2) &
                  + w(i,j,k+1,ivz)*sk(i,j,k,2) - w(i,j,k-1,ivz)*sk(i,j,k-1,2)
             wwz = w(i+1,j,k,ivz)*si(i,j,k,3) - w(i-1,j,k,ivz)*si(i-1,j,k,3) &
                  + w(i,j+1,k,ivz)*sj(i,j,k,3) - w(i,j-1,k,ivz)*sj(i,j-1,k,3) &
                  + w(i,j,k+1,ivz)*sk(i,j,k,3) - w(i,j,k-1,ivz)*sk(i,j,k-1,3)

             ! Compute the components of the stress tensor.
             ! The combination of the current scaling of the velocity
             ! gradients (2*vol) and the definition of the stress tensor,
             ! leads to the factor 1/(4*vol).

             fact = fourth/vol(i,j,k)

             ! -- Calcs for strain --
             sxx = two*fact*uux
             syy = two*fact*vvy
             szz = two*fact*wwz

             sxy = fact*(uuy + vvx)
             sxz = fact*(uuz + wwx)
             syz = fact*(vvz + wwy)

             ! Compute 2/3 * divergence of velocity squared

             div2 = f23*(sxx+syy+szz)**2

             ! Compute strain production term

             strainMag2 = two*(sxy**2 + sxz**2 + syz**2) &
                  +           sxx**2 + syy**2 + szz**2

             ! -- Calcs for vorticity --

             ! Compute the three components of the vorticity vector.
             ! Substract the part coming from the rotating frame.

             vortx = two*fact*(wwy - vvz) - two*omegax
             vorty = two*fact*(uuz - wwx) - two*omegay
             vortz = two*fact*(vvx - uuy) - two*omegaz

             if (turbProd == strain) then
                sqrtProd = sqrt(max(two*strainMag2-div2, eps))
             else
                sqrtProd = sqrt(vortx**2 + vorty**2 + vortz**2)
             end if

             ! Compute the laminar kinematic viscosity, the inverse of
             ! wall distance squared, the ratio chi (ratio of nuTilde
             ! and nu) and the functions fv1 and fv2. The latter corrects
             ! the production term near a viscous wall.

             nu       = rlv(i,j,k)/w(i,j,k,irho)
             dist2Inv = one/(d2Wall(i,j,k)**2)
             chi      = w(i,j,k,itu1)/nu
             chi2     = chi*chi
             chi3     = chi*chi2
             fv1      = chi3/(chi3+cv13)
             fv2      = one - chi/(one + chi*fv1)

             ! The function ft2, which is designed to keep a laminar
             ! solution laminar. When running in fully turbulent mode
             ! this function should be set to 0.0.

             ft2 = zero
             if (useft2SA) then
                ft2 = rsaCt3*exp(-rsaCt4*chi2)
             end if

             ! Correct the production term to account for the influence
             ! of the wall.

             sst = sqrtProd + w(i,j,k,itu1)*fv2*kar2Inv*dist2Inv

             ! Add rotation term (useRotationSA defined in inputParams.F90)

             if (useRotationSA) then
                sst = sst + rsaCrot*min(zero,sqrt(two*strainMag2))
             end if

             ! Make sure that this term remains positive
             ! (the function fv2 is negative between chi = 1 and 18.4,
             ! which can cause sst to go negative, which is undesirable).

             sst = max(sst,xminn)

             ! Compute the function fw. The argument rr is cut off at 10
             ! to avoid numerical problems. This is ok, because the
             ! asymptotical value of fw is then already reached.

             rr     = w(i,j,k,itu1)*kar2Inv*dist2Inv/sst
             rr     = min(rr,10.0_realType)
             gg     = rr + rsaCw2*(rr**6 - rr)
             gg6    = gg**6
             termFw = ((one + cw36)/(gg6 + cw36))**sixth
             fwSa   = gg*termFw

             ! Compute the source term; some terms are saved for the
             ! linearization. The source term is stored in dvt.

             term1 = rsaCb1*(one-ft2)*sqrtProd*term1Fact
             term2 = dist2Inv*(kar2Inv*rsaCb1*((one-ft2)*fv2 + ft2) &
                  -           rsaCw1*fwSa)

             dw(i, j, k, itu1) = dw(i, j, k, itu1) + (term1 + term2*w(i,j,k,itu1))*w(i,j,k,itu1)

          enddo
       enddo
    enddo
  end subroutine saSource

  subroutine saViscous
    ! ---------------------------------------------
    !                    SA Viscous Term
    ! ---------------------------------------------

    use constants
    use sa, only : cv13, kar2Inv, cw36, cb3Inv
    use paramTurb
    implicit none

    ! Variables for sa Viscous
    real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
    real(kind=realType) :: xa, ya, za, ttm, ttp, cnud, cam, cap
    real(kind=realType) :: nutm, nutp, num, nup, cdm, cdp
    real(kind=realType) :: c1m, c1p, c10, b1, c1, d1, qs, nu
    integer(Kind=intType) :: i, j, k

    ! Set model constants
    cv13    = rsaCv1**3
    kar2Inv = one/(rsaK**2)
    cw36    = rsaCw3**6
    cb3Inv  = one/rsaCb3

    !
    !       Viscous terms in k-direction.
    !
    do k=2, kl
       do j=2, jl
          do i=2, il

             ! Compute the metrics in zeta-direction, i.e. along the
             ! line k = constant.

             voli  = one/vol(i,j,k)
             volmi = two/(vol(i,j,k) + vol(i,j,k-1))
             volpi = two/(vol(i,j,k) + vol(i,j,k+1))

             xm = sk(i,j,k-1,1)*volmi
             ym = sk(i,j,k-1,2)*volmi
             zm = sk(i,j,k-1,3)*volmi
             xp = sk(i,j,k,  1)*volpi
             yp = sk(i,j,k,  2)*volpi
             zp = sk(i,j,k,  3)*volpi

             xa  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))*voli
             ya  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))*voli
             za  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))*voli
             ttm = xm*xa + ym*ya + zm*za
             ttp = xp*xa + yp*ya + zp*za

             ! Computation of the viscous terms in zeta-direction; note
             ! that cross-derivatives are neglected, i.e. the mesh is
             ! assumed to be orthogonal.
             ! Furthermore, the grad(nu)**2 has been rewritten as
             ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
             ! The second derivative in zeta-direction is constructed as
             ! the central difference of the first order derivatives, i.e.
             ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
             ! In this way the metric can be taken into account.

             ! Compute the diffusion coefficients multiplying the nodes
             ! k+1, k and k-1 in the second derivative. Make sure that
             ! these coefficients are nonnegative.

             cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
             cam  =  ttm*cnud
             cap  =  ttp*cnud

             nutm = half*(w(i,j,k-1,itu1) + w(i,j,k,itu1))
             nutp = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
             nu   = rlv(i,j,k)/w(i,j,k,irho)
             num  = half*(rlv(i,j,k-1)/w(i,j,k-1,irho) + nu)
             nup  = half*(rlv(i,j,k+1)/w(i,j,k+1,irho) + nu)
             cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
             cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

             c1m = max(cdm+cam, zero)
             c1p = max(cdp+cap, zero)
             c10 = c1m + c1p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, c1 and d1.

             dw(i,j,k,itu1) = dw(i,j,k,itu1)      + c1m*w(i,j,k-1,itu1) &
                  - c10*w(i,j,k,itu1) + c1p*w(i,j,k+1,itu1)
          end do
       enddo
    enddo
    !
    !       Viscous terms in j-direction.
    !
    do k=2, kl
       do j=2, jl
          do i=2, il

             ! Compute the metrics in eta-direction, i.e. along the
             ! line j = constant.

             voli  = one/vol(i,j,k)
             volmi = two/(vol(i,j,k) + vol(i,j-1,k))
             volpi = two/(vol(i,j,k) + vol(i,j+1,k))

             xm = sj(i,j-1,k,1)*volmi
             ym = sj(i,j-1,k,2)*volmi
             zm = sj(i,j-1,k,3)*volmi
             xp = sj(i,j,  k,1)*volpi
             yp = sj(i,j,  k,2)*volpi
             zp = sj(i,j,  k,3)*volpi

             xa  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))*voli
             ya  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))*voli
             za  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))*voli
             ttm = xm*xa + ym*ya + zm*za
             ttp = xp*xa + yp*ya + zp*za

             ! Computation of the viscous terms in eta-direction; note
             ! that cross-derivatives are neglected, i.e. the mesh is
             ! assumed to be orthogonal.
             ! Furthermore, the grad(nu)**2 has been rewritten as
             ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
             ! The second derivative in eta-direction is constructed as
             ! the central difference of the first order derivatives, i.e.
             ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
             ! In this way the metric can be taken into account.

             ! Compute the diffusion coefficients multiplying the nodes
             ! j+1, j and j-1 in the second derivative. Make sure that
             ! these coefficients are nonnegative.

             cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
             cam  =  ttm*cnud
             cap  =  ttp*cnud

             nutm = half*(w(i,j-1,k,itu1) + w(i,j,k,itu1))
             nutp = half*(w(i,j+1,k,itu1) + w(i,j,k,itu1))
             nu   = rlv(i,j,k)/w(i,j,k,irho)
             num  = half*(rlv(i,j-1,k)/w(i,j-1,k,irho) + nu)
             nup  = half*(rlv(i,j+1,k)/w(i,j+1,k,irho) + nu)
             cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
             cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

             c1m = max(cdm+cam, zero)
             c1p = max(cdp+cap, zero)
             c10 = c1m + c1p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, c1 and d1.

             dw(i,j,k,itu1) = dw(i,j,k,itu1)      + c1m*w(i,j-1,k,itu1) &
                  - c10*w(i,j,k,itu1) + c1p*w(i,j+1,k,itu1)

          enddo
       enddo
    enddo
    !
    !       Viscous terms in i-direction.
    !
    do k=2, kl
       do j=2, jl
          do i=2, il

             ! Compute the metrics in xi-direction, i.e. along the
             ! line i = constant.

             voli  = one/vol(i,j,k)
             volmi = two/(vol(i,j,k) + vol(i-1,j,k))
             volpi = two/(vol(i,j,k) + vol(i+1,j,k))

             xm = si(i-1,j,k,1)*volmi
             ym = si(i-1,j,k,2)*volmi
             zm = si(i-1,j,k,3)*volmi
             xp = si(i,  j,k,1)*volpi
             yp = si(i,  j,k,2)*volpi
             zp = si(i,  j,k,3)*volpi

             xa  = half*(si(i,j,k,1) + si(i-1,j,k,1))*voli
             ya  = half*(si(i,j,k,2) + si(i-1,j,k,2))*voli
             za  = half*(si(i,j,k,3) + si(i-1,j,k,3))*voli
             ttm = xm*xa + ym*ya + zm*za
             ttp = xp*xa + yp*ya + zp*za

             ! Computation of the viscous terms in xi-direction; note
             ! that cross-derivatives are neglected, i.e. the mesh is
             ! assumed to be orthogonal.
             ! Furthermore, the grad(nu)**2 has been rewritten as
             ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
             ! The second derivative in xi-direction is constructed as
             ! the central difference of the first order derivatives, i.e.
             ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
             ! In this way the metric can be taken into account.

             ! Compute the diffusion coefficients multiplying the nodes
             ! i+1, i and i-1 in the second derivative. Make sure that
             ! these coefficients are nonnegative.

             cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
             cam  =  ttm*cnud
             cap  =  ttp*cnud

             nutm = half*(w(i-1,j,k,itu1) + w(i,j,k,itu1))
             nutp = half*(w(i+1,j,k,itu1) + w(i,j,k,itu1))
             nu   = rlv(i,j,k)/w(i,j,k,irho)
             num  = half*(rlv(i-1,j,k)/w(i-1,j,k,irho) + nu)
             nup  = half*(rlv(i+1,j,k)/w(i+1,j,k,irho) + nu)
             cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
             cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

             c1m = max(cdm+cam, zero)
             c1p = max(cdp+cap, zero)
             c10 = c1m + c1p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, c1 and d1.

             dw(i,j,k,itu1) = dw(i,j,k,itu1)      + c1m*w(i-1,j,k,itu1) &
                  - c10*w(i,j,k,itu1) + c1p*w(i+1,j,k,itu1)
          enddo
       enddo
    enddo
  end subroutine saViscous

  subroutine saAdvection
    ! ---------------------------------------------
    !                SA Advection
    ! ---------------------------------------------
    use constants
    use inputDiscretization, only : orderTurb
    use iteration, only : groundlevel
    use turbMod, only : secondOrd
    implicit none

    ! Variables for sa Advection
    real(kind=realType) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk, qs
    real(kind=realType) :: voli, xa, ya, za
    integer(kind=intType), parameter :: nAdv=1
    integer(kind=intType) :: offset, i, j, k, ii, jj

    ! Determine whether or not a second order discretization for the
    ! advective terms must be used.
    secondOrd = .false.
    if(groundLevel == 1_intType .and. &
         orderTurb == secondOrder) secondOrd = .true.

    offset=itu1-1
    do k=2, kl
       do j=2, jl
          do i=2, il

             ! Compute the grid velocity if present.
             ! It is taken as the average of k and k-1,

             voli = half/vol(i,j,k)
             qs = (sFaceK(i,j,k) + sFaceK(i,j,k-1))*voli

             ! Compute the normal velocity, where the normal direction
             ! is taken as the average of faces k and k-1.

             xa = (sk(i,j,k,1) + sk(i,j,k-1,1))*voli
             ya = (sk(i,j,k,2) + sk(i,j,k-1,2))*voli
             za = (sk(i,j,k,3) + sk(i,j,k-1,3))*voli

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

             ! Determine the situation we are having here, i.e. positive
             ! or negative normal velocity.

             velKdir: if(uu > zero) then

                ! Velocity has a component in positive k-direction.
                ! Loop over the number of advection equations.

                do ii=1,nAdv

                   ! Set the value of jj such that it corresponds to the
                   ! turbulent entry in w.

                   jj = ii + offset

                   ! Check whether a first or a second order discretization
                   ! must be used.

                   if( secondOrd ) then

                      ! Second order; store the three differences for the
                      ! discretization of the derivative in k-direction.

                      dwtm1 = w(i,j,k-1,jj) - w(i,j,k-2,jj)
                      dwt   = w(i,j,k,  jj) - w(i,j,k-1,jj)
                      dwtp1 = w(i,j,k+1,jj) - w(i,j,k,  jj)

                      ! Construct the derivative in this cell center. This
                      ! is the first order upwind derivative with two
                      ! nonlinear corrections.

                      dwtk = dwt

                      if(dwt*dwtp1 > zero) then
                         if(abs(dwt) < abs(dwtp1)) then
                            dwtk = dwtk + half*dwt
                         else
                            dwtk = dwtk + half*dwtp1
                         endif
                      endif

                      if(dwt*dwtm1 > zero) then
                         if(abs(dwt) < abs(dwtm1)) then
                            dwtk = dwtk - half*dwt
                         else
                            dwtk = dwtk - half*dwtm1
                         endif
                      endif

                   else

                      ! 1st order upwind scheme.

                      dwtk = w(i,j,k,jj) - w(i,j,k-1,jj)

                   endif

                   ! Update the residual. The convective term must be
                   ! substracted, because it appears on the other side of
                   ! the equation as the source and viscous terms.

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtk
                enddo

             else velKdir

                ! Velocity has a component in negative k-direction.
                ! Loop over the number of advection equations
                do ii=1,nAdv

                   ! Set the value of jj such that it corresponds to the
                   ! turbulent entry in w.

                   jj = ii + offset

                   ! Check whether a first or a second order discretization
                   ! must be used.

                   if( secondOrd ) then

                      ! Store the three differences for the discretization of
                      ! the derivative in k-direction.

                      dwtm1 = w(i,j,k,  jj) - w(i,j,k-1,jj)
                      dwt   = w(i,j,k+1,jj) - w(i,j,k,  jj)
                      dwtp1 = w(i,j,k+2,jj) - w(i,j,k+1,jj)

                      ! Construct the derivative in this cell center. This is
                      ! the first order upwind derivative with two nonlinear
                      ! corrections.

                      dwtk = dwt

                      if(dwt*dwtp1 > zero) then
                         if(abs(dwt) < abs(dwtp1)) then
                            dwtk = dwtk - half*dwt
                         else
                            dwtk = dwtk - half*dwtp1
                         endif
                      endif

                      if(dwt*dwtm1 > zero) then
                         if(abs(dwt) < abs(dwtm1)) then
                            dwtk = dwtk + half*dwt
                         else
                            dwtk = dwtk + half*dwtm1
                         endif
                      endif

                   else

                      ! 1st order upwind scheme.

                      dwtk = w(i,j,k+1,jj) - w(i,j,k,jj)

                   endif

                   ! Update the residual. The convective term must be
                   ! substracted, because it appears on the other side
                   ! of the equation as the source and viscous terms.

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtk
                end do
             endif velKdir
          enddo
       enddo
    enddo

    !
    !       Upwind discretization of the convective term in j (eta)
    !       direction. Either the 1st order upwind or the second order
    !       fully upwind interpolation scheme, kappa = -1, is used in
    !       combination with the minmod limiter.
    !       The possible grid velocity must be taken into account.
    !
    do k=2, kl
       do j=2, jl
          do i=2, il


             ! Compute the grid velocity if present.
             ! It is taken as the average of j and j-1,

             voli = half/vol(i,j,k)
             qs = (sFaceJ(i,j,k) + sFaceJ(i,j-1,k))*voli

             ! Compute the normal velocity, where the normal direction
             ! is taken as the average of faces j and j-1.

             xa = (sj(i,j,k,1) + sj(i,j-1,k,1))*voli
             ya = (sj(i,j,k,2) + sj(i,j-1,k,2))*voli
             za = (sj(i,j,k,3) + sj(i,j-1,k,3))*voli

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

             ! Determine the situation we are having here, i.e. positive
             ! or negative normal velocity.

             velJdir: if(uu > zero) then

                ! Velocity has a component in positive j-direction.
                ! Loop over the number of advection equations.
                do ii=1,nAdv

                   ! Set the value of jj such that it corresponds to the
                   ! turbulent entry in w.

                   jj = ii + offset

                   ! Check whether a first or a second order discretization
                   ! must be used.

                   if( secondOrd ) then

                      ! Second order; store the three differences for the
                      ! discretization of the derivative in j-direction.

                      dwtm1 = w(i,j-1,k,jj) - w(i,j-2,k,jj)
                      dwt   = w(i,j,  k,jj) - w(i,j-1,k,jj)
                      dwtp1 = w(i,j+1,k,jj) - w(i,j,  k,jj)

                      ! Construct the derivative in this cell center. This is
                      ! the first order upwind derivative with two nonlinear
                      ! corrections.

                      dwtj = dwt

                      if(dwt*dwtp1 > zero) then
                         if(abs(dwt) < abs(dwtp1)) then
                            dwtj = dwtj + half*dwt
                         else
                            dwtj = dwtj + half*dwtp1
                         endif
                      endif

                      if(dwt*dwtm1 > zero) then
                         if(abs(dwt) < abs(dwtm1)) then
                            dwtj = dwtj - half*dwt
                         else
                            dwtj = dwtj - half*dwtm1
                         endif
                      endif

                   else

                      ! 1st order upwind scheme.

                      dwtj = w(i,j,k,jj) - w(i,j-1,k,jj)

                   endif

                   ! Update the residual. The convective term must be
                   ! substracted, because it appears on the other side of
                   ! the equation as the source and viscous terms.

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtj
                enddo

             else velJdir

                ! Velocity has a component in negative j-direction.
                ! Loop over the number of advection equations.
                do ii=1,nAdv

                   ! Set the value of jj such that it corresponds to the
                   ! turbulent entry in w.

                   jj = ii + offset

                   ! Check whether a first or a second order discretization
                   ! must be used.

                   if( secondOrd ) then

                      ! Store the three differences for the discretization of
                      ! the derivative in j-direction.

                      dwtm1 = w(i,j,  k,jj) - w(i,j-1,k,jj)
                      dwt   = w(i,j+1,k,jj) - w(i,j,  k,jj)
                      dwtp1 = w(i,j+2,k,jj) - w(i,j+1,k,jj)

                      ! Construct the derivative in this cell center. This is
                      ! the first order upwind derivative with two nonlinear
                      ! corrections.

                      dwtj = dwt

                      if(dwt*dwtp1 > zero) then
                         if(abs(dwt) < abs(dwtp1)) then
                            dwtj = dwtj - half*dwt
                         else
                            dwtj = dwtj - half*dwtp1
                         endif
                      endif

                      if(dwt*dwtm1 > zero) then
                         if(abs(dwt) < abs(dwtm1)) then
                            dwtj = dwtj + half*dwt
                         else
                            dwtj = dwtj + half*dwtm1
                         endif
                      endif

                   else

                      ! 1st order upwind scheme.

                      dwtj = w(i,j+1,k,jj) - w(i,j,k,jj)

                   endif

                   ! Update the residual. The convective term must be
                   ! substracted, because it appears on the other side
                   ! of the equation as the source and viscous terms.

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtj
                enddo
             endif velJdir
          enddo
       enddo
    enddo
    !
    !       Upwind discretization of the convective term in i (xi)
    !       direction. Either the 1st order upwind or the second order
    !       fully upwind interpolation scheme, kappa = -1, is used in
    !       combination with the minmod limiter.
    !       The possible grid velocity must be taken into account.
    !
    qs = zero
    do k=2, kl
       do j=2, jl
          do i=2, il
             ! Compute the grid velocity if present.
             ! It is taken as the average of i and i-1,

             voli = half/vol(i,j,k)
             qs = (sFaceI(i,j,k) + sFaceI(i-1,j,k))*voli

             ! Compute the normal velocity, where the normal direction
             ! is taken as the average of faces i and i-1.

             xa = (si(i,j,k,1) + si(i-1,j,k,1))*voli
             ya = (si(i,j,k,2) + si(i-1,j,k,2))*voli
             za = (si(i,j,k,3) + si(i-1,j,k,3))*voli

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

             ! Determine the situation we are having here, i.e. positive
             ! or negative normal velocity.

             velIdir: if(uu > zero) then

                ! Velocity has a component in positive i-direction.
                ! Loop over the number of advection equations.
                do ii=1,nAdv

                   ! Set the value of jj such that it corresponds to the
                   ! turbulent entry in w.

                   jj = ii + offset

                   ! Check whether a first or a second order discretization
                   ! must be used.

                   if( secondOrd ) then

                      ! Second order; store the three differences for the
                      ! discretization of the derivative in i-direction.

                      dwtm1 = w(i-1,j,k,jj) - w(i-2,j,k,jj)
                      dwt   = w(i,  j,k,jj) - w(i-1,j,k,jj)
                      dwtp1 = w(i+1,j,k,jj) - w(i,  j,k,jj)

                      ! Construct the derivative in this cell center. This is
                      ! the first order upwind derivative with two nonlinear
                      ! corrections.

                      dwti = dwt

                      if(dwt*dwtp1 > zero) then
                         if(abs(dwt) < abs(dwtp1)) then
                            dwti = dwti + half*dwt
                         else
                            dwti = dwti + half*dwtp1
                         endif
                      endif

                      if(dwt*dwtm1 > zero) then
                         if(abs(dwt) < abs(dwtm1)) then
                            dwti = dwti - half*dwt
                         else
                            dwti = dwti - half*dwtm1
                         endif
                      endif

                   else

                      ! 1st order upwind scheme.

                      dwti = w(i,j,k,jj) - w(i-1,j,k,jj)

                   endif

                   ! Update the residual. The convective term must be
                   ! substracted, because it appears on the other side of
                   ! the equation as the source and viscous terms.

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwti
                enddo

             else velIdir

                ! Velocity has a component in negative i-direction.
                ! Loop over the number of advection equations.
                do ii=1,nAdv

                   ! Set the value of jj such that it corresponds to the
                   ! turbulent entry in w.

                   jj = ii + offset

                   ! Check whether a first or a second order discretization
                   ! must be used.

                   if( secondOrd ) then

                      ! Second order; store the three differences for the
                      ! discretization of the derivative in i-direction.

                      dwtm1 = w(i,  j,k,jj) - w(i-1,j,k,jj)
                      dwt   = w(i+1,j,k,jj) - w(i,  j,k,jj)
                      dwtp1 = w(i+2,j,k,jj) - w(i+1,j,k,jj)

                      ! Construct the derivative in this cell center. This is
                      ! the first order upwind derivative with two nonlinear
                      ! corrections.

                      dwti = dwt

                      if(dwt*dwtp1 > zero) then
                         if(abs(dwt) < abs(dwtp1)) then
                            dwti = dwti - half*dwt
                         else
                            dwti = dwti - half*dwtp1
                         endif
                      endif

                      if(dwt*dwtm1 > zero) then
                         if(abs(dwt) < abs(dwtm1)) then
                            dwti = dwti + half*dwt
                         else
                            dwti = dwti + half*dwtm1
                         endif
                      endif

                   else

                      ! 1st order upwind scheme.

                      dwti = w(i+1,j,k,jj) - w(i,j,k,jj)

                   endif

                   ! Update the residual. The convective term must be
                   ! substracted, because it appears on the other side
                   ! of the equation as the source and viscous terms.

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwti

                   ! Update the central jacobian. First the term which is
                   ! always present, i.e. -uu.
                enddo

             endif velIdir
          enddo
       enddo
    enddo
  end subroutine saAdvection

  subroutine saResScale

    !
    !  Multiply the residual by the volume and store this in dw; this
    ! * is done for monitoring reasons only. The multiplication with the
    ! * volume is present to be consistent with the flow residuals; also
    !  the negative value is taken, again to be consistent with the
    ! * flow equations. Also multiply by iblank so that no updates occur
    !  in holes or the overset boundary.
    use constants
    implicit none

    ! Local variables
    integer(kind=intType) :: i,j,k,ii
    real(kind=realType) :: rblank

    do k=2, kl
       do j=2, jl
          do i=2, il
             rblank = max(real(iblank(i,j,k), realType), zero)
             dw(i,j,k,itu1) = -volRef(i,j,k)*dw(i,j,k,itu1)*rblank
          enddo
       enddo
    enddo

  end subroutine saResScale

  subroutine timeStep(updateDtl)
    ! ---------------------------------------------
    !               Spectral Radius
    ! ---------------------------------------------

    use constants
    use blockPointers, only : sectionID
    use flowvarRefState, only : pInfCorr, rhoInf, gammaInf, viscous, timeRef
    use inputPhysics, only : equationMode
    use inputDiscretization, only : adis
    use section, only : sections
    use inputTimeSpectral, only : nTimeIntervalsSpectral

    implicit none

    ! Input
    logical, intent(in), optional :: updateDtl

    ! Local parameters.
    real(kind=realType), parameter :: b = 2.0_realType

    ! Variables for spectral Radius
    real(kind=realType) :: plim, rlim, clim2
    real(kind=realType) :: cc2, qsi, qsj, qsk, sx, sy, sz, rmu
    real(kind=realType) :: ri, rj, rk, rij, rjk, rki
    real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
    real(kind=realType) :: sFace, tmp, uux, uuy, uuz
    logical :: doScaling, updateDt
    integer(kind=intType) :: i, j, k

    updateDt = .False.
    if (present(updateDtl)) then
       updateDt = .True.
    end if

    ! Set the value of plim. To be fully consistent this must have
    ! the dimension of a pressure. Therefore a fraction of pInfCorr
    ! is used. Idem for rlim; compute clim2 as well.

    plim  = 0.001_realType*pInfCorr
    rlim  = 0.001_realType*rhoInf
    clim2 = 0.000001_realType*gammaInf*pInfCorr/rhoInf
    doScaling = .True.

    ! Initialize sFace to zero. This value will be used if the
    ! block is not moving.

    sFace = zero
    !
    !           Inviscid contribution, depending on the preconditioner.
    !           Compute the cell centered values of the spectral radii.
    !
    ! Note:  DON'T change the ranges for i. It will mess up dtl.
    ! we don't copy the spectral-radii and dtl to keep the
    ! code simple, therefore this loop needs the full single
    ! halo range.
    do k=1,ke
       do j=1,je
          do i=1,ie

             ! Compute the velocities and speed of sound squared.

             uux  = w(i,j,k,ivx)
             uuy  = w(i,j,k,ivy)
             uuz  = w(i,j,k,ivz)
             cc2 = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
             cc2 = max(cc2,clim2)

             ! Set the dot product of the grid velocity and the
             ! normal in i-direction for a moving face. To avoid
             ! a number of multiplications by 0.5 simply the sum
             ! is taken.

             sFace = sFaceI(i-1,j,k) + sFaceI(i,j,k)

             ! Spectral radius in i-direction.

             sx = si(i-1,j,k,1) + si(i,j,k,1)
             sy = si(i-1,j,k,2) + si(i,j,k,2)
             sz = si(i-1,j,k,3) + si(i,j,k,3)

             qsi = uux*sx + uuy*sy + uuz*sz - sFace

             ri = half*(abs(qsi) &
                  +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

             ! The grid velocity in j-direction.
             sFace = sFaceJ(i,j-1,k) + sFaceJ(i,j,k)

             ! Spectral radius in j-direction.

             sx = sj(i,j-1,k,1) + sj(i,j,k,1)
             sy = sj(i,j-1,k,2) + sj(i,j,k,2)
             sz = sj(i,j-1,k,3) + sj(i,j,k,3)

             qsj = uux*sx + uuy*sy + uuz*sz - sFace

             rj = half*(abs(qsj) &
                  +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

             ! The grid velocity in k-direction.
             sFace = sFaceK(i,j,k-1) + sFaceK(i,j,k)

             ! Spectral radius in k-direction.

             sx = sk(i,j,k-1,1) + sk(i,j,k,1)
             sy = sk(i,j,k-1,2) + sk(i,j,k,2)
             sz = sk(i,j,k-1,3) + sk(i,j,k,3)

             qsk = uux*sx + uuy*sy + uuz*sz - sFace

             rk = half*(abs(qsk) &
                  +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

             ! Store in tdl if required
             if (updateDt) then
                dtl(i,j,k) = ri + rj + rk
             end if

             ! Avoid division by zero by clipping radi, radJ and
             ! radK.

             ri = max(ri, eps)
             rj = max(rj, eps)
             rk = max(rk, eps)

             ! Compute the scaling in the three coordinate
             ! directions.

             rij = (ri/rj)**adis
             rjk = (rj/rk)**adis
             rki = (rk/ri)**adis

             ! Create the scaled versions of the aspect ratios.
             ! Note that the multiplication is done with radi, radJ
             ! and radK, such that the influence of the clipping
             ! is negligible.

             radi(i,j,k) = ri*(one + one/rij + rki)
             radJ(i,j,k) = rj*(one + one/rjk + rij)
             radK(i,j,k) = rk*(one + one/rki + rjk)
          end do
       enddo
    enddo


    ! The rest is only necessary if the timeStep needs to be computed
    if (updateDt) then

       viscousTerm: if( viscous ) then

          ! Loop over the owned cell centers.

          do k=2,kl
             do j=2,jl
                do i=2,il

                   ! Compute the effective viscosity coefficient. The
                   ! factor 0.5 is a combination of two things. In the
                   ! standard central discretization of a second
                   ! derivative there is a factor 2 multiplying the
                   ! central node. However in the code below not the
                   ! average but the sum of the left and the right face
                   ! is taken and squared. This leads to a factor 4.
                   ! Combining both effects leads to 0.5. Furthermore,
                   ! it is divided by the volume and density to obtain
                   ! the correct dimensions and multiplied by the
                   ! non-dimensional factor factVis.

                   rmu = rlv(i,j,k)
                   rmu = rmu + rev(i,j,k)
                   rmu = half*rmu/(w(i,j,k,irho)*vol(i,j,k))

                   ! Add the viscous contribution in i-direction to the
                   ! (inverse) of the time step.

                   sx = si(i,j,k,1) + si(i-1,j,k,1)
                   sy = si(i,j,k,2) + si(i-1,j,k,2)
                   sz = si(i,j,k,3) + si(i-1,j,k,3)

                   vsi        = rmu*(sx*sx + sy*sy + sz*sz)
                   dtl(i,j,k) = dtl(i,j,k) + vsi

                   ! Add the viscous contribution in j-direction to the
                   ! (inverse) of the time step.

                   sx = sj(i,j,k,1) + sj(i,j-1,k,1)
                   sy = sj(i,j,k,2) + sj(i,j-1,k,2)
                   sz = sj(i,j,k,3) + sj(i,j-1,k,3)

                   vsj        = rmu*(sx*sx + sy*sy + sz*sz)
                   dtl(i,j,k) = dtl(i,j,k) + vsj

                   ! Add the viscous contribution in k-direction to the
                   ! (inverse) of the time step.

                   sx = sk(i,j,k,1) + sk(i,j,k-1,1)
                   sy = sk(i,j,k,2) + sk(i,j,k-1,2)
                   sz = sk(i,j,k,3) + sk(i,j,k-1,3)

                   vsk        = rmu*(sx*sx + sy*sy + sz*sz)
                   dtl(i,j,k) = dtl(i,j,k) + vsk

                enddo
             enddo
          enddo
       endif viscousTerm


       ! For the spectral mode an additional term term must be
       ! taken into account, which corresponds to the contribution
       ! of the highest frequency.

       if(equationMode == timeSpectral) then

          tmp = nTimeIntervalsSpectral*pi*timeRef &
               / sections(sectionID)%timePeriod

          ! Loop over the owned cell centers and add the term.

          do k=2,kl
             do j=2,jl
                do i=2,il
                   dtl(i,j,k) = dtl(i,j,k) + tmp*vol(i,j,k)
                enddo
             enddo
          enddo

       endif

       ! Currently the inverse of dt/vol is stored in dtl. Invert
       ! this value such that the time step per unit cfl number is
       ! stored and correct in cases of high gradients.

       do k=2,kl
          do j=2,jl
             do i=2,il
                dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                     /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
                dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                     /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
                dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                     /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
                rfl = one/(one + b*(dpi  +dpj  +dpk))

                dtl(i,j,k) = rfl/dtl(i,j,k)
             enddo
          enddo
       enddo
    end if


  end subroutine timeStep

  subroutine inviscidCentralFlux

    ! ---------------------------------------------
    !               Inviscid central flux
    ! ---------------------------------------------
    use constants
    use blockPointers, only : blockIsMoving, nBkGlobal
    use flowVarRefState, only : timeRef
    use cgnsGrid, only : cgnsDoms
    use inputPhysics, only : equationMode
    implicit none

    ! Variables for inviscid central flux
    real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
    real(kind=realType) :: pa, vnp, vnm, fs, sFace
    integer(kind=intType) :: i, j, k
    real(kind=realType) :: wwx, wwy, wwz, rvol

    do k=2, kl
       do j=2, jl
          do i=1, il

             ! Set the dot product of the grid velocity and the
             ! normal in i-direction for a moving face.

             sFace = sFaceI(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i+1,j,k,ivx)*sI(i,j,k,1) &
                  + w(i+1,j,k,ivy)*sI(i,j,k,2) &
                  + w(i+1,j,k,ivz)*sI(i,j,k,3)
             vnm = w(i,  j,k,ivx)*sI(i,j,k,1) &
                  + w(i,  j,k,ivy)*sI(i,j,k,2) &
                  + w(i,  j,k,ivz)*sI(i,j,k,3)
             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half
             if(porI(i,j,k) == noFlux)    porFlux = zero
             if(porI(i,j,k) == boundFlux) then
                porVel = zero
                vnp    = sFace
                vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities relative to the grid for
             ! the face as well as the mass fluxes.

             qsp = (vnp -sFace)*porVel
             qsm = (vnm -sFace)*porVel

             rqsp = qsp*w(i+1,j,k,irho)
             rqsm = qsm*w(i,  j,k,irho)

             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i+1,j,k) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i+1,j,k. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i+1,j,k,irho) = dw(i+1,j,k,irho) - fs
             dw(i,  j,k,irho) = dw(i,  j,k,irho) + fs

             fs = rqsp*w(i+1,j,k,ivx) + rqsm*w(i,j,k,ivx) &
                  + pa*sI(i,j,k,1)
             dw(i+1,j,k,imx) = dw(i+1,j,k,imx) - fs
             dw(i,  j,k,imx) = dw(i,  j,k,imx) + fs

             fs = rqsp*w(i+1,j,k,ivy) + rqsm*w(i,j,k,ivy) &
                  + pa*sI(i,j,k,2)
             dw(i+1,j,k,imy) = dw(i+1,j,k,imy) - fs
             dw(i,  j,k,imy) = dw(i,  j,k,imy) + fs

             fs = rqsp*w(i+1,j,k,ivz) + rqsm*w(i,j,k,ivz) &
                  + pa*sI(i,j,k,3)
             dw(i+1,j,k,imz) = dw(i+1,j,k,imz) - fs
             dw(i,  j,k,imz) = dw(i,  j,k,imz) + fs

             fs = qsp*w(i+1,j,k,irhoE) + qsm*w(i,j,k,irhoE) &
                  + porFlux*(vnp*p(i+1,j,k) + vnm*p(i,j,k))
             dw(i+1,j,k,irhoE) = dw(i+1,j,k,irhoE) - fs
             dw(i,  j,k,irhoE) = dw(i,  j,k,irhoE) + fs
          enddo
       enddo
    enddo

    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Set the dot product of the grid velocity and the
             ! normal in j-direction for a moving face.

             sFace = sFaceJ(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i,j+1,k,ivx)*sJ(i,j,k,1) &
                  + w(i,j+1,k,ivy)*sJ(i,j,k,2) &
                  + w(i,j+1,k,ivz)*sJ(i,j,k,3)
             vnm = w(i,j,  k,ivx)*sJ(i,j,k,1) &
                  + w(i,j,  k,ivy)*sJ(i,j,k,2) &
                  + w(i,j,  k,ivz)*sJ(i,j,k,3)

             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half
             if(porJ(i,j,k) == noFlux)    porFlux = zero
             if(porJ(i,j,k) == boundFlux) then
                porVel = zero
                vnp    = sFace
                vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities for the face as well as the
             ! mass fluxes.

             qsp = (vnp - sFace)*porVel
             qsm = (vnm - sFace)*porVel

             rqsp = qsp*w(i,j+1,k,irho)
             rqsm = qsm*w(i,j,  k,irho)

             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i,j+1,k) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i,j+1,k. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i,j+1,k,irho) = dw(i,j+1,k,irho) - fs
             dw(i,j,  k,irho) = dw(i,j,  k,irho) + fs

             fs = rqsp*w(i,j+1,k,ivx) + rqsm*w(i,j,k,ivx) &
                  + pa*sJ(i,j,k,1)
             dw(i,j+1,k,imx) = dw(i,j+1,k,imx) - fs
             dw(i,j,  k,imx) = dw(i,j,  k,imx) + fs

             fs = rqsp*w(i,j+1,k,ivy) + rqsm*w(i,j,k,ivy) &
                  + pa*sJ(i,j,k,2)
             dw(i,j+1,k,imy) = dw(i,j+1,k,imy) - fs
             dw(i,j,  k,imy) = dw(i,j,  k,imy) + fs

             fs = rqsp*w(i,j+1,k,ivz) + rqsm*w(i,j,k,ivz) &
                  + pa*sJ(i,j,k,3)
             dw(i,j+1,k,imz) = dw(i,j+1,k,imz) - fs
             dw(i,j,  k,imz) = dw(i,j,  k,imz) + fs

             fs = qsp*w(i,j+1,k,irhoE) + qsm*w(i,j,k,irhoE) &
                  + porFlux*(vnp*p(i,j+1,k) + vnm*p(i,j,k))
             dw(i,j+1,k,irhoE) = dw(i,j+1,k,irhoE) - fs
             dw(i,j,  k,irhoE) = dw(i,j,  k,irhoE) + fs
          enddo
       enddo
    enddo

    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Set the dot product of the grid velocity and the
             ! normal in k-direction for a moving face.

             sFace = sFaceK(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i,j,k+1,ivx)*sK(i,j,k,1) &
                  + w(i,j,k+1,ivy)*sK(i,j,k,2) &
                  + w(i,j,k+1,ivz)*sK(i,j,k,3)
             vnm = w(i,j,k,  ivx)*sK(i,j,k,1) &
                  + w(i,j,k,  ivy)*sK(i,j,k,2) &
                  + w(i,j,k,  ivz)*sK(i,j,k,3)

             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! block boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half

             if(porK(i,j,k) == noFlux)    porFlux = zero
             if(porK(i,j,k) == boundFlux) then
                porVel = zero
                vnp    = sFace
                vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities for the face as well as the
             ! mass fluxes.

             qsp = (vnp - sFace)*porVel
             qsm = (vnm - sFace)*porVel

             rqsp = qsp*w(i,j,k+1,irho)
             rqsm = qsm*w(i,j,k,  irho)

             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i,j,k+1) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i,j,k+1. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i,j,k+1,irho) = dw(i,j,k+1,irho) - fs
             dw(i,j,k,  irho) = dw(i,j,k,  irho) + fs

             fs = rqsp*w(i,j,k+1,ivx) + rqsm*w(i,j,k,ivx) &
                  + pa*sK(i,j,k,1)
             dw(i,j,k+1,imx) = dw(i,j,k+1,imx) - fs
             dw(i,j,k,  imx) = dw(i,j,k,  imx) + fs

             fs = rqsp*w(i,j,k+1,ivy) + rqsm*w(i,j,k,ivy) &
                  + pa*sK(i,j,k,2)
             dw(i,j,k+1,imy) = dw(i,j,k+1,imy) - fs
             dw(i,j,k,  imy) = dw(i,j,k,  imy) + fs

             fs = rqsp*w(i,j,k+1,ivz) + rqsm*w(i,j,k,ivz) &
                  + pa*sK(i,j,k,3)
             dw(i,j,k+1,imz) = dw(i,j,k+1,imz) - fs
             dw(i,j,k,  imz) = dw(i,j,k,  imz) + fs

             fs = qsp*w(i,j,k+1,irhoE) + qsm*w(i,j,k,irhoE) &
                  + porFlux*(vnp*p(i,j,k+1) + vnm*p(i,j,k))
             dw(i,j,k+1,irhoE) = dw(i,j,k+1,irhoE) - fs
             dw(i,j,k,  irhoE) = dw(i,j,k,  irhoE) + fs

          enddo
       enddo
    enddo

    rotation: if(blockIsMoving .and. equationMode == steady) then

       ! Compute the three nonDimensional angular velocities.

       wwx = timeRef*cgnsDoms(nbkGlobal)%rotRate(1)
       wwy = timeRef*cgnsDoms(nbkGlobal)%rotRate(2)
       wwz = timeRef*cgnsDoms(nbkGlobal)%rotRate(3)

       ! Loop over the internal cells of this block to compute the
       ! rotational terms for the momentum equations.
       do k=2, kl
          do j=2, jl
             do i=2, il
                rvol = w(i, j, k, irho)*vol(i, j, k)
                dw(i,j,k,imx) = dw(i,j,k,imx) &
                     + rvol*(wwy*w(i,j,k,ivz) - wwz*w(i,j,k,ivy))
                dw(i,j,k,imy) = dw(i,j,k,imy) &
                     + rvol*(wwz*w(i,j,k,ivx) - wwx*w(i,j,k,ivz))
                dw(i,j,k,imz) = dw(i,j,k,imz) &
                     + rvol*(wwx*w(i,j,k,ivy) - wwy*w(i,j,k,ivx))
             enddo
          end do
       end do
    endif rotation

  end subroutine inviscidCentralFlux

  subroutine inviscidDissFluxMatrix
    !
    !       inviscidDissFluxMatrix computes the matrix artificial
    !       dissipation term. Instead of the spectral radius, as used in
    !       the scalar dissipation scheme, the absolute value of the flux
    !       jacobian is used. This leads to a less diffusive and
    !       consequently more accurate scheme. It is assumed that the
    !       pointers in blockPointers already point to the correct block.
    !
    use constants
    use flowVarRefState, only : pInfCorr
    use inputDiscretization, only: vis2, vis4
    use inputPhysics, only : equations
    use iteration, only : rFil
    use utils, only : getCorrectForK, myDim
    implicit none
    !
    !      Local parameters.
    !
    real(kind=realType), parameter :: dpMax        = 0.25_realType
    real(kind=realType), parameter :: epsAcoustic  = 0.25_realType
    real(kind=realType), parameter :: epsShear     = 0.025_realType
    real(kind=realType), parameter :: omega        = 0.5_realType
    real(kind=realType), parameter :: oneMinOmega = one - omega
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ind, ii

    real(kind=realType) :: plim, sface
    real(kind=realType) :: sfil, fis2, fis4
    real(kind=realType) :: gammaAvg, gm1, ovgm1, gm53
    real(kind=realType) :: ppor, rrad, dis2, dis4
    real(kind=realType) :: dp1, dp2, tmp, fs
    real(kind=realType) :: ddw1, ddw2, ddw3, ddw4, ddw5, ddw6
    real(kind=realType) :: dr, dru, drv, drw, dre, drk, sx, sy, sz
    real(kind=realType) :: uAvg, vAvg, wAvg, a2Avg, aAvg, hAvg
    real(kind=realType) :: alphaAvg, unAvg, ovaAvg, ova2Avg
    real(kind=realType) :: kAvg, lam1, lam2, lam3, area
    real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
    logical :: correctForK

    ! Set the value of plim. To be fully consistent this must have
    ! the dimension of a pressure. Therefore a fraction of pInfCorr
    ! is used.

    plim = 0.001_realType*pInfCorr

    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.

    correctForK = getCorrectForK()

    ! Initialize sface to zero. This value will be used if the
    ! block is not moving.

    sface = zero

    ! Set a couple of constants for the scheme.

    fis2 = rFil*vis2
    fis4 = rFil*vis4
    sfil = one - rFil

    ! Initialize the dissipative residual to a certain times,
    ! possibly zero, the previously stored value.

    fw = sfil*fw

    ! Compute the pressure sensor for each cell, in each direction:
    do k=1,ke
       do j=1,je
          do i=singleHaloStart,ie
             dss(i,j,k,1) =abs((p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k))        &
                  /     (omega*(p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k)) &
                  +      oneMinOmega*(abs(p(i+1,j,k) - p(i,j,k))      &
                  +                   abs(p(i,j,k) - p(i-1,j,k))) + plim))


             dss(i,j,k,2) =abs((p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k))        &
                  /     (omega*(p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k)) &
                  +      oneMinOmega*(abs(p(i,j+1,k) - p(i,j,k))      &
                  +                   abs(p(i,j,k) - p(i,j-1,k))) + plim))

             dss(i,j,k,3) =  abs((p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1))        &
                  /     (omega*(p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1)) &
                  +      oneMinOmega*(abs(p(i,j,k+1) - p(i,j,k))      &
                  +                   abs(p(i,j,k) - p(i,j,k-1))) + plim))
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the i-direction.
    !
    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = one
             dis2 = ppor*fis2*min(dpMax, max(dss(i,j,k,1), dss(i+1,j,k,1)))
             dis4 = dim(ppor*fis4, dis2)

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw1 = w(i+1,j,k,irho) - w(i,j,k,irho)
             dr  = dis2*ddw1 &
                  - dis4*(w(i+2,j,k,irho) - w(i-1,j,k,irho) - three*ddw1)

             ddw2 = w(i+1,j,k,irho)*w(i+1,j,k,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw2                             &
                  - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,ivx) &
                  -       w(i-1,j,k,irho)*w(i-1,j,k,ivx) - three*ddw2)

             ddw3 = w(i+1,j,k,irho)*w(i+1,j,k,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw3                             &
                  - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,ivy) &
                  -       w(i-1,j,k,irho)*w(i-1,j,k,ivy) - three*ddw3)

             ddw4 = w(i+1,j,k,irho)*w(i+1,j,k,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw4                             &
                  - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,ivz) &
                  -       w(i-1,j,k,irho)*w(i-1,j,k,ivz) - three*ddw4)

             ddw5 = w(i+1,j,k,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw5 &
                  - dis4*(w(i+2,j,k,irhoE) - w(i-1,j,k,irhoE) - three*ddw5)

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.
             drk   = zero
             kAvg = zero

             if( correctForK ) then
                ddw6 = w(i+1,j,k,irho)*w(i+1,j,k,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw6                              &
                     - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,itu1) &
                     -       w(i-1,j,k,irho)*w(i-1,j,k,itu1) - three*ddw6)

                kAvg = half*(w(i,j,k,itu1) + w(i+1,j,k,itu1))
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i+1,j,k) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i+1,j,k,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i+1,j,k,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i+1,j,k,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i+1,j,k)*p(i+1,j,k)/w(i+1,j,k,irho) &
                  +       gamma(i,  j,k)*p(i,  j,k)/w(i,  j,k,irho))

             area = sqrt(si(i,j,k,1)**2 + si(i,j,k,2)**2 + si(i,j,k,3)**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = si(i,j,k,1)*tmp
             sy   = si(i,j,k,2)*tmp
             sz   = si(i,j,k,3)*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             sface = sFaceI(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = max(lam1,epsAcoustic*rrad)*area
             lam2 = max(lam2,epsAcoustic*rrad)*area
             lam3 = max(lam3,epsShear*rrad)*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do
    !
    !       Dissipative fluxes in the j-direction.
    !
    do k=2,kl
       do j=1,jl
          do i=2,il


             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = one

             dis2 = ppor*fis2*min(dpMax, max(dss(i,j,k,2), dss(i,j+1,k,2)))
             dis4 = dim(ppor*fis4, dis2)

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw1 = w(i,j+1,k,irho) - w(i,j,k,irho)
             dr  = dis2*ddw1 &
                  - dis4*(w(i,j+2,k,irho) - w(i,j-1,k,irho) - three*ddw1)

             ddw2 = w(i,j+1,k,irho)*w(i,j+1,k,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw2                             &
                  - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,ivx) &
                  -       w(i,j-1,k,irho)*w(i,j-1,k,ivx) - three*ddw2)

             ddw3 = w(i,j+1,k,irho)*w(i,j+1,k,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw3 &
                  - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,ivy) &
                  -       w(i,j-1,k,irho)*w(i,j-1,k,ivy) - three*ddw3)

             ddw4 = w(i,j+1,k,irho)*w(i,j+1,k,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw4 &
                  - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,ivz) &
                  -       w(i,j-1,k,irho)*w(i,j-1,k,ivz) - three*ddw4)

             ddw5 = w(i,j+1,k,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw5 &
                  - dis4*(w(i,j+2,k,irhoE) - w(i,j-1,k,irhoE) - three*ddw5)

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.
             drk   = zero
             kAvg = zero

             if( correctForK ) then
                ddw6 = w(i,j+1,k,irho)*w(i,j+1,k,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw6                              &
                     - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,itu1) &
                     -       w(i,j-1,k,irho)*w(i,j-1,k,itu1) - three*ddw6)

                kAvg = half*(w(i,j,k,itu1) + w(i,j+1,k,itu1))
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i,j+1,k) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i,j+1,k,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i,j+1,k,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i,j+1,k,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i,j+1,k)*p(i,j+1,k)/w(i,j+1,k,irho) &
                  +       gamma(i,j,  k)*p(i,j,  k)/w(i,j,  k,irho))

             area = sqrt(sj(i,j,k,1)**2 + sj(i,j,k,2)**2 + sj(i,j,k,3)**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = sj(i,j,k,1)*tmp
             sy   = sj(i,j,k,2)*tmp
             sz   = sj(i,j,k,3)*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             sface = sFaceJ(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = max(lam1,epsAcoustic*rrad)*area
             lam2 = max(lam2,epsAcoustic*rrad)*area
             lam3 = max(lam3,epsShear*rrad)*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do
    !
    !       Dissipative fluxes in the k-direction.
    !
    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = one

             dis2 = ppor*fis2*min(dpMax, max(dss(i,j,k,3), dss(i,j,k+1,3)))
             dis4 = dim(ppor*fis4, dis2)

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw1 = w(i,j,k+1,irho) - w(i,j,k,irho)
             dr  = dis2*ddw1 &
                  - dis4*(w(i,j,k+2,irho) - w(i,j,k-1,irho) - three*ddw1)

             ddw2 = w(i,j,k+1,irho)*w(i,j,k+1,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw2                             &
                  - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,ivx) &
                  -       w(i,j,k-1,irho)*w(i,j,k-1,ivx) - three*ddw2)

             ddw3 = w(i,j,k+1,irho)*w(i,j,k+1,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw3 &
                  - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,ivy) &
                  -       w(i,j,k-1,irho)*w(i,j,k-1,ivy) - three*ddw3)

             ddw4 = w(i,j,k+1,irho)*w(i,j,k+1,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw4 &
                  - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,ivz) &
                  -       w(i,j,k-1,irho)*w(i,j,k-1,ivz) - three*ddw4)

             ddw5 = w(i,j,k+1,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw5 &
                  - dis4*(w(i,j,k+2,irhoE) - w(i,j,k-1,irhoE) - three*ddw5)

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.
             drk   = zero
             kAvg = zero

             if( correctForK ) then
                ddw6 = w(i,j,k+1,irho)*w(i,j,k+1,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw6                              &
                     - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,itu1) &
                     -       w(i,j,k-1,irho)*w(i,j,k-1,itu1) - three*ddw6)

                kAvg = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i,j,k+1) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i,j,k+1,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i,j,k+1,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i,j,k+1,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i,j,k+1)*p(i,j,k+1)/w(i,j,k+1,irho) &
                  +       gamma(i,j,k)  *p(i,j,k)  /w(i,j,k,  irho))

             area = sqrt(sk(i,j,k,1)**2 + sk(i,j,k,2)**2 + sk(i,j,k,3)**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = sk(i,j,k,1)*tmp
             sy   = sk(i,j,k,2)*tmp
             sz   = sk(i,j,k,3)*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             sface = sFaceK(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = max(lam1,epsAcoustic*rrad)*area
             lam2 = max(lam2,epsAcoustic*rrad)*area
             lam3 = max(lam3,epsShear*rrad)*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do

  end subroutine inviscidDissFluxMatrix


  subroutine inviscidDissFluxScalar
    ! ---------------------------------------------
    !             Inviscid Diss Flux Scalar
    ! ---------------------------------------------

    use constants
    use flowVarRefState, only : pInfCorr
    use inputDiscretization, only: vis2, vis4
    use inputPhysics, only : equations
    use iteration, only : rFil
    use flowVarRefState, only : gammaInf, pInfCorr, rhoInf
    implicit none

    ! Variables for inviscid diss flux scalar
    real(kind=realType), parameter :: dssMax = 0.25_realType
    real(kind=realType) :: sslim, rhoi
    real(kind=realType) :: sfil, fis2, fis4
    real(kind=realType) :: ppor, rrad, dis2, dis4, fs
    real(kind=realType) :: ddw1,ddw2,ddw3,ddw4,ddw5
    integer(kind=intType) :: i, j, k

    ! Determine the variables used to compute the switch.
    ! For the inviscid case this is the pressure; for the viscous
    ! case it is the entropy.

    select case (equations)
    case (EulerEquations)

       ! Inviscid case. Pressure switch is based on the pressure.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of pressure and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr

       ! Copy the pressure in ss. Only need the entries used in the
       ! discretization, i.e. not including the corner halo's, but we'll
       ! just copy all anyway.

       ss = P
       !===============================================================

    case (NSEquations, RANSEquations)

       ! Viscous case. Pressure switch is based on the entropy.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of entropy and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

       ! Store the entropy in ss. See above.
       do k=0, kb
          do j=0, jb
             do i=doubleHaloStart, ib
                ss(i,j,k) = p(i,j,k)/(w(i,j,k,irho)**gamma(i,j,k))
             end do
          end do
       end do
    end select

    ! Compute the pressure sensor for each cell, in each direction:
    do k=1,ke
       do j=1,je
          do i=singleHaloStart,ie
             dss(i,j,k,1) = abs((ss(i+1,j,k) - two*ss(i,j,k) + ss(i-1,j,k)) &
                  /     (ss(i+1,j,k) + two*ss(i,j,k) + ss(i-1,j,k) + sslim))

             dss(i,j,k,2) = abs((ss(i,j+1,k) - two*ss(i,j,k) + ss(i,j-1,k)) &
                  /     (ss(i,j+1,k) + two*ss(i,j,k) + ss(i,j-1,k) + sslim))

             dss(i,j,k,3) = abs((ss(i,j,k+1) - two*ss(i,j,k) + ss(i,j,k-1)) &
                  /     (ss(i,j,k+1) + two*ss(i,j,k) + ss(i,j,k-1) + sslim))
          end do
       end do
    end do

    ! Set a couple of constants for the scheme.

    fis2 = rFil*vis2
    fis4 = rFil*vis4
    sfil = one - rFil

    ! Initialize the dissipative residual to a certain times,
    ! possibly zero, the previously stored value. Owned cells
    ! only, because the halo values do not matter.

    fw = sfil*fw
    !
    !       Dissipative fluxes in the i-direction.
    !
    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,1), dss(i+1,j,k,1)))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw1 = w(i+1,j,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw1 &
                  - dis4*(w(i+2,j,k,irho) - w(i-1,j,k,irho) - three*ddw1)

             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw2 = w(i+1,j,k,ivx)*w(i+1,j,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw2 &
                  - dis4*(w(i+2,j,k,ivx)*w(i+2,j,k,irho) - w(i-1,j,k,ivx)*w(i-1,j,k,irho) - three*ddw2)

             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw3 = w(i+1,j,k,ivy)*w(i+1,j,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw3 &
                  - dis4*(w(i+2,j,k,ivy)*w(i+2,j,k,irho) - w(i-1,j,k,ivy)*w(i-1,j,k,irho) - three*ddw3)

             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw4 = w(i+1,j,k,ivz)*w(i+1,j,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw4 &
                  - dis4*(w(i+2,j,k,ivz)*w(i+2,j,k,irho) - w(i-1,j,k,ivz)*w(i-1,j,k,irho) - three*ddw4)

             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw5 = (w(i+1,j,k,irhoE) + P(i+1,j,K))- (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw5 &
                  - dis4*((w(i+2,j,k,irhoE) + P(i+2,j,k)) - (w(i-1,j,k,irhoE) + P(i-1,j,k)) - three*ddw5)

             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the j-direction.
    !
    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,2),dss(i,j+1,k,2)))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw1 = w(i,j+1,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw1 &
                  - dis4*(w(i,j+2,k,irho) - w(i,j-1,k,irho) - three*ddw1)

             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw2 = w(i,j+1,k,ivx)*w(i,j+1,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw2 &
                  - dis4*(w(i,j+2,k,ivx)*w(i,j+2,k,irho) - w(i,j-1,k,ivx)*w(i,j-1,k,irho) - three*ddw2)

             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw3 = w(i,j+1,k,ivy)*w(i,j+1,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw3 &
                  - dis4*(w(i,j+2,k,ivy)*w(i,j+2,k,irho) - w(i,j-1,k,ivy)*w(i,j-1,k,irho) - three*ddw3)

             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw4 = w(i,j+1,k,ivz)*w(i,j+1,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw4 &
                  - dis4*(w(i,j+2,k,ivz)*w(i,j+2,k,irho) - w(i,j-1,k,ivz)*w(i,j-1,k,irho) - three*ddw4)

             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw5 = (w(i,j+1,k,irhoE) + P(i,j+1,k)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw5 &
                  - dis4*((w(i,j+2,k,irhoE) + P(i,j+2,k)) - (w(i,j-1,k,irhoE) + P(i,j-1,k)) - three*ddw5)

             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the k-direction.
    !
    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,3), dss(i,j,k+1,3)))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw1 = w(i,j,k+1,irho) - w(i,j,k,irho)
             fs  = dis2*ddw1 &
                  - dis4*(w(i,j,k+2,irho) - w(i,j,k-1,irho) - three*ddw1)

             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw2 = w(i,j,k+1,ivx)*w(i,j,k+1,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw2 &
                  - dis4*(w(i,j,k+2,ivx)*w(i,j,k+2,irho) - w(i,j,k-1,ivx)*w(i,j,k-1,irho) - three*ddw2)

             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw3 = w(i,j,k+1,ivy)*w(i,j,k+1,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw3 &
                  - dis4*(w(i,j,k+2,ivy)*w(i,j,k+2,irho) - w(i,j,k-1,ivy)*w(i,j,k-1,irho) - three*ddw3)

             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw4 = w(i,j,k+1,ivz)*w(i,j,k+1,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw4 &
                  - dis4*(w(i,j,k+2,ivz)*w(i,j,k+2,irho) - w(i,j,k-1,ivz)*w(i,j,k-1,irho) - three*ddw4)

             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw5 = (w(i,j,k+1,irhoE) + P(i,j,k+1)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw5 &
                  - dis4*((w(i,j,k+2,irhoE) + P(i,j,k+2)) - (w(i,j,k-1,irhoE) + P(i,j,k-1)) - three*ddw5)

             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
  end subroutine inviscidDissFluxScalar

  subroutine inviscidUpwindFlux(fineGrid)
    !
    !       inviscidUpwindFlux computes the artificial dissipation part of
    !       the Euler fluxes by means of an approximate solution of the 1D
    !       Riemann problem on the face. For first order schemes,
    !       fineGrid == .false., the states in the cells are assumed to
    !       be constant; for the second order schemes on the fine grid a
    !       nonlinear reconstruction of the left and right state is done
    !       for which several options exist.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block.
    !
    use constants
    use flowVarRefState, only : kPresent, nw, nwf, rgas, tref
    use blockPointers, only : rotMatrixI, rotMatrixJ, rotMatrixK
    use inputDiscretization, only: limiter, precond, riemann, &
         riemannCoarse, orderTurb, kappaCoef
    use inputPhysics, only : equations
    use iteration, only : rFil, currentLevel, groundLevel
    use utils, only : getCorrectForK, terminate
    use flowUtils, only : eTot
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: fineGrid
    !
    !      Local variables.
    !
    integer(kind=porType) :: por

    integer(kind=intType) :: nwInt
    integer(kind=intType) :: i, j, k, ind
    integer(kind=intType) :: limUsed, riemannUsed

    real(kind=realType) :: sx, sy, sz, omk, opk, sFil, gammaFace
    real(kind=realType) :: factMinmod, sFace

    real(kind=realType), dimension(nw)  :: left, right
    real(kind=realType), dimension(nw)  :: du1, du2, du3
    real(kind=realType), dimension(nwf) :: flux

    logical :: firstOrderK, correctForK, rotationalPeriodic

    ! Check if the formulation for rotational periodic problems
    ! must be used.

    if( associated(rotMatrixI) ) then
       rotationalPeriodic = .true.
    else
       rotationalPeriodic = .false.
    endif

    ! Initialize the dissipative residual to a certain times,
    ! possibly zero, the previously stored value. Owned cells
    ! only, because the halo values do not matter.

    sFil = one - rFil

    do k=2,kl
       do j=2,jl
          do i=2,il
             fw(i,j,k,irho)  = sFil*fw(i,j,k,irho)
             fw(i,j,k,imx)   = sFil*fw(i,j,k,imx)
             fw(i,j,k,imy)   = sFil*fw(i,j,k,imy)
             fw(i,j,k,imz)   = sFil*fw(i,j,k,imz)
             fw(i,j,k,irhoE) = sFil*fw(i,j,k,irhoE)
          enddo
       enddo
    enddo

    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.
    correctForK = getCorrectForK()

    ! Compute the factor used in the minmod limiter.

    factMinmod = (three-kappaCoef) &
         / max(1.e-10_realType, one-kappaCoef)

    ! Determine the limiter scheme to be used. On the fine grid the
    ! user specified scheme is used; on the coarse grid a first order
    ! scheme is computed.

    limUsed = firstOrder
    if( fineGrid ) limUsed = limiter

    ! Determine the riemann solver which must be used.

    riemannUsed = riemannCoarse
    if( fineGrid ) riemannUsed = riemann

    ! Store 1-kappa and 1+kappa a bit easier and multiply it by 0.25.

    omk = fourth*(one - kappaCoef)
    opk = fourth*(one + kappaCoef)

    ! Set the number of variables to be interpolated depending
    ! whether or not a k-equation is present. If a k-equation is
    ! present also set the logical firstOrderK. This indicates
    ! whether or not only a first order approximation is to be used
    ! for the turbulent kinetic energy.

    if( correctForK ) then
       if(orderTurb == firstOrder) then
          nwInt       = nwf
          firstOrderK = .true.
       else
          nwInt       = itu1
          firstOrderK = .false.
       endif
    else
       nwInt       = nwf
       firstOrderK = .false.
    endif
    !
    !       Flux computation. A distinction is made between first and
    !       second order schemes to avoid the overhead for the first order
    !       scheme.
    !
    orderTest: if(limUsed == firstOrder) then
       !
       !         First order reconstruction. The states in the cells are
       !         constant. The left and right states are constructed easily.
       !
       ! Fluxes in the i-direction.

       do k=2,kl
          do j=2,jl
             do i=1,il

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = si(i,j,k,1); sy = si(i,j,k,2); sz = si(i,j,k,3)
                por = porI(i,j,k)
                sFace = sFaceI(i,j,k)

                ! Determine the left and right state.

                left(irho)  = w(i,j,k,irho)
                left(ivx)   = w(i,j,k,ivx)
                left(ivy)   = w(i,j,k,ivy)
                left(ivz)   = w(i,j,k,ivz)
                left(irhoE) = p(i,j,k)
                if( correctForK ) left(itu1) = w(i,j,k,itu1)

                right(irho)  = w(i+1,j,k,irho)
                right(ivx)   = w(i+1,j,k,ivx)
                right(ivy)   = w(i+1,j,k,ivy)
                right(ivz)   = w(i+1,j,k,ivz)
                right(irhoE) = p(i+1,j,k)
                if( correctForK ) right(itu1) = w(i+1,j,k,itu1)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i+1,j,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i+1,j,k,irho)  = fw(i+1,j,k,irho)  - flux(irho)
                fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   - flux(imx)
                fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   - flux(imy)
                fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   - flux(imz)
                fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) - flux(irhoE)

             enddo
          enddo
       enddo

       ! Fluxes in j-direction.

       do k=2,kl
          do j=1,jl
             do i=2,il

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sj(i,j,k,1); sy = sj(i,j,k,2); sz = sj(i,j,k,3)
                por = porJ(i,j,k)
                sFace = sFaceJ(i,j,k)

                ! Determine the left and right state.

                left(irho)  = w(i,j,k,irho)
                left(ivx)   = w(i,j,k,ivx)
                left(ivy)   = w(i,j,k,ivy)
                left(ivz)   = w(i,j,k,ivz)
                left(irhoE) = p(i,j,k)
                if( correctForK ) left(itu1) = w(i,j,k,itu1)

                right(irho)  = w(i,j+1,k,irho)
                right(ivx)   = w(i,j+1,k,ivx)
                right(ivy)   = w(i,j+1,k,ivy)
                right(ivz)   = w(i,j+1,k,ivz)
                right(irhoE) = p(i,j+1,k)
                if( correctForK ) right(itu1) = w(i,j+1,k,itu1)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j+1,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j+1,k,irho)  = fw(i,j+1,k,irho)  - flux(irho)
                fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   - flux(imx)
                fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   - flux(imy)
                fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   - flux(imz)
                fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) - flux(irhoE)
             enddo
          enddo
       enddo

       ! Fluxes in k-direction.

       do k=1,kl
          do j=2,jl
             do i=2,il

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sk(i,j,k,1); sy = sk(i,j,k,2); sz = sk(i,j,k,3)
                por = porK(i,j,k)
                sFace = sFaceK(i,j,k)

                ! Determine the left and right state.

                left(irho)  = w(i,j,k,irho)
                left(ivx)   = w(i,j,k,ivx)
                left(ivy)   = w(i,j,k,ivy)
                left(ivz)   = w(i,j,k,ivz)
                left(irhoE) = p(i,j,k)
                if( correctForK ) left(itu1) = w(i,j,k,itu1)

                right(irho)  = w(i,j,k+1,irho)
                right(ivx)   = w(i,j,k+1,ivx)
                right(ivy)   = w(i,j,k+1,ivy)
                right(ivz)   = w(i,j,k+1,ivz)
                right(irhoE) = p(i,j,k+1)
                if( correctForK ) right(itu1) = w(i,j,k+1,itu1)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j,k+1))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j,k+1,irho)  = fw(i,j,k+1,irho)  - flux(irho)
                fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   - flux(imx)
                fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   - flux(imy)
                fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   - flux(imz)
                fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) - flux(irhoE)

             enddo
          enddo
       enddo

       !      ==================================================================

    else orderTest

       !      ==================================================================
       !
       !         Second order reconstruction of the left and right state.
       !         The three differences used in the, possibly nonlinear,
       !         interpolation are constructed here; the actual left and
       !         right states, or at least the differences from the first
       !         order interpolation, are computed in the subroutine
       !         leftRightState.
       !
       ! Fluxes in the i-direction.

       do k=2,kl
          do j=2,jl
             do i=1,il

                ! Store the three differences used in the interpolation
                ! in du1, du2, du3.

                du1(irho) = w(i,  j,k,irho) - w(i-1,j,k,irho)
                du2(irho) = w(i+1,j,k,irho) - w(i,  j,k,irho)
                du3(irho) = w(i+2,j,k,irho) - w(i+1,j,k,irho)

                du1(ivx) = w(i,  j,k,ivx) - w(i-1,j,k,ivx)
                du2(ivx) = w(i+1,j,k,ivx) - w(i,  j,k,ivx)
                du3(ivx) = w(i+2,j,k,ivx) - w(i+1,j,k,ivx)

                du1(ivy) = w(i,  j,k,ivy) - w(i-1,j,k,ivy)
                du2(ivy) = w(i+1,j,k,ivy) - w(i,  j,k,ivy)
                du3(ivy) = w(i+2,j,k,ivy) - w(i+1,j,k,ivy)

                du1(ivz) = w(i,  j,k,ivz) - w(i-1,j,k,ivz)
                du2(ivz) = w(i+1,j,k,ivz) - w(i,  j,k,ivz)
                du3(ivz) = w(i+2,j,k,ivz) - w(i+1,j,k,ivz)

                du1(irhoE) = p(i,  j,k) - p(i-1,j,k)
                du2(irhoE) = p(i+1,j,k) - p(i,  j,k)
                du3(irhoE) = p(i+2,j,k) - p(i+1,j,k)

                if( correctForK ) then
                   du1(itu1) = w(i,  j,k,itu1) - w(i-1,j,k,itu1)
                   du2(itu1) = w(i+1,j,k,itu1) - w(i,  j,k,itu1)
                   du3(itu1) = w(i+2,j,k,itu1) - w(i+1,j,k,itu1)
                endif

                ! Compute the differences from the first order scheme.

                call leftRightState(du1, du2, du3, rotMatrixI, &
                     left, right)

                ! Add the first order part to the currently stored
                ! differences, such that the correct state vector
                ! is stored.

                left(irho)  = left(irho)  + w(i,j,k,irho)
                left(ivx)   = left(ivx)   + w(i,j,k,ivx)
                left(ivy)   = left(ivy)   + w(i,j,k,ivy)
                left(ivz)   = left(ivz)   + w(i,j,k,ivz)
                left(irhoE) = left(irhoE) + p(i,j,k)

                right(irho)  = right(irho)  + w(i+1,j,k,irho)
                right(ivx)   = right(ivx)   + w(i+1,j,k,ivx)
                right(ivy)   = right(ivy)   + w(i+1,j,k,ivy)
                right(ivz)   = right(ivz)   + w(i+1,j,k,ivz)
                right(irhoE) = right(irhoE) + p(i+1,j,k)

                if( correctForK ) then
                   left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                   right(itu1) = right(itu1) + w(i+1,j,k,itu1)
                endif

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = si(i,j,k,1); sy = si(i,j,k,2); sz = si(i,j,k,3)
                por = porI(i,j,k)
                sFace = sFaceI(i,j,k)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i+1,j,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i+1,j,k,irho)  = fw(i+1,j,k,irho)  - flux(irho)
                fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   - flux(imx)
                fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   - flux(imy)
                fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   - flux(imz)
                fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) - flux(irhoE)

             enddo
          enddo
       enddo

       ! Fluxes in the j-direction.

       do k=2,kl
          do j=1,jl
             do i=2,il

                ! Store the three differences used in the interpolation
                ! in du1, du2, du3.

                du1(irho) = w(i,  j,k,irho) - w(i,j-1,k,irho)
                du2(irho) = w(i,j+1,k,irho) - w(i,  j,k,irho)
                du3(irho) = w(i,j+2,k,irho) - w(i,j+1,k,irho)

                du1(ivx) = w(i,  j,k,ivx) - w(i,j-1,k,ivx)
                du2(ivx) = w(i,j+1,k,ivx) - w(i,  j,k,ivx)
                du3(ivx) = w(i,j+2,k,ivx) - w(i,j+1,k,ivx)

                du1(ivy) = w(i,  j,k,ivy) - w(i,j-1,k,ivy)
                du2(ivy) = w(i,j+1,k,ivy) - w(i,  j,k,ivy)
                du3(ivy) = w(i,j+2,k,ivy) - w(i,j+1,k,ivy)

                du1(ivz) = w(i,  j,k,ivz) - w(i,j-1,k,ivz)
                du2(ivz) = w(i,j+1,k,ivz) - w(i,  j,k,ivz)
                du3(ivz) = w(i,j+2,k,ivz) - w(i,j+1,k,ivz)

                du1(irhoE) = p(i,  j,k) - p(i,j-1,k)
                du2(irhoE) = p(i,j+1,k) - p(i,  j,k)
                du3(irhoE) = p(i,j+2,k) - p(i,j+1,k)

                if( correctForK ) then
                   du1(itu1) = w(i,  j,k,itu1) - w(i,j-1,k,itu1)
                   du2(itu1) = w(i,j+1,k,itu1) - w(i,  j,k,itu1)
                   du3(itu1) = w(i,j+2,k,itu1) - w(i,j+1,k,itu1)
                endif

                ! Compute the differences from the first order scheme.

                call leftRightState(du1, du2, du3, rotMatrixJ, &
                     left, right)

                ! Add the first order part to the currently stored
                ! differences, such that the correct state vector
                ! is stored.

                left(irho)  = left(irho)  + w(i,j,k,irho)
                left(ivx)   = left(ivx)   + w(i,j,k,ivx)
                left(ivy)   = left(ivy)   + w(i,j,k,ivy)
                left(ivz)   = left(ivz)   + w(i,j,k,ivz)
                left(irhoE) = left(irhoE) + p(i,j,k)

                right(irho)  = right(irho)  + w(i,j+1,k,irho)
                right(ivx)   = right(ivx)   + w(i,j+1,k,ivx)
                right(ivy)   = right(ivy)   + w(i,j+1,k,ivy)
                right(ivz)   = right(ivz)   + w(i,j+1,k,ivz)
                right(irhoE) = right(irhoE) + p(i,j+1,k)

                if( correctForK ) then
                   left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                   right(itu1) = right(itu1) + w(i,j+1,k,itu1)
                endif

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sj(i,j,k,1); sy = sj(i,j,k,2); sz = sj(i,j,k,3)
                por = porJ(i,j,k)
                sFace = sFaceJ(i,j,k)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j+1,k))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j+1,k,irho)  = fw(i,j+1,k,irho)  - flux(irho)
                fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   - flux(imx)
                fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   - flux(imy)
                fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   - flux(imz)
                fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) - flux(irhoE)
             enddo
          enddo
       enddo

       ! Fluxes in the k-direction.

       do k=1,kl
          do j=2,jl
             do i=2,il

                ! Store the three differences used in the interpolation
                ! in du1, du2, du3.

                du1(irho) = w(i,  j,k,irho) - w(i,j,k-1,irho)
                du2(irho) = w(i,j,k+1,irho) - w(i,  j,k,irho)
                du3(irho) = w(i,j,k+2,irho) - w(i,j,k+1,irho)

                du1(ivx) = w(i,  j,k,ivx) - w(i,j,k-1,ivx)
                du2(ivx) = w(i,j,k+1,ivx) - w(i,  j,k,ivx)
                du3(ivx) = w(i,j,k+2,ivx) - w(i,j,k+1,ivx)

                du1(ivy) = w(i,  j,k,ivy) - w(i,j,k-1,ivy)
                du2(ivy) = w(i,j,k+1,ivy) - w(i,  j,k,ivy)
                du3(ivy) = w(i,j,k+2,ivy) - w(i,j,k+1,ivy)

                du1(ivz) = w(i,  j,k,ivz) - w(i,j,k-1,ivz)
                du2(ivz) = w(i,j,k+1,ivz) - w(i,  j,k,ivz)
                du3(ivz) = w(i,j,k+2,ivz) - w(i,j,k+1,ivz)

                du1(irhoE) = p(i,  j,k) - p(i,j,k-1)
                du2(irhoE) = p(i,j,k+1) - p(i,  j,k)
                du3(irhoE) = p(i,j,k+2) - p(i,j,k+1)

                if( correctForK ) then
                   du1(itu1) = w(i,  j,k,itu1) - w(i,j,k-1,itu1)
                   du2(itu1) = w(i,j,k+1,itu1) - w(i,  j,k,itu1)
                   du3(itu1) = w(i,j,k+2,itu1) - w(i,j,k+1,itu1)
                endif

                ! Compute the differences from the first order scheme.

                call leftRightState(du1, du2, du3, rotMatrixK, &
                     left, right)

                ! Add the first order part to the currently stored
                ! differences, such that the correct state vector
                ! is stored.

                left(irho)  = left(irho)  + w(i,j,k,irho)
                left(ivx)   = left(ivx)   + w(i,j,k,ivx)
                left(ivy)   = left(ivy)   + w(i,j,k,ivy)
                left(ivz)   = left(ivz)   + w(i,j,k,ivz)
                left(irhoE) = left(irhoE) + p(i,j,k)

                right(irho)  = right(irho)  + w(i,j,k+1,irho)
                right(ivx)   = right(ivx)   + w(i,j,k+1,ivx)
                right(ivy)   = right(ivy)   + w(i,j,k+1,ivy)
                right(ivz)   = right(ivz)   + w(i,j,k+1,ivz)
                right(irhoE) = right(irhoE) + p(i,j,k+1)

                if( correctForK ) then
                   left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                   right(itu1) = right(itu1) + w(i,j,k+1,itu1)
                endif

                ! Store the normal vector, the porosity and the
                ! mesh velocity if present.

                sx = sk(i,j,k,1); sy = sk(i,j,k,2); sz = sk(i,j,k,3)
                por = porK(i,j,k)
                sFace = sFaceK(i,j,k)

                ! Compute the value of gamma on the face. Take an
                ! arithmetic average of the two states.

                gammaFace = half*(gamma(i,j,k) + gamma(i,j,k+1))

                ! Compute the dissipative flux across the interface.

                call riemannFlux(left, right, flux)

                ! And scatter it to the left and right.

                fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
                fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
                fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
                fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

                fw(i,j,k+1,irho)  = fw(i,j,k+1,irho)  - flux(irho)
                fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   - flux(imx)
                fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   - flux(imy)
                fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   - flux(imz)
                fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) - flux(irhoE)

             enddo
          enddo
       enddo

    endif orderTest

    !      ==================================================================

  contains

    subroutine leftRightState(du1, du2, du3, rotMatrix, left, right)
      !
      !         leftRightState computes the differences in the left and
      !         right state compared to the first order interpolation. For a
      !         monotonic second order discretization the interpolations
      !         need to be nonlinear. The linear second order scheme can be
      !         stable (depending on the value of kappa), but it will have
      !         oscillations near discontinuities.
      !
      implicit none
      !
      !        Local parameter.
      !
      real(kind=realType), parameter :: epsLim = 1.e-10_realType
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(:), intent(inout) :: du1, du2, du3
      real(kind=realType), dimension(:), intent(out)   :: left, right

      real(kind=realType), dimension(:,:,:,:,:), pointer :: rotMatrix
      !
      !        Local variables.
      !
      integer(kind=intType) :: l

      real(kind=realType) :: rl1, rl2, rr1, rr2, tmp, dvx, dvy, dvz

      real(kind=realType), dimension(3,3) :: rot

      ! Check if the velocity components should be transformed to
      ! the cylindrical frame.

      if( rotationalPeriodic ) then

         ! Store the rotation matrix a bit easier. Note that the i,j,k
         ! come from the main subroutine.

         rot(1,1) = rotMatrix(i+ii-2, j+jj-2,k+kk-2, 1,1)
         rot(1,2) = rotMatrix(i+ii-2, j+jj-2,k+kk-2, 1,2)
         rot(1,3) = rotMatrix(i+ii-2, j+jj-2,k+kk-2, 1,3)

         rot(2,1) = rotMatrix(i+ii-2, j+jj-2, k+kk-2, 2,1)
         rot(2,2) = rotMatrix(i+ii-2, j+jj-2, k+kk-2, 2,2)
         rot(2,3) = rotMatrix(i+ii-2, j+jj-2, k+kk-2, 2,3)

         rot(3,1) = rotMatrix(i+ii-2, j+jj-2, k+kk-2, 3,1)
         rot(3,2) = rotMatrix(i+ii-2, j+jj-2, k+kk-2, 3,2)
         rot(3,3) = rotMatrix(i+ii-2, j+jj-2, k+kk-2, 3,3)

         ! Apply the transformation to the velocity components
         ! of du1, du2 and du3.

         dvx = du1(ivx); dvy = du1(ivy); dvz = du1(ivz)
         du1(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
         du1(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
         du1(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

         dvx = du2(ivx); dvy = du2(ivy); dvz = du2(ivz)
         du2(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
         du2(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
         du2(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

         dvx = du3(ivx); dvy = du3(ivy); dvz = du3(ivz)
         du3(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
         du3(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
         du3(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

      endif

      ! Determine the limiter used.

      select case (limUsed)

      case (noLimiter)

         ! Linear interpolation; no limiter.
         ! Loop over the number of variables to be interpolated.

         do l=1,nwInt
            left(l)  =  omk*du1(l) + opk*du2(l)
            right(l) = -omk*du3(l) - opk*du2(l)
         enddo

         !          ==============================================================

      case (vanAlbeda)

         ! Nonlinear interpolation using the van albeda limiter.
         ! Loop over the number of variables to be interpolated.

         do l=1,nwInt

            ! Compute the limiter argument rl1, rl2, rr1 and rr2.
            ! Note the cut off to 0.0.

            tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
            rl1 = max(zero, &
                 du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
            rl2 = max(zero,du1(l)*tmp)

            rr1 = max(zero,du3(l)*tmp)
            rr2 = max(zero, &
                 du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))

            ! Compute the corresponding limiter values.

            rl1 = rl1*(rl1 + one)/(rl1*rl1 + one)
            rl2 = rl2*(rl2 + one)/(rl2*rl2 + one)
            rr1 = rr1*(rr1 + one)/(rr1*rr1 + one)
            rr2 = rr2*(rr2 + one)/(rr2*rr2 + one)

            ! Compute the nonlinear corrections to the first order
            ! scheme.

            left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
            right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)

         enddo

         !          ==============================================================

      case (minmod)

         ! Nonlinear interpolation using the minmod limiter.
         ! Loop over the number of variables to be interpolated.

         do l=1,nwInt

            ! Compute the limiter argument rl1, rl2, rr1 and rr2.
            ! Note the cut off to 0.0.

            tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
            rl1 = max(zero, &
                 du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
            rl2 = max(zero,du1(l)*tmp)

            rr1 = max(zero,du3(l)*tmp)
            rr2 = max(zero, &
                 du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))

            ! Compute the corresponding limiter values.

            rl1 = min(one, factMinmod*rl1)
            rl2 = min(one, factMinmod*rl2)
            rr1 = min(one, factMinmod*rr1)
            rr2 = min(one, factMinmod*rr2)

            ! Compute the nonlinear corrections to the first order
            ! scheme.

            left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
            right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)

         enddo

      end select

      ! In case only a first order scheme must be used for the
      ! turbulent transport equations, set the correction for the
      ! turbulent kinetic energy to 0.

      if( firstOrderK ) then
         left(itu1)  = zero
         right(itu1) = zero
      endif

      ! For rotational periodic problems transform the velocity
      ! differences back to Cartesian again. Note that now the
      ! transpose of the rotation matrix must be used.

      if( rotationalPeriodic ) then

         ! Left state.

         dvx = left(ivx); dvy = left(ivy); dvz = left(ivz)
         left(ivx) = rot(1,1)*dvx + rot(2,1)*dvy + rot(3,1)*dvz
         left(ivy) = rot(1,2)*dvx + rot(2,2)*dvy + rot(3,2)*dvz
         left(ivz) = rot(1,3)*dvx + rot(2,3)*dvy + rot(3,3)*dvz

         ! Right state.

         dvx = right(ivx); dvy = right(ivy); dvz = right(ivz)
         right(ivx) = rot(1,1)*dvx + rot(2,1)*dvy + rot(3,1)*dvz
         right(ivy) = rot(1,2)*dvx + rot(2,2)*dvy + rot(3,2)*dvz
         right(ivz) = rot(1,3)*dvx + rot(2,3)*dvy + rot(3,3)*dvz

      endif

    end subroutine leftRightState

    !        ================================================================

    subroutine riemannFlux(left, right, flux)
      !
      !         riemannFlux computes the flux for the given face and left
      !         and right states.
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(*), intent(in)  :: left, right
      real(kind=realType), dimension(*), intent(out) :: flux
      !
      !        Local variables.
      !
      real(kind=realType) :: porFlux, rFace
      real(kind=realType) :: Etl, Etr, z1l, z1r, tmp
      real(kind=realType) :: dr, dru, drv, drw, drE, drk
      real(kind=realType) :: rAvg, uAvg, vAvg, wAvg, hAvg, kAvg
      real(kind=realType) :: alphaAvg, a2Avg, aAvg, unAvg
      real(kind=realType) :: ovaAvg, ova2Avg, area, Eta
      real(kind=realType) :: gm1, gm53
      real(kind=realType) :: lam1, lam2, lam3
      real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
      real(kind=realType), dimension(2) :: ktmp

      ! Set the porosity for the flux. The default value, 0.5*rFil, is
      ! a scaling factor where an rFil != 1 is taken into account.

      porFlux = half*rFil
      if(por == noFlux .or. por == boundFlux)    porFlux = zero

      ! Abbreviate some expressions in which gamma occurs.

      gm1  = gammaFace - one
      gm53 = gammaFace - five*third

      ! Determine which riemann solver must be solved.

      select case (riemannUsed)

      case (Roe)

         ! Determine the preconditioner used.

         select case (precond)

         case (noPrecond)

            ! No preconditioner used. Use the Roe scheme of the
            ! standard equations.

            ! Compute the square root of the left and right densities
            ! and the inverse of the sum.

            z1l = sqrt(left(irho))
            z1r = sqrt(right(irho))
            tmp = one/(z1l + z1r)

            ! Compute some variables depending whether or not a
            ! k-equation is present.

            if( correctForK ) then

               ! Store the left and right kinetic energy in ktmp,
               ! which is needed to compute the total energy.

               ktmp(1) = left(itu1)
               ktmp(2) = right(itu1)

               ! Store the difference of the turbulent kinetic energy
               ! per unit volume, i.e. the conserved variable.

               drk = right(irho)*right(itu1) - left(irho)*left(itu1)

               ! Compute the average turbulent energy per unit mass
               ! using Roe averages.

               kAvg = tmp*(z1l*left(itu1) + z1r*right(itu1))

            else

               ! Set the difference of the turbulent kinetic energy
               ! per unit volume and the averaged kinetic energy per
               ! unit mass to zero.

               drk  = 0.0
               kAvg = 0.0

            endif

            ! Compute the total energy of the left and right state.
            call etot(left(irho), left(ivx), left(ivy), left(ivz), &
                 left(irhoe), ktmp(1), Etl, correctForK)

            call etot(right(irho), right(ivx), right(ivy), right(ivz), &
                 right(irhoe), ktmp(2), Etr, correctForK)

            ! Compute the difference of the conservative mean
            ! flow variables.

            dr  = right(irho) - left(irho)
            dru = right(irho)*right(ivx) - left(irho)*left(ivx)
            drv = right(irho)*right(ivy) - left(irho)*left(ivy)
            drw = right(irho)*right(ivz) - left(irho)*left(ivz)
            drE = Etr - Etl

            ! Compute the Roe average variables, which can be
            ! computed directly from the average Roe vector.

            rAvg = fourth*(z1r + z1l)**2
            uAvg = tmp*(z1l*left(ivx) + z1r*right(ivx))
            vAvg = tmp*(z1l*left(ivy) + z1r*right(ivy))
            wAvg = tmp*(z1l*left(ivz) + z1r*right(ivz))
            hAvg = tmp*((Etl+left(irhoE)) /z1l &
                 +      (Etr+right(irhoE))/z1r)

            ! Compute the unit vector and store the area of the
            ! normal. Also compute the unit normal velocity of the face.

            area  = sqrt(sx**2 + sy**2 + sz**2)
            tmp   = one/max(1.e-25_realType,area)
            sx    = sx*tmp
            sy    = sy*tmp
            sz    = sz*tmp
            rFace = sFace*tmp

            ! Compute some dependent variables at the Roe
            ! average state.

            alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
            a2Avg    = abs(gm1*(hAvg - alphaAvg) - gm53*kAvg)
            aAvg     = sqrt(a2Avg)
            unAvg    = uAvg*sx + vAvg*sy + wAvg*sz

            ovaAvg  = one/aAvg
            ova2Avg = one/a2Avg

            ! Set for a boundary the normal velocity to rFace, the
            ! normal velocity of the boundary.

            if(por == boundFlux) unAvg = rFace

            ! Compute the coefficient eta for the entropy correction.
            ! At the moment a 1D entropy correction is used, which
            ! removes expansion shocks. Although it also reduces the
            ! carbuncle phenomenon, it does not remove it completely.
            ! In other to do that a multi-dimensional entropy fix is
            ! needed, see Sanders et. al, JCP, vol. 145, 1998,
            ! pp. 511 - 537. Although relatively easy to implement,
            ! an efficient implementation requires the storage of
            ! all the left and right states, which is rather
            ! expensive in terms of memory.

            eta = half*(abs((left(ivx) - right(ivx))*sx        &
                 +           (left(ivy) - right(ivy))*sy        &
                 +           (left(ivz) - right(ivz))*sz)       &
                 +       abs(sqrt(gammaFace*left(irhoE)/left(irho)) &
                 -           sqrt(gammaFace*right(irhoE)/right(irho))))

            ! Compute the absolute values of the three eigenvalues.

            lam1 = abs(unAvg - rFace + aAvg)
            lam2 = abs(unAvg - rFace - aAvg)
            lam3 = abs(unAvg - rFace)

            ! Apply the entropy correction to the eigenvalues.

            tmp = two*eta
            if(lam1 < tmp) lam1 = eta + fourth*lam1*lam1/eta
            if(lam2 < tmp) lam2 = eta + fourth*lam2*lam2/eta
            if(lam3 < tmp) lam3 = eta + fourth*lam3*lam3/eta

            ! Multiply the eigenvalues by the area to obtain
            ! the correct values for the dissipation term.

            lam1 = lam1*area
            lam2 = lam2*area
            lam3 = lam3*area

            ! Some abbreviations, which occur quite often in the
            ! dissipation terms.

            abv1 = half*(lam1 + lam2)
            abv2 = half*(lam1 - lam2)
            abv3 = abv1 - lam3

            abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                 -      wAvg*drw + drE) - gm53*drk
            abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

            abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
            abv7 = abv2*abv4*ovaAvg  + abv3*abv5

            ! Compute the dissipation term, -|a| (wr - wl), which is
            ! multiplied by porFlux. Note that porFlux is either
            ! 0.0 or 0.5*rFil.

            flux(irho)  = -porFlux*(lam3*dr  + abv6)
            flux(imx)   = -porFlux*(lam3*dru + uAvg*abv6 &
                 + sx*abv7)
            flux(imy)   = -porFlux*(lam3*drv + vAvg*abv6 &
                 + sy*abv7)
            flux(imz)   = -porFlux*(lam3*drw + wAvg*abv6 &
                 + sz*abv7)
            flux(irhoE) = -porFlux*(lam3*drE + hAvg*abv6 &
                 + unAvg*abv7)

            !          tmp = max(lam1,lam2,lam3)

            !          flux(irho)  = -porFlux*(tmp*dr)
            !          flux(imx)   = -porFlux*(tmp*dru)
            !          flux(imy)   = -porFlux*(tmp*drv)
            !          flux(imz)   = -porFlux*(tmp*drw)
            !          flux(irhoE) = -porFlux*(tmp*drE)

         case (Turkel)
            call terminate(&
                 "riemannFlux",&
                 "Turkel preconditioner not implemented yet")

         case (ChoiMerkle)
            call terminate("riemannFlux",&
                 "choi merkle preconditioner not implemented yet")

         end select

      case (vanLeer)
         call terminate("riemannFlux", "van leer fvs not implemented yet")

      case (ausmdv)
         call terminate("riemannFlux","ausmdv fvs not implemented yet")

      end select

    end subroutine riemannFlux

  end subroutine inviscidUpwindFlux


  subroutine inviscidDissFluxScalarApprox
    ! ---------------------------------------------
    !             Inviscid Diss Flux Scalar
    ! ---------------------------------------------

    use constants
    use flowVarRefState, only : pInfCorr
    use inputDiscretization, only: vis2, vis4, sigma
    use inputPhysics, only : equations
    use iteration, only : rFil
    use flowVarRefState, only : gammaInf, pInfCorr, rhoInf
    implicit none

    ! Variables for inviscid diss flux scalar
    real(kind=realType), parameter :: dssMax = 0.25_realType
    real(kind=realType) :: sslim, rhoi
    real(kind=realType) :: sfil, fis2, fis4
    real(kind=realType) :: ppor, rrad, dis2, dis4, fs
    real(kind=realType) :: ddw
    integer(kind=intType) :: i, j, k
    select case (equations)
    case (EulerEquations)

       ! Inviscid case. Pressure switch is based on the pressure.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of pressure and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr

       !===============================================================

    case (NSEquations, RANSEquations)

       ! Viscous case. Pressure switch is based on the entropy.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of entropy and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

    end select

    ! Compute the pressure sensor for each cell, in each direction:
    do k=1,ke
       do j=1,je
          do i=singleHaloStart,ie
             dss(i,j,k,1) = abs((ss(i+1,j,k) - two*ss(i,j,k) + ss(i-1,j,k)) &
                  /     (ss(i+1,j,k) + two*ss(i,j,k) + ss(i-1,j,k) + sslim))

             dss(i,j,k,2) = abs((ss(i,j+1,k) - two*ss(i,j,k) + ss(i,j-1,k)) &
                  /     (ss(i,j+1,k) + two*ss(i,j,k) + ss(i,j-1,k) + sslim))

             dss(i,j,k,3) = abs((ss(i,j,k+1) - two*ss(i,j,k) + ss(i,j,k-1)) &
                  /     (ss(i,j,k+1) + two*ss(i,j,k) + ss(i,j,k-1) + sslim))
          end do
       end do
    end do

    ! Set a couple of constants for the scheme.
    fis2 = vis2
    fis4 = vis4
    !
    !       Dissipative fluxes in the i-direction.
    !
    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,1), dss(i+1,j,k,1))) + sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i+1,j,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i+1,j,k,ivx)*w(i+1,j,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i+1,j,k,ivy)*w(i+1,j,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.
             ddw = w(i+1,j,k,ivz)*w(i+1,j,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.
             ddw = (w(i+1,j,k,irhoE) + P(i+1,j,K))- (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw

             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do
    !
    !       Dissipative fluxes in the j-direction.
    !
    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,2),dss(i,j+1,k,2))) +sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j+1,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i,j+1,k,ivx)*w(i,j+1,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j+1,k,ivy)*w(i,j+1,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j+1,k,ivz)*w(i,j+1,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = (w(i,j+1,k,irhoE) + P(i,j+1,k)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw

             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the k-direction.
    !
    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,3), dss(i,j,k+1,3))) + sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j,k+1,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i,j,k+1,ivx)*w(i,j,k+1,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j,k+1,ivy)*w(i,j,k+1,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j,k+1,ivz)*w(i,j,k+1,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.
             ddw = (w(i,j,k+1,irhoE) + P(i,j,k+1)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw

             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
  end subroutine inviscidDissFluxScalarApprox

  subroutine inviscidDissFluxMatrixApprox
    !
    !       inviscidDissFluxMatrix computes the matrix artificial
    !       dissipation term. Instead of the spectral radius, as used in
    !       the scalar dissipation scheme, the absolute value of the flux
    !       jacobian is used. This leads to a less diffusive and
    !       consequently more accurate scheme. It is assumed that the
    !       pointers in blockPointers already point to the correct block.
    !
    use constants
    use flowVarRefState, only : pInfCorr
    use inputDiscretization, only: vis2, vis4, sigma
    use inputPhysics, only : equations
    use iteration, only : rFil
    use utils, only : getCorrectForK, myDim
    implicit none
    !
    !      Local parameters.
    !
    real(kind=realType), parameter :: dpMax        = 0.25_realType
    real(kind=realType), parameter :: epsAcoustic  = 0.25_realType
    real(kind=realType), parameter :: epsShear     = 0.025_realType
    real(kind=realType), parameter :: omega        = 0.5_realType
    real(kind=realType), parameter :: oneMinOmega = one - omega
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ind, ii

    real(kind=realType) :: plim, sface
    real(kind=realType) :: sfil, fis2, fis4
    real(kind=realType) :: gammaAvg, gm1, ovgm1, gm53
    real(kind=realType) :: ppor, rrad, dis2, dis4
    real(kind=realType) :: dp1, dp2, tmp, fs
    real(kind=realType) :: ddw, ddw6
    real(kind=realType) :: dr, dru, drv, drw, dre, drk, sx, sy, sz
    real(kind=realType) :: uAvg, vAvg, wAvg, a2Avg, aAvg, hAvg
    real(kind=realType) :: alphaAvg, unAvg, ovaAvg, ova2Avg
    real(kind=realType) :: kAvg, lam1, lam2, lam3, area
    real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7
    logical :: correctForK

    ! Set the value of plim. To be fully consistent this must have
    ! the dimension of a pressure. Therefore a fraction of pInfCorr
    ! is used.

    plim = 0.001_realType*pInfCorr

    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.

    correctForK = getCorrectForK()

    ! Initialize sface to zero. This value will be used if the
    ! block is not moving.

    sface = zero

    ! Set a couple of constants for the scheme.

    fis2 = rFil*vis2
    fis4 = rFil*vis4
    sfil = one - rFil

    ! Initialize the dissipative residual to a certain times,
    ! possibly zero, the previously stored value.

    fw = sfil*fw

    ! Compute the pressure sensor for each cell, in each direction:
    do k=1,ke
       do j=1,je
          do i=singleHaloStart,ie
             dss(i,j,k,1) =abs((ss(i+1,j,k) - two*ss(i,j,k) + ss(i-1,j,k))        &
                  /     (omega*(ss(i+1,j,k) + two*ss(i,j,k) + ss(i-1,j,k)) &
                  +      oneMinOmega*(abs(ss(i+1,j,k) - ss(i,j,k))      &
                  +                   abs(ss(i,j,k) - ss(i-1,j,k))) + plim))


             dss(i,j,k,2) =abs((ss(i,j+1,k) - two*ss(i,j,k) + ss(i,j-1,k))        &
                  /     (omega*(ss(i,j+1,k) + two*ss(i,j,k) + ss(i,j-1,k)) &
                  +      oneMinOmega*(abs(ss(i,j+1,k) - ss(i,j,k))      &
                  +                   abs(ss(i,j,k) - ss(i,j-1,k))) + plim))

             dss(i,j,k,3) =  abs((ss(i,j,k+1) - two*ss(i,j,k) + ss(i,j,k-1))        &
                  /     (omega*(ss(i,j,k+1) + two*ss(i,j,k) + ss(i,j,k-1)) &
                  +      oneMinOmega*(abs(ss(i,j,k+1) - ss(i,j,k))      &
                  +                   abs(ss(i,j,k) - ss(i,j,k-1))) + plim))
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the i-direction.
    !
    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = one

             dis2 = fis2*ppor*min(dpMax,max(dss(i,j,k,1),dss(i+1,j,k,1)))&
                    +sigma*fis4*ppor

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw = w(i+1,j,k,irho) - w(i,j,k,irho)
             dr  = dis2*ddw

             ddw = w(i+1,j,k,irho)*w(i+1,j,k,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw

             ddw = w(i+1,j,k,irho)*w(i+1,j,k,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw

             ddw = w(i+1,j,k,irho)*w(i+1,j,k,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw

             ddw = w(i+1,j,k,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.
             drk   = zero
             kAvg = zero

             if( correctForK ) then
                ddw6 = w(i+1,j,k,irho)*w(i+1,j,k,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw6                              &
                     - dis4*(w(i+2,j,k,irho)*w(i+2,j,k,itu1) &
                     -       w(i-1,j,k,irho)*w(i-1,j,k,itu1) - three*ddw6)

                kAvg = half*(w(i,j,k,itu1) + w(i+1,j,k,itu1))
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i+1,j,k) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i+1,j,k,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i+1,j,k,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i+1,j,k,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i+1,j,k)*p(i+1,j,k)/w(i+1,j,k,irho) &
                  +       gamma(i,  j,k)*p(i,  j,k)/w(i,  j,k,irho))

             area = sqrt(si(i,j,k,1)**2 + si(i,j,k,2)**2 + si(i,j,k,3)**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = si(i,j,k,1)*tmp
             sy   = si(i,j,k,2)*tmp
             sz   = si(i,j,k,3)*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             sface = sFaceI(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = max(lam1,epsAcoustic*rrad)*area
             lam2 = max(lam2,epsAcoustic*rrad)*area
             lam3 = max(lam3,epsShear*rrad)*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do
    !
    !       Dissipative fluxes in the j-direction.
    !
    do k=2,kl
       do j=1,jl
          do i=2,il


             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = one

             dis2 = fis2*ppor*min(dpMax,max(dss(i,j,k,2),dss(i,j+1,k,2)))&
                    +sigma*fis4*ppor

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw = w(i,j+1,k,irho) - w(i,j,k,irho)
             dr  = dis2*ddw

             ddw = w(i,j+1,k,irho)*w(i,j+1,k,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw

             ddw = w(i,j+1,k,irho)*w(i,j+1,k,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw

             ddw = w(i,j+1,k,irho)*w(i,j+1,k,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw

             ddw = w(i,j+1,k,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.
             drk   = zero
             kAvg = zero

             if( correctForK ) then
                ddw6 = w(i,j+1,k,irho)*w(i,j+1,k,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw6                              &
                     - dis4*(w(i,j+2,k,irho)*w(i,j+2,k,itu1) &
                     -       w(i,j-1,k,irho)*w(i,j-1,k,itu1) - three*ddw6)

                kAvg = half*(w(i,j,k,itu1) + w(i,j+1,k,itu1))
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i,j+1,k) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i,j+1,k,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i,j+1,k,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i,j+1,k,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i,j+1,k)*p(i,j+1,k)/w(i,j+1,k,irho) &
                  +       gamma(i,j,  k)*p(i,j,  k)/w(i,j,  k,irho))

             area = sqrt(sj(i,j,k,1)**2 + sj(i,j,k,2)**2 + sj(i,j,k,3)**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = sj(i,j,k,1)*tmp
             sy   = sj(i,j,k,2)*tmp
             sz   = sj(i,j,k,3)*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             sface = sFaceJ(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = max(lam1,epsAcoustic*rrad)*area
             lam2 = max(lam2,epsAcoustic*rrad)*area
             lam3 = max(lam3,epsShear*rrad)*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do
    !
    !       Dissipative fluxes in the k-direction.
    !
    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = one

             dis2 = fis2*ppor*min(dpMax,max(dss(i,j,k,3),dss(i,j,k+1,3)))&
                    +sigma*fis4*ppor

             ! Construct the vector of the first and third differences
             ! multiplied by the appropriate constants.

             ddw = w(i,j,k+1,irho) - w(i,j,k,irho)
             dr  = dis2*ddw

             ddw = w(i,j,k+1,irho)*w(i,j,k+1,ivx) &
                  - w(i,j,k,irho)*w(i,j,k,ivx)
             dru = dis2*ddw

             ddw = w(i,j,k+1,irho)*w(i,j,k+1,ivy) &
                  - w(i,j,k,irho)*w(i,j,k,ivy)
             drv = dis2*ddw

             ddw = w(i,j,k+1,irho)*w(i,j,k+1,ivz) &
                  - w(i,j,k,irho)*w(i,j,k,ivz)
             drw = dis2*ddw

             ddw = w(i,j,k+1,irhoE) - w(i,j,k,irhoE)
             dre = dis2*ddw

             ! In case a k-equation is present, compute the difference
             ! of rhok and store the average value of k. If not present,
             ! set both these values to zero, such that later on no
             ! decision needs to be made anymore.
             drk   = zero
             kAvg = zero

             if( correctForK ) then
                ddw6 = w(i,j,k+1,irho)*w(i,j,k+1,itu1) &
                     - w(i,j,k,irho)*w(i,j,k,itu1)
                drk = dis2*ddw6                              &
                     - dis4*(w(i,j,k+2,irho)*w(i,j,k+2,itu1) &
                     -       w(i,j,k-1,irho)*w(i,j,k-1,itu1) - three*ddw6)

                kAvg = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
             endif

             ! Compute the average value of gamma and compute some
             ! expressions in which it occurs.

             gammaAvg = half*(gamma(i,j,k+1) + gamma(i,j,k))
             gm1      = gammaAvg - one
             ovgm1    = one/gm1
             gm53     = gammaAvg - five*third

             ! Compute the average state at the interface.

             uAvg  = half*(w(i,j,k+1,ivx) + w(i,j,k,ivx))
             vAvg  = half*(w(i,j,k+1,ivy) + w(i,j,k,ivy))
             wAvg  = half*(w(i,j,k+1,ivz) + w(i,j,k,ivz))
             a2Avg = half*(gamma(i,j,k+1)*p(i,j,k+1)/w(i,j,k+1,irho) &
                  +       gamma(i,j,k)  *p(i,j,k)  /w(i,j,k,  irho))

             area = sqrt(sk(i,j,k,1)**2 + sk(i,j,k,2)**2 + sk(i,j,k,3)**2)
             tmp  = one/max(1.e-25_realType,area)
             sx   = sk(i,j,k,1)*tmp
             sy   = sk(i,j,k,2)*tmp
             sz   = sk(i,j,k,3)*tmp

             alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
             hAvg     = alphaAvg + ovgm1*(a2Avg - gm53*kAvg)
             aAvg     = sqrt(a2Avg)
             unAvg    = uAvg*sx + vAvg*sy + wAvg*sz
             ovaAvg   = one/aAvg
             ova2Avg  = one/a2Avg

             ! The mesh velocity if the face is moving. It must be
             ! divided by the area to obtain a true velocity.

             sface = sFaceK(i,j,k)*tmp

             ! Compute the absolute values of the three eigenvalues
             ! and make sure they don't become zero by cutting them
             ! off to a certain minimum.

             lam1 = abs(unAvg - sface + aAvg)
             lam2 = abs(unAvg - sface - aAvg)
             lam3 = abs(unAvg - sface)

             rrad = lam3 + aAvg

             ! Multiply the eigenvalues by the area to obtain
             ! the correct values for the dissipation term.

             lam1 = max(lam1,epsAcoustic*rrad)*area
             lam2 = max(lam2,epsAcoustic*rrad)*area
             lam3 = max(lam3,epsShear*rrad)*area

             ! Some abbreviations, which occur quite often in the
             ! dissipation terms.

             abv1 = half*(lam1 + lam2)
             abv2 = half*(lam1 - lam2)
             abv3 = abv1 - lam3

             abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                  -      wAvg*drw + dre) - gm53*drk
             abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

             abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
             abv7 = abv2*abv4*ovaAvg  + abv3*abv5

             ! Compute and scatter the dissipative flux.
             ! Density.

             fs               = lam3*dr  + abv6
             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             fs              = lam3*dru + uAvg*abv6 + sx*abv7
             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             fs              = lam3*drv + vAvg*abv6 + sy*abv7
             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             fs              = lam3*drw + wAvg*abv6 + sz*abv7
             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             fs                = lam3*dre + hAvg*abv6 + unAvg*abv7
             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do

  end subroutine inviscidDissFluxMatrixApprox


  subroutine computeSpeedOfSoundSquared
    ! ---------------------------------------------
    !           Compute the speed of sound squared
    ! ---------------------------------------------
    use constants
    use utils, only : getCorrectForK
    implicit none

    ! Variables for speed of sound
    logical :: correctForK
    real(kind=realType) :: pp
    real(kind=realType), parameter :: twoThird = two*third
    integer(kind=intType) :: i, j, k

    ! Determine if we need to correct for K
    correctForK = getCorrectForK()

    if (correctForK) then
       do k=1,ke
          do j=1,je
             do i=singleHaloStart,ie
                pp = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
                aa(i,j,k) = gamma(i,j,k)*pp/w(i,j,k,irho)
             enddo
          enddo
       enddo
    else
       do k=1,ke
          do j=1,je
             do i=singleHaloStart,ie
                aa(i,j,k) = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
             enddo
          enddo
       enddo
    end if
  end subroutine computeSpeedOfSoundSquared

  subroutine allNodalGradients

    ! ---------------------------------------------
    !           Compute nodal gradients
    ! ---------------------------------------------
    use constants
    implicit none

    ! Variables for nodal gradients
    real(kind=realType) :: a2, oVol, uBar, vBar, wBar, sx, sy, sz
    integer(kind=intType) :: i, j, k

    ! Zero just the required part of the nodal gradients since the
    ! first value may be useful.
    ux(nodeStart:, :, :) = zero
    uy(nodeStart:, :, :) = zero
    uz(nodeStart:, :, :) = zero

    vx(nodeStart:, :, :) = zero
    vy(nodeStart:, :, :) = zero
    vz(nodeStart:, :, :) = zero

    wx(nodeStart:, :, :) = zero
    wy(nodeStart:, :, :) = zero
    wz(nodeStart:, :, :) = zero

    qx(nodeStart:, :, :) = zero
    qy(nodeStart:, :, :) = zero
    qz(nodeStart:, :, :) = zero

    ! First part. Contribution in the k-direction.
    ! The contribution is scattered to both the left and right node
    ! in k-direction.

    do k=1, ke
       do j=1, jl
          do i=nodeStart, il

             ! Compute 8 times the average normal for this part of
             ! the control volume. The factor 8 is taken care of later
             ! on when the division by the volume takes place.

             sx = sk(i,j,k-1,  1) + sk(i+1,j,k-1,  1) &
                  + sk(i,j+1,k-1,1) + sk(i+1,j+1,k-1,1) &
                  + sk(i,j,  k,  1) + sk(i+1,j,  k,  1) &
                  + sk(i,j+1,k  ,1) + sk(i+1,j+1,k  ,1)
             sy = sk(i,j,k-1,  2) + sk(i+1,j,k-1,  2) &
                  + sk(i,j+1,k-1,2) + sk(i+1,j+1,k-1,2) &
                  + sk(i,j,  k,  2) + sk(i+1,j,  k,  2) &
                  + sk(i,j+1,k  ,2) + sk(i+1,j+1,k  ,2)
             sz = sk(i,j,k-1,  3) + sk(i+1,j,k-1,  3) &
                  + sk(i,j+1,k-1,3) + sk(i+1,j+1,k-1,3) &
                  + sk(i,j,  k,  3) + sk(i+1,j,  k,  3) &
                  + sk(i,j+1,k  ,3) + sk(i+1,j+1,k  ,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,  k,ivx) + w(i+1,j,  k,ivx) &
                  +         w(i,j+1,k,ivx) + w(i+1,j+1,k,ivx))
             vbar = fourth*(w(i,j,  k,ivy) + w(i+1,j,  k,ivy) &
                  +         w(i,j+1,k,ivy) + w(i+1,j+1,k,ivy))
             wbar = fourth*(w(i,j,  k,ivz) + w(i+1,j,  k,ivz) &
                  +         w(i,j+1,k,ivz) + w(i+1,j+1,k,ivz))

             a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j+1,k) + aa(i+1,j+1,k))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(k > 1) then
                ux(i,j,k-1) = ux(i,j,k-1) + ubar*sx
                uy(i,j,k-1) = uy(i,j,k-1) + ubar*sy
                uz(i,j,k-1) = uz(i,j,k-1) + ubar*sz

                vx(i,j,k-1) = vx(i,j,k-1) + vbar*sx
                vy(i,j,k-1) = vy(i,j,k-1) + vbar*sy
                vz(i,j,k-1) = vz(i,j,k-1) + vbar*sz

                wx(i,j,k-1) = wx(i,j,k-1) + wbar*sx
                wy(i,j,k-1) = wy(i,j,k-1) + wbar*sy
                wz(i,j,k-1) = wz(i,j,k-1) + wbar*sz

                qx(i,j,k-1) = qx(i,j,k-1) - a2*sx
                qy(i,j,k-1) = qy(i,j,k-1) - a2*sy
                qz(i,j,k-1) = qz(i,j,k-1) - a2*sz
             endif

             if(k < ke) then
                ux(i,j,k) = ux(i,j,k) - ubar*sx
                uy(i,j,k) = uy(i,j,k) - ubar*sy
                uz(i,j,k) = uz(i,j,k) - ubar*sz

                vx(i,j,k) = vx(i,j,k) - vbar*sx
                vy(i,j,k) = vy(i,j,k) - vbar*sy
                vz(i,j,k) = vz(i,j,k) - vbar*sz

                wx(i,j,k) = wx(i,j,k) - wbar*sx
                wy(i,j,k) = wy(i,j,k) - wbar*sy
                wz(i,j,k) = wz(i,j,k) - wbar*sz

                qx(i,j,k) = qx(i,j,k) + a2*sx
                qy(i,j,k) = qy(i,j,k) + a2*sy
                qz(i,j,k) = qz(i,j,k) + a2*sz
             endif
          end do
       enddo
    enddo


    ! Second part. Contribution in the j-direction.
    ! The contribution is scattered to both the left and right node
    ! in j-direction.

    do k=1, kl
       do j=1, je
          do i=nodeStart, il

             ! Compute 8 times the average normal for this part of
             ! the control volume. The factor 8 is taken care of later
             ! on when the division by the volume takes place.

             sx = sj(i,j-1,k,  1) + sj(i+1,j-1,k,  1) &
                  + sj(i,j-1,k+1,1) + sj(i+1,j-1,k+1,1) &
                  + sj(i,j,  k,  1) + sj(i+1,j,  k,  1) &
                  + sj(i,j,  k+1,1) + sj(i+1,j,  k+1,1)
             sy = sj(i,j-1,k,  2) + sj(i+1,j-1,k,  2) &
                  + sj(i,j-1,k+1,2) + sj(i+1,j-1,k+1,2) &
                  + sj(i,j,  k,  2) + sj(i+1,j,  k,  2) &
                  + sj(i,j,  k+1,2) + sj(i+1,j,  k+1,2)
             sz = sj(i,j-1,k,  3) + sj(i+1,j-1,k,  3) &
                  + sj(i,j-1,k+1,3) + sj(i+1,j-1,k+1,3) &
                  + sj(i,j,  k,  3) + sj(i+1,j,  k,  3) &
                  + sj(i,j,  k+1,3) + sj(i+1,j,  k+1,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,k,  ivx) + w(i+1,j,k,  ivx) &
                  +         w(i,j,k+1,ivx) + w(i+1,j,k+1,ivx))
             vbar = fourth*(w(i,j,k,  ivy) + w(i+1,j,k,  ivy) &
                  +         w(i,j,k+1,ivy) + w(i+1,j,k+1,ivy))
             wbar = fourth*(w(i,j,k,  ivz) + w(i+1,j,k,  ivz) &
                  +         w(i,j,k+1,ivz) + w(i+1,j,k+1,ivz))

             a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j,k+1) + aa(i+1,j,k+1))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(j > 1) then
                ux(i,j-1,k) = ux(i,j-1,k) + ubar*sx
                uy(i,j-1,k) = uy(i,j-1,k) + ubar*sy
                uz(i,j-1,k) = uz(i,j-1,k) + ubar*sz

                vx(i,j-1,k) = vx(i,j-1,k) + vbar*sx
                vy(i,j-1,k) = vy(i,j-1,k) + vbar*sy
                vz(i,j-1,k) = vz(i,j-1,k) + vbar*sz

                wx(i,j-1,k) = wx(i,j-1,k) + wbar*sx
                wy(i,j-1,k) = wy(i,j-1,k) + wbar*sy
                wz(i,j-1,k) = wz(i,j-1,k) + wbar*sz

                qx(i,j-1,k) = qx(i,j-1,k) - a2*sx
                qy(i,j-1,k) = qy(i,j-1,k) - a2*sy
                qz(i,j-1,k) = qz(i,j-1,k) - a2*sz
             endif

             if(j < je) then
                ux(i,j,k) = ux(i,j,k) - ubar*sx
                uy(i,j,k) = uy(i,j,k) - ubar*sy
                uz(i,j,k) = uz(i,j,k) - ubar*sz

                vx(i,j,k) = vx(i,j,k) - vbar*sx
                vy(i,j,k) = vy(i,j,k) - vbar*sy
                vz(i,j,k) = vz(i,j,k) - vbar*sz

                wx(i,j,k) = wx(i,j,k) - wbar*sx
                wy(i,j,k) = wy(i,j,k) - wbar*sy
                wz(i,j,k) = wz(i,j,k) - wbar*sz

                qx(i,j,k) = qx(i,j,k) + a2*sx
                qy(i,j,k) = qy(i,j,k) + a2*sy
                qz(i,j,k) = qz(i,j,k) + a2*sz
             endif
          end do
       enddo
    enddo
    !
    ! Third part. Contribution in the i-direction.
    ! The contribution is scattered to both the left and right node
    ! in i-direction.
    !
    do k=1,kl
       do j=1,jl
          do i=nodeStart,ie

             ! Compute 8 times the average normal for this part of
             ! the control volume. The factor 8 is taken care of later
             ! on when the division by the volume takes place.

             sx = si(i-1,j,k,  1) + si(i-1,j+1,k,  1) &
                  + si(i-1,j,k+1,1) + si(i-1,j+1,k+1,1) &
                  + si(i,  j,k,  1) + si(i,  j+1,k,  1) &
                  + si(i,  j,k+1,1) + si(i,  j+1,k+1,1)
             sy = si(i-1,j,k,  2) + si(i-1,j+1,k,  2) &
                  + si(i-1,j,k+1,2) + si(i-1,j+1,k+1,2) &
                  + si(i,  j,k,  2) + si(i,  j+1,k,  2) &
                  + si(i,  j,k+1,2) + si(i,  j+1,k+1,2)
             sz = si(i-1,j,k,  3) + si(i-1,j+1,k,  3) &
                  + si(i-1,j,k+1,3) + si(i-1,j+1,k+1,3) &
                  + si(i,  j,k,  3) + si(i,  j+1,k,  3) &
                  + si(i,  j,k+1,3) + si(i,  j+1,k+1,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,k,  ivx) + w(i,j+1,k,  ivx) &
                  +         w(i,j,k+1,ivx) + w(i,j+1,k+1,ivx))
             vbar = fourth*(w(i,j,k,  ivy) + w(i,j+1,k,  ivy) &
                  +         w(i,j,k+1,ivy) + w(i,j+1,k+1,ivy))
             wbar = fourth*(w(i,j,k,  ivz) + w(i,j+1,k,  ivz) &
                  +         w(i,j,k+1,ivz) + w(i,j+1,k+1,ivz))

             a2 = fourth*(aa(i,j,k) + aa(i,j+1,k) + aa(i,j,k+1) + aa(i,j+1,k+1))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(i > nodeStart) then
                ux(i-1,j,k) = ux(i-1,j,k) + ubar*sx
                uy(i-1,j,k) = uy(i-1,j,k) + ubar*sy
                uz(i-1,j,k) = uz(i-1,j,k) + ubar*sz

                vx(i-1,j,k) = vx(i-1,j,k) + vbar*sx
                vy(i-1,j,k) = vy(i-1,j,k) + vbar*sy
                vz(i-1,j,k) = vz(i-1,j,k) + vbar*sz

                wx(i-1,j,k) = wx(i-1,j,k) + wbar*sx
                wy(i-1,j,k) = wy(i-1,j,k) + wbar*sy
                wz(i-1,j,k) = wz(i-1,j,k) + wbar*sz

                qx(i-1,j,k) = qx(i-1,j,k) - a2*sx
                qy(i-1,j,k) = qy(i-1,j,k) - a2*sy
                qz(i-1,j,k) = qz(i-1,j,k) - a2*sz
             endif

             if(i < ie) then
                ux(i,j,k) = ux(i,j,k) - ubar*sx
                uy(i,j,k) = uy(i,j,k) - ubar*sy
                uz(i,j,k) = uz(i,j,k) - ubar*sz

                vx(i,j,k) = vx(i,j,k) - vbar*sx
                vy(i,j,k) = vy(i,j,k) - vbar*sy
                vz(i,j,k) = vz(i,j,k) - vbar*sz

                wx(i,j,k) = wx(i,j,k) - wbar*sx
                wy(i,j,k) = wy(i,j,k) - wbar*sy
                wz(i,j,k) = wz(i,j,k) - wbar*sz

                qx(i,j,k) = qx(i,j,k) + a2*sx
                qy(i,j,k) = qy(i,j,k) + a2*sy
                qz(i,j,k) = qz(i,j,k) + a2*sz
             endif
          enddo
       enddo
    enddo

    ! Divide by 8 times the volume to obtain the correct gradients.

    do k=1,kl
       do j=1,jl
          do i=nodeStart,il

             ! Compute the inverse of 8 times the volume for this node.

             oVol = one/(vol(i,  j,  k) + vol(i,  j,  k+1) &
                  +      vol(i+1,j,  k) + vol(i+1,j,  k+1) &
                  +      vol(i,  j+1,k) + vol(i,  j+1,k+1) &
                  +      vol(i+1,j+1,k) + vol(i+1,j+1,k+1))

             ! Compute the correct velocity gradients and "unit" heat
             ! fluxes. The velocity gradients are stored in ux, etc.

             ux(i,j,k) = ux(i,j,k)*oVol
             uy(i,j,k) = uy(i,j,k)*oVol
             uz(i,j,k) = uz(i,j,k)*oVol

             vx(i,j,k) = vx(i,j,k)*oVol
             vy(i,j,k) = vy(i,j,k)*oVol
             vz(i,j,k) = vz(i,j,k)*oVol

             wx(i,j,k) = wx(i,j,k)*oVol
             wy(i,j,k) = wy(i,j,k)*oVol
             wz(i,j,k) = wz(i,j,k)*oVol

             qx(i,j,k) = qx(i,j,k)*oVol
             qy(i,j,k) = qy(i,j,k)*oVol
             qz(i,j,k) = qz(i,j,k)*oVol
          end do
       enddo
    enddo
  end subroutine allNodalGradients

  subroutine viscousFlux(storeWallTensor)
    ! ---------------------------------------------
    !                   Viscous Flux
    ! ---------------------------------------------

    use constants
    use inputPhysics, only : useQCR, prandtl, prandtlturb
    use flowvarRefState, only : eddyModel
    use iteration, only : rFil
    use blockPointers, only : bil => il, bjl=>jl, bkl=>kl, &
         viscIminPointer, viscImaxPointer, viscSubFace, &
         viscJminPointer, viscJmaxPointer, &
         viscKminPointer, viscKmaxPointer
    implicit none

    ! Input
    logical, intent(in), optional :: storeWallTensor

    ! Variables for viscous flux
    real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
    real(kind=realType) :: gm1, factLamHeat, factTurbHeat
    real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
    real(kind=realType) :: q_x, q_y, q_z
    real(kind=realType) :: corr, ssx, ssy, ssz, fracDiv, snrm
    real(kind=realType) :: tauxx, tauyy, tauzz
    real(kind=realType) :: tauxy, tauxz, tauyz
    real(kind=realType) :: tauxxS, tauyyS, tauzzS
    real(kind=realType) :: tauxyS, tauxzS, tauyzS
    real(kind=realType) :: ubar, vbar, wbar
    real(kind=realType) :: exx, eyy, ezz
    real(kind=realType) :: exy, exz, eyz
    real(kind=realType) :: Wxx, Wyy, Wzz
    real(kind=realType) :: Wxy, Wxz, Wyz, Wyx, Wzx, Wzy
    real(kind=realType) :: den, Ccr1
    real(kind=realType) :: fmx, fmy, fmz, frhoE, fact
    integer(kind=intType) :: i, j, k, io, jo, ko
    real(kind=realType), parameter :: xminn = 1.e-10_realType
    real(kind=realType), parameter :: twoThird = two*third
    real(kind=realType), dimension(9, 2:max(il,jl), 2: max(jl,kl), 2) :: tmpStore

    logical :: storeWall

    storeWall = .False.
    if (present(storeWallTensor)) then
       storeWall = .True.
    end if

    ! Set QCR parameters
    Ccr1 = 0.3_realType
    rFilv = rFil

    ! The diagonals of the vorticity tensor components are always zero
    Wxx = zero
    Wyy = zero
    Wzz = zero
    !
    !         viscous fluxes in the k-direction.
    !
    mue = zero
    do k=1,kl
       do j=2,jl
          do i=2,il


             ! Set the value of the porosity. If not zero, it is set
             ! to average the eddy-viscosity and to take the factor
             ! rFilv into account.

             por = half*rFilv
             if(porK(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i,j,k+1))
             mue = por*(rev(i,j,k) + rev(i,j,k+1))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) + gamma(i,j,k+1)) - one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the gradients at the face by averaging the four
             ! nodal values.

             u_x = fourth*(ux(i-1,j-1,k) + ux(i,j-1,k) &
                  +         ux(i-1,j,  k) + ux(i,j,  k))
             u_y = fourth*(uy(i-1,j-1,k) + uy(i,j-1,k) &
                  +         uy(i-1,j,  k) + uy(i,j,  k))
             u_z = fourth*(uz(i-1,j-1,k) + uz(i,j-1,k) &
                  +         uz(i-1,j,  k) + uz(i,j,  k))

             v_x = fourth*(vx(i-1,j-1,k) + vx(i,j-1,k) &
                  +         vx(i-1,j,  k) + vx(i,j,  k))
             v_y = fourth*(vy(i-1,j-1,k) + vy(i,j-1,k) &
                  +         vy(i-1,j,  k) + vy(i,j,  k))
             v_z = fourth*(vz(i-1,j-1,k) + vz(i,j-1,k) &
                  +         vz(i-1,j,  k) + vz(i,j,  k))

             w_x = fourth*(wx(i-1,j-1,k) + wx(i,j-1,k) &
                  +         wx(i-1,j,  k) + wx(i,j,  k))
             w_y = fourth*(wy(i-1,j-1,k) + wy(i,j-1,k) &
                  +         wy(i-1,j,  k) + wy(i,j,  k))
             w_z = fourth*(wz(i-1,j-1,k) + wz(i,j-1,k) &
                  +         wz(i-1,j,  k) + wz(i,j,  k))

             q_x = fourth*(qx(i-1,j-1,k) + qx(i,j-1,k) &
                  +         qx(i-1,j,  k) + qx(i,j,  k))
             q_y = fourth*(qy(i-1,j-1,k) + qy(i,j-1,k) &
                  +         qy(i-1,j,  k) + qy(i,j,  k))
             q_z = fourth*(qz(i-1,j-1,k) + qz(i,j-1,k) &
                  +         qz(i-1,j,  k) + qz(i,j,  k))


             ! The gradients in the normal direction are corrected, such
             ! that no averaging takes places here.
             ! First determine the vector in the direction from the
             ! cell center k to cell center k+1.

             ssx = eighth*(x(i-1,j-1,k+1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i-1,j,  k+1,1) - x(i-1,j,  k-1,1) &
                  +         x(i,  j-1,k+1,1) - x(i,  j-1,k-1,1) &
                  +         x(i,  j,  k+1,1) - x(i,  j,  k-1,1))
             ssy = eighth*(x(i-1,j-1,k+1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i-1,j,  k+1,2) - x(i-1,j,  k-1,2) &
                  +         x(i,  j-1,k+1,2) - x(i,  j-1,k-1,2) &
                  +         x(i,  j,  k+1,2) - x(i,  j,  k-1,2))
             ssz = eighth*(x(i-1,j-1,k+1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i-1,j,  k+1,3) - x(i-1,j,  k-1,3) &
                  +         x(i,  j-1,k+1,3) - x(i,  j-1,k-1,3) &
                  +         x(i,  j,  k+1,3) - x(i,  j,  k-1,3))

             ! Determine the length of this vector and create the
             ! unit normal.

             snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = snrm*ssx
             ssy = snrm*ssy
             ssz = snrm*ssz

             ! Correct the gradients.

             corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                  - (w(i,j,k+1,ivx) - w(i,j,k,ivx))*snrm
             u_x  = u_x - corr*ssx
             u_y  = u_y - corr*ssy
             u_z  = u_z - corr*ssz

             corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                  - (w(i,j,k+1,ivy) - w(i,j,k,ivy))*snrm
             v_x  = v_x - corr*ssx
             v_y  = v_y - corr*ssy
             v_z  = v_z - corr*ssz

             corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                  - (w(i,j,k+1,ivz) - w(i,j,k,ivz))*snrm
             w_x  = w_x - corr*ssx
             w_y  = w_y - corr*ssy
             w_z  = w_z - corr*ssz

             corr = q_x*ssx + q_y*ssy + q_z*ssz &
                  + (aa(i,j,k+1) - aa(i,j,k))*snrm
             q_x  = q_x - corr*ssx
             q_y  = q_y - corr*ssy
             q_z  = q_z - corr*ssz

             ! Compute the stress tensor and the heat flux vector.
             ! We remove the viscosity from the stress tensor (tau)
             ! to define tauS since we still need to separate between
             ! laminar and turbulent stress for QCR.
             ! Therefore, laminar tau = mue*tauS, turbulent
             ! tau = mue*tauS, and total tau = mut*tauS.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxxS = two*u_x - fracDiv
             tauyyS = two*v_y - fracDiv
             tauzzS = two*w_z - fracDiv

             tauxyS = u_y + v_x
             tauxzS = u_z + w_x
             tauyzS = v_z + w_y

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Add QCR corrections if necessary
             if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
                     v_x*v_x + v_y*v_y + v_z*v_z + &
                     w_x*w_x + w_y*w_y + w_z*w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue*Ccr1/den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact*(Wxy*tauxyS + Wxz*tauxzS)*two
                eyy = fact*(Wyx*tauxyS + Wyz*tauyzS)*two
                ezz = fact*(Wzx*tauxzS + Wzy*tauyzS)*two

                exy = fact*(Wxy*tauyyS + Wxz*tauyzS + &
                     Wyx*tauxxS + Wyz*tauxzS)
                exz = fact*(Wxy*tauyzS + Wxz*tauzzS + &
                     Wzx*tauxxS + Wzy*tauxyS)
                eyz = fact*(Wyx*tauxzS + Wyz*tauzzS + &
                     Wzx*tauxyS + Wzy*tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut*tauxxS - exx
                tauyy = mut*tauyyS - eyy
                tauzz = mut*tauzzS - ezz
                tauxy = mut*tauxyS - exy
                tauxz = mut*tauxzS - exz
                tauyz = mut*tauyzS - eyz

             else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut*tauxxS
                tauyy = mut*tauyyS
                tauzz = mut*tauzzS
                tauxy = mut*tauxyS
                tauxz = mut*tauxzS
                tauyz = mut*tauyzS

             end if

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i,j,k+1,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i,j,k+1,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i,j,k+1,ivz))

             ! Compute the viscous fluxes for this k-face.

             fmx   = tauxx*sk(i,j,k,1) + tauxy*sk(i,j,k,2) &
                  + tauxz*sk(i,j,k,3)
             fmy   = tauxy*sk(i,j,k,1) + tauyy*sk(i,j,k,2) &
                  + tauyz*sk(i,j,k,3)
             fmz   = tauxz*sk(i,j,k,1) + tauyz*sk(i,j,k,2) &
                  + tauzz*sk(i,j,k,3)
             frhoE =         (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1)
             frhoE = frhoE + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2)
             frhoE = frhoE + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3)
             frhoE = frhoE -  q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

             ! Update the residuals of cell k and k+1.

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

             fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   + fmx
             fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   + fmy
             fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   + fmz
             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + frhoE

             ! Temporarily store the shear stress and heat flux, even
             ! if we won't need it. This can still vectorize

             if (k == 1) then
                tmpStore(1, i, j, 1) = tauxx
                tmpStore(2, i, j, 1) = tauyy
                tmpStore(3, i, j, 1) = tauzz
                tmpStore(4, i, j, 1) = tauxy
                tmpStore(5, i, j, 1) = tauxz
                tmpStore(6, i, j, 1) = tauyz

                tmpStore(7, i, j, 1) = q_x
                tmpStore(8, i, j, 1) = q_y
                tmpStore(9, i, j, 1) = q_z
             end if

             if (k == kl) then
                tmpStore(1, i, j, 2) = tauxx
                tmpStore(2, i, j, 2) = tauyy
                tmpStore(3, i, j, 2) = tauzz
                tmpStore(4, i, j, 2) = tauxy
                tmpStore(5, i, j, 2) = tauxz
                tmpStore(6, i, j, 2) = tauyz

                tmpStore(7, i, j, 2) = q_x
                tmpStore(8, i, j, 2) = q_y
                tmpStore(9, i, j, 2) = q_z
             end if

          end do
       enddo
    end do

    ! Save into the subface if necessary
    if (storeWall) then

       origKMin: if (kk-1 == 1) then
          do j=2, jl
             do i=2, il
                io = i + ii - 2
                jo = j + jj - 2

                if (viscKminPointer(io, jo) > 0) then
                   viscSubface(viscKminPointer(io, jo))%tau(io, jo, :) = tmpStore(1:6, i, j, 1)
                   viscSubface(viscKminPointer(io, jo))%q(io, jo, :) = tmpStore(7:9, i, j, 1)
                endif
             end do
          end do
       end if origKMin

       origKMax: if (kk + nz - 1 == bkl) then
          do j=2, jl
             do i=2, il
                io = i + ii - 2
                jo = j + jj - 2
                if(viscKmaxPointer(io, jo) > 0) then
                   viscSubface(viscKmaxPointer(io, jo))%tau(io, jo, :) = tmpStore(1:6, i, j, 2)
                   viscSubface(viscKmaxPointer(io, jo))%q(io, jo, :) = tmpStore(7:9, i, j, 2)
                endif
             end do
          end do
       end if origKMax
    end if
    !
    !         Viscous fluxes in the j-direction.
    !
    do k=2,kl
       do j=1,jl
          do i=2,il


             ! Set the value of the porosity. If not zero, it is set
             ! to average the eddy-viscosity and to take the factor
             ! rFilv into account.

             por = half*rFilv
             if(porJ(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i,j+1,k))
             mue = por*(rev(i,j,k) + rev(i,j+1,k))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) + gamma(i,j+1,k)) - one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the gradients at the face by averaging the four
             ! nodal values.

             u_x = fourth*(ux(i-1,j,k-1) + ux(i,j,k-1) &
                  +        ux(i-1,j,k  ) + ux(i,j,k  ))
             u_y = fourth*(uy(i-1,j,k-1) + uy(i,j,k-1) &
                  +        uy(i-1,j,k  ) + uy(i,j,k  ))
             u_z = fourth*(uz(i-1,j,k-1) + uz(i,j,k-1) &
                  +        uz(i-1,j,k  ) + uz(i,j,k  ))

             v_x = fourth*(vx(i-1,j,k-1) + vx(i,j,k-1) &
                  +        vx(i-1,j,k  ) + vx(i,j,k  ))
             v_y = fourth*(vy(i-1,j,k-1) + vy(i,j,k-1) &
                  +        vy(i-1,j,k  ) + vy(i,j,k  ))
             v_z = fourth*(vz(i-1,j,k-1) + vz(i,j,k-1) &
                  +        vz(i-1,j,k  ) + vz(i,j,k  ))

             w_x = fourth*(wx(i-1,j,k-1) + wx(i,j,k-1) &
                  +        wx(i-1,j,k  ) + wx(i,j,k  ))
             w_y = fourth*(wy(i-1,j,k-1) + wy(i,j,k-1) &
                  +        wy(i-1,j,k  ) + wy(i,j,k  ))
             w_z = fourth*(wz(i-1,j,k-1) + wz(i,j,k-1) &
                  +        wz(i-1,j,k  ) + wz(i,j,k  ))

             q_x = fourth*(qx(i-1,j,k-1) + qx(i,j,k-1) &
                  +        qx(i-1,j,k  ) + qx(i,j,k  ))
             q_y = fourth*(qy(i-1,j,k-1) + qy(i,j,k-1) &
                  +        qy(i-1,j,k  ) + qy(i,j,k  ))
             q_z = fourth*(qz(i-1,j,k-1) + qz(i,j,k-1) &
                  +        qz(i-1,j,k  ) + qz(i,j,k  ))

             ! The gradients in the normal direction are corrected, such
             ! that no averaging takes places here.
             ! First determine the vector in the direction from the
             ! cell center j to cell center j+1.

             ssx = eighth*(x(i-1,j+1,k-1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i-1,j+1,k,  1) - x(i-1,j-1,k,  1) &
                  +         x(i,  j+1,k-1,1) - x(i,  j-1,k-1,1) &
                  +         x(i,  j+1,k,  1) - x(i,  j-1,k,  1))
             ssy = eighth*(x(i-1,j+1,k-1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i-1,j+1,k,  2) - x(i-1,j-1,k,  2) &
                  +         x(i,  j+1,k-1,2) - x(i,  j-1,k-1,2) &
                  +         x(i,  j+1,k,  2) - x(i,  j-1,k,  2))
             ssz = eighth*(x(i-1,j+1,k-1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i-1,j+1,k,  3) - x(i-1,j-1,k,  3) &
                  +         x(i,  j+1,k-1,3) - x(i,  j-1,k-1,3) &
                  +         x(i,  j+1,k,  3) - x(i,  j-1,k,  3))

             ! Determine the length of this vector and create the
             ! unit normal.

             snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = snrm*ssx
             ssy = snrm*ssy
             ssz = snrm*ssz

             ! Correct the gradients.

             corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                  - (w(i,j+1,k,ivx) - w(i,j,k,ivx))*snrm
             u_x  = u_x - corr*ssx
             u_y  = u_y - corr*ssy
             u_z  = u_z - corr*ssz

             corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                  - (w(i,j+1,k,ivy) - w(i,j,k,ivy))*snrm
             v_x  = v_x - corr*ssx
             v_y  = v_y - corr*ssy
             v_z  = v_z - corr*ssz

             corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                  - (w(i,j+1,k,ivz) - w(i,j,k,ivz))*snrm
             w_x  = w_x - corr*ssx
             w_y  = w_y - corr*ssy
             w_z  = w_z - corr*ssz

             corr = q_x*ssx + q_y*ssy + q_z*ssz &
                  + (aa(i,j+1,k) - aa(i,j,k))*snrm
             q_x  = q_x - corr*ssx
             q_y  = q_y - corr*ssy
             q_z  = q_z - corr*ssz

             ! Compute the stress tensor and the heat flux vector.
             ! We remove the viscosity from the stress tensor (tau)
             ! to define tauS since we still need to separate between
             ! laminar and turbulent stress for QCR.
             ! Therefore, laminar tau = mue*tauS, turbulent
             ! tau = mue*tauS, and total tau = mut*tauS.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxxS = two*u_x - fracDiv
             tauyyS = two*v_y - fracDiv
             tauzzS = two*w_z - fracDiv

             tauxyS = u_y + v_x
             tauxzS = u_z + w_x
             tauyzS = v_z + w_y

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Add QCR corrections if necessary
             if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
                     v_x*v_x + v_y*v_y + v_z*v_z + &
                     w_x*w_x + w_y*w_y + w_z*w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue*Ccr1/den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact*(Wxy*tauxyS + Wxz*tauxzS)*two
                eyy = fact*(Wyx*tauxyS + Wyz*tauyzS)*two
                ezz = fact*(Wzx*tauxzS + Wzy*tauyzS)*two

                exy = fact*(Wxy*tauyyS + Wxz*tauyzS + &
                     Wyx*tauxxS + Wyz*tauxzS)
                exz = fact*(Wxy*tauyzS + Wxz*tauzzS + &
                     Wzx*tauxxS + Wzy*tauxyS)
                eyz = fact*(Wyx*tauxzS + Wyz*tauzzS + &
                     Wzx*tauxyS + Wzy*tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut*tauxxS - exx
                tauyy = mut*tauyyS - eyy
                tauzz = mut*tauzzS - ezz
                tauxy = mut*tauxyS - exy
                tauxz = mut*tauxzS - exz
                tauyz = mut*tauyzS - eyz

             else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut*tauxxS
                tauyy = mut*tauyyS
                tauzz = mut*tauzzS
                tauxy = mut*tauxyS
                tauxz = mut*tauxzS
                tauyz = mut*tauyzS

             end if

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i,j+1,k,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i,j+1,k,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i,j+1,k,ivz))

             ! Compute the viscous fluxes for this j-face.

             fmx   = tauxx*sj(i,j,k,1) + tauxy*sj(i,j,k,2) &
                  + tauxz*sj(i,j,k,3)
             fmy   = tauxy*sj(i,j,k,1) + tauyy*sj(i,j,k,2) &
                  + tauyz*sj(i,j,k,3)
             fmz   = tauxz*sj(i,j,k,1) + tauyz*sj(i,j,k,2) &
                  + tauzz*sj(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sj(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sj(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sj(i,j,k,3) &
                  - q_x*sj(i,j,k,1) - q_y*sj(i,j,k,2) - q_z*sj(i,j,k,3)

             ! Update the residuals of cell j and j+1.

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

             fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   + fmx
             fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   + fmy
             fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   + fmz
             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + frhoE

             ! Temporarily store the shear stress and heat flux, even
             ! if we won't need it. This can still vectorize

             if (j == 1) then
                tmpStore(1, i, k, 1) = tauxx
                tmpStore(2, i, k, 1) = tauyy
                tmpStore(3, i, k, 1) = tauzz
                tmpStore(4, i, k, 1) = tauxy
                tmpStore(5, i, k, 1) = tauxz
                tmpStore(6, i, k, 1) = tauyz

                tmpStore(7, i, k, 1) = q_x
                tmpStore(8, i, k, 1) = q_y
                tmpStore(9, i, k, 1) = q_z
             end if

             if (j == jl) then
                tmpStore(1, i, k, 2) = tauxx
                tmpStore(2, i, k, 2) = tauyy
                tmpStore(3, i, k, 2) = tauzz
                tmpStore(4, i, k, 2) = tauxy
                tmpStore(5, i, k, 2) = tauxz
                tmpStore(6, i, k, 2) = tauyz

                tmpStore(7, i, k, 2) = q_x
                tmpStore(8, i, k, 2) = q_y
                tmpStore(9, i, k, 2) = q_z
             end if
          enddo
       enddo
    enddo
    ! Save into the subface if necessary
    if (storeWall) then
       origJMin: if (jj-1 == 1) then
          do k=2, kl
             do i=2, il
                io = i + ii - 2
                ko = k + kk - 2

                if (viscJminPointer(io, ko) > 0) then
                   viscSubface(viscJminPointer(io, ko))%tau(io, ko, :) = tmpStore(1:6, i, k, 1)
                   viscSubface(viscJminPointer(io, ko))%q(io, ko, :) = tmpStore(7:9, i, k, 1)
                endif
             end do
          end do
       end if origJMin

       origJMax: if (jj + ny - 1 == bjl) then
          do k=2, kl
             do i=2, il
                io = i + ii - 2
                ko = k + kk - 2
                if(viscJmaxPointer(io, ko) > 0) then
                   viscSubface(viscJmaxPointer(io, ko))%tau(io, ko, :) = tmpStore(1:6, i, k, 2)
                   viscSubface(viscJmaxPointer(io, ko))%q(io, ko, :) = tmpStore(7:9, i, k, 2)
                endif
             end do
          end do
       end if origJMax
    end if
    !
    !         Viscous fluxes in the i-direction.
    !
    do k=2, kl
       do j=2, jl
          do i=1, il
             ! Set the value of the porosity. If not zero, it is set
             ! to average the eddy-viscosity and to take the factor
             ! rFilv into account.

             por = half*rFilv
             if(porI(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i+1,j,k))
             mue = por*(rev(i,j,k) + rev(i+1,j,k))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) + gamma(i+1,j,k)) - one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the gradients at the face by averaging the four
             ! nodal values.

             u_x = fourth*(ux(i,j-1,k-1) + ux(i,j,k-1) &
                  +        ux(i,j-1,k  ) + ux(i,j,k  ))
             u_y = fourth*(uy(i,j-1,k-1) + uy(i,j,k-1) &
                  +        uy(i,j-1,k  ) + uy(i,j,k  ))
             u_z = fourth*(uz(i,j-1,k-1) + uz(i,j,k-1) &
                  +        uz(i,j-1,k  ) + uz(i,j,k  ))

             v_x = fourth*(vx(i,j-1,k-1) + vx(i,j,k-1) &
                  +        vx(i,j-1,k  ) + vx(i,j,k  ))
             v_y = fourth*(vy(i,j-1,k-1) + vy(i,j,k-1) &
                  +        vy(i,j-1,k  ) + vy(i,j,k  ))
             v_z = fourth*(vz(i,j-1,k-1) + vz(i,j,k-1) &
                  +        vz(i,j-1,k  ) + vz(i,j,k  ))

             w_x = fourth*(wx(i,j-1,k-1) + wx(i,j,k-1) &
                  +        wx(i,j-1,k  ) + wx(i,j,k  ))
             w_y = fourth*(wy(i,j-1,k-1) + wy(i,j,k-1) &
                  +        wy(i,j-1,k  ) + wy(i,j,k  ))
             w_z = fourth*(wz(i,j-1,k-1) + wz(i,j,k-1) &
                  +        wz(i,j-1,k  ) + wz(i,j,k  ))

             q_x = fourth*(qx(i,j-1,k-1) + qx(i,j,k-1) &
                  +        qx(i,j-1,k  ) + qx(i,j,k  ))
             q_y = fourth*(qy(i,j-1,k-1) + qy(i,j,k-1) &
                  +        qy(i,j-1,k  ) + qy(i,j,k  ))
             q_z = fourth*(qz(i,j-1,k-1) + qz(i,j,k-1) &
                  +        qz(i,j-1,k  ) + qz(i,j,k  ))

             ! The gradients in the normal direction are corrected, such
             ! that no averaging takes places here.
             ! First determine the vector in the direction from the
             ! cell center i to cell center i+1.

             ssx = eighth*(x(i+1,j-1,k-1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i+1,j-1,k,  1) - x(i-1,j-1,k,  1) &
                  +         x(i+1,j,  k-1,1) - x(i-1,j,  k-1,1) &
                  +         x(i+1,j,  k,  1) - x(i-1,j,  k,  1))
             ssy = eighth*(x(i+1,j-1,k-1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i+1,j-1,k,  2) - x(i-1,j-1,k,  2) &
                  +         x(i+1,j,  k-1,2) - x(i-1,j,  k-1,2) &
                  +         x(i+1,j,  k,  2) - x(i-1,j,  k,  2))
             ssz = eighth*(x(i+1,j-1,k-1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i+1,j-1,k,  3) - x(i-1,j-1,k,  3) &
                  +         x(i+1,j,  k-1,3) - x(i-1,j,  k-1,3) &
                  +         x(i+1,j,  k,  3) - x(i-1,j,  k,  3))

             ! Determine the length of this vector and create the
             ! unit normal.

             snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = snrm*ssx
             ssy = snrm*ssy
             ssz = snrm*ssz

             ! Correct the gradients.

             corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                  - (w(i+1,j,k,ivx) - w(i,j,k,ivx))*snrm
             u_x  = u_x - corr*ssx
             u_y  = u_y - corr*ssy
             u_z  = u_z - corr*ssz

             corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                  - (w(i+1,j,k,ivy) - w(i,j,k,ivy))*snrm
             v_x  = v_x - corr*ssx
             v_y  = v_y - corr*ssy
             v_z  = v_z - corr*ssz

             corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                  - (w(i+1,j,k,ivz) - w(i,j,k,ivz))*snrm
             w_x  = w_x - corr*ssx
             w_y  = w_y - corr*ssy
             w_z  = w_z - corr*ssz

             corr = q_x*ssx + q_y*ssy + q_z*ssz &
                  + (aa(i+1,j,k) - aa(i,j,k))*snrm
             q_x  = q_x - corr*ssx
             q_y  = q_y - corr*ssy
             q_z  = q_z - corr*ssz

             ! Compute the stress tensor and the heat flux vector.
             ! We remove the viscosity from the stress tensor (tau)
             ! to define tauS since we still need to separate between
             ! laminar and turbulent stress for QCR.
             ! Therefore, laminar tau = mue*tauS, turbulent
             ! tau = mue*tauS, and total tau = mut*tauS.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxxS = two*u_x - fracDiv
             tauyyS = two*v_y - fracDiv
             tauzzS = two*w_z - fracDiv

             tauxyS = u_y + v_x
             tauxzS = u_z + w_x
             tauyzS = v_z + w_y

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Add QCR corrections if necessary
             if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
                     v_x*v_x + v_y*v_y + v_z*v_z + &
                     w_x*w_x + w_y*w_y + w_z*w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue*Ccr1/den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact*(Wxy*tauxyS + Wxz*tauxzS)*two
                eyy = fact*(Wyx*tauxyS + Wyz*tauyzS)*two
                ezz = fact*(Wzx*tauxzS + Wzy*tauyzS)*two

                exy = fact*(Wxy*tauyyS + Wxz*tauyzS + &
                     Wyx*tauxxS + Wyz*tauxzS)
                exz = fact*(Wxy*tauyzS + Wxz*tauzzS + &
                     Wzx*tauxxS + Wzy*tauxyS)
                eyz = fact*(Wyx*tauxzS + Wyz*tauzzS + &
                     Wzx*tauxyS + Wzy*tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut*tauxxS - exx
                tauyy = mut*tauyyS - eyy
                tauzz = mut*tauzzS - ezz
                tauxy = mut*tauxyS - exy
                tauxz = mut*tauxzS - exz
                tauyz = mut*tauyzS - eyz

             else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut*tauxxS
                tauyy = mut*tauyyS
                tauzz = mut*tauzzS
                tauxy = mut*tauxyS
                tauxz = mut*tauxzS
                tauyz = mut*tauyzS

             end if

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i+1,j,k,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i+1,j,k,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i+1,j,k,ivz))

             ! Compute the viscous fluxes for this i-face.

             fmx   = tauxx*si(i,j,k,1) + tauxy*si(i,j,k,2) &
                  + tauxz*si(i,j,k,3)
             fmy   = tauxy*si(i,j,k,1) + tauyy*si(i,j,k,2) &
                  + tauyz*si(i,j,k,3)
             fmz   = tauxz*si(i,j,k,1) + tauyz*si(i,j,k,2) &
                  + tauzz*si(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*si(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*si(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*si(i,j,k,3) &
                  - q_x*si(i,j,k,1) - q_y*si(i,j,k,2) - q_z*si(i,j,k,3)

             ! Update the residuals of cell i and i+1.
             fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   + fmx
             fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   + fmy
             fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   + fmz
             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + frhoE

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE


             ! Temporarily store the shear stress and heat flux, even
             ! if we won't need it. This can still vectorize

             if (i == 1) then
                tmpStore(1, j, k, 1) = tauxx
                tmpStore(2, j, k, 1) = tauyy
                tmpStore(3, j, k, 1) = tauzz
                tmpStore(4, j, k, 1) = tauxy
                tmpStore(5, j, k, 1) = tauxz
                tmpStore(6, j, k, 1) = tauyz

                tmpStore(7, j, k, 1) = q_x
                tmpStore(8, j, k, 1) = q_y
                tmpStore(9, j, k, 1) = q_z
             end if

             if (i == il) then
                tmpStore(1, j, k, 2) = tauxx
                tmpStore(2, j, k, 2) = tauyy
                tmpStore(3, j, k, 2) = tauzz
                tmpStore(4, j, k, 2) = tauxy
                tmpStore(5, j, k, 2) = tauxz
                tmpStore(6, j, k, 2) = tauyz

                tmpStore(7, j, k, 2) = q_x
                tmpStore(8, j, k, 2) = q_y
                tmpStore(9, j, k, 2) = q_z
             end if
          enddo
       enddo
    enddo
    ! Save into the subface if necessary
    if (storeWall) then
       origIMin: if (ii-1 == 1) then
          do k=2, kl
             do j=2, jl
                jo = j + jj - 2
                ko = k + kk - 2

                if (viscIminPointer(jo, ko) > 0) then
                   viscSubface(viscIminPointer(jo, ko))%tau(jo, ko, :) = tmpStore(1:6, j, k, 1)
                   viscSubface(viscIminPointer(jo, ko))%q(jo, ko, :) = tmpStore(7:9, j, k, 1)
                endif
             end do
          end do
       end if origIMin

       origIMax: if (ii + nx - 1 == bil) then
          do k=2, kl
             do j=2, jl
                jo = j + jj - 2
                ko = k + kk - 2
                if(viscImaxPointer(jo, ko) > 0) then
                   viscSubface(viscImaxPointer(jo, ko))%tau(jo, ko, :) = tmpStore(1:6, j, k, 2)
                   viscSubface(viscImaxPointer(jo, ko))%q(jo, ko, :) = tmpStore(7:9, j, k, 2)
                endif
             end do
          end do
       end if origIMax
    end if
  end subroutine viscousFlux

  subroutine viscousFluxApprox

    use constants
    use flowVarRefState
    use inputPhysics
    use iteration
    implicit none
    !
    !      Local parameter.
    !
    real(kind=realType), parameter :: twoThird = two*third
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k
    integer(kind=intType) :: ii, jj, kk

    real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
    real(kind=realType) :: gm1, factLamHeat, factTurbHeat
    real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
    real(kind=realType) :: q_x, q_y, q_z, ubar, vbar, wbar
    real(kind=realType) :: corr, ssx, ssy, ssz, ss, fracDiv
    real(kind=realType) :: tauxx, tauyy, tauzz
    real(kind=realType) :: tauxy, tauxz, tauyz
    real(kind=realType) :: fmx, fmy, fmz, frhoE
    real(kind=realType) :: dd
    logical :: correctForK

    mue = zero
    rFilv = rFil

    ! Viscous fluxes in the I-direction

    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the vector from the center of cell i to cell i+1
             ssx = eighth*(x(i+1,j-1,k-1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i+1,j-1,k,  1) - x(i-1,j-1,k,  1) &
                  +         x(i+1,j,  k-1,1) - x(i-1,j,  k-1,1) &
                  +         x(i+1,j,  k,  1) - x(i-1,j,  k,  1))
             ssy = eighth*(x(i+1,j-1,k-1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i+1,j-1,k,  2) - x(i-1,j-1,k,  2) &
                  +         x(i+1,j,  k-1,2) - x(i-1,j,  k-1,2) &
                  +         x(i+1,j,  k,  2) - x(i-1,j,  k,  2))
             ssz = eighth*(x(i+1,j-1,k-1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i+1,j-1,k,  3) - x(i-1,j-1,k,  3) &
                  +         x(i+1,j,  k-1,3) - x(i-1,j,  k-1,3) &
                  +         x(i+1,j,  k,  3) - x(i-1,j,  k,  3))

             ! And determine one/ length of vector squared
             ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = ss*ssx
             ssy = ss*ssy
             ssz = ss*ssz

             ! Now compute each gradient
             dd = w(i+1, j, k, ivx)-w(i, j, k, ivx)
             u_x = dd*ssx
             u_y = dd*ssy
             u_z = dd*ssz

             dd = w(i+1, j, k, ivy)-w(i, j, k, ivy)
             v_x = dd*ssx
             v_y = dd*ssy
             v_z = dd*ssz

             dd = w(i+1, j, k, ivz)-w(i, j, k, ivz)
             w_x = dd*ssx
             w_y = dd*ssy
             w_z = dd*ssz

             dd = aa(i+1, j, k)-aa(i, j, k)
             q_x = -dd*ssx
             q_y = -dd*ssy
             q_z = -dd*ssz

             por = half*rFilv
             if(porI(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i+1,j,k))
             mue = por*(rev(i,j,k) + rev(i+1,j,k))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) +gamma(i+1,j,k))- one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the stress tensor and the heat flux vector.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxx = mut*(two*u_x - fracDiv)
             tauyy = mut*(two*v_y - fracDiv)
             tauzz = mut*(two*w_z - fracDiv)

             tauxy = mut*(u_y + v_x)
             tauxz = mut*(u_z + w_x)
             tauyz = mut*(v_z + w_y)

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i+1,j,k,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i+1,j,k,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i+1,j,k,ivz))

             ! Compute the viscous fluxes for this i-face.

             fmx   = tauxx*si(i,j,k,1) + tauxy*si(i,j,k,2) + tauxz*si(i,j,k,3)
             fmy   = tauxy*si(i,j,k,1) + tauyy*si(i,j,k,2) + tauyz*si(i,j,k,3)
             fmz   = tauxz*si(i,j,k,1) + tauyz*si(i,j,k,2) + tauzz*si(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*si(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*si(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*si(i,j,k,3) &
                  - q_x*si(i,j,k,1) - q_y*si(i,j,k,2) - q_z*si(i,j,k,3)

             ! Update the residuals of cell i and i+1.
             fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   + fmx
             fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   + fmy
             fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   + fmz
             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + frhoE

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE


          end do
       end do
    end do

    ! Viscous fluxes in the J-direction

    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Compute the vector from the center of cell j to cell j+1
             ssx = eighth*(x(i-1,j+1,k-1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i-1,j+1,k,  1) - x(i-1,j-1,k,  1) &
                  +         x(i,  j+1,k-1,1) - x(i,  j-1,k-1,1) &
                  +         x(i,  j+1,k,  1) - x(i,  j-1,k,  1))
             ssy = eighth*(x(i-1,j+1,k-1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i-1,j+1,k,  2) - x(i-1,j-1,k,  2) &
                  +         x(i,  j+1,k-1,2) - x(i,  j-1,k-1,2) &
                  +         x(i,  j+1,k,  2) - x(i,  j-1,k,  2))
             ssz = eighth*(x(i-1,j+1,k-1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i-1,j+1,k,  3) - x(i-1,j-1,k,  3) &
                  +         x(i,  j+1,k-1,3) - x(i,  j-1,k-1,3) &
                  +         x(i,  j+1,k,  3) - x(i,  j-1,k,  3))

             ! And determine one/ length of vector squared
             ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = ss*ssx
             ssy = ss*ssy
             ssz = ss*ssz

             ! Now compute each gradient
             dd = w(i, j+1, k, ivx)-w(i, j, k, ivx)
             u_x = dd*ssx
             u_y = dd*ssy
             u_z = dd*ssz

             dd = w(i, j+1, k, ivy)-w(i, j, k, ivy)
             v_x = dd*ssx
             v_y = dd*ssy
             v_z = dd*ssz

             dd = w(i, j+1, k, ivz)-w(i, j, k, ivz)
             w_x = dd*ssx
             w_y = dd*ssy
             w_z = dd*ssz

             dd = aa(i, j+1, k)-aa(i, j, k)
             q_x = -dd*ssx
             q_y = -dd*ssy
             q_z = -dd*ssz

             por = half*rFilv
             if(porJ(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i,j+1,k))
             mue = por*(rev(i,j,k) + rev(i,j+1,k))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) + gamma(i,j+1,k)) - one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the stress tensor and the heat flux vector.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxx = mut*(two*u_x - fracDiv)
             tauyy = mut*(two*v_y - fracDiv)
             tauzz = mut*(two*w_z - fracDiv)

             tauxy = mut*(u_y + v_x)
             tauxz = mut*(u_z + w_x)
             tauyz = mut*(v_z + w_y)

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i,j+1,k,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i,j+1,k,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i,j+1,k,ivz))

             ! Compute the viscous fluxes for this j-face.

             fmx   = tauxx*sj(i,j,k,1) + tauxy*sj(i,j,k,2) + tauxz*sj(i,j,k,3)
             fmy   = tauxy*sj(i,j,k,1) + tauyy*sj(i,j,k,2) + tauyz*sj(i,j,k,3)
             fmz   = tauxz*sj(i,j,k,1) + tauyz*sj(i,j,k,2) + tauzz*sj(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sj(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sj(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sj(i,j,k,3) &
                  - q_x*sj(i,j,k,1) - q_y*sj(i,j,k,2) - q_z*sj(i,j,k,3)

             ! Update the residuals of cell j and j+1.

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

             fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   + fmx
             fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   + fmy
             fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   + fmz
             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + frhoE

          end do
       end do
    end do

    ! Viscous fluxes in the K-direction

    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the vector from the center of cell k to cell k+1
             ssx = eighth*(x(i-1,j-1,k+1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i-1,j,  k+1,1) - x(i-1,j,  k-1,1) &
                  +         x(i,  j-1,k+1,1) - x(i,  j-1,k-1,1) &
                  +         x(i,  j,  k+1,1) - x(i,  j,  k-1,1))
             ssy = eighth*(x(i-1,j-1,k+1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i-1,j,  k+1,2) - x(i-1,j,  k-1,2) &
                  +         x(i,  j-1,k+1,2) - x(i,  j-1,k-1,2) &
                  +         x(i,  j,  k+1,2) - x(i,  j,  k-1,2))
             ssz = eighth*(x(i-1,j-1,k+1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i-1,j,  k+1,3) - x(i-1,j,  k-1,3) &
                  +         x(i,  j-1,k+1,3) - x(i,  j-1,k-1,3) &
                  +         x(i,  j,  k+1,3) - x(i,  j,  k-1,3))
             ! And determine one/ length of vector squared
             ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = ss*ssx
             ssy = ss*ssy
             ssz = ss*ssz

             ! Now compute each gradient
             dd = w(i, j, k+1, ivx)-w(i, j, k, ivx)
             u_x = dd*ssx
             u_y = dd*ssy
             u_z = dd*ssz

             dd = w(i, j, k+1, ivy)-w(i, j, k, ivy)
             v_x = dd*ssx
             v_y = dd*ssy
             v_z = dd*ssz

             dd = w(i, j, k+1, ivz)-w(i, j, k, ivz)
             w_x = dd*ssx
             w_y = dd*ssy
             w_z = dd*ssz

             dd = aa(i, j, k+1)-aa(i, j, k)
             q_x = -dd*ssx
             q_y = -dd*ssy
             q_z = -dd*ssz

             por = half*rFilv
             if(porK(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i,j,k+1))
             mue = por*(rev(i,j,k) + rev(i,j,k+1))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) + gamma(i,j,k+1)) - one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the stress tensor and the heat flux vector.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxx = mut*(two*u_x - fracDiv)
             tauyy = mut*(two*v_y - fracDiv)
             tauzz = mut*(two*w_z - fracDiv)

             tauxy = mut*(u_y + v_x)
             tauxz = mut*(u_z + w_x)
             tauyz = mut*(v_z + w_y)

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i,j,k+1,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i,j,k+1,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i,j,k+1,ivz))

             ! Compute the viscous fluxes for this j-face.

             fmx   = tauxx*sk(i,j,k,1) + tauxy*sk(i,j,k,2) + tauxz*sk(i,j,k,3)
             fmy   = tauxy*sk(i,j,k,1) + tauyy*sk(i,j,k,2) + tauyz*sk(i,j,k,3)
             fmz   = tauxz*sk(i,j,k,1) + tauyz*sk(i,j,k,2) + tauzz*sk(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3) &
                  - q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

             ! Update the residuals of cell j and j+1.

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

             fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   + fmx
             fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   + fmy
             fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   + fmz
             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + frhoE

          end do
       end do
    end do
  end subroutine viscousFluxApprox

  subroutine sumDwandFw

    ! ---------------------------------------------
    !                 Sum dw and fw/res scale
    ! ---------------------------------------------
    use constants
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use inputIteration, only : turbResScale
    implicit none

    ! Variables for final summing
    integer(kind=intType) :: nTurb, i, j, k, l
    real(kind=realType) :: oVol, rBlank

    nTurb = nt2-nt1+1
    do l=1,nwf
       do k=2, kl
          do j=2, jl
             do i=2, il
                rblank = max(real(iblank(i,j,k), realType), zero)
                dw(i, j, k, l) = (dw(i, j, k, l) + fw(i, j, k, l))*rBlank
             end do
          end do
       end do
    end do
  end subroutine sumDwandFw

  subroutine resScale

    use constants
    use flowVarRefState, only : nwf, nt1, nt2
    use inputIteration, only : turbResScale
    implicit none

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii, nTurb
    real(kind=realType) :: ovol

    ! Divide through by the reference volume
    nTurb = nt2-nt1+1
    do k=2,kl
       do j=2,jl
          do i=2,il

             oVol = one/volRef(i,j,k)
             dw(i, j, k, 1:nwf) = dw(i,j, k, 1:nwf)* ovol
             dw(i, j, k, nt1:nt2) = dw(i, j, k, nt1:nt2) * ovol * turbResScale(1:nTurb)
          enddo
       enddo
    enddo

  end subroutine resScale

end module blockette
