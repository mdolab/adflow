!
!      ******************************************************************
!      *                                                                *
!      * File:          initFlow.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-06-2003                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initFlow
!
!      ******************************************************************
!      *                                                                *
!      * initFlow sets the prescribed boundary data, allocates the      *
!      * memory for and initializes the flow variables. In case a       *
!      * restart is performed the owned variables are read from the     *
!      * previous solution file(s).                                     *
!      *                                                                *
!      ******************************************************************
!
       use bleedFlows
       use block
       use communication
       use inputIO
       use inputTimeSpectral
       use IOModule
       use iteration
       use restartMod
       use bcdata, only : setSupersonicInletFreeStream, setbcdataFineGrid,&
            setBCDataCoarseGrid, nonDimBoundData, setInletFreeStreamTurb
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: sps, level, nLevels

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of multigrid levels.

       nLevels = ubound(flowDoms,2)

       ! Allocate the memory for the prescribed boundary data at the
       ! boundary faces and determine the data for both the fine grid.

       groundLevel = 1

       call setBCDataFineGrid(.true.)

       ! As some boundary conditions can be treated in multiple ways,
       ! some memory allocated must be released again.

       call releaseExtraMemBcs

       ! Determine the reference state.
       call referenceState

       ! Determine the prescribed data on the coarse grid levels
       ! by interpolation.

       call setBCDataCoarseGrid

       ! Non-dimensionalize the boundary data.

       call nonDimBoundData

       ! Set the primitive variables to the free stream values for
       ! the supersonic inflow faces for which this data has not
       ! been prescribed.

       call setSupersonicInletFreeStream

       ! Set the turbulent quantities to the free stream values for
       ! the inflow faces for which this data has not been prescribed.
       ! These can be both subsonic and supersonic inflow faces.

       call setInletFreestreamTurb

       ! Determine for the time spectral mode the matrices for the
       ! time derivatives.
       call timeSpectralMatrices

       ! Loop over the number of spectral solutions to allocate
       ! the memory for the w-variables and p on the fine grid.
       do sps=1,nTimeIntervalsSpectral
         call allocMemFlovarPart1(sps, 1_intType)
       enddo

       ! Allocate the memory for the solution variables on the coarse
       ! grid levels and the memory for the dependent flow variables,
       ! residuals, etc, on all multigrid levels.
       do sps=1,nTimeIntervalsSpectral
         call allocMemFlovarPart2(sps, 1_intType)

         do level=2,nLevels
           call allocMemFlovarPart1(sps, level)
           call allocMemFlovarPart2(sps, level)
         enddo
       enddo

       ! Initialize free stream field
       call initFlowfield

       ! Initialize the prescribed boundary data for domain interfaces.
       ! This is just to avoid weird behavior in initDepvarAndHalos due
       ! to a wrong initialization. The actual values should be
       ! overwritten by the coupler.

       call initBCDataDomainInterfaces

       ! Initialize the dependent flow variables and the halo values.

       call initDepvarAndHalos(halosRead)
  
     end subroutine initFlow


  subroutine initBCDataDomainInterfaces
    !
    !      ******************************************************************
    !      *                                                                *
    !      * initBCDataDomainInterfaces initializes the prescribed boundary *
    !      * data for domain interfaces. It is just an initialization such  *
    !      * the initial halo computations behave normally. During the      *
    !      * actual computation this data should be renewed constantly by   *
    !      * the coupler.                                                   *
    !      *                                                                *
    !      ******************************************************************
    !
    use blockPointers
    use flowVarRefState
    use inputIteration
    use inputTimeSpectral
    use BCRoutines
    use utils, only : setPointers
    use flowUtils, only : computePtot, computeTtot
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, nn, mm, l, sps

    real(kind=realType) :: rho, vvx, vvy, vvz, pres, vtotInv

    ! Loop over the number of spectral solutions and blocks of
    ! the multigrid start level.

    spectralLoop: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, mgStartlevel, sps)

          ! Loop over the boundary condition subfaces of this block.

          bocos: do mm=1,nBocos

             ! Check for interface boundary condition types.

             select case (BCType(mm))

             case (DomainInterfaceAll)

                ! All data must be prescribed. Nullify the pointers
                ! and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the state vector to set the prescribed state.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                      BCData(mm)%rho(i,j)  = ww2(i,j,irho)
                      BCData(mm)%velx(i,j) = ww2(i,j,ivx)
                      BCData(mm)%vely(i,j) = ww2(i,j,ivy)
                      BCData(mm)%velz(i,j) = ww2(i,j,ivz)
                      BCData(mm)%ps(i,j)   = pp2(i,j)

                      do l=nt1,nt2
                         BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                      enddo
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceRhoUVW)

                ! Density, velocity and turbulent variables are
                ! prescribed. Nullify the pointers and set them to the
                ! correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the state vector to set the prescribed state.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                      BCData(mm)%rho(i,j)  = ww2(i,j,irho)
                      BCData(mm)%velx(i,j) = ww2(i,j,ivx)
                      BCData(mm)%vely(i,j) = ww2(i,j,ivy)
                      BCData(mm)%velz(i,j) = ww2(i,j,ivz)

                      do l=nt1,nt2
                         BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                      enddo
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceP)

                ! Pressure must be prescribed. Nullify the pointers
                ! and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the pressure to set the prescribed value.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                      BCData(mm)%ps(i,j) = pp2(i,j)
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceRho)

                ! Density must be prescribed. Nullify the pointers
                ! and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the density to set the prescribed value.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                      BCData(mm)%rho(i,j) = ww2(i,j,irho)
                   enddo
                enddo

                !=========================================================

             case (DomainInterfaceTotal)

                ! Total conditions must be prescribed. Nullify the
                ! pointers and set them to the correct subface.

                call setBCPointers(mm, .False.)

                ! Loop over the generic subface and simply extrapolate
                ! the total conditions and the turbulence variables
                ! to set the prescribed value.

                do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                      ! Store the variables a bit easier.

                      rho  = ww2(i,j,irho)
                      vvx   = ww2(i,j,ivx)
                      vvy   = ww2(i,j,ivy)
                      vvz   = ww2(i,j,ivz)
                      pres = pp2(i,j)

                      ! Compute the total pressure, total temperature
                      ! and total entahlpy.

                      call computePtot(rho, vvx, vvy, vvz, pres, &
                           BCData(mm)%ptInlet(i,j))
                      call computeTtot(rho, vvx, vvy, vvz, pres, &
                           BCData(mm)%ttInlet(i,j))
                      BCData(mm)%htInlet(i,j) = (ww2(i,j,irhoE) + pres) &
                           / rho

                      ! Determine the velocity direction.

                      vtotInv = one/max(eps,sqrt(vvx*vvx + vvy*vvy + vvz*vvz))
                      BCData(mm)%flowXdirInlet(i,j) = vvx*vtotInv
                      BCData(mm)%flowYdirInlet(i,j) = vvy*vtotInv
                      BCData(mm)%flowZdirInlet(i,j) = vvz*vtotInv

                      ! Simply extrapolate the turbulence variables.

                      do l=nt1,nt2
                         BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                      enddo
                   enddo
                enddo

             end select

          enddo bocos
       enddo domains
    enddo spectralLoop

  end subroutine initBCDataDomainInterfaces
