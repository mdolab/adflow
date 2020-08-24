module turbAPI

contains
  subroutine turbSolveDDADI
    !
    !       turbSolveDDADI solves the turbulent transport equations
    !       separately, i.e. the mean flow variables are kept constant
    !       and the turbulent variables are updated.
    !
    use constants
    use blockPointers, only : nDom
    use flowVarRefState
    use inputDiscretization
    use inputIteration
    use inputPhysics
    use iteration
    use turbMod
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use sa
    use kw
    use kt
    use SST
    use vf
    use haloExchange, only : whalo2
    use utils, only : setPointers
    use turbUtils, only : unsteadyTurbSpectral
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: iter, sps, nn

    ! Loop over the number of iterations for the turbulence.

    do iter=1,nSubIterTurb

       ! Compute the quantities for certain turbulence models that
       ! need to be communicated between blocks.

       if (turbModel == menterSST) then
          call f1SST
       end if

       ! Compute the time derivative for the time spectral mode.
       select case(turbModel)
       case(spalartAllmaras)
          call unsteadyTurbSpectral(itu1,itu1)
       case (komegaWilcox, komegaModified, menterSST, ktau)
          call unsteadyTurbSpectral(itu1,itu2)
       case (v2f)
          call unsteadyTurbSpectral(itu1,itu3)
       end select

       ! Loop over the number of spectral solutions.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

          ! Loop over the number of blocks.

          domains: do nn=1,nDom

             ! setPointers for this block:
             call setPointers(nn, currentLevel, sps)

             ! Now call the selected turbulence model
             select case (turbModel)

             case (spalartAllmaras)
                call sa_block(.false.)

             case (komegaWilcox, komegaModified)
                call kw_block(.false.)

             case (menterSST)
                call SST_block(.false.)

             case (ktau)
                call kt_block(.false.)

             case (v2f)
                call vf_block(.false.)

             end select

          end do domains
       end do spectralLoop

       ! Exchange the halo data. As it is guaranteed that we are on the
       ! finest mesh, exchange both layers of halo's.

       call whalo2(groundLevel, nt1, nt2, .false., .false., .true.)

    enddo

  end subroutine turbSolveDDADI

  subroutine turbResidual
    !
    !       turbResidual computes the residual of the residual of the
    !       turbulent transport equations on the current multigrid level.
    !
    use constants
    use blockPointers, only : nDom
    use inputDiscretization
    use inputPhysics
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration
    use turbMod
    use sa
    use kt
    use kw
    use SST
    use vf
    use utils, only : setPointers
    implicit none

    integer(kind=intType) :: nn, sps

    ! Compute the quantities for certain turbulence models that
    ! need to be communicated between blocks.

    if (turbModel == menterSST) then
       call f1SST
    end if

    ! Loop over the number of spectral solutions.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Loop over the number of blocks.

       domains: do nn=1,nDom

          ! setPointers for this block:
          call setPointers(nn, currentLevel, sps)

          ! Now call the selected turbulence model
          select case (turbModel)

          case (spalartAllmaras)
             call sa_block(.True.)

          case (komegaWilcox, komegaModified)
             call kw_block(.True.)

          case (menterSST)
             call SST_block(.True.)

          case (ktau)
             call kt_block(.True.)

          case (v2f)
             call vf_block(.True.)

          end select

       end do domains
    end do spectralLoop

  end subroutine turbResidual


end module turbAPI
