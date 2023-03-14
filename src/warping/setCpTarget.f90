subroutine setCpTargets(cptarget, npts, sps_in)
    use constants
    use blockPointers
    use flowVarRefState
    use inputTimeSpectral
    use communication
    use inputPhysics
    use utils, only: setPointers
    implicit none
    !
    !      Arguments.
    !
    integer(kind=intType), intent(in) :: npts, sps_in
    real(kind=realType), intent(in) :: cptarget(npts)
    !
    !      Local variables.
    !
    integer(kind=intType) :: mm, nn, i, j, ii, sps
    !integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

    sps = sps_in

    ii = 0
    domains: do nn = 1, nDom
        call setPointers(nn, 1_intType, sps)

        ! Loop over the number of viscous boundary subfaces of this block.
        ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
        ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
        bocos: do mm = 1, nBocos
            if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
                do j = BCData(mm)%jnBeg, BCData(mm)%jnEnd
                    do i = BCData(mm)%inBeg, BCData(mm)%inEnd
                        ii = ii + 1
                        BCData(mm)%CpTarget(i, j) = cptarget(ii)
                    end do
                end do
            end if

        end do bocos
    end do domains

end subroutine setCpTargets

subroutine getCpTargets(cptarget, npts, sps)

    use constants
    use blockPointers, only: nDom, nBocos, BCData, BCType
!  use flowVarRefState, only : TRef
    use utils, only: setPointers
    implicit none

    ! Input Variables
    integer(kind=intType), intent(in) :: npts, sps
    real(kind=realType), intent(out) :: cptarget(npts)

    ! Local Variables
    integer(kind=intType) :: mm, nn, i, j, ii
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

    ii = 0
    domains: do nn = 1, nDom
        call setPointers(nn, 1_intType, sps)
        ! Loop over the number of viscous boundary subfaces of this block.
        bocos: do mm = 1, nBocos
            wall: if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic .or. &
                      BCType(mm) == NSWallIsoThermal) then
                do j = BCData(mm)%jnBeg, BCData(mm)%jnEnd
                    do i = BCData(mm)%inBeg, BCData(mm)%inEnd
                        ii = ii + 1
                        CpTarget(ii) = BCData(mm)%CpTarget(i, j)
                    end do
                end do
            end if wall
        end do bocos
    end do domains
end subroutine getCpTargets
