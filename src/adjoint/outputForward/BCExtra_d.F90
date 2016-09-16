module BCExtra_d

  use constants
  !use BCRoutines_d
  implicit none
  save

contains

  subroutine applyAllBC_block_d(secondHalo)

    ! Apply BC's for a single block

    use constants
    use BCRoutines_d
    use blockPointers, only : bcType, nBocos, nViscBocos, bcdatad
    use flowVarRefState, only : kPresent
    use inputDiscretization, only : precond
    use iteration, only : currentLevel, groundLevel
    use BCPointers_d
    use utils, only : getCorrectForK, setBCPointers_d
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo

    ! Local variables.
    logical :: correctForK
    integer(kind=intType) :: nn
    !
    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.
    correctForK = getCorrectForK()

    ! Apply all the boundary conditions. The order is important!  Only
    ! some of them have been AD'ed

    ! ------------------------------------
    !  Symmetry Boundary Condition 
    ! ------------------------------------
    do nn=1, nBocos
       if (bcType(nn) == symm) then 
          call setBCPointers_d(nn, .False.)
          call bcSymm1stHalo_d(nn)
       end if
    end do

    if (secondHalo) then 
       do nn=1, nBocos
          if (bcType(nn) == symm) then 
             call setBCPointers_d(nn, .False.)
             call bcSymm2ndHalo_d(nn)
          end if
       end do
    end if

    ! ! ------------------------------------
    ! !  Symmetry Polar Boundary Condition
    ! ! ------------------------------------
    ! !$AD II-LOOP
    ! do nn=1, nBocos
    !    if (BCType(nn) == symmPolar) then
    !       call setBCPointers_d(nn, .True.)
    !       call bcSymmPolar1stHalo_d(nn)
    !    end if
    ! end do

    ! if (secondHalo) then
    !    !$AD II-LOOP
    !    do nn=1, nBocos
    !       if (BCType(nn) == symmPolar) then
    !          call setBCPointers_d(nn, .True.)
    !          call bcSymmPolar2ndHalo_d(nn)
    !       end if
    !    end do
    ! end if

    ! ------------------------------------
    !  adibatic Wall Boundary Condition 
    ! ------------------------------------
    do nn=1, nViscBocos
       if (bcType(nn) == NSWallAdiabatic) then 
          call setBCPointers_d(nn, .False.)
          call bcNSWallAdiabatic_d(nn, secondHalo, correctForK)
       end if
    end do

    ! ------------------------------------
    !  Isothermal Wall Boundary Condition 
    ! ------------------------------------
    do nn=1, nViscBocos
       if (bcType(nn) == NSWallIsoThermal) then 
          call setBCPointers_d(nn, .False.)
          call bcNSWallIsothermal_d(nn, secondHalo, correctForK)
       end if
    end do

    ! ------------------------------------
    !  Farfield Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == farField) then
          call setBCPointers_d(nn, .False.)
          call bcFarField_d(nn, secondHalo, correctForK)
       end if
    end do

    ! ------------------------------------
    !  Euler Wall Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == EulerWall) then
          call setBCPointers_d(nn, .True.)
          call bcEulerWall_d(nn, secondHalo, correctForK)
       end if
    end do

  end subroutine applyAllBC_block_d

end module BCExtra_d
