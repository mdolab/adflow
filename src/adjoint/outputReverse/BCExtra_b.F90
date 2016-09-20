module BCExtra_b

  use constants
  implicit none
  save
contains 

  subroutine applyAllBC_block_b(secondHalo)

    ! Apply BC's for a single block
    use constants
    use bcroutines_b ! All of them
    use blockPointers, only : nBocos, nViscBocos, BCType
    use utils, only : setBCPointers_d, getCorrectForK
    implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo

    ! Local variables.
    logical :: correctForK
    integer(kind=intType) ::mm
    !
    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.
    correctForK = getCorrectForK()

    ! ------------------------------------
    !  Farfield Boundary Condition 
    ! ------------------------------------
    do mm=1,nBocos
       if (bcType(mm) == farField) then
          call setBCPointers_d(mm, .False.)
          call bcFarField_b(mm, secondHalo, correctForK)
       end if
    end do

    do mm=1,nBocos
       if (bcType(mm) == eulerWall) then
          call setBCPointers_d(mm, .False.)
          call bcEulerWall_b(mm, secondHalo, correctForK)
       end if
    end do

    ! ------------------------------------
    !  Adibatic Wall Boundary Condition 
    ! ------------------------------------
    DO mm=1,nviscbocos
       IF (bctype(mm) .EQ. nswalladiabatic) THEN
          CALL setBCPointers_d(mm, .false.)
          CALL BCNSWALLADIABATIC_B(mm, secondhalo, correctfork)
       END IF
    END DO
   
    ! ------------------------------------
    !  Symmetry Boundary Condition 
    ! ------------------------------------
    if (secondHalo) then 
       DO mm=1,nbocos
          IF (bctype(mm) .EQ. symm) THEN
             CALL setBCPointers_d(mm, .false.)
             CALL BCSYMM2ndhalo_b(mm)
          END IF
       END DO
    END if

    if (secondHalo) then 
       DO mm=1,nbocos
          IF (bctype(mm) .EQ. symm) THEN
             CALL setBCPointers_d(mm, .false.)
             CALL BCSYMM1sthalo_b(mm)
          END IF
       END DO
    END if
  end subroutine applyAllBC_block_b
end module BCExtra_b
