module petscCompat
! Include PETSc version macros
#include <petsc/finclude/petsc.h>
    ! Module for PETSc compatibility functions
    ! This module provides compatibility for various subroutines in PETSc
    ! This is because of the changes in the API between different versions of PETSc
    ! This module is used to ensure that the code works with different versions of PETSc
    use petscvec
    use petscmat
    use petscksp
    use petscdm
    use petscdmda
    use precision
    implicit none
    save

contains

    subroutine VecGetArrayCompat(v, array, ierr)
        ! PETSc compatibility routine to get array from vector
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecGetArrayF90(v, array, ierr)
#else
        call VecGetArray(v, array, ierr)
#endif

    end subroutine VecGetArrayCompat

    subroutine VecRestoreArrayCompat(v, array, ierr)
        ! PETSc compatibility routine to restore array
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecRestoreArrayF90(v, array, ierr)
#else
        call VecRestoreArray(v, array, ierr)
#endif

    end subroutine VecRestoreArrayCompat

    subroutine VecGetArrayReadCompat(v, array, ierr)
        ! PETSc compatibility routine to get read-only array from vector
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecGetArrayReadF90(v, array, ierr)
#else
        call VecGetArrayRead(v, array, ierr)
#endif

    end subroutine VecGetArrayReadCompat

    subroutine VecRestoreArrayReadCompat(v, array, ierr)
        ! PETSc compatibility routine to restore array
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecRestoreArrayReadF90(v, array, ierr)
#else
        call VecRestoreArrayRead(v, array, ierr)
#endif

    end subroutine VecRestoreArrayReadCompat

    subroutine DMDAVecGetArrayCompat(da, v, array, ierr)
        ! PETSc compatibility routine to get array from vector
        implicit none
        DM :: da
        Vec :: v
        real(kind=realType), dimension(:, :, :), pointer :: array
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call DMDAVecGetArrayF90(da, v, array, ierr)
#else
        call DMDAVecGetArray(da, v, array, ierr)
#endif

    end subroutine DMDAVecGetArrayCompat

    subroutine DMDAVecRestoreArrayCompat(da, v, array, ierr)
        ! PETSc compatibility routine to restore array
        implicit none
        DM :: da
        Vec :: v
        real(kind=realType), dimension(:, :, :), pointer :: array
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call DMDAVecRestoreArrayF90(da, v, array, ierr)
#else
        call DMDAVecRestoreArray(da, v, array, ierr)
#endif

    end subroutine DMDAVecRestoreArrayCompat

    subroutine VecCreateMPIWithArrayCompat(comm, bs, n, nGlobal, array, vv, ierr)
        implicit none
        integer, intent(in) :: comm
        integer(kind=intType), intent(in) :: bs, n, nGlobal
        real(kind=realType), intent(in), optional, dimension(:) :: array
        Vec, intent(out) :: vv
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,22,0)
        if (present(array)) then
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, array, vv, ierr)
        else
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, PETSC_NULL_SCALAR_ARRAY, vv, ierr)
        end if
#else
        if (present(array)) then
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, array, vv, ierr)
        else
            call VecCreateMPIWithArray(comm, bs, n, nGlobal, PETSC_NULL_SCALAR, vv, ierr)
        end if

#endif

    end subroutine VecCreateMPIWithArrayCompat

    subroutine PCASMGetSubKSPCompat(pc, nlocal, first, subKSP, ierr)
        use utils, only: EChk
        implicit none
        PC, intent(in) :: pc
        integer(kind=intType), intent(out) :: nlocal, first
        KSP, pointer, intent(out) :: subKSP(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,23,0)
        ! PETSc >= 3.23, just call directly
        call PCASMGetSubKSP(pc, nlocal, first, subKSP, ierr)
#else
        ! Older PETSc, first query to get the size, then allocate, then call again
        call PCASMGetSubKSP(pc, nlocal, first, PETSC_NULL_KSP, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (nlocal > 0) then
            allocate (subKSP(nlocal))
            call PCASMGetSubKSP(pc, nlocal, first, subKSP, ierr)
        else
            nullify (subKSP)
        end if
#endif

    end subroutine PCASMGetSubKSPCompat

    subroutine PCASMRestoreSubKSPCompat(pc, nlocal, first, subKSP, ierr)

        implicit none
        PC, intent(in) :: pc
        integer(kind=intType), intent(in) :: nlocal, first
        KSP, pointer, intent(inout) :: subKSP(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,23,0)
#if PETSC_VERSION_GE(3,23,0) && PETSC_VERSION_LT(3,24,0)
        ! PETSc 3.23 specific handling
        ! Note on PETSc 3.23 : It appears that PETSc does not have PCASMRestoreSubKSP, this is added in 3.24.
        ierr = 0
#else
        ! PETSc >= 3.24, just call directly
        call PCASMRestoreSubKSP(pc, nlocal, first, subKSP, ierr)
#endif
#else
        ! Older PETSc, check if pointer is associated, deallocate
        if (associated(subKSP)) deallocate (subKSP)
        ierr = 0
#endif

    end subroutine PCASMRestoreSubKSPCompat

    subroutine KSPSetResidualHistoryCompat(ksp, a, na, reset, ierr)
        implicit none
        KSP, intent(in) :: ksp
        real(kind=alwaysRealType), dimension(:), intent(in) :: a
        PetscCount, intent(in) :: na
        logical, intent(in) :: reset
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,22,0)
        call KSPSetResidualHistory(ksp, a, int(na, kind=intType), reset, ierr)
#else
        call KSPSetResidualHistory(ksp, a, na, reset, ierr)
#endif

    end subroutine KSPSetResidualHistoryCompat


end module petscCompat
