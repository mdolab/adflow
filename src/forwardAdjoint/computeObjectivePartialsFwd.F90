Subroutine computeObjectivePartialsFwd(costFunction)

  !     ******************************************************************
  !     *                                                                *
  !     * Assemble the objective partials from the data precomputed      *
  !     * from the forward mode computations. 
  !     *                                                                *
  !     * costFunction: The index of the cost function to use.           *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointVars
  use ADjointPETSc
  use blockPointers   
  use communication  
  use costFunctions
  use inputTimeSpectral
  use inputPhysics
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: costFunction

  ! Working Variables
  integer(kind=intType) :: i, ierr, sps
  real(kind=realtype) :: val
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, forceb
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: moment, momentb
  real(kind=realType) :: alpha, alphab, beta, betab
  real(kind=realType) :: objValue, objValueb
  integer(kind=intType) :: liftIndex, idim
  ! Compute the requiquired sensitivity of the objective with respect
  ! to the forces, moments and extra variables.
#ifndef USE_COMPLEX       
  do sps=1, nTimeIntervalsSpectral
     call getSolution(sps)
     force(1, sps) = functionValue(costFuncForceX)
     force(2, sps) = functionValue(costFuncForceY)
     force(3, sps) = functionValue(costFuncForceZ)
     moment(1, sps) = functionValue(costFuncMomX)
     moment(2, sps) = functionValue(costFuncMomY)
     moment(3, sps) = functionValue(costFuncMomZ)
  end do

  objValueb = one
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  call getCostFunction_b(costFunction, force, forceb, moment, momentb, &
       alpha, alphab, beta, betab, liftIndex, objValue, objValueb)

  !******************************************! 
  !               dIdw                       ! 
  !******************************************! 

  ! Zero Entries and multiply through by reverse-mode derivatives
  call VecZeroEntries(dJdw, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  do sps=1, nTimeIntervalsSpectral
     do i=1, 3
        call VecAXPY(dJdw, forceb(i, sps), FMw(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAXPY(dJdw, momentb(i, sps), FMw(i+3, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  ! Assemble dJdw
  call VecAssemblyBegin(dJdw, ierr)
  call EChk(ierr, __FILE__, __LINE__)
 call VecAssemblyEnd(dJdw, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !******************************************!
  !          dIdx CALCULATIONS               !
  !******************************************!

  ! Zero Entries and multiply by reverse-mode derivatives
  call VecZeroEntries(dJdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do sps=1, nTimeIntervalsSpectral
     do i = 1, 3
        call VecAXPY(dJdx, forceb(i, sps), FMx(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAXPY(dJdx, momentb(i, sps), FMx(i+3, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  ! Assemble dJdx
  call VecAssemblyBegin(dJdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(dJdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !******************************************!
  !          dIda CALCULATIONS               !
  !******************************************!
  if (nDesignExtra > 0) then
     dIda = zero

     if (myid == 0) then
        if (nDesignAoA >= 0) then 
           dIda(nDesignAoA + 1) = dIda(nDesignAoA+1) + alphab
        end if
        if (nDesignSSA >= 0) then 
           dIda(nDesignSSA + 1) = dIda(nDesignSSA+1) + betab
        end if
        if (nDesignMach > 0) then 
           dIda(nDesignMach + 1) = dIda(nDesignMach+1) + machb + machcoefb
        end if
        
        if (nDesignMachGrid > 0) then 
           dIda(nDesignMachGrid + 1) = dIda(nDesignMachGrid+1) + machgridb + machcoefb
        end if

        if (nDesignLengthRef > 0) then
           dIda(nDesignLengthRef + 1) = dIda(nDesignLengthRef+1) + lengthrefb
        end if

        ! Explict dependence on pointRef....Only on one proc since
        ! they are the same across all procs, and the sum would nProc
        ! too large if they were on all procs. 
        if (nDesignPointRefX > 0) then
           dIda(nDesignPointRefX + 1) = dIda(nDesignPointRefX + 1) + pointrefb(1)
        end if

        if (nDesignPointRefY > 0) then
           dIda(nDesignPointRefY + 1) = dIda(nDesignPointRefY + 1) + pointrefb(2)
        end if

        if (nDesignPointRefZ > 0) then
           dIda(nDesignPointRefZ + 1) = dIda(nDesignPointRefZ + 1) + pointrefb(3)
        end if

     end if

     
     ! These three are a little different; The derivative wrt to the
     ! forces and moments for each spectral instance are computed when
     ! the extra residual matrix is computed. This is the only
     ! sensivitiy for pointRef. We also know the objective derivative
     ! wrt forces and moment so we chain-rule them together. Note that
     ! there is no dependence of 'force' on pointRef so it is not
     ! included here. Also these derivatives DO need to be summed over
     ! all procs
     
     ! add missing dependence of mach
     if (nDesignMach > 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignMach+1) = dIda(nDesignMach+1) + &
                   dFMdExtra(idim, nDesignmach+1, sps)*forceb(idim, sps) + &
                   dFMdExtra(idim+3, nDesignmach+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignPointRefX > 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignPointRefX+1) = dIda(nDesignPointRefX+1) + &
                   dFMdExtra(3+idim, nDesignPointRefX+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignPointRefY > 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignPointRefY+1) = dIda(nDesignPointRefY+1) + &
                   dFMdExtra(3+idim, nDesignPointRefY+1, sps)*momentb(idim, sps)
           end do
        end do
     end if
     if (nDesignPointRefZ > 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignPointRefZ+1) = dIda(nDesignPointRefZ+1) + &
                   dFMdExtra(3+idim, nDesignPointRefZ+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignLengthRef > 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignLengthRef+1) = dIda(nDesignLengthRef+1) + &
                   dFMdExtra(idim, nDesignLengthRef+1, sps)*forceb(idim, sps)
              dIda(nDesignLengthRef+1) = dIda(nDesignLengthRef+1) + &
                   dFMdExtra(idim+3, nDesignLengthRef+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

  end if
#else
  print *,'Cost Function routines are not complexified'
  stop
#endif

end subroutine computeObjectivePartialsFwd

! Add two functions to return dIdw and dIdx. dIda is available
! directly in python in dIda in the adjointVars module. 

subroutine getdIdw(output, nstate)

#ifndef USE_NO_PETSC	
  ! #define PETSC_AVOID_MPIF_
  ! #include "finclude/petscdef.h"

  use ADjointPETSc, only : dJdw
  use constants

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in):: nstate
  real(kind=realType), dimension(nstate), intent(inout) :: output
  real(kind=realType), pointer :: dJdw_pointer(:)

  ! Local Variables
  integer(kind=intType) :: i, ierr

  ! Copy out adjoint vector:
  call VecGetArrayF90(dJdw, dJdw_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do a straight copy:
  do i=1, nstate
     output(i) = dJdw_pointer(i)
  end do

  call VecRestoreArrayF90(dJdw, dJdw_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

end subroutine getdIdw

subroutine getdIdx(output, ndof)

#ifndef USE_NO_PETSC	
  ! #define PETSC_AVOID_MPIF_
  ! #include "finclude/petscdef.h"

  use ADjointPETSc, only : dJdx
  use constants

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in):: ndof
  real(kind=realType), dimension(ndof), intent(inout) :: output
  real(kind=realType), pointer :: dJdx_pointer(:)

  ! Local Variables
  integer(kind=intType) :: i, ierr

  ! Copy out adjoint vector:
  call VecGetArrayF90(dJdx, dJdx_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do a straight copy:
  do i=1, ndof
     output(i) = dJdx_pointer(i)
  end do
  
  call VecRestoreArrayF90(dJdx, dJdx_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

end subroutine getdIdx

subroutine zeroObjPartials(stateSetup, spatialSetup)
#ifndef USE_NO_PETSC
  use precision 
  use ADjointVars ! include costFunctions
  use ADjointPETSc

  integer(kind=intType) :: ierr
  logical, intent(in) :: stateSetup, spatialSetup

  if (stateSetup) then
     call VecZeroEntries(dJdw, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  if (spatialSetup) then
     call VecZeroEntries(dJdx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if
#endif
end subroutine zeroObjPartials
