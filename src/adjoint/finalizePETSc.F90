subroutine finalizePETSc
  !
  !      Finalize PETSc by calling the appropriate routine              
  !      PetscFinalize provided in the PETSc library. This              
  !      automatically calls MPI_Finalize().                            
  !
  use ADjointPETSc, only : PETScIerr
  implicit none

#ifndef USE_NO_PETSC

  call PetscFinalize(PETScIerr)

#endif
end subroutine finalizePETSc
