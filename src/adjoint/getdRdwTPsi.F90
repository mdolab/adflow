subroutine getdRdwTPsi(ndof,dRdwTPsi)
#ifndef USE_NO_PETSC

  !
  use communication
  use ADjointPETSc, only: dRdwT,psi,wVec
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(out)  :: dRdwtPsi(ndof)

  integer(kind=intType) :: ierr

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,dRdwTPsi,wVec,ierr)
  call MatMult(dRdwT,psi,wVec,ierr)
  call VecDestroy(wVec,ierr)

#endif
end subroutine getdRdwTPsi
