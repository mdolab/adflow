
subroutine getdRdwTPsi(ndof,dRdwTPsi)
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

end subroutine getdRdwTPsi

subroutine getdRdwTPsi2(ndof,dRdwTPsi)
  !
  use communication
  use ADjointPETSc, only: dRdwT,psi,wVec
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(out)  :: dRdwtPsi(ndof)

  integer(kind=intType) :: ierr,i

  real(kind=realType), dimension(2) :: time
  call cpu_time(time(1))
  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,dRdwTPsi,wVec,ierr)
  do i=1,1000


     call MatMult(dRdwT,psi,wVec,ierr)


  end do
  call VecDestroy(wVec,ierr)
  call cpu_time(time(2))
  print *,'Time for 10 matvecs:',time(2)-time(1)


end subroutine getdRdwTPsi2
