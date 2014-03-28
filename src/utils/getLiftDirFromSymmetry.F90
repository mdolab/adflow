
subroutine getLiftDirFromSymmetry(liftDir)

! The purpose of this function is to determine what coordinate
! direction the mirror plane is in. It does NOT handle multiple mirror
! planes. It is used just to determine what the lift direction is. 

  use constants
  use bcTypes
  use blockPointers
  use communication
  implicit none

  ! Output
  integer(kind=intType), intent(out) :: liftDir
  integer(kind=intType), dimension(3) :: sym_local, sym

  ! Working
  integer(kind=intType) :: nn, i_index(1), mm, ierr, jerr
  real(kind=realType), dimension(:, :, :), pointer :: xx
  real(kind=realType) :: cp(3), v1(3), v2(3)
  ! Loop over each block and each subFace

  sym_local = 0_intType
  sym       = 0_intType
  liftDir   = 0_intType
  do nn=1, nDom
     call setPointers(nn, 1, 1)
     do mm=1,nBocos
        if (bcType(mm) == symm) then

           select case (BCFaceID(mm))
           case (iMin)
              xx => x(1, :, :, :)
           case (iMax)
              xx => x(il, :, :, :)
           case (jMin)
              xx => x(:, 1, :, :)
           case (jMax)
              xx => x(:, jl, :, :)
           case (kMin)
              xx => x(:, :, 1, :)
           case (kMax)
              xx => x(:, :, kl, :)
           end select

           ! Take the cross product
           v1(:) = xx(bcData(mm)%inEnd, bcData(mm)%jnEnd, :) - &
                   xx(bcData(mm)%inBeg, bcData(mm)%jnBeg, :)
           v2(:) = xx(bcData(mm)%inBeg, bcData(mm)%jnEnd, :) - &
                   xx(bcData(mm)%inEnd, bcData(mm)%jnBeg, :)

           ! Cross Product
           cp(1) = (v1(2)*v2(3) - v1(3)*v2(2))
           cp(2) = (v1(3)*v2(1) - v1(1)*v2(3))
           cp(3) = (v1(1)*v2(2) - v1(2)*v2(1))  

           ! Only interesed in abs values
           cp = abs(cp)

           ! Location, ie coordiante direction of dominate direction
           i_index = maxloc(real(cp))

           sym_local(i_index(1)) = 1_intType
        end if
     end do
  end do

  ! Now we have a bunch of sym_locals, mpi_allreduce them and SUM

  call MPI_Allreduce (sym_local, sym, 3, sumb_integer, &
       MPI_SUM, sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now we should make sure that only ONE of the values is
  ! non-zero. If more than one value is zero, it means we have
  ! multiple symmetry planes which we can't support.
  if (sym(1) == 0 .and. sym(2) == 0 .and. sym(3) == 0) then
     ! Pass - no sym, can't determine lift dir:
  else if(sym(1) .ne. 0 .and. sym(2) == 0 .and. sym(3) == 0) then
     ! Pass - x dir can't be symmetry
  else if(sym(1) == 0 .and. sym(2) .ne. 0 .and. sym(3) == 0) then
     liftDir = 3
  else if(sym(1) == 0 .and. sym(2) == 0 .and. sym(3) .ne. 0) then
     liftDir = 2
  else
     ! Multiple orientations...can't do anything
  end if

end subroutine getLiftDirFromSymmetry
