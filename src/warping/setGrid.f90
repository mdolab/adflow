!
! ***********************************
! *  File: setGrid.f90
! *  Author: Gaetan Kenway
! *  Started: 05-29-2010
! *  Modified: 05-29-2010
! ***********************************

subroutine setGrid(grid,ndof)

  ! The purpose of this routine is to set the grid dof as returned by
  ! the external warping. 

  use precision 
  use blockPointers
  use communication
  use warpingPetsc
  implicit none

  integer(kind=intType),intent(in) :: ndof
  real(kind=realType) ,intent(in) :: grid(ndof)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,counter,idim

  ! This is very straight forward...loop over all domains and set all elements
  counter = 0
  do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     do k=1,kl
        do j=1,jl
           do i=1,il
              do idim=1,3
                 counter = counter + 1
                 X(i,j,k,idim) = grid(counter)
              end do
           end do
        end do
     end do
  end do
  call xhalo(1_intType)
end subroutine setGrid
