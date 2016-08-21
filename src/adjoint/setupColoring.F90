! This routine setups "coloring" stencils for various sizes

subroutine setup_PC_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ib, jb, kb
  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k

  call setPointers(nn, level, 1) ! Just to get the correct sizes

  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           flowDoms(nn, 1, 1)%color(i, j, k) = &
                mod(i + 5*j + 4*k, 7) + 1
        end do
     end do
  end do
  
  nColor = 7

end subroutine setup_PC_coloring

subroutine setup_dRdw_euler_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ib, jb, kb

  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k

  call setPointers(nn, level, 1) ! Just to get the correct sizes
  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           flowDoms(nn, 1, 1)%color(i, j, k) = &
                mod( i + 3*j + 4*k , 13) + 1
        end do
     end do
  end do
  
  nColor = 13

end subroutine setup_dRdw_euler_coloring

subroutine setup_dRdw_visc_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ib, jb, kb

  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k

  call setPointers(nn, level, 1) ! Just to get the correct sizes
  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           flowDoms(nn, 1, 1)%color(i,j,k) = &
                mod( i + 19*j + 11*k ,35) + 1
        end do
     end do
  end do
  
  nColor = 35

end subroutine setup_dRdw_visc_coloring


! -------------------------------------------------------------
!                   Debugging Color Colorings
! -------------------------------------------------------------

subroutine setup_3x3x3_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ib, jb, kb

  implicit none

  ! This is a dense 3x3x3 cube for debugging only
  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k, modi, modj, modk

  call setPointers(nn, level, 1) ! Just to get the correct sizes

  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           modi = mod(i, 3)
           modj = mod(j, 3)
           modk = mod(k, 3)

           flowDoms(nn, 1, 1)%color(i, j, k) = modi + 3*modj + 9*modk + 1

        end do
     end do
  end do
  
  nColor = 27
end subroutine setup_3x3x3_coloring

subroutine setup_4x4x4_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ib, jb, kb

  implicit none

  ! This is a dense 4x4x4 cube for debugging drdx only
  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k, modi, modj, modk

  call setPointers(nn, level, 1) ! Just to get the correct sizes

  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           modi = mod(i, 4)
           modj = mod(j, 4)
           modk = mod(k, 4)

           flowDoms(nn, 1, 1)%color(i, j, k) = modi + 4*modj + 16*modk + 1

        end do
     end do
  end do
  
  nColor = 64
end subroutine setup_4x4x4_coloring

subroutine setup_5x5x5_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ib, jb, kb

  implicit none

  ! This is a dense 5x5x5 cube for debugging only
  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k, modi, modj, modk

  call setPointers(nn, level, 1) ! Just to get the correct sizes

  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           modi = mod(i, 5)
           modj = mod(j, 5)
           modk = mod(k, 5)

           flowDoms(nn, 1, 1)%color(i, j, k) = modi + 5*modj + 25*modk + 1

        end do
     end do
  end do
  
  nColor = 125
end subroutine setup_5x5x5_coloring

subroutine setup_BF_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ib, jb, kb

  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k

  ! This is a REALLY brute force coloring for debugging

  call setPointers(nn, level, 1) ! Just to get the correct sizes

  do k=0, kb
     do j=0, jb
        do i=0, ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)

           flowDoms(nn, 1, 1)%color(i, j, k) = i + j*(ib+1) + k*((ib+1)*(jb+1)) + 1
        end do
     end do
  end do
  
  nColor = (ib+1)*(jb+1)*(kb+1)
end subroutine setup_BF_coloring

subroutine setup_BF_node_coloring(nn, level, nColor)

  use constants
  use blockPointers, only : flowDoms, ie, je, ke

  implicit none

  ! This is a REALLY brute force coloring for debugging
  ! Input parameters
  integer(kind=intType), intent(in) :: nn, level

  ! Output parameters
  integer(kind=intTYpe), intent(out) :: nColor

  ! Working 
  integer(kind=intType) :: i, j, k

  call setPointers(nn, level, 1) ! Just to get the correct sizes

  do k=0, ke
     do j=0, je
        do i=0, ie
           ! Add the extra one for 1-based numbering (as opposed to zero-based)

           flowDoms(nn, 1, 1)%color(i, j, k) = i + j*(ie+1) + k*((ie+1)*(je+1)) + 1

        end do
     end do
  end do
  
  nColor = (ie+1)*(je+1)*(ke+1)

end subroutine setup_BF_Node_coloring
