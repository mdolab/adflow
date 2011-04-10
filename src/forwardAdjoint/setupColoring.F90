! This routine setups "coloring" stencils for various sizes

subroutine setup_PC_coloring(nn,nColor)
  use blockPointers
  implicit none

  ! We assume that setPointers has already been called for this block, nn

  integer(kind=intType) :: i,j,k,nn
  integer(kind=intType) :: nColor

  call setPointers(nn,1,1) ! Just to get the correct sizes

  do k=0,kb
     do j=0,jb
        do i=0,ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)

           flowDomsd(1)%color(i,j,k) = &
                mod( (i+1) + 5*(j+1) + 4*(k+1),7) + 1

        end do
     end do
  end do
  
  nColor = 7

end subroutine setup_PC_coloring

subroutine setup_dRdw_euler_coloring(nn,nColor)
  use blockPointers
  implicit none

  ! We assume that setPointers has already been called for this block, nn


  integer(kind=intType) :: i,j,k,nn
  integer(kind=intType) :: nColor

  call setPointers(nn,1,1) ! Just to get the correct sizes
  do k=0,kb
     do j=0,jb
        do i=0,ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)


           flowDomsd(1)%color(i,j,k) = &
                mod( i + 14*j + 4*k ,17) + 1

        end do
     end do
  end do
  
  nColor = 17

end subroutine setup_dRdw_euler_coloring

subroutine setup_dRdw_visc_coloring(nn,nColor)
  use blockPointers
  implicit none

  ! We assume that setPointers has already been called for this block, nn


  integer(kind=intType) :: i,j,k,nn
  integer(kind=intType) :: nColor

  call setPointers(nn,1,1) ! Just to get the correct sizes
  do k=0,kb
     do j=0,jb
        do i=0,ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)

           flowDomsd(1)%color(i,j,k) = 0
           print * ,'This coloring is not implemented yet'
           stop

        end do
     end do
  end do
  
  nColor = 13

end subroutine setup_dRdw_visc_coloring

! -------------------------------------------------------------
!                   Debugging Color Colorings
! -------------------------------------------------------------

subroutine setup_3x3x3_coloring(nn,nColor)
  use blockPointers
  implicit none

  ! This is a dense 3x3x3 cube for debugging only

  integer(kind=intType) :: i,j,k,nn
  integer(kind=intType) :: nColor,modi,modj,modk

  call setPointers(nn,1,1) ! Just to get the correct sizes

  do k=0,kb
     do j=0,jb
        do i=0,ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           modi = mod(i,3)
           modj = mod(j,3)
           modk = mod(k,3)

           flowDomsd(1)%color(i,j,k) = modi + 3*modj + 9*modk + 1

        end do
     end do
  end do
  
  nColor = 27
end subroutine setup_3x3x3_coloring

subroutine setup_5x5x5_coloring(nn,nColor)
  use blockPointers
  implicit none

  ! This is a dense 5x5x5 cube for debugging only

  integer(kind=intType) :: i,j,k,nn
  integer(kind=intType) :: nColor,modi,modj,modk

  call setPointers(nn,1,1) ! Just to get the correct sizes

  do k=0,kb
     do j=0,jb
        do i=0,ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)
           modi = mod(i,5)
           modj = mod(j,5)
           modk = mod(k,5)

           flowDomsd(1)%color(i,j,k) = modi + 5*modj + 25*modk + 1

        end do
     end do
  end do
  
  nColor = 125
end subroutine setup_5x5x5_coloring

subroutine setup_BF_coloring(nn,nColor)
  use blockPointers
  implicit none

  ! This is a REALLY brute force coloring for debugging

  integer(kind=intType) :: i,j,k,nn
  integer(kind=intType) :: nColor,modi,modj,modk

  call setPointers(nn,1,1) ! Just to get the correct sizes

  do k=0,kb
     do j=0,jb
        do i=0,ib
           ! Add the extra one for 1-based numbering (as opposed to zero-based)

           flowDomsd(1)%color(i,j,k) = i + j*(ib+1) + k*((ib+1)*(jb+1)) + 1

        end do
     end do
  end do
  
  nColor = (ib+1)*(jb+1)*(kb+1)
end subroutine setup_BF_coloring
