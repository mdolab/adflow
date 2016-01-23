subroutine irregularCellCorrection(level, sps)

  use blockPointers

  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps
  
  ! Working
  integer(kind=intType) :: i, j, k, nn

  do nn=1, nDom
     call setPointers(nn, level, sps)
     
     do k=2, kl
        do j=2, jl
           do i=2, il
              if (fringes(i, j, k)%isDonor .and. &
                   fringes(i, j, k)%donorProc /= -1) then 
                 ! Clear this fringe
                 call emptyFringe(fringes(i, j, k))
              end if
           end do
        end do
     end do
  end do
  
  
end subroutine irregularCellCorrection
