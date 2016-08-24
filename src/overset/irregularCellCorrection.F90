subroutine irregularCellCorrection(level, sps)

  use blockPointers
  use utils, only : setPointers
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps
  
  ! Working
  integer(kind=intType) :: i, j, k, nn
  logical :: isDonor
  do nn=1, nDom
     call setPointers(nn, level, sps)
     
     do k=2, kl
        do j=2, jl
           do i=2, il
              if (isDonor(fringes(i, j, k)%status) .and. &
                   fringes(i, j, k)%donorProc /= -1) then 
                 ! Clear this fringe
                 call emptyFringe(fringes(i, j, k))
              end if
           end do
        end do
     end do
  end do
  
  
end subroutine irregularCellCorrection
