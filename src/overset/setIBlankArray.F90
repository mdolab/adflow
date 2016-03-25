subroutine setIblankArray(level, sps)

  use blockPointers
  use communication
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps
  
  ! Working
  integer(kind=intType) :: i, j, k, nn

  do nn=1, nDom
     call setPointers(nn, level, sps)
     iBlank = 1
     do k=2, kl
        do j=2, jl
           do i=2, il
              if (fringes(i, j, k)%donorProc /= -1) then
                 iblank(i, j, k) = -1
              end if
              if (fringes(i, j, k)%isHole) then 
                 iBlank(i, j, k) = 0
              end if


              if (fringes(i, j, k)%isFlooded) then 
                 iBlank(i, j, k) = -2
              end if

              if (fringes(i, j, k)%isFloodSeed) then 
                 iBlank(i, j, k) = -3
              end if

           end do
        end do
     end do
  end do

  ! Update the iblank info. 
  
  call exchangeIblanks(level, sps, commPatternCell_2nd, internalCell_2nd)
end subroutine setIblankArray
