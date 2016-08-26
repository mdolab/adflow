!
subroutine binSearchPocketEdgeType(arr, search, nn, searchInd)

!  binSearchPocketEdgeType does binary searche for            
!  pocketEdgeType 'search' edge and returns index 'searchInd' 
!  where 'search' lies in arr.                                

  use overset
  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: nn
  integer(kind=intType), intent(out) :: searchInd

  type(pocketEdge), intent(in) :: search
  type(pocketEdge), dimension(*), intent(in) :: arr

  ! Local variables
  integer(kind=intType) :: first, last, middle

  first = 1
  last = nn

  middle = (first+last)/2

  do while (first <= last)
     if (arr(middle) < search) then
        first = middle + 1
     else if (arr(middle) == search) then
        searchInd = middle
        exit
     else
        last = middle -1
     end if

     middle = (first+last)/2
  end do !while

  if (first > last) then
     print*, ' binSearchPocketEdgeType fails for Edge with nodes ',&
             search%n1, search%n2
     STOP 
  end if
end subroutine binSearchPocketEdgeType
