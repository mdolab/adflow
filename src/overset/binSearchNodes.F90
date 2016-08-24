!
subroutine binSearchNodes(arr, searchNode, nn, searchInd)

! **************************************************************
! * binSearchNodes does binary search for a node 'searchNode'  *
! * in arr(1:nn) and returns index 'searchInd' where           *
! * 'searchNode' lies in arr. searchInd = -1 if not found.     *
! **************************************************************

  use overset
  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: nn, searchNode
  integer(kind=intType), intent(out) :: searchInd
  integer(kind=intType), intent(in) :: arr(nn)

  ! Local variables
  integer(kind=intType) :: first, last, middle

  first = 1
  last = nn

  middle = (first+last)/2

  do while (first <= last)
     if (arr(middle) < searchNode) then
        first = middle + 1
     else if (arr(middle) == searchNode) then
        searchInd = middle
        exit
     else
        last = middle -1
     end if

     middle = (first+last)/2
  end do !while

  if (first > last) then
     searchInd = -1
     print*, ' binSearchNode fails for searchNode ',searchNode
     STOP 
  end if
end subroutine binSearchNodes
