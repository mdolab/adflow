
subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)

  ! Given a list of N points (pts) in three space, with possible
  ! duplicates, (to within tol) return a list of the nUnqiue
  ! uniquePoints of points and a link array of length N, that points
  ! into the unique list

  use precision 
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: N
  real(kind=realType), intent(in), dimension(:, :) :: pts
  real(kind=realType), intent(in) :: tol

  ! Output Parametres
  real(kind=realType), intent(out), dimension(:,:) :: uniquePts
  integer(kind=intType), intent(out), dimension(:) :: link
  integer(kind=intType), intent(out) :: nUnique

  ! Working Parameters
  real(kind=realType), allocatable, dimension(:) :: dists, tmp
  integer(kind=intType), allocatable, dimension(:) :: ind
  integer(kind=intType) :: i, j, nTmp, link_counter, ii
  logical cont, cont2
  integer(kind=intType), parameter :: dupMax = 50
  integer(kind=intType) :: tmpInd(dupMax), subLink(dupMax), nSubUnique
  real(kind=realType) :: subPts(3, dupMax), subUniquePts(3, dupMax)
  interface
     subroutine pointReduceBruteForce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none

       real(kind=realType), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), dimension(:, :) :: uniquePts
       integer(kind=intType), dimension(:) :: link
       integer(kind=intType) :: nUnique
     end subroutine pointReduceBruteForce

  end interface

  if (N == 0) then
     nUnique = 0
     return
  end If

  ! Allocate dists, and the ind pointer
  allocate(dists(N), tmp(N), ind(N))

  ! Compute distances of all points from the first
  do i=1, N
     dists(i) = sqrt((pts(1,i)-pts(1,1))**2 + (pts(2,i)-pts(2,1))**2 + (pts(3, i)-pts(3,1))**2)
     tmp(i) = dists(i)
     ind(i) = i
  end do

  ! Do an argsort on the distances
  call ArgQsort(tmp, N, ind)

  i = 1
  cont = .True.
  link_counter = 0
  nUnique = 0

  do while(cont)
     cont2 = .True.

     j = i
     nTmp = 0
     do while(cont2)
        if (abs(dists(ind(i))-dists(ind(j))) < tol) then
           nTmp = nTmp + 1
           if (nTmp > dupMax) then
              ! Need to make dupMax larger!
              call ECHK(-99, __FILE__, __LINE__)
           end if

           tmpInd(nTmp) = ind(j)
           j = j + 1
           if (j == N+1) then ! Overrun check
              cont2 = .False.
           end if
        else
           cont2 = .False.
        end if
     end do

     ! Copy the points that have the same distance into subPts. Note
     ! these may NOT be the same, since two points can have the same
     ! distnace, but not be co-incident (ie (1,0,0), (-1,0,0))
     do ii=1,nTmp
        subPts(:, ii) = pts(:, tmpInd(ii))
     end do

     ! Brute Force Search them 
     call pointReduceBruteForce(subPts, nTmp, tol, subUniquePts, subLink, nSubUnique)

     do ii=1,nSubUnique
        nUnique = nUnique + 1
        uniquePts(:, nUnique) = subUniquePts(:,ii)
     end do

     do ii=1,nTmp
        link(tmpInd(ii)) = subLink(ii) + link_counter
     end do

     link_counter = link_counter +  maxval(subLink) 

     i = j - 1 + 1
     if (i == N+1) then
        cont = .False.
     end if
  end do
  deallocate(dists, tmp, ind)

end subroutine pointReduce

subroutine pointReduceBruteForce(pts, N, tol, uniquePts, link, nUnique)

  ! Given a list of N points (pts) in three space, with possible
  ! duplicates, (to within tol) return a list of the nUnqiue
  ! uniquePoints of points and a link array of length N, that points
  ! into the unique list

  use precision 
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: N
  real(kind=realType), intent(in), dimension(:, :) :: pts
  real(kind=realType), intent(in) :: tol

  ! Output Parametres
  real(kind=realType), intent(out), dimension(:,:) :: uniquePts
  integer(kind=intType), intent(out), dimension(:) :: link
  integer(kind=intType), intent(out) :: nUnique

  ! Working parameters
  integer(kind=intType) :: i, j
  real(kind=realType) :: dist
  logical :: found_it

  ! First point is *always* unique
  uniquePts(:, 1) = pts(:, 1)
  link(:) = 0
  link(1) = 1
  nUnique = 1

  ! Loop over remainder of points
  do i=2,N
     found_it = .False.
     ! Loop over found unique points (this is where the n^2 comes from)

     uniqueLoop: do j=1,nUnique
        dist = sqrt((pts(1, i)-uniquePts(1, j))**2 + &
             (pts(2, i) - uniquePts(2, j))**2 + &
             (pts(3, i) - uniquePts(3, j))**2)

        ! If pt(i) is within tolerance of the jth unique point, set
        ! link for ith point to the jth unique element and break loop
        if (dist < tol) then
           link(i) = j
           found_it = .True. 
           exit uniqueLoop
        end if
     end do uniqueLoop

     ! If we never found the point, it is unique so add it it to the
     ! unique list list:
     if (.not. found_it) then
        nUnique = nUnique + 1
        uniquePts(:, nUnique) = pts(:, i)
        link(i) = j
     end if
  end do

  ! Agorithm observations: It would appear this algorithm runs in
  ! linear time if all the points are the same, and quadratic time if
  ! all the points are unique, and somewhere inbetween
  ! otherwise. Since this is called from pointReduce, it is most
  ! likely that there are going to be more duplicates than the unique
  ! points. 

end subroutine pointReduceBruteForce

recursive subroutine ArgQSort(A, nA, ind)

  ! Do an ArgQuickSort. Adapted from
  ! http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#FPr. Thw
  ! index 'ind' is initialzed inside the algorithm

  use constants
  implicit none

  ! DUMMY ARGUMENTS
  integer(kind=intType), intent(in) :: nA
  real(kind=realType), dimension(nA), intent(inout) :: A
  integer(kind=intType), dimension(nA), intent(inout) :: ind

  ! LOCAL VARIABLES
  integer(kind=intType) :: left, right, itemp, i
  real(kind=realType) :: random, pivot, temp
  integer(kind=intType) :: marker

  if (nA > 1) then

     ! random pivot (not best performance, but avoids worst-case)
     ! call random_number(random)

     ! Just take mid point...complex code doesn't like random_number
     random = half
     i = int(random*real(nA-1))+1
   
     pivot = A(i) 
     left = 0
     right = nA + 1

     do while (left < right)
        right = right - 1
        do while (A(right) > pivot)
           right = right - 1
        end do
        left = left + 1
        do while (A(left) < pivot)
           left = left + 1
        end do
        if (left < right) then
           ! Swap value 
           temp = A(left)
           A(left) = A(right)
           A(right) = temp
           ! ! And swap index
           iTemp = ind(left)
           ind(left) = ind(right)
           ind(right) = itemp

        end if
     end do

     if (left == right) then
        marker = left + 1
     else
        marker = left
     end if

     call argQSort(A(:marker-1), marker-1, ind(:marker-1))
     call argQSort(A(marker:), nA-marker+1, ind(marker:))

  end if

end subroutine ArgQSort
