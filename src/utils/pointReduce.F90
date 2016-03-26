subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)

  ! Given a list of N points (pts) in three space, with possible
  ! duplicates, (to within tol) return a list of the nUnique
  ! uniquePoints of points and a link array of length N, that points
  ! into the unique list
  use precision
  use kdtree2_module
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: N
  real(kind=realType), intent(in), dimension(3, N) :: pts
  real(kind=realType), intent(in) :: tol

  ! Output Parametres
  real(kind=realType), intent(out), dimension(3, N) :: uniquePts
  integer(kind=intType), intent(out), dimension(N) :: link
  integer(kind=intType), intent(out) :: nUnique

  ! Working paramters
  type(kdtree2), pointer :: mytree
  real(kind=realType) :: tol2, timeb, timea
  integer(kind=intType) :: nFound, i, j, nAlloc
  type(kdtree2_result), allocatable, dimension(:) :: results

  if (N==0) then 
     nUnique = 0
     return 
  end if

  ! We will use the KD_tree to do most of the heavy lifting here:
 
  mytree => kdtree2_create(pts, sort=.True.)

  ! KD tree works with the square of the tolerance
  tol2 = tol**2

  ! Unlikely we'll have more than 20 points same, but there is a
  ! safetly check anwyay.
  nalloc = 20
  allocate(results(nalloc))

  link = 0
  nUnique = 0

  ! Loop over all nodes
  do i=1, N
     if (link(i) == 0) then 
        call kdtree2_r_nearest(mytree, pts(:, i), tol2, nFound, nAlloc, results)

        ! Expand if necesary and re-run
        if (nfound > nalloc) then 
           deallocate(results)
           nalloc = nfound
           allocate(results(nalloc))
           call kdtree2_r_nearest(mytree, pts(:, i), tol2, nFound, nAlloc, results)
        end if

        if (nFound == 1) then 
           ! This one is easy, it is already a unique node
           nUnique = nUnique + 1
           link(i) = nUnique
           uniquePts(:, nUnique) = pts(:, i)
        else
           if (link(i) == 0) then 
              ! This node hasn't been assigned yet:
              nUnique = nUnique + 1
              uniquePts(:, nUnique) = pts(:, i)

              do j=1, nFound
                 link(results(j)%idx) = nUnique
              end do
           end if
        end if
     end if
  end do

  ! Done with the tree and the result vector
  call kdtree2_destroy(mytree)
  deallocate(results)

end subroutine pointReduce
