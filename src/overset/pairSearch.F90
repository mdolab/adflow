!
!      ******************************************************************
!      *                                                                *
!      * pairSearch searches block 'A' for potential donors in block 'B'*
!      *                                                                *
!      ******************************************************************

subroutine pairSearch(A, B)

  use constants
  use overset
  use inputOverset
  use adtLocalSearch

  implicit none

  ! Input/Ouput
  type(oversetBlock), intent(inout) :: A, B

  ! Working Varaibles
  integer(kind=intType) :: i, j, k, ii, jj, kk, iii, jjj, kkk, l, mm, mmm, tmp
  integer(kind=intType) :: nCoor, nHexa, nInterpol, elemID, nalloc
  logical :: invalidDonors
  real(kind=realType) :: uvw(4)
  integer(Kind=intType) ::intInfo(3), jjADT, localSearches
  real(kind=realType), dimension(:, :), allocatable ::  arrDonor

  real(kind=realType) :: donorQual
  real(kind=realType) :: uvwQuadratic(3), xElem(3, -1:1, -1:1, -1:1)

  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: BB
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew

  ! Now search the coordinates from *A* using the
  ! containment search

  localSearches = 0
  if (oversetInterpolation == linear) then 

     allocate(arrDonor(1, B%ie*B%je*B%ke))
     mm = 0
     do k=1,B%ke
        do j=1, B%je
           do i=1, B%ie
              mm = mm + 1
              arrDonor(1, mm) = B%qualDonor(i, j, k)
           end do
        end do
     end do
     nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 
  else
     ! For the quadratic we dont...we just take the index of the
     ! primal volume it found and use that for the volume comparison. 
     nInterpol = 0
  end if

  nAlloc = ubound(ADTs, 1)
  do jjAdt=1,nAlloc
     if(B%adtName == ADTs(jjAdt)%adtID) exit
  enddo
  
  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(20), frontLeaves(25), frontLeavesNew(25))

  ! Search the cells one at a time:
    do k=2, A%kl
     do j=2, A%jl
        do i=2, A%il


           ! We can make a short cut: If the min cell volume of 'B' is
           ! larger than the cell volume i'm looking for, we can skip
           ! the search *completely*, since it there cannot possibly
           ! be a donor for this cell. This can significantly reduce
           ! the total cost. The reason this is so important is that
           ! typically background meshes (especially the largest one)
           ! fully contain all the other meshes. This means that those
           ! cells will *always* find "potential" donors on the
           ! background mesh, which means it does a costly newton
           ! search. However, the vast majority of cells in the
           ! nearfield meshes have cell volumes smaller than the
           ! smallest volume on the background mesh. Therefore, there
           ! no possible way that a potential donor can be found and a
           ! very simple check (done below) can elminate the current
           ! cell. Initial experiments indicate 30-40% in total search
           ! time. Of course, if the cell in question is a
           ! forcedReciver, we can't make this short-cut.
           
           if (B%minVolume > A%qualRecv(i, j, k) .and. &
                .not. A%forceRecv(i, j, k)) then 
              cycle
           end if

           ! Otherwise, we do have to do search. Note that this uses a
           ! lower-level local API that doesn't deal with any of the
           ! parallel ADT stuff. 

           call containmentTreeSearchSinglePoint(jjADT, A%xSearch(:, i,j, k), &
                intInfo, uvw, arrDonor, nInterpol, BB, frontLeaves, frontLeavesNew)

           ! For reporting purposes keep track of the number of
           ! searches 
           localSearches = localSearches + 1

           elemFound: if (.not. intInfo(1) < 0) then 
    
              ! This is where it splits again:
              if (oversetInterpolation == linear) then 
                 
                 donorQual = uvw(4)
                 
                 ! Two conditions are check: IF the donorQual is
                 ! better than my quality, we accept the donor. If we
                 ! have a forced receiver that does not yet have donor
                 ! assigned to it, we also accept it, since it may be
                 ! only one we get. 

                 ! Unwind the element indices for the donor
                 ! Remember we have (il, jl, kl) elements in the dual
                 ! mesh. 
                 elemID = intInfo(3) 
                 tmp = elemID - 1
                 ii = mod(tmp, B%il) + 1
                 jj = mod(tmp/B%il, B%jl) + 1
                 kk = tmp/(B%il*B%jl) + 1

                 ! We need to check if any of the 8 dual nodes to make
                 ! sure that they are valid donors

                 invalidDonors = .False.
                 do iii=ii, ii+1
                    do jjj=jj, jj+1
                       do kkk=kk, kk+1
                          if (B%invalidDonor(iii, jjj, kkk)) then 
                             invalidDonors = .True.
                          end if
                       end do
                    end do
                 end do

                 if ((.not. invalidDonors) .and. ( &
                      donorQual < A%qualRecv(i, j, k) .or.  &
                      A%forceRecv(i, j, K) .and. &
                      A%donors(i, j, k)%donorProcID ==-1)) then 

                    ! Save the necessary donor infor in my (A's) donor
                    ! array. If there was an exsting donor there, and
                    ! we found better one, we just overwrite what was
                    ! there. 
                    A%donors(i, j, k)%donorProcID = B%procID
                    A%donors(i, j, k)%donorBlockID = B%globalBlockID
                    A%donors(i, j, k)%frac = uvw(1:3)
                    A%donors(i, j, k)%ind = (/ii, jj, kk/)
                    
                    ! Donor is better than me: Flag myself as a fringe.
                    A%iBlank(i, j, k) = -1
                    A%recvStatus(i, j, k) = .True. 

                    ! Set my receiver quality as what we found from
                    ! the donor. This way any potential new donor has
                    ! to beat the current best. 
                    A%qualRecv(i, j, k) = donorQual
                    
                    ! Notet that we don't modify any information on
                    ! block "B" since this could be a remote block. We
                    ! store all information we need on the receiver
                    ! side. We'll deal with the communication
                    ! elsewhere.



                    ! Save the fractions. 
                    ! A%donorFrac(:, A%nFringe) = uvw(:, mm)

                    ! ! Flag donor cell(s) as donor
                    ! B%donorStatus(ii  , jj  , kk  ) = B%donorStatus(ii  , jj  , kk  ) + 1
                    ! B%donorStatus(ii+1, jj  , kk  ) = B%donorStatus(ii+1, jj  , kk  ) + 1
                    ! B%donorStatus(ii  , jj+1, kk  ) = B%donorStatus(ii  , jj+1, kk  ) + 1
                    ! B%donorStatus(ii+1, jj+1, kk  ) = B%donorStatus(ii+1, jj+1, kk  ) + 1 
                    ! B%donorStatus(ii  , jj  , kk+1) = B%donorStatus(ii  , jj  , kk+1) + 1
                    ! B%donorStatus(ii+1, jj  , kk+1) = B%donorStatus(ii+1, jj  , kk+1) + 1
                    ! B%donorStatus(ii  , jj+1, kk+1) = B%donorStatus(ii  , jj+1, kk+1) + 1
                    ! B%donorStatus(ii+1, jj+1, kk+1) = B%donorStatus(ii+1, jj+1, kk+1) + 1
                    
                    ! Getting the global set of indices is slightly
                    ! tricker than it sounds. We frist have to
                    ! reconstruct the i,j,k indices of mm (which is in
                    ! block B ordering). From this "cell" we can get the
                    ! i,j,k coordinates of the nodes, which are actually
                    ! the cells of the primal mesh. We can then index
                    ! into the globalCell array
                    
                    ! ! that unwinds our index. Now do the 2x2x2 loop and reassemble the index
                    ! A%donorIndices(1, A%nFringe) = B%globalCell(ii  , jj  , kk  )
                    ! A%donorIndices(2, A%nFringe) = B%globalCell(ii+1, jj  , kk  )
                    ! A%donorIndices(3, A%nFringe) = B%globalCell(ii  , jj+1, kk  )
                    ! A%donorIndices(4, A%nFringe) = B%globalCell(ii+1, jj+1, kk  )
                    ! A%donorIndices(5, A%nFringe) = B%globalCell(ii  , jj  , kk+1)
                    ! A%donorIndices(6, A%nFringe) = B%globalCell(ii+1, jj  , kk+1)
                    ! A%donorIndices(7, A%nFringe) = B%globalCell(ii  , jj+1, kk+1)
                    ! A%donorIndices(8, A%nFringe) = B%globalCell(ii+1, jj+1, kk+1)
                    ! do l=1,8
                    !    if (A%donorIndices(l, A%nFringe) < 0) then 
                    !       A%donorIndices(l, A%nFringe) = 0
                    !    end if
                    ! end do
                 end if
              else
                 ! Quadratic interpolation 

                 ! Unwind the element indices for the donor
                 ! Remember we have (nx, ny, nz) elements in the primal mesh
                 ! mesh. 
                 ! elemID = elementID(mm)
                 ! tmp = elemID - 1
                 ! ii = mod(tmp, B%nx) + 1
                 ! jj = mod(tmp/B%nx, B%ny) + 1
                 ! kk = tmp/(B%nx*B%ny) + 1

                 ! ! The extra +1 is due that qualDonor has first level
                 ! ! halos and the ii, jj, kk indices above start at 1
                 ! ! for the first real cell. 

                 ! if (B%qualDonor(ii+1, jj+1, kk+1) < overlapFactor*A%qualRecv(i, j, k)) then 
                   
                    
                 !    ! Increment the number of fringes on this block
                 !    A%nFringe = A%nFringe + 1
                    
                 !    ! Save the indices of this fringe:
                 !    A%fringeIndices(1, A%nFringe) = i
                 !    A%fringeIndices(2, A%nFringe) = j
                 !    A%fringeIndices(3, A%nFringe) = k
                    
                 !    ! Donor is better than me: Flag myself as a fringe.
                 !    A%iBlank(i, j, k) = -1
                    
                 !    ! Set my quality as a receiver as my donor
                 !    A%qualRecv(i, j, k) = donorQual
                    
                 !    ! Unwind the element indices for the donor
                 !    ! Remember we have (nx, ny, nz) elements in the primal mesh
                 !    tmp = elemID - 1
                 !    ii = mod(tmp, B%nx) + 2
                 !    jj = mod(tmp/B%nx, B%ny) + 2
                 !    kk = tmp/(B%nx*B%ny) + 2
                    
                 !    ! For the quadratic interpolation we're not done
                 !    ! yet: We have only determined what primal cell
                 !    ! the point lies in. We will now form a quadratic
                 !    ! (3x3x3) element from the dual mesh around this
                 !    ! cell. Since the newton search was done on the
                 !    ! primal cell, this will give us a (very) good
                 !    ! starting point for this secondary Newton Raphson
                 !    ! iteration. 



                 !    ! ||====================================||
                 !    ! ||             ||           ||        ||
                 !    ! ||     +------ ||-----+-----||----+   ||
                 !    ! ||     |       ||     |     ||    |   ||
                 !    ! ||====================================||
                 !    ! ||     |       ||     |     ||    |   ||
                 !    ! ||     +       ||     +     ||    +   ||
                 !    ! ||     |       ||     |  x  ||    |   ||
                 !    ! ||====================================||
                 !    ! ||     |       ||     |     ||    |   ||
                 !    ! ||     +------ ||-----+-----||----+   ||
                 !    ! ||             ||           ||        ||
                 !    ! ||=============||===========||========||                    
                 !    ! 
                 !    ! Double liens are the primal mesh. Single lines
                 !    ! are the dual mesh formed from the cell
                 !    ! centers. We now now the cell that the 'x' is in
                 !    ! due to the ADT search above. Now we want to know
                 !    ! what the (u,v,w) coordinates of the point 'x' is
                 !    ! in the quadratic element formed by the thin
                 !    ! lines. The uvw cordiantes of 'x' in the center
                 !    ! element can be resued. 

                 !    l = 0
                 !    do kkk=-1,1
                 !       do jjj=-1,1
                 !          do iii=-1, 1

                 !             ! Reconstitute the xDual index
                 !             mmm = B%ie*B%je * (kk + kkk - 1) + B%ie * (jj + jjj - 1) + ii + iii

                 !             xElem(:, iii, jjj, kkk) = B%xDual(:, mmm)

                 !             l = l + 1
                 !             ! We know what the indices will be so set
                 !             ! them now while we're at it.
                 !             A%donorIndices(l, A%nFringe) = B%globalCell(ii+iii, jj+jjj, kk+kkk)

                 !             if (A%donorIndices(l, A%nFringe) < 0) then 
                 !                A%donorIndices(l, A%nFringe) = 0
                 !             end if
                 !          end do
                 !       end do
                 !    end do

                 !    ! Estimate our uvw. Since we use finite elements,
                 !    ! uvw will be in range [-1, 1]. Since uvw was the
                 !    ! range 0,1 the effect of the larger mesh cancels
                 !    ! out the former 0,1 range giving just an offset
                 !    ! of 1/2
                 !    uvwQuadratic = uvw(:, mm) - half
                    
                 !    call computeQuadraticWeights(A%xSearch(:, mm), xElem, uvwQuadratic)

                 !    ! Now that we know the uvw coefficients in the
                 !    ! quadratic element we can set the weights and
                 !    ! we're good to go. 
                    
                 !    ! Save the fractions. 
                 !    A%donorFrac(:, A%nFringe) = uvwQuadratic
                             
              ! end if
              end if
           end if elemFound
        end do
     end do
  end do

  ! Ditch temporary memory
  if (nInterpol == 1) then 
     deallocate(arrDonor)
  end if

  totalSearches = totalSearches + localSearches
  searchCosts = searchCosts + localSearches * log(one*(B%ie * B%je * B%ke))
  deallocate(BB, frontLeaves, frontLeavesNew)
  
end subroutine pairSearch
