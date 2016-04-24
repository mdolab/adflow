subroutine fringeSearch(oBlock, oFringe, bWall, fWall)

  use constants
  use overset
  use inputOverset
  use adtLocalSearch
  implicit none

  type(oversetBlock), intent(inout) :: oBlock
  type(oversetFringe), intent(inout) :: oFringe
  type(oversetWall), intent(inout) :: bWall, fWall
  integer(kind=intType) :: idom, jdom

  ! Working Varaibles
  integer(kind=intType) :: nInterpol, elemID, nalloc, intInfo(3), i, ii, jj, kk, j, nn
  integer(kind=intTYpe) :: iii, jjj, kkk, n, myind,  nx, ny, nz, myindex
  logical :: invalid, failed
  real(kind=realType) :: uvw(5), donorQual, xx(4), pt(3)
  real(kind=realType), dimension(:, :), allocatable :: offset
  real(kind=realType) :: oneMinusU, oneMinusV, oneMinusW, weight(8)
  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: BB
  type(adtBBoxTargetType), dimension(:), pointer :: BB2
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew

  nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(20), BB2(20), frontLeaves(25), frontLeavesNew(25), stack(100))

  ! Number of fringes we have:
  n = size(oFringe%x, 2)

  ! Offset vector:
  allocate(offset(3, n))
  offset = zero

  ! Determine if we have a wall-wall overlap:
  if (bWall%nNodes /= 0 .and. fWall%nNodes /= 0) then 
     call surfaceCorrection(oBlock, oFringe, bWall, fWall, offset, n)
  end if

  ! Search the cells one at a time:
  do i=1, n

     ! Cells with negative quality are those inside the body. Don't
     ! bother to search for them since don't need donors. 
     if (oFringe%quality(i) < zero) then 
        cycle
     end if

     ! We can take a little short cut here: If we know that the min
     ! volume of the oBlock is LARGER than the cell we're searching
     ! for, it would never be selcted. However, for flood seeding, we
     ! want to know what is regardless so there is a sepcial attribute
     ! of fringe to make sure we get an interpolant back, if there is
     ! one. Note that the check is for whether to do the operation or
     ! not, so we check that the minimum volume in the oBlock is LESS
     ! than the volume of the coordinate we're searching for.
     shortCut: if (  oBlock%minVol < oFringe%quality(i) .or.&
          (oFringe%isWall(i)>0)) then 

        ! Compute the potentailly offset point to search for. 
        xx(1:3) = oFringe%x(:, i) + offset(:, i)

        call containmentTreeSearchSinglePoint(oBlock%ADT, xx, intInfo, uvw, &
             oBlock%qualDonor, nInterpol, BB, frontLeaves, frontLeavesNew, failed)
        
        if (intInfo(1) >= 0 .and. failed) then 
           ! we "found" a point but it is garbage. Do the failsafe search
           xx(4) = large
           call minDistanceTreeSearchSinglePoint(oBlock%ADT, xx, intInfo, uvw, &
             oBlock%qualDonor, nInterpol, BB2, frontLeaves, frontLeavesNew)
           
        end if
        
        elemFound: if (intInfo(1) >= 0) then 

           ! Check if our solution is actually any good by evaluating the
           ! interpolated point. 


           ! Donor and block and index information for this donor. 
           donorQual = uvw(4)
           elemID = intInfo(3) - 1 ! Make it zero based
           ii = mod(elemID, oBlock%il) + 1
           jj = mod(elemID/oBlock%il, oBlock%jl) + 1
           kk = elemID/(oBlock%il*oBlock%jl) + 1

           ! If we found a donor and our fringe is a wall, we will
           ! record if only if it is also not "near" another wall (on
           ! the oBlock)

           if ((oFringe%isWall(i) > 0) .and.&
                .not. (oBlock%nearWall(ii, jj, kk) == 1)) then 
              nLocalWallFringe = nLocalWallFringe + 1
              
              ! Check if we still have enough room in our
              ! localWallFringe array and realloc if necessary
              ! necessary
              
              if (nLocalWallFringe > size(localWallFringes)) then 
                 nn = size(localWallFringes)
                 tmpFringePtr => localWallFringes            ! Pointer to existing data
                 allocate(localWallFringes(2*nn))            ! Allocate new space
                 localWallFringes(1:nn) = tmpFringePtr(1:nn) ! Copy exsitng values
                 deallocate(tmpFringePtr)                    ! Free original memory 
              end if
              
              ! Now record the information
              localWallFringes(nLocalWallFringe)%donorProc = oBlock%proc
              localWallFringes(nLocalWallFringe)%donorBlock = oBlock%block

               if (uvw(1) < half) then 
                  localWallFringes(nLocalWallFringe)%dI = ii
               else
                  localWallFringes(nLocalWallFringe)%dI = ii +1 
               end if

               if (uvw(2) < half) then 
                  localWallFringes(nLocalWallFringe)%dJ = jj
               else
                  localWallFringes(nLocalWallFringe)%dJ = jj +1 
               end if

               if (uvw(3) < half) then 
                  localWallFringes(nLocalWallFringe)%dK = kk
               else
                  localWallFringes(nLocalWallFringe)%dK = kk +1 
               end if
                 
           end if

           ! This check is for the actual donors. The '<' is really
           ! the 'implicit hole' method where we take the best quality
           ! donor.
         
           if ( donorQual < oFringe%quality(i)) then 
              invalid = .False.
              do kkk=0,1
                 do jjj=0,1
                    do iii=0,1
                       if (oBlock%invalidDonor(ii+iii, jj+jjj, kk+kkk)==1) then 
                          invalid = .True.
                       end if
                    end do
                 end do
              end do

              if (.not. invalid) then 

                 oFringe%donorProc(i) = oBlock%proc
                 oFringe%donorBlock(i) = oBlock%block
                 oFringe%dI(i) = ii
                 oFringe%dJ(i) = jj
                 oFringe%dK(i) = kk
                 oFringe%donorFrac(:, i) = uvw(1:3)

                 ! This is important: We need to reset the quality of
                 ! this fringe to the donorQuality. This is necessary
                 ! since we may be searching these same fringes local
                 ! against another oBlock and we will only accept a
                 ! different donor if it is better than the best one
                 ! so far.
                 oFringe%quality(i) = donorQual

                 ! Save the global indices as well. 
                 oFringe%gInd(1, i) = oBlock%globalCell(ii  , jj  , kk  )
                 oFringe%gInd(2, i) = oBlock%globalCell(ii+1, jj  , kk  )
                 oFringe%gInd(3, i) = oBlock%globalCell(ii  , jj+1, kk  )
                 oFringe%gInd(4, i) = oBlock%globalCell(ii+1, jj+1, kk  )
                 oFringe%gInd(5, i) = oBlock%globalCell(ii  , jj  , kk+1)
                 oFringe%gInd(6, i) = oBlock%globalCell(ii+1, jj  , kk+1)
                 oFringe%gInd(7, i) = oBlock%globalCell(ii  , jj+1, kk+1)
                 oFringe%gInd(8, i) = oBlock%globalCell(ii+1, jj+1, kk+1)
              end if
           end if
        end if elemFound
     end if shortCut
  end do
  deallocate(offset, BB, BB2, frontLeaves, frontLeavesNew, stack)

end subroutine fringeSearch
