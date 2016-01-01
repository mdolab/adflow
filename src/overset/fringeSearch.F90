
subroutine fringeSearch(oBlock, fringes, n)

  use constants
  use overset
  use inputOverset
  use adtLocalSearch
  implicit none

  type(oversetBlock), intent(inout) :: oBlock
  type(fringeType), dimension(n), intent(inout) :: fringes
  integer(kind=intType), intent(in) :: n

  ! Working Varaibles
  integer(kind=intType) :: nInterpol, elemID, nalloc, intInfo(3), i, ii, jj, kk, j, nn
  integer(kind=intTYpe) :: iii, jjj, kkk
  logical :: invalid
  real(kind=realType) :: uvw(4), donorQual

  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: BB
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew

  nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(20), frontLeaves(25), frontLeavesNew(25))

  ! Search the cells one at a time:
  do i=1, n

     ! We can take a little short cut here: If we know that the min
     ! volume of the oBlock is LARGER than the cell we're searching
     ! for, it would never be selcted. However, for flood seeding, we
     ! want to know what is regardless so there is a sepcial attribute
     ! of fringe to make sure we get an interpolant back, if there is
     ! one. Note that the check is for whether to do the operation or
     ! not, so we check that the minimum volume in the oBlock is LESS
     ! than the volume of the coordinate we're searching for.

     shortCut: if (  oBlock%minVol < fringes(i)%quality .or. fringes(i)%isWall) then 

        call containmentTreeSearchSinglePoint(oBlock%ADT, &
             fringes(i)%x, intInfo, uvw, oBlock%qualDonor, &
             nInterpol, BB, frontLeaves, frontLeavesNew)

        elemFound: if (.not. intInfo(1) < 0) then 

           ! Donor and block and index information for this donor. 
           donorQual = uvw(4)
           elemID = intInfo(3) - 1 ! Make it zero based
           ii = mod(elemID, oBlock%il) + 1
           jj = mod(elemID/oBlock%il, oBlock%jl) + 1
           kk = elemID/(oBlock%il*oBlock%jl) + 1
              
           ! If we found a donor and our fringe is a wall, we will
           ! record if only if it is also not "near" another wall (on
           ! the oBlock)
           if (fringes(i)%isWall .and. .not. oBlock%nearWall(ii, jj, kk) == 1) then 
              
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
              localWallFringes(nLocalWallFringe)%dI = ii 
              localWallFringes(nLocalWallFringe)%dJ = jj
              localWallFringes(nLocalWallFringe)%dK = kk

           end  if

           ! This is for the actual donors. The '<' is really the
           ! 'implicit hole' method where we take the best quality donor. 
           if ( donorQual < fringes(i)%quality) then

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
                 fringes(i)%donorProc = oBlock%proc
                 fringes(i)%donorBlock = oBlock%block
                 fringes(i)%dI = ii 
                 fringes(i)%dJ = jj
                 fringes(i)%dK = kk
                 
                 fringes(i)%donorFrac = uvw(1:3)
                 fringes(i)%isCompute = .False.

                 ! This is important: We need to reset the quality of
                 ! this fringe to the donorQuality. This is necessary
                 ! since we may be searching these same fringes local
                 ! against another oBlock and we will only accept a
                 ! different donor if it is better than the best one
                 ! so far.
                 fringes(i)%quality = donorQual

                 
                 ! Save the global indices as well. 
                 fringes(i)%gInd(1) = oBlock%globalCell(ii  , jj  , kk  )
                 fringes(i)%gInd(2) = oBlock%globalCell(ii+1, jj  , kk  )
                 fringes(i)%gInd(3) = oBlock%globalCell(ii  , jj+1, kk  )
                 fringes(i)%gInd(4) = oBlock%globalCell(ii+1, jj+1, kk  )
                 fringes(i)%gInd(5) = oBlock%globalCell(ii  , jj  , kk+1)
                 fringes(i)%gInd(6) = oBlock%globalCell(ii+1, jj  , kk+1)
                 fringes(i)%gInd(7) = oBlock%globalCell(ii  , jj+1, kk+1)
                 fringes(i)%gInd(8) = oBlock%globalCell(ii+1, jj+1, kk+1)
              end if
           end if
        end if elemFound
     end if shortCut
  end do

end subroutine fringeSearch
