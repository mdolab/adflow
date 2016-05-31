subroutine initBCDataiBlank(level, sps)

  use blockPointers
  use BCTypes
  use communication
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, sps

  ! Local variables
  integer(kind=intType) :: mm, nn, i, j, k, iBeg, iEnd, jBeg, jEnd
  logical :: side(4), isWallType

  integer(kind=intType), dimension(:, :), pointer :: ibp, gcp, frx
  integer(kind=intType), dimension(:, :, :), pointer :: forcedRecv
  integer(kind=intType), dimension(:, :), allocatable :: toFlip

  ! This routine initializes the surface cell iblank based on the
  ! volume iblank. It is not a straight copy since we a little
  ! preprocessing 

  ! This is a little trickier than it seems. The reason is that we
  ! will allow a limited number of interpolated cells to be used
  ! directly in the integration provided the meet certain criteria.
  
  ! Consider the following
  !  +------+--------+-------+
  !  |ib=1  |  ib=1  | ib= 1 |
  !  |      |        |       |
           !  |      |        |       |
  !  +------+========+-------+
  !  |ib=1  || ib=-1 || ib=1 |
  !  |      ||       ||      |
  !  |      ||       ||      |
  !==+======+--------+=======+==
  !  |ib=-1 |  ib=-1 | ib=-1 |
  !  |      |        |       |
  !  |      |        |       |
  !  +------+--------+-------+
  !
  ! The boundary between real/interpolated cells is marked by double
  ! lines. For zipper mesh purposes, it is generally going to be
  ! better to treat the center cell, as a regular force integration
  ! cell (ie surface iblank=1). The criteria for selection of these
  ! cells is:
  
  ! 1. The cell must not have been a forced receiver (ie at overset outer
  ! bound)
  ! 2. Any pair *opposite* sides of the cell must be compute cells.
  
  ! This criterial allows one-cell wide 'slits' to be pre-eliminated.


  domainLoop: do nn=1, nDom
     call setPointers(nn, level, sps)
     allocate(forcedRecv(1:ie, 1:je, 1:ke))
     call flagForcedReceivers(forcedRecv)

     bocoLoop: do mm=1, nBocos
        wallType: if (isWallType(BCType(mm))) then
           select case (BCFaceID(mm))
           case (iMin)
              ibp => iblank(2, :, :)
              gcp => globalCell(2, :, :)
              frx => forcedRecv(2, :, :)
           case (iMax)
              ibp => iblank(il, :, :)
              gcp => globalCell(il, :, :)
              frx => forcedRecv(il, :, :)
           case (jMin)
              ibp => iblank(:, 2, :)
              gcp => globalCell(:, 2, :)
              frx => forcedRecv(:, 2, :)
           case (jMax)
              ibp => iblank(:, jl, :)
              gcp => globalCell(:, jl, :)
              frx => forcedRecv(:, jl, :)
           case (kMin)
              ibp => iblank(:, :, 2)
              gcp => globalCell(:, :, 2)
              frx => forcedRecv(:, :, 2)
           case (kMax)
              ibp => iblank(:, :, kl)
              gcp => globalCell(:, :, kl)
              frx => forcedRecv(:, :, kl)
           end select

           ! -------------------------------------------------
           ! Step 1: Set the (haloed) cell iBlanks directly from
           ! the volume iBlanks
           ! -------------------------------------------------
           jBeg = BCData(mm)%jnBeg+1 ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg+1 ; iEnd = BCData(mm)%inEnd
        
           ! Just set the cell iblank directly from the cell iblank
           ! above it. Remember the +1 in ibp is for the pointer
           ! offset. These ranges *ALWAYS* give 1 level of halos
           do j=jBeg-1, jEnd+1
              do i=iBeg-1, iEnd+1
                 BCData(mm)%iBlank(i,j) = ibp(i+1, j+1)
              end do
           end do

           ! Make bounds a little easier to read. Owned cells only
           ! from now on.
     
           ! -------------------------------------------------
           ! Step 2: Slit elimination
           ! -------------------------------------------------

           ! Now we loop back through the cells again. For
           ! interpolated cells with iblank=-1 we see if it satifies
           ! the criteria above. If so we flag it wil "toFlip" =
           ! 1. Note that we can't set a particular iblank directly
           ! since that could cause an "avalance" effect with the
           ! later iterations using the updated iblank from a previous
           ! iteration.
           allocate(toFlip(iBeg:iEnd, jBeg:jEnd))
           toFlip = 0
           do j=jBeg, jEnd
              do i=iBeg, iEnd

                 ! We *might* add it if the interpolated cell is
                 ! touching two real cell on opposite sides.
                 if (BCData(mm)%iBlank(i, j) == -1 .and. validCell(i, j)) then

                    ! Reset the side flag
                    side = .False.

                    if (validCell(i-1, j) .and. BCData(mm)%iBlank(i-1, j) == 1) then
                       side(1) = .True.
                    end if

                    if (validCell(i+1, j) .and. BCData(mm)%iBlank(i+1, j) == 1) then
                       side(2) = .True.
                    end if

                    if (validCell(i, j-1) .and. BCData(mm)%iBlank(i, j-1) ==1 ) then
                       side(3) = .True.
                    end if

                    if (validCell(i, j+1) .and. BCData(mm)%iBlank(i, j+1) == 1) then
                       side(4) = .True.
                    end if

                    if ((side(1) .and. side(2)) .or. (side(3) .and. side(4)))  then
                       toFlip(i,j) = 1
                    end if
                 end if
              end do
           end do

           ! Now just set the cell surface iblank to 1 if we
           ! determined above we need to flip  the cell

           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 if (toFlip(i, j) == 1) then
                    BCData(mm)%iBlank(i, j) = 1
                 end if
              end do
           end do
           deallocate(toFlip)
        end if wallType
     end do bocoLoop
     deallocate(forcedRecv)
  end do domainLoop

contains
  function validCell(i, j)
    implicit none
    integer(kind=intType), intent(in) :: i, j
    logical :: validCell

    ! for our purposes here, a valid cell is one that:
    ! 1. Is not a boundary halo. ie has globalCell >= 0
    ! 2. It is not a force receiver.

    validCell = .False.
    if (gcp(i+1, j+1) >= 0 .and. frx(i, j) == 0) then
       validCell = .True.
    end if
  end function validCell
end subroutine initBCDataiBlank
