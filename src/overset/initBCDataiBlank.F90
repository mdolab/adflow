!
! Initialize BCData cell data. 
!
subroutine initBCDataiBlank(level, sps)

  use blockPointers
  use BCTypes

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, sps

  ! Local variables
  integer(kind=intType) :: iStart, jStart, kStart, kEnd
  integer(kind=intType) :: mm, nn, i, j, k, ibc, jbc
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  integer(kind=intType) :: ins, ine, jns, jne, ics, ice, jcs, jce 
  logical :: side(4), isWallType

  integer(kind=intType), dimension(:, :), pointer :: ibp, gcp, frx
  integer(kind=intType), dimension(:, :, :), pointer :: tmp
  integer(kind=intType), dimension(:, :), allocatable :: toFlip

  ! This is a little trickier than it seems. The reason is that we
  ! will allow a limited number of interpolated cells to be used
  ! directly in the integratio provided the meet certain criteria. 

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
  ! 2. Any pair *oppose* sides of the cell must be compute cells. 

  ! This criterial allows one-cell wide 'slits' to be pre-eliminated. 

  domainLoop: do nn=1, nDom
     call setPointers(nn, level, sps)

     allocate(tmp(1:ie, 1:je, 1:ke))
     call flagForcedReceivers(tmp)

     bocoLoop: do mm=1, nBocos
        wallType: if (isWallType(BCType(mm))) then 

           select case (BCFaceID(mm))
           case (iMin)
              ibp => iblank(2, :, :)
              gcp => globalCell(2, :, :)
              frx => tmp(2, :, :)
           case (iMax)
              ibp => iblank(il, :, :)
              gcp => globalCell(il, :, :)
              frx => tmp(il, :, :)
           case (jMin)
              ibp => iblank(:, 2, :)
              gcp => globalCell(:, 2, :)
              frx => tmp(:, 2, :)
           case (jMax)
              ibp => iblank(:, jl, :)
              gcp => globalCell(:, jl, :)
              frx => tmp(:, jl, :)
           case (kMin)
              ibp => iblank(:, :, 2)
              gcp => globalCell(:, :, 2)
              frx => tmp(:, :, 2)
           case (kMax)
              ibp => iblank(:, :, kl)
              gcp => globalCell(:, :, kl)
              frx => tmp(:, :, kl)
           end select
           
           ! Make bounds a little easier to read. 1-level halos
           jBeg = BCData(mm)%jcBeg ; jEnd = BCData(mm)%jcEnd
           iBeg = BCData(mm)%icBeg ; iEnd = BCData(mm)%icEnd

           ! Just set the cell iblank directly from the cell iblank
           ! above it. Remember the +1 in ibp is for the pointer offset. 
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 BCData(mm)%iBlank(i,j) = ibp(i+1, j+1)
              end do
           end do

           ! Make bounds a little easier to read. Owned cells only
           ! from now on.
           jBeg = BCData(mm)%jnBeg+1 ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg+1 ; iEnd = BCData(mm)%inEnd

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
     deallocate(tmp)
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

