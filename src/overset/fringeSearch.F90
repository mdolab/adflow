subroutine fringeSearch(oBlock, oFringe)

  use constants
  use block, only : fringeType
  use overset, only : oversetBlock, oversetFringe, localWallFringes, &
       tmpFringePtr, nLocalWallFringe
  use inputOverset, onlY : overlapFactor, oversetProjTol
  use adtLocalSearch, only :  mindistancetreesearchsinglepoint, &
       containmenttreesearchsinglepoint
  use adtData, only : adtBBoxTargetType
  use adtUtils, only : stack
  use utils, only : mynorm2
  use oversetUtilities, only : fracToWeights2, addToFringeList, tic, toc, windIndex
  implicit none

  type(oversetBlock), intent(inout) :: oBlock
  type(oversetFringe), intent(inout) :: oFringe

  integer(kind=intType) :: idom, jdom

  ! Working Varaibles
  integer(kind=intType) :: nInterpol, elemID, nalloc, intInfo(3), intInfo2(3)
  integer(kind=intType) :: i, ii, jj, kk, j, nn, myI, myJ, myK
  integer(kind=intTYpe) :: iii, jjj, kkk, n, myind,  nx, ny, nz, myindex
  logical :: invalid, failed
  real(kind=realType) :: uu, vv, ww, err1, err2
  real(kind=realType) :: uvw(5), uvw2(5), xx(4), pt(3), xcheck(3)
  real(kind=realType), dimension(:, :), allocatable :: offset
  real(kind=realType) :: oneMinusU, oneMinusV, oneMinusW, weight(8)
  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: BB
  type(adtBBoxTargetType), dimension(:), pointer :: BB2
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew
  type(fringeType) :: fringe
  nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(20), BB2(20), frontLeaves(25), frontLeavesNew(25), stack(100))

  ! Number of fringes we have:
  n = size(oFringe%x, 2)

  ! Offset vector:
  allocate(offset(3, n))
  offset = zero
  call tic(iSurfaceCorrection)
  call surfaceCorrection(oBlock, oFringe, offset, n)
  call toc(iSurfaceCorrection)

  call tic(iDonorSearch)
  ! Search the cells one at a time:
  do i=1, n

     ! Compute the potentailly offset point to search for. 
     xx(1:3) = oFringe%x(:, i) + offset(:, i)

     call containmentTreeSearchSinglePoint(oBlock%ADT, xx, intInfo, uvw, &
          oBlock%qualDonor, nInterpol, BB, frontLeaves, frontLeavesNew, failed)
     
     if (intInfo(1) >= 0) then 
        call fracToWeights2(uvw(1:3), weight)
        xcheck = zero
        do j=1,8
           xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
        end do
        
        if (mynorm2(xcheck - xx(1:3)) > oversetProjTol) then 
           failed = .True.
        end if
     end if
     
     if (intInfo(1) >= 0 .and. failed) then 
        ! we "found" a point but it is garbage. Do the failsafe search
        xx(4) = large
        call minDistanceTreeSearchSinglePoint(oBlock%ADT, xx, intInfo, uvw, &
             oBlock%qualDonor, nInterpol, BB2, frontLeaves, frontLeavesNew)
        
        ! Check this one:
        call fracToWeights2(uvw(1:3), weight)
        xcheck = zero
        do j=1,8
           xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
        end do
        
        ! Since this is the last line of defence, relax the tolerance a bit
        if (mynorm2(xcheck - xx(1:3)) > 100*oversetProjTol) then 
           ! This fringe has not found a donor
           intInfo(1) = -1
        else
           ! This one has now passed.
           
           ! Important! uvw(4) is the distance squared for this search
           ! not
           uvw(4) = uvw(5)
        end if
        
     end if
     
     elemFound: if (intInfo(1) >= 0) then 

        ! Recompute the i,j,k indices on the donor
        elemID = intInfo(3) - 1 ! Make it zero based for the modding. 
        ii = mod(elemID, oBlock%il) + 1
        jj = mod(elemID/oBlock%il, oBlock%jl) + 1
        kk = elemID/(oBlock%il*oBlock%jl) + 1
        
        ! Now record the information onto the fringe
        fringe%donorProc = oBlock%proc
        fringe%donorBlock= oBlock%block
        fringe%dIndex    = windIndex(ii, jj, kk, oBlock%il, oBlock%jl, oBlock%kl)
        fringe%donorFrac = uvw(1:3)
        fringe%quality   = uvw(4)
                
        ! Also save the information about where it came from,
        ! we need this to combine everything together at the end. 
        fringe%myBlock = oFringe%block

        myI = mod((i-1), oFringe%nx) + 2
        myJ = mod((i-1)/oFringe%nx, oFringe%ny) + 2
        myK = (i-1)/(oFringe%nx*oFringe%ny) + 2
        fringe%myIndex = windIndex(myI, myJ, myK, oFringe%il, oFringe%jl, oFringe%kl)
       
        ! Store the donor in the big flat list if it isn't invalid
        invalid = .False.
        do kkk=0,1
           do jjj=0,1
              do iii=0,1
                 if (oBlock%invalidDonor(ii+iii, jj+jjj, kk+kkk)  .ne. 0) then 
                    invalid = .True.
                 end if
              end do
           end do
        end do
             
        if (.not. invalid) then 
           call addToFringeList(oFringe%fringes, oFringe%nDonor, fringe)
        end if

        ! Save the fringe to the wallList. Note that we have to do
        ! this *after* the actual fringeList becuase we may modify the
        ! dI, dJ, dK here.

        if ((oFringe%isWall(i) > 0) .and. .not. (oBlock%nearWall(ii, jj, kk) == 1)) then 
 
           ! Here we have to recompute the i,j,k indices since we may
           ! need to modify them based on the frac. 
           
           if (uvw(1) >= half) then 
              ii = ii + 1
           end if
           
           if (uvw(2) >= half) then 
              jj =jj + 1
           end if
           
           if (uvw(3) >= half) then 
              kk =kk + 1
           end if

           ! Recompute the full index for the wall. 
           fringe%dIndex = windIndex(ii, jj, kk, oBlock%il, oBlock%jl, oBlock%kl)

           call addToFringeList(localWallFringes, nLocalWallFringe, fringe)
        end if
     end if elemFound
  end do
  deallocate(offset, BB, BB2, frontLeaves, frontLeavesNew, stack)
  call toc(iDonorSearch)
end subroutine fringeSearch
