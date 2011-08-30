subroutine statePreAllocation(onProc,offProc,wSize,stencil,N_stencil)

  ! This is a generic function that determines the correct
  ! pre-allocation for on and off processor parts. It take in a
  ! "stencil" definition (look at modules/stencil.f90 for the
  ! definitions) and uses this to determine on and off proc values. 

  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    
  use BCTypes

  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(in)  :: N_stencil
  integer(kind=intType),intent(in)  :: stencil(N_stencil,3)
  integer(kind=intType),intent(out) :: onProc(wSize),offProc(wSize)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,sps,ii,jj
  integer(kind=intType) :: cell(3)
  integer(kind=intTYpe) :: onAdd,offAdd
  ii = 0
  ! Set the onProc values for each cell to the number of "OFF" time
  ! spectral instances. The "on" spectral instances are accounted for
  ! in the stencil

  onProc(:) = nTimeIntervalsSpectral-1
  offProc(:) = 0_intType 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 
                 ! Increment ii ONLY for each each movement of center cell
                 ii = ii + 1

                 ! Loop over the cells in the provided stencil:
                 do jj=1,N_stencil
                                  
                    ! Determine the cell we are dealing with 
                    cell = (/i,j,k/) + stencil(jj,:)

                    ! Fully inside:
                    if (cell(1) >=2 .and. cell(1) <= il .and. &
                        cell(2) >=2 .and. cell(2) <= jl .and. &
                        cell(3) >=2 .and. cell(3) <= kl) then
                       
                       onProc(ii) = onProc(ii) + 1

                    else ! BC or B2B

                       onAdd = 0
                       offAdd = 0
                  
                       ! Basically what we need to determine if
                       ! "cell_to_check" is on this proc or an off proc

                       ! Low I check
                       if (cell(1) < 2) then
                          call checkCell(iMin,j,k,onAdd,offAdd,1)
                       end if

                       ! High I Check
                       if (cell(1) > il) then
                          call checkCell(iMax,j,k,onAdd,offAdd,1)
                       end if

                       ! Low J check
                       if (cell(2) < 2) then
                          call checkCell(jMin,i,k,onAdd,offAdd,1)
                       end if
                    
                       ! High J Check
                       if (cell(2) > jl) then
                          call checkCell(jMax,i,k,onAdd,offAdd,1)
                       end if
                    
                       ! Low K check
                       if (cell(3) < 2) then
                          call checkCell(kMin,i,j,onAdd,offAdd,1)
                       end if

                       ! High K check
                       if (cell(3) > kl) then
                          call checkCell(kMax,i,j,onAdd,offAdd,1)
                       end if

                       if (offAdd >= 1) then
                          offProc(ii) = offProc(ii) + 1
                       else if (onAdd > 0) then
                          onProc(ii) = onProc(ii) + 1
                       end if
                          
                    end if

                 end do ! Stencil Loop
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop
  
end subroutine statePreAllocation

subroutine drdxPreAllocation(onProc,offProc,wSize)

  ! Get a good estimate of the number of non zero rows for the
  ! on-diagonal and off-diagonal portions of the matrix
  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    
  use BCTypes

  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(out) :: onProc(wSize),offProc(wSize)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,l,sps,ii

  onProc(:) = 8*3+8*3*(nTimeIntervalsSpectral-1) ! ALWAYS have the center cell ON-PROCESSOR
  offProc(:) = 0_intType 

  ii = 0 

  ! This is for the "Regular" drdx calculation. i.e. xadjb
 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 do l=1,nw
                    ii = ii + 1
                    if (i-1 < 2) then 
                       call checkCell(iMin,j,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (i+1 > il) then 
                       call checkCell(iMax,j,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (j-1 < 2) then 
                       call checkCell(jMin,i,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (j+1 > jl) then 
                       call checkCell(jMax,i,k,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (k-1 < 2) then 
                       call checkCell(kMin,i,j,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if

                    if (k+1 > kl) then 
                       call checkCell(kMax,i,j,onProc(ii),offProc(ii),12)
                    else
                       onProc(ii) = onProc(ii) + 12
                    end if
                 end do ! l loop
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop

  ! However, drdx is more complex since we ALSO have
  ! xblockcorners. These however, only show up for the cells that are
  ! along a symmetry plane. Lets try to estimate those. 

  ! THIS MAY NOT WORK FOR SPS CASE!!! 
  ii = 0
 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    ii = ii + 1
                    call checkCellSym(i,j,k,onProc(ii),24)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine drdxPreAllocation

subroutine checkCell(iface,i,j,onProc,offProc,addVal)

  use blockPointers
  use BCTypes
  use communication
  implicit None

  ! Subroutine Arguments

  integer(kind=intType), intent(in) :: iface,i,j
  integer(kind=intType), intent(inout) :: onProc,offProc
  integer(kind=intType), intent(in) :: addVal

  !local Variables
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd,mm,ll

  ! It is assumed blockPointers are already set for this block

  ! Basically what we want to do is take the cell defined by index i
  ! and j on face defined by iface and determine:
  ! 1. What block subface it is on.
  ! 2. Determine if this is a boundary condition or a block-match
  ! 3. If its a block-match is the connecting block on- or off-processor


  ! If its a BC condition, nothing is added onProc and offProc
  ! If its a block-match on-proc, addVal is added to onProc
  ! If its a block-match off-proc, addVal is added to offProc

  n1to1Loop: do mm=1,nsubFace

     ! Store the correct index for this subface, i.e. add the
     ! offset from the boundary subfaces.

     if (BCFaceID(mm) == iface) then ! Check Face
        if (iface == iMin .or. iface == iMax) then
           iBeg = jcBeg(mm) ; iEnd = jcEnd(mm)
           jBeg = kcBeg(mm) ; jEnd = kcEnd(mm)
        else if(iface == jMin .or. iface == jMax) then
           iBeg = icBeg(mm) ; iEnd = icEnd(mm)
           jBeg = kcBeg(mm) ; jEnd = kcEnd(mm)
        else
           iBeg = icBeg(mm) ; iEnd = icEnd(mm)
           jBeg = jcBeg(mm) ; jEnd = jcEnd(mm)
        end if

        ! Check to make sure cell is on this (possible) sub-face
        if (i>=iBeg .and. i<=iEnd .and. j>=jBeg .and. j<= jEnd) then

           if (neighproc(mm) == myid) then
              onProc = onProc + addVal
           else
              if (neighproc(mm) >=0) then
                 offProc = offProc + addVal
              end if
           end if
        end if

     end if
  end do n1to1Loop
end subroutine checkCell

subroutine checkCellSym(i,j,k,onProc,addVal)

  use blockPointers
  use BCTypes
  use communication
  implicit None

  ! Subroutine Arguments

  integer(kind=intType), intent(in) :: i,j,k
  integer(kind=intType), intent(inout) :: onProc
  integer(kind=intType), intent(in) :: addVal

  !local Variables
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd,mm,ll

  n1to1Loop: do mm=1,nsubFace
     if (BCType(mm) == Symm) then 
        ! Is cell i,j,k "on" the symmetry plane
        select case (BCFaceID(mm))
        case (iMin)
           if (i == 2) then
              onProc = onProc + addVal
           end if
        case (iMax)
           if (i == il) then
              onProc = onProc + addVal
           end if
        case (jMin)
           if (j == 2) then
              onProc = onProc + addVal
           end if
        case (jMax)
           if (j == jl) then
              onProc = onProc + addVal
           end if
        case (kMin)
           if (k == 2) then
              onProc = onProc + addVal
           end if
        case (kMax)
           if (k == kl) then
              onProc = onProc + addVal
           end if
        end select
     end if
  end do n1to1Loop
end subroutine checkCellSym

