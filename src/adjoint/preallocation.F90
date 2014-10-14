subroutine statePreAllocation(onProc, offProc, wSize, stencil, N_stencil, &
     level)

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
  integer(kind=intType), intent(in)  :: wSize
  integer(kind=intType), intent(in)  :: N_stencil
  integer(kind=intType), intent(in)  :: stencil(N_stencil, 3)
  integer(kind=intType), intent(out) :: onProc(wSize), offProc(wSize)
  integer(kind=intType), intent(in)  :: level

  ! Local Variables
  integer(kind=intType) :: nn, i, j, k, sps, ii, jj, iii, jjj, kkk
  integer(kind=intType) :: iRowStart, iRowEnd
  
  ! Zero the cell movement counter
  ii = 0
  
  ! Set the onProc values for each cell to the number of "OFF" time
  ! spectral instances. The "on" spectral instances are accounted for
  ! in the stencil  
  onProc(:) = nTimeIntervalsSpectral-1
  offProc(:) = 0_intType 

  ! Determine the range of onProc in dRdwT
  iRowStart = flowDoms(1, 1, 1)%globalCell(2,2,2)
  call SetPointers(nDom, 1, nTimeIntervalsSpectral)
  iRowEnd   = flowDoms(nDom, 1, nTimeIntervalsSpectral)%globalCell(il, jl, kl)

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        ! Loop over each Cell
        do k=2, kl
           do j=2, jl
              do i=2, il 

                 ! Increment ii ONLY for each each movement of center cell
                 ii = ii + 1

                 ! Loop over the cells in the provided stencil:
                 do jj=1, N_stencil
                    
                    ! Determine the cell we are dealing with 
                    iii = stencil(jj, 1)
                    jjj = stencil(jj, 2)
                    kkk = stencil(jj, 3)
                    
                    ! Check if it is onProc
                    if (globalCell(i + iii, j + jjj, k + kkk) >= 0) then ! Real cell
                       if (globalCell(i + iii, j + jjj, k+kkk) >= irowStart .and. &
                            globalCell(i + iii, j + jjj, k+kkk) <= irowEnd) then
                          
                          ! Increase onProc
                          onProc(ii) = onProc(ii) + 1
                       else
                          ! Increase offProc
                          offProc(ii) = offProc(ii) + 1
                       end if
                    end if
                 end do ! Stencil Loop
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop

end subroutine statePreAllocation

subroutine drdxPreAllocation(onProc, offProc, xSize, level)

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
  use stencils
  implicit none

  ! Subroutine Arguments
  integer(kind=intType), intent(in)  :: xSize
  integer(kind=intType), intent(out) :: onProc(xSize), offProc(xSize)
  integer(kind=intType), intent(in)  :: level

  ! Local Variables
  integer(kind=intType) :: nn, sps, ii, jj, mm, iii, jjj, kkk
  integer(kind=intType) :: inode, jnode, knode, iDim, irowStart, irowEnd

  integer(kind=intType), dimension(:, :), pointer :: stencil
  integer(kind=intType) :: n_stencil
  integer(kind=intType) :: ijk(3), cells_on_face, nState
  logical :: faces(6)
  logical :: is_corner, is_a_corner

  if (frozenTurbulence) then
     nState = nwf
  else
     nState = nw
  end if

  onProc(:) = 8_intType*nState*(nTimeIntervalsSpectral-1)
  offProc(:) = 0_intType 

  ii = 0  
  call initialize_stencils
  
  if (viscous) then
     stencil => visc_drdx_stencil
     n_stencil = N_visc_drdx
  else
     stencil => euler_drdx_stencil
     n_stencil = N_euler_drdx
  endif

  ! Determine the range of onProc in dRdx
  iRowStart = flowDoms(1, 1, 1)%globalCell(2,2,2)
  call SetPointers(nDom, 1, nTimeIntervalsSpectral)
  iRowEnd   = flowDoms(nDom, 1, nTimeIntervalsSpectral)%globalCell(il, jl, kl)

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        ! Loop over each Cell
        do kNode=1, kl
           do jNode=1, jl
              do iNode=1, il 
                 do iDim = 1, 3
                    
                    ! Increment ii ONLY for each each movement of center cell
                    ii = ii + 1
                    
                    ! Loop over the cells in the provided stencil:
                    do jj=1, n_stencil
                       
                       ! Determine the cell we are dealing with 
                       iii = stencil(jj, 1)
                       jjj = stencil(jj, 2)
                       kkk = stencil(jj, 3)
                       
                       ! Check if it is onProc
                       if (globalCell(iNode + iii, jNode + jjj, kNode + kkk) >= 0) then ! Real cell
                          if (globalCell(iNode + iii, jNode + jjj, kNode+kkk) >= irowStart .and. &
                               globalCell(iNode + iii, jNode + jjj, kNode+kkk) <= irowEnd) then
                          
                             ! Increase onProc
                             onProc(ii) = onProc(ii) + nState
                          else
                             ! Increase offProc
                             offProc(ii) = offProc(ii) + nState
                          end if
                       end if
                    end do ! Stencil Loop
                 end do ! Dim loop
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop

  ! However, drdx is more complex since we ALSO have
  ! xblockcorners. These however, only show up for the cells that are
  ! along a symmetry plane. Lets try to estimate those. The current
  ! implementation seems to over estimate the number of non-zeros and
  ! this doesn't appear necessary. However, this is so tiny, we leave
  ! it in anyway, just in case. 

  ! THIS MAY NOT WORK FOR SPS CASE!!! 
  ii = 0
  ! This is for the "Regular" drdx calculation.
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        ! Loop over each Node
        do kNode=1, kl
           do jNode=1, jl
              do iNode=1, il 
                 do iDim=1, 3
                    ii = ii + 1 ! Continuous counter

                    ijk = (/iNode, jNode, kNode/)
                    ! Determine if this is a corner:
                    
                    is_a_corner = IS_CORNER(ijk)

                    ! Determine what faces this corner is apart of
                    if (is_a_corner) then
                       call ON_WHICH_FACE(ijk, faces)

                       ! Determine if these faces are Symmetry faces
                       
                       symLoop: do mm=1, nsubFace
                          if (BCType(mm) == Symm) then 
                             if (faces(BCFaceID(mm))) then
                                cells_on_face = (icEnd(mm)-icBeg(mm) + 1) * &
                                     (jcEnd(mm)-jcBeg(mm) + 1) * &
                                     (kcEnd(mm)-kcBeg(mm) + 1)

                                onProc(ii) = onProc(ii) + nState*cells_on_face
                             end if
                          end if
                       end do symLoop
                    end if ! Corner Check
                 end do ! idim loop
              end do ! i Node Loop
           end do !j Node Loop
        end do !k Node Loop
     end do ! sps Loop
  end do! Domain Loop
end subroutine drdxPreAllocation
