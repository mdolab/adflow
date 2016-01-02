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
  call setPointers(nDom, 1, nTimeIntervalsSpectral)
  iRowEnd   = flowDoms(nDom, 1, nTimeIntervalsSpectral)%globalCell(il, jl, kl)

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        ! Loop over each real cell
        do k=2, kl
           do j=2, jl
              do i=2, il 

                 ! Increment ii ONLY for each each movement of center cell
                 ii = ii + 1

                 ! Loop over the cells in the provided stencil:
                 do jj=1, N_stencil
                    
                    ! Determine the cell we are dealing with 
                    iii = stencil(jj, 1) + i
                    jjj = stencil(jj, 2) + j 
                    kkk = stencil(jj, 3) + k 
                    
                    ! Check if the cell in question is a fringe or not:
                    if (iblank(iii, jjj, kkk) == 1) then 

                       ! Check if it is onProc
                       if (globalCell(iii, jjj, kkk) >= 0) then ! Real cell

                          if (globalCell(iii, jjj, kkk) >= irowStart .and. &
                               globalCell(iii, jjj, kkk) <= irowEnd) then
                             
                             ! Cell is on processor
                             onProc(ii) = onProc(ii) + 1
                          else
                             ! Cell is off processor
                             offProc(ii) = offProc(ii) + 1
                          end if
                       end if

                    else if (iblank(iii, jjj, kkk) == -1) then 
                       ! Fringe cell. We won't actually be putting the
                       ! values at this location, but we will looop
                       ! here over the donors for this fringe cell and
                       ! add non-zero entries as necessary. We still
                       ! have to do the same onProc/offProc check as
                       ! above since while the donor for a fringe cell
                       ! will be from another block that block could
                       ! be on-processor or off-processor. 
                       onProc(ii) = onProc(ii) + 1

                    else if (iblank(iii, jjj, kkk) == 0) then 
                       ! blanked cell. We just have the single
                       ! identity on the diagonal block. This will of
                       ! course always be onProc. 
                       onProc(ii) = onProc(ii) + 1
                    end if
                       


                 end do ! Stencil Loop
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop

end subroutine statePreAllocation

