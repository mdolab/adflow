subroutine statePreAllocation(onProc, offProc, wSize, stencil, N_stencil, &
     level)

  ! This is a generic function that determines the correct
  ! pre-allocation for on and off processor parts. It take in a
  ! "stencil" definition (look at modules/stencil.f90 for the
  ! definitions) and uses this to determine on and off proc values. 

  use constants
  use blockPointers, only : nDom, il, jl, kl, globalCell, flowDoms
  use inputTimeSpectral , only : nTimeIntervalsSpectral

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

