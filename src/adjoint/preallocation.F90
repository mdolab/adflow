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
  integer(kind=intType) :: nn, i, j, k, sps, ii, jj, kk, iii, jjj, kkk, n, m, gc
  integer(kind=intType) :: iRowStart, iRowEnd
  integer(kind=intType), dimension((N_stencil-1)*8) :: cellBuffer, dummy

  logical :: overset
  
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

                 ! Increment the running ii counter ONLY for each each
                 ! movement of center cell
                 ii = ii + 1

                 ! Reset the running tally of the number of neighbours
                 n = 0

                 blankedTest: if (iblank(i, j, k) == 1) then 
                    
                    ! Short-cut flag for cells without interpolated
                    ! cells in it's stencil
                    overset = .False.

                    ! Loop over the cells in the provided stencil:
                    do jj=1, N_stencil
                    
                       ! Determine the cell we are dealing with 
                       iii = stencil(jj, 1) + i
                       jjj = stencil(jj, 2) + j 
                       kkk = stencil(jj, 3) + k 
                    
                       ! Index of the cell we are dealing with. Make
                       ! code easier to read
                       gc = globalCell(iii, jjj, kkk)

                       ! Check if the cell in question is a fringe or not:
                       if (iblank(iii, jjj, kkk) == 1) then 
                          ! regular cell, add to our list, if it is
                          ! not a boundary
                          if (gc >= 0) then 
                             n = n + 1
                             cellBuffer(n) = gc
                          end if

                       else if (iblank(iii, jjj, kkk) == -1) then 
                          ! Fringe cell. What we do here is loop over
                          ! the donors for this cell and add any
                          ! entries that are real cells
                          overset = .True.
                          do kk=1,8
                             gc = fringes(iii, jjj, kkk)%gInd(kk)
                             if (gc >= 0) then 
                                n = n + 1
                                cellBuffer(n) = gc
                             end if
                          end do
                       end if
                    end do

                    ! We have now added 'n' cells to our buffer. For
                    ! the overset interpolation case, it is possible
                    ! (actually highly likely) that the same donor
                    ! cells are used in multiple fringes. To avoid
                    ! allocating more space than necessary, we
                    ! unique-ify the values, producing 'm' unique
                    ! values. If overset wasn't present, we can be
                    ! sure that m=n and we simply don't do the unique
                    ! operation. 
                    
                    if (overset) then 
                       call unique(cellBuffer, n, m, dummy)
                    else
                       m = n
                    end if
                    
                    ! Now we loop over the total number of
                    ! (unique) neighbours we have and assign them
                    ! to either an on-proc or an off-proc entry:
                    do jj=1, m
                       gc = cellBuffer(jj)
                     
                       if (gc >= irowStart .and. gc <= iRowEnd) then 
                          onProc(ii) = onProc(ii) + 1
                       else
                          offProc(ii) = offProc(ii) + 1
                       end if
                    end do
                 else
                    ! Blanked and interpolated cells only need a single
                    ! non-zero per row for the identity on the diagonal.
                    onProc(ii) = onProc(ii) + 1
                 end if blankedTest
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! sps loop
  end do ! Domain Loop

end subroutine statePreAllocation
