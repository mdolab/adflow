
subroutine statePreAllocation(onProc, offProc, wSize, stencil, N_stencil, &
     level, transposed)

  ! This is a generic function that determines the correct
  ! pre-allocation for on and off processor parts of the TRANSPOSED
  ! matrix. With overset, it is quite tricky to determine the
  ! transpose sparsity structure exactly, so we use an alternative
  ! approach. We proceed to determine the non-zeros of the untranposed
  ! matrix, but instead of assigning a non-zero the row we're looping
  ! over, we assign it to the column, which will become a row in the
  ! tranposed matrix. Since this requires communication we use a petsc
  ! vector for doing off processor values. This is not strictly
  ! correct since we will be using the real values as floats, but
  ! since the number of non-zeros per row is always going to be
  ! bounded, we don't have to worry about the integer/floating point
  ! conversions. 


  use blockPointers
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Subroutine Arguments
  integer(kind=intType), intent(in)  :: wSize
  integer(kind=intType), intent(in)  :: N_stencil
  integer(kind=intType), intent(in)  :: stencil(N_stencil, 3)
  integer(kind=intType), intent(out) :: onProc(wSize), offProc(wSize)
  integer(kind=intType), intent(in)  :: level
  logical, intent(in) :: transposed

  ! Local Variables
  integer(kind=intType) :: nn, i, j, k, sps, ii, jj, kk, iii, jjj, kkk, n, m, gc
  integer(kind=intType) :: iRowStart, iRowEnd, ierr
  integer(kind=intType), dimension((N_stencil-1)*8) :: cellBuffer, dummy
  Vec offProcVec
  logical :: overset
  real(kind=realType), pointer :: tmpPointer(:)

  
  call vecCreateMPI(sumb_comm_world, wSize, PETSC_DETERMINE, offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  

  ! Zero the cell movement counter
  ii = 0
  
  ! Set the onProc values for each cell to the number of "OFF" time
  ! spectral instances. The "on" spectral instances are accounted for
  ! in the stencil  
  onProc(:) = nTimeIntervalsSpectral-1
  offProc(:) = 0
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

                    ! -------------------- Non-transposed code ----------------
                    if (.not. transposed) then 
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
                       ! -------------------- Ttransposed code ----------------

                       ! Now we ALSO loop over the total number of
                       ! (unique) neighbours. However, instead of
                       ! adding to the non-zeros to the on/offproc for
                       ! row 'ii', we add them to the column index
                       ! which will be the row index for the
                       ! transposed matrix.
                       do jj=1, m
                          gc = cellBuffer(jj)
                          
                          if (gc >= irowStart .and. gc <= iRowEnd) then 
                             ! On processor values can be dealt with
                             ! directly since the diagonal part is square. 
                             onProc(gc-iRowStart + 1) = onProc(gc-iRowStart+1)  +1
                          else
                             ! The offproc values need to be sent to
                             ! the other processors and summed.
                             call VecSetValue(offProcVec, gc, real(1), ADD_VALUES, ierr)
                             call EChk(ierr, __FILE__, __LINE__)
                          end if
                       end do
                    end if
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

  ! Assemble the offproc vector. This doesn't take any time for the
  ! non-transposed operation.
  call VecAssemblyBegin(offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyEnd(offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  if (transposed) then 
     ! Pull the local vector out and convert it back to integers. 
     call VecGetArrayF90(offProcVec, tmpPointer, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     do i=1,wSize
        offProc(i) = int(tmpPointer(i) + half) ! Make sure, say 14.99999 is 15. 
     end do

     call VecRestoreArrayF90(offProcVec, tmpPointer, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  ! Done with the temporary offProcVec
  call vecDestroy(offProcVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine statePreAllocation

