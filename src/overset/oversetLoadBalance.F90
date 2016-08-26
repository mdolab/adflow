!
!       oversetLoadBalance determine the deistributation of donor and  
!       receiver blocks that will result in approximate even load      
!       balancing. The sparse matrix structrue of the overla is        
!       provided. This computation runs on all processors.             

subroutine oversetLoadBalance(overlap)

  use constants
  use communication
  use overset
  implicit none

  ! Input/Output
  type(CSRMatrix), intent(inout) :: overlap

  ! Working paramters
  integer(kind=intType) :: curRow, jj, jj1, iProc, iRow
  real(kind=realType) :: evenCost, potentialSum, targetCost
  real(Kind=realType) :: totalSearch, totalBuild

  real(kind=realType), dimension(0:nProc-1) :: procCosts
  real(kind=realType), dimension(0:nProc) :: cumProcCosts
  real(kind=realType), dimension(overlap%nRow) :: buildCost
  real(kind=realType), parameter :: tol=0.1_realType
!  real(kind=realType), parameter :: K=10_realType
  logical, dimension(overlap%nnz) :: blockTaken
  logical :: increment

  ! Pointers to make code a litte easier to read
  integer(kind=intType), pointer, dimension(:) :: rowPtr, assignedProc
  real(kind=realType), pointer, dimension(:) :: data

  ! Set the couple  of pointers
  rowPtr => overlap%rowPtr
  assignedProc => overlap%assignedProc
  data => overlap%data

  ! Determine the total search cost:
  totalSearch = sum(overlap%data)
 
  ! Target amount of work for each processor
  evenCost = totalSearch / nProc

  ! Initialize the taken processor to False
  blockTaken = .False. 

  ! Initialzie assignedProc to -1 since there could be entries we can
  ! ignore. 
  assignedProc(:) = -1
  procCosts = zero
  cumProcCosts(0) = zero

  ! Initialize the starting point
  jj = 1
  iProc = 0 

  ! Find the first row with non-zeros
  curRow = 1
  do while(rowPtr(curRow+1)-rowPtr(curRow) == 0) 
     curRow = curRow + 1
  end do
   
  masterLoop: do while (curRow <= overlap%nRow .and. iProc <= nProc) 

     ! Normally we increment
     increment = .True. 
     
     ! This is our current target cost.
     targetCost = evenCost*(iProc + 1)

     ! It is still possible that data(jj) is zero. That's ok...we'll
     ! explictly ignore them. 
     if (data(jj) /= zero .and. .not. (blockTaken(jj))) then 
     
        if (procCosts(iProc) == 0 .or. iProc == nProc-1) then 
           ! Must be added
           procCosts(iProc) = procCosts(iProc) + data(jj)
           blockTaken(jj) = .True.
           assignedProc(jj) = iProc

        else

           ! There is already something in there. See what the
           !  potential sum will be:
           potentialSum = cumProcCosts(iProc) + procCosts(iProc) + data(jj)

           if (potentialSum < targetCost - tol*evenCost) then 
              !  We are not close to our limit yet so just add it normally
              procCosts(iProc) = procCosts(iProc) + data(jj)
              blockTaken(jj) = .True.
              assignedProc(jj) = iProc

           else if (potentialSum >= targetCost - tol*evenCost  .and. &
                    potentialSum <= targetCost + tol*evenCost) then 
                
              ! This one looks perfect. Call it a day...add it and
              !  move on to the next proc

              procCosts(iProc) = procCosts(iProc) + data(jj)
              blockTaken(jj) = .True.
              assignedProc(jj) = iProc

              ! Processor can be incremented
              cumProcCosts(iProc+1) = cumProcCosts(iProc) + procCosts(iProc)
              iProc = iProc + 1
           else
              ! This means potentialSum > targetCost + tol*evenCost

              ! This is somewhat bad news...this may be *horrendly*
              ! load balanced. The algorithm dictates we *MUST*
              ! finish this proc no matter what before we go back to
              ! the outer loop. Essentially we know jj is bad,
              ! instead scan over the rest of the row and see if we
              ! can add something else that is decent. 
              increment = .False.

              restOfRow: do jj1=jj+1, rowPtr(curRow+1)-1
                    
                 potentialSum = cumProcCosts(iProc) + procCosts(iProc) + data(jj1)

                 if (data(jj1) /= zero .and. .not. (blockTaken(jj1))) then 

                    if (potentialSum < targetCost - tol*evenCost) then 
                       !Huh...that one fit in without going
                       ! over....add it and kep going in the loop

                       procCosts(iProc) = procCosts(iProc) + data(jj1)
                       blockTaken(jj1) = .True.
                       assignedProc(jj1) = iProc

                    else if (potentialSum >= targetCost - tol*evenCost  .and. &
                         potentialSum <= targetCost + tol*evenCost) then 

                       ! This one fit in perfectly. 
                       procCosts(iProc) = procCosts(iProc) + data(jj1)
                       blockTaken(jj1) = .True.
                       assignedProc(jj1) = iProc

                       ! No need to keep going
                       exit  restOfRow

                    end if
                 end if
              end do restOfRow

              ! Well, the loop finished, we may or may not have
              ! added something. If so great...if not, oh well. We
              ! just keep going to the next proc. That's the greedy
              ! algorithm for you. 

              ! Processor can be incremented
              cumProcCosts(iProc+1) = cumProcCosts(iProc) + procCosts(iProc)
              iProc = iProc + 1
           end if
        end if
     end if

     ! Move 1 in jj, until we reach the end and wrap around.
     if (increment) then 
        jj = jj + 1

        ! Switch to the next row:
        if (jj == rowPtr(curRow+1)) then 
           
           ! This is really tricky...we know we're at the end of the
           ! row, but we have to SKIP OVER THE EMPTY rows, or else the
           ! algorithm will crap out. Keep incrementing the curRow
           ! until we get a row with something in it. Make sure we
           ! don't go out the end, so check again nRow

           findNextNonZeroRow: do while(jj == rowPtr(curRow+1))
              curRow = curRow + 1
              if (curRow > overlap%nRow) then 
                 exit findNextNonZeroRow
              end if
           end do findNextNonZeroRow
        end if
     end if
  end do masterLoop

end subroutine oversetLoadBalance

