subroutine exchangeStatus(level, sps, commPattern, internal)
  !
  !       ExchangeIsCompute exchanges the isCompute flag for the 1 to 1  
  !       connectivity for the given level and sps instance.             
  !
  use blockPointers
  use communication
  use utils, only : setPointers
  use haloExchange, only : whalo1to1intgeneric

  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, sps

  type(commType),          dimension(*), intent(in) :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal
  integer(kind=intType) :: i, j, k, nn

  ! We can't put a pointer directly to status since it is inside of
  ! fringes. So we allocate the first intComm pointer and copy in the
  ! value, exchange and then copy back out. 
  domainLoop:do nn=1, nDom
     call setPointers(nn, level, sps)
     allocate(flowDoms(nn, level, sps)%intCommVars(1)%var(1:ib+1, 1:jb+1, 1:kb+1))
     do k=0, kb
        do j=0, jb
           do i=0, ib
              flowDoms(nn, level, sps)%intCommVars(1)%var(i+1, j+1, k+1) = &
                   flowDoms(nn, level, sps)%fringes(i, j, k)%status
           end do
        end do
     end do
  end do domainLoop
  
 ! Run the generic integer exchange
 call wHalo1to1IntGeneric(1, level, sps, commPattern, internal)
 
  domainLoop2:do nn=1, nDom
     call setPointers(nn, level, sps)
     do k=0, kb
        do j=0, jb
           do i=0, ib
              flowDoms(nn, level, sps)%fringes(i, j, k)%status = & 
                   flowDoms(nn, level, sps)%intCommVars(1)%var(i+1, j+1, k+1)
           end do
        end do
     end do
     deallocate(flowDoms(nn, level, sps)%intCommVars(1)%var)
  end do domainLoop2
end subroutine exchangeStatus
