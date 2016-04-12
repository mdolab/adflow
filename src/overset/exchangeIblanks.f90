subroutine exchangeIblanks(level, sps, commPattern, internal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * ExchangeIblank exchanges the 1 to 1 internal halo's for the    *
  !      * given level and sps instance.                                  *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use communication
  use inputTimeSPectral
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, sps

  type(commType),          dimension(*), intent(in) :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal

 integer(kind=intType) :: nn

 domainLoop:do nn=1, nDom
    flowDoms(nn, level, sps)%intCommVars(1)%var => &
         flowDoms(nn, level, sps)%iblank(:, :, :)
 end do domainLoop

 ! Run the generic integer exchange
 call wHalo1to1IntGeneric(1, level, sps, commPattern, internal)
 
end subroutine exchangeIblanks
