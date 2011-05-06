!
!      ******************************************************************
!      *                                                                *
!      * File:          setPointers_offTSInstance.f90                   *
!      *                                                                *
!      ******************************************************************
!
subroutine setPointersOffTSInstance(nn,sps,sps2)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * setPointers_offTSInstance calls the normal setPointers but also*
  !      * sets w_offTimeInstance and vol_offTimeInstance which are       *
  !      * required for the forward mode AD calculations                  *
  !      *                                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  implicit none
  !
  !      Subroutine arguments
  !
  integer(kind=intType), intent(in) :: nn,sps,sps2
  
  call setPointers(nn,1,sps)

  w_offTimeInstance  => flowDoms(nn,1,sps2)%w
  vol_offTimeInstance  => flowDoms(nn,1,sps2)%vol

end subroutine setPointersOffTSInstance
