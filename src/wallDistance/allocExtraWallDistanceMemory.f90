!
!      ******************************************************************
!      *                                                                *
!      * File:          allocExtraWallDistanceMemory.f90                *
!      * Author:        Gaetan K.W. Kenway                              *
!      * Starting date: 04-12-2012                                      *
!      * Last modified: 04-12-2012                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine allocExtraWallDistanceMemory(level, sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subrouine allocates the extra memory required for the fast*
  !      * wall distance calculations. 
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use viscSurface
  implicit none

  !      Subroutine arguments.
  integer(kind=intType) :: level, sps

  !      Local variables.
  integer(kind=intType) :: nn, ierr, size

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  ! Check if elemID and uv are allocated, and allocate if not
  do nn=1,nDom
     call setPointers(nn,level,sps)
     if (.not. associated(flowDoms(nn,level,sps)%elemID)) then
        allocate(flowDoms(nn,level,sps)%elemID(2:il,2:jl,2:kl),stat=ierr)
        allocate(flowDoms(nn,level,sps)%uv(2,2:il,2:jl,2:kl),stat=ierr)
     end if
  end do
  
  if (allocated(unique_face_info(level,sps)%id)) then
     deallocate(unique_face_info(level,sps)%id)
  end if
  
  size = unique_face_info(level,sps)%n
  allocate(unique_face_info(level,sps)%id(size))
  
end subroutine allocExtraWallDistanceMemory
