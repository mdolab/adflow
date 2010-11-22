!
!     ******************************************************************
!     *                                                                *
!     * File:          releaseMemADjoint.f90                           *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 07-20-2006                                      *
!     * Last modified: 02-07-2007                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine releaseMemADjoint(level,sps)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Release all memory allocated for the adjoint solver by         *
  !     * "allocMemADjoint".                                             *
  !     *                                                                *
  !     * This routine was based on:                                     *
  !     * /utils/releaseMemory.f90                                       *
  !     *                                                                *
  !     ******************************************************************
  !
  use BCTypes
  use blockPointers
  implicit none
  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, sps
  !
  !     Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: nn, mm

  logical :: deallocationFailure
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  ! Loop over the number of blocks.

  domainLoop: do nn=1,nDom

     ! Initialize deallocationFailure to .false.

     deallocationFailure = .false.

     ! Deallocate the memory for the global node numbering, 
     ! including halos, if they exist.

     if( associated(flowDoms(nn,level,sps)%globalNode) ) then
        deallocate(flowDoms(nn,level,sps)%globalNode, stat=ierr)
        if(ierr /= 0) deallocationFailure = .true.
     endif


     ! Deallocate the memory for the global Cell numbering, 
     ! including halos, if they exist.

     if( associated(flowDoms(nn,level,sps)%globalCell) ) then
        deallocate(flowDoms(nn,level,sps)%globalCell, stat=ierr)
        if(ierr /= 0) deallocationFailure = .true.
     endif

     ! Deallocate the memory for the adjoint solution, 
     ! including halos, if they exist.

     if( associated(flowDoms(nn,level,sps)%psiAdj) ) then
        deallocate(flowDoms(nn,level,sps)%psiAdj, stat=ierr)
        if(ierr /= 0) deallocationFailure = .true.
     endif


     ! Nullify the pointers, such that no attempt is made to
     ! release the memory again.

     nullify(flowDoms(nn,level,sps)%globalNode)
     nullify(flowDoms(nn,level,sps)%psiAdj)
     nullify(flowDoms(nn,level,sps)%globalCell)

     ! Check for errors in the deallocation.

     if( deallocationFailure ) &
          call terminate("releaseMemADjoint", &
          "Something went wrong when deallocating memory")

  enddo domainLoop

end subroutine releaseMemADjoint
