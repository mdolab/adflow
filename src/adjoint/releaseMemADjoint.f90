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
subroutine releaseMemADjoint()
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
  use ADjointPETSc
  use blockPointers
  use inputTimeSpectral
  implicit none
  !
  !

  !
  !     Local variables.
  !
  integer :: ierr, level,sps

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


  level = 1_intType

  domainLoop: do nn=1,nDom

     do sps=1,nTimeIntervalsSpectral
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

     ! Nullify the pointers, such that no attempt is made to
     ! release the memory again.

     nullify(flowDoms(nn,level,sps)%globalNode)
     nullify(flowDoms(nn,level,sps)%globalCell)

     ! Check for errors in the deallocation.

     if( deallocationFailure ) &
          call terminate("releaseMemADjoint", &
          "Something went wrong when deallocating memory")
  end do
enddo domainLoop
  
  ! Destroy the two empty vectors:
  call vecDestroy(fVec1,PETScIerr)
  call vecDestroy(fVec2,PETScIerr)

end subroutine releaseMemADjoint
