!
!      ******************************************************************
!      *                                                                *
!      * File:          determineNSlices.F90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-16-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineNSlices(nSlices, slideID, commSlide)
!
!      ******************************************************************
!      *                                                                *
!      * determineNSlices determines the number of periodic slices      *
!      * present in a complete rotation for the given part of the       *
!      * sliding mesh interface slideID.                                *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(inout) :: nSlices
       integer(kind=intType), intent(in)    :: slideID
       integer,               intent(in)    :: commSlide
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nn, mm
!
!      ******************************************************************
!      *                                                                *
!      * Begin executation.                                             *
!      *                                                                *
!      ******************************************************************
!
       ! Find the first local block that participates to the given
       ! slide id.

       domains: do nn=1,nDom
         do mm=1,flowDoms(nn,1,1)%nBocos
           if(flowDoms(nn,1,1)%BCType(mm)   == slidingInterface .and. &
              flowDoms(nn,1,1)%groupNum(mm) == slideID) exit domains
         enddo
       enddo domains

       ! Determine the number of slices if a local domain was found
       ! that participates in this part of the interface.

       if(nn <= nDom) then

         ! Determine the corresponding section and set nSlices to the
         ! number of slices of that section.

         nn = flowDoms(nn,1,1)%sectionID
         nSlices = sections(nn)%nSlices

       endif

       ! Make sure that this value is known on all processors that
       ! participate in this sliding interface.

       nn = nSlices
       call mpi_allreduce(nn, nSlices, 1, sumb_integer, mpi_max, &
                          commSlide, ierr)

       end subroutine determineNSlices
