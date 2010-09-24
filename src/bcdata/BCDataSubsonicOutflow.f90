!
!      ******************************************************************
!      *                                                                *
!      * File:          BCDataSubsonicOutflow.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-24-2004                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine BCDataSubsonicOutflow(boco)
!
!      ******************************************************************
!      *                                                                *
!      * BCDataSubsonicOutflow tries to extract the static pressure     *
!      * for the currently active boundary face, which is a subsonic    *
!      * outflow boundary.                                              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType) :: boco
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j

       real(kind=realType) :: mult, trans

       character(len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the buffer bcVarArray, which is used
       ! for the interpolation and set the cgns names.

       nbcVar = 1
       allocate(bcVarArray(iBeg:iEnd,jBeg:jEnd,nbcVar), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("BCDataSubsonicOutflow", &
                        "Memory allocation failure for bcVarArray")

       bcVarNames(1) = cgnsPressure

       ! Try to determine the static pressure from the data set.

       call extractFromDataSet(BCFaceID(boco))

       ! Write an error message and terminate if it was not
       ! possible to determine the static pressure.

       if(.not. bcVarPresent(1)) then

         write(errorMessage,100)            &
               trim(cgnsDoms(nbkGlobal)%zonename), &
               trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
 100     format("Zone ",a,", boundary subface ",a, &
                ": Static pressure not specified for subsonic outlet")

         call terminate("BCDataSubsonicOutflow", errorMessage)

       endif

       ! Convert to SI-units and store the pressure in ps.

       call siPressure(mass(1), length(1), time(1), mult, trans)

       do j=jBeg,jEnd
         do i=iBeg,iEnd
           BCData(boco)%ps(i,j) = mult*bcVarArray(i,j,1) + trans
         enddo
       enddo

       ! Release the memory of the bcVarArray.

       deallocate(bcVarArray, stat=ierr)
       if(ierr /= 0)                               &
         call terminate("BCDataSubsonicOutflow", &
                        "Deallocation failure for bcVarArray")

       end subroutine BCDataSubsonicOutflow
