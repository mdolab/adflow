!
!      ******************************************************************
!      *                                                                *
!      * File:          BCDataIsothermalWall.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-24-2004                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine BCDataIsothermalWall(boco)
!
!      ******************************************************************
!      *                                                                *
!      * BCDataIsothermalWall tries to extract the wall temperature     *
!      * for the currently active boundary face, which is an isothermal *
!      * viscous wall.                                                  *
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
       if(ierr /= 0)                            &
         call terminate("BCDataIsothermalWall", &
                        "Memory allocation failure for bcVarArray")

       bcVarNames(1) = cgnsTemp

       ! Try to determine the temperature from the data set.

       call extractFromDataSet(BCFaceID(boco))

       ! Write an error message and terminate if it was not
       ! possible to determine the temperature.

       if(.not. bcVarPresent(1)) then

         write(errorMessage,100)                    &
               trim(cgnsDoms(nbkGlobal)%zonename), &
               trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
 100     format("Zone ",a,", boundary subface ",a, &
                ": Wall temperature not specified for isothermal wall")

         call terminate("BCDataIsothermalWall", errorMessage)

       endif

       ! Convert to si-units and store the temperature in TNS_Wall.

       call siTemperature(temp(1), mult, trans)

       do j=jBeg,jEnd
         do i=iBeg,iEnd
           BCData(boco)%TNS_Wall(i,j) = mult*bcVarArray(i,j,1) + trans
         enddo
       enddo

       ! Release the memory of the bcVarArray.

       deallocate(bcVarArray, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("BCDataIsothermalWall", &
                        "Deallocation failure for bcVarArray")

       end subroutine BCDataIsothermalWall
