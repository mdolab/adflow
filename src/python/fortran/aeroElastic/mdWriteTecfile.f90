!
!      ******************************************************************
!      *                                                                *
!      * File:          mdWriteTecfile.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-27-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdWriteTecfile(famID, tecFile)
!
!      ******************************************************************
!      *                                                                *
!      * mdWriteTecfile writes for the given family an ascii tecplot    *
!      * of the variable currently stored in mdSurfVal.                 *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use constants
       use mdData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: famID
       character(len=*),      intent(in) :: tecFile
!
!      Local parameter.
!
       integer, parameter :: writeUnit = 31
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: startInd,   endInd
       integer(kind=intType) :: startPatch, endPatch
       integer(kind=intType) :: modFamID
       integer(kind=intType) :: mm, nn, ii, jj, kk, i, j, k

       character(len=2*maxStringLen) :: string
       character(len=7) :: int1String, int2String, int3String
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Perform a check to see if this routine is called correctly.
       ! If not, terminate the program.

       if(famID == 0 .and. cgnsNfamilies > 0) then
         if(myID == 0)                      &
           call terminate("mdWriteTecfile", &
                          "Family ID 0 is only allowed when no family &
                          &info is present in the grid")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Return directly if myID /= 0. Only one processor has to do
       ! the writing of the file.

       if(myID /= 0) return

       ! Determine the starting and ending indices in the arrays for
       ! the coordinates and the solution variables as well as the
       ! starting and ending indices for the patches.

       modFamID = max(famID, 1_intType)

       startInd = mdNsurfNodes(0,    modFamID) + 1
       endInd   = mdNsurfNodes(nProc,modFamID)

       startPatch = mdNsurfPatches(0,    modFamID) + 1
       endPatch   = mdNsurfPatches(nProc,modFamID)

       ! Open the file for writing.

       open(unit=writeUnit, file=tecFile, status="unknown", &
            action="write", iostat=ierr)
       if(ierr /= 0) then

         write(string,*) "Tecplot file", trim(tecFile), &
                         " could not be opened for writing"
         call terminate("mdWriteTecfile", string)

       endif

       ! Write the variable definition and initialize the counter
       ! mm to the starting value in the arrays.

       write(writeUnit,100)
       mm = startInd

       ! Loop over the number of patches.

       patchLoop: do nn=startPatch,endPatch

         ! Store the nodal dimensions of the patch a bit easier.

         ii = mdPatchDimensions(1,nn)
         jj = mdPatchDimensions(2,nn)
         kk = mdPatchDimensions(3,nn)

         ! Write the zone title.

         write(int1String,"(i7)") nn
         int1String = adjustl(int1String)
         write(writeUnit,101) trim(int1String)

         ! Write the zone header.

         write(int1String,"(i7)") ii
         write(int2String,"(i7)") jj
         write(int3String,"(i7)") kk

         int1String = adjustl(int1String)
         int2String = adjustl(int2String)
         int3String = adjustl(int3String)

         write(writeUnit,102) trim(int1String), trim(int2String), &
                              trim(int3String)
         write(writeUnit,103)

         ! Loop over the nodes of this patch and write its data.

         do k=1,kk
           do j=1,jj
             do i=1,ii
               write(writeUnit,104) mdSurfXx(1,mm), mdSurfXx(2,mm), &
                                    mdSurfXx(3,mm), mdSurfVal(mm)
               mm = mm + 1
             enddo
           enddo
         enddo

       enddo patchLoop

       ! Close the file.

       close(unit=writeUnit)

       ! Format statements.

 100   format('variables = "x" "y" "z" "varName"')
 101   format('zone T="Zone', a, '"')
 102   format('i=',a, ', j=', a, ', k=', a, ', ZONETYPE=Ordered')
 103   format('datapacking=point')
 104   format(4(1x,e15.8))

       end subroutine mdWriteTecfile
