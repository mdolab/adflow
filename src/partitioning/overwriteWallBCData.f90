!
!       File:          overwriteFamilyData.f90                         
!       Author:        Edwin van der Weide                             
!       Starting date: 02-11-2004                                      
!       Last modified: 10-29-2007                                      
!
       subroutine overwriteWallBCData(nZone,cgnsNBocosP)
!
!       OverwriteFamilyData overwrites family data specified in the    
!       cgns file with the data specified in the parameter file. This  
!       can be the rotation information or boundary condition          
!       information of the families. The corresponding values          
!       specified in the grid file, if any, are ignored later on.      
!
       use cgnsGrid
       use communication
       use inputIO
       use killSignals
       implicit none

       integer :: nZone, cgnsNBocosP
!
!      Local parameter.
!
       integer, parameter :: readUnit = 32
!
!      Local variables.
!
       integer :: ios, pos, ierr

       integer(kind=intType) :: nn

       character(len=2*maxStringLen) :: errorMessage
       character(len=maxStringLen)   :: keyword, value
       character(len=512)            :: string
       character(len=15)             :: familyTest


       if (fromPython) return

       ! Open the parameter file for reading. This should normally be
       ! no problem, because the input parameters have already been read
       ! from it. But a check does not really hurt.

       open(unit=readUnit, file=paramFile, status="old", &
            action="read", iostat=ios)

       ! Print an error message if the file could not be opened.
       ! The message is only printed by processor 0, while the others
       ! wait to get killed.

       if(ios /= 0) then

         if(myID == 0) then
           write(errorMessage,*) "Parameter file ", trim(paramFile), &
                                 " not found anymore."
           call returnFail("overwriteWallBCData", errorMessage)
         endif

         call mpi_barrier(SUmb_comm_world, ierr)

       endif

       ! Loop to read the data

       dataLoop: do

         ! Read a string from the file. In case the end of the file
         ! has been reached, exit the loop.

         read(unit=readUnit, fmt="(a512)", iostat=ios) string
         if(ios /= 0) exit dataLoop

         ! Replace all the tab and return characters by spaces.

         call replaceTabsAndReturns(string)

         ! Get rid of the leading and trailing spaces in string.

         string = adjustl(string)
         string = trim(string)

         ! In case this is an empty string or if the first character
         ! is a comment sign, continue with the next line.

         if(len_trim(string) == 0) cycle
         if(string(:1) == "#") cycle

         ! Find a possible comment sign somewhere in the string.
         ! If present the info following the comment sign is ignored.

         pos = index(string, "#")
         if(pos > 0) then
           string = string(:pos-1)
           string = trim(string)
         endif

         ! Search for the : in the string. If not present, continue
         ! with the next line.

         pos = index(string, ":")
         if(pos == 0) cycle

         ! Create the strings keyword and value and get rid of the
         ! leading and trailing spaces of keyword. As this operation has
         ! already been applied for string, only a trim needs to be done.

         keyword = string(:pos-1)
         keyword = trim(keyword)

         value = string(pos+1:)

         ! Extract the substring needed for the family test and
         ! remove that part from the string keyword.

         familyTest = keyword(:14)
         call convertToLowerCase(familyTest)

         keyword = keyword(16:)
         keyword = adjustl(keyword)
         keyword = trim(keyword)

         ! Check whether the first 15 characters of keyword are
         ! "rotating family" or "boundary family". If so call the
         ! corresponding routine to extract that info.

         
         if(familyTest == "wall bc family") &
              call overwriteWallBCInfo(keyword, value,nZone,cgnsNBocosP)

       enddo dataLoop

       ! Close the file

       close(unit=readUnit)

     end subroutine overwriteWallBCData


