!
!      ******************************************************************
!      *                                                                *
!      * File:          overwriteFamilyData.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-11-2004                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine overwriteWallBCData(nZone,cgnsNBocosP)
!
!      ******************************************************************
!      *                                                                *
!      * OverwriteFamilyData overwrites family data specified in the    *
!      * cgns file with the data specified in the parameter file. This  *
!      * can be the rotation information or boundary condition          *
!      * information of the families. The corresponding values          *
!      * specified in the grid file, if any, are ignored later on.      *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use inputIO
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

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
           call terminate("overwriteFamilyData", errorMessage)
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


!========================================================================

       subroutine overwriteWallBCInfo(famName, value, nZone,cgnsNBocosP)
!
!      ******************************************************************
!      *                                                                *
!      * overwriteWallBCInfo overwrites the (possible) wall boundary    *
!      * conditions specified in the cgns file                          *
!      *  eran-cbd routine                                              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use su_cgns
       implicit none
!
!      Subroutine arguments.
!
       integer :: nZone, nn,cgnsNBocosP
       character(len=*), intent(inout) :: value
       character(len=*), intent(in)    :: famName

!      Local variables
!
       integer :: ierr, pos, ios

       integer(kind=intType) :: i=0, j ,nWallBc, nFamLocated,mm
       integer(kind=intType),dimension(:),allocatable ::  indFamLoc
       
       character(len=512)            :: string


!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

       ! Search the name in sortedFamName and check if it was found.

!
! Test for the possibility to overwrite the contribution
           ! to the forces on a per family basis.

       nWallBc = 0
       nFamLocated = 0
       allocate(indFamLoc(cgnsNBocosP))
       bocos: do nn=1,cgnsNBocosP
         
          select case (cgnsDoms(nZone)%bocoInfo(nn)%BCTypeCGNS)

            case(BCWall,BCWallInviscid,BCWallViscous,&
                 BCWallViscousHeatFlux,BCWallViscousIsothermal)
               nWallBc = nWallBc + 1
               i=cgnsDoms(nZone)%bocoInfo(nn)%idWBC

               if(trim(famName) == trim(WallBCNames(i)))then
                  nFamLocated = nFamLocated + 1
                  indFamLoc(nFamLocated) = nn
               end if
            end select
             
       enddo bocos

       wallBCfound : if(nFamLocated > 0)then


          ! Get rid of the leading and trailing spaces in value and
          ! store it in string

          value  = adjustl(value)
          value  = trim(value)
          string = value

         ! Convert value to lower case and check if a colon is present.

          call convertToLowerCase(value)
          pos = index(value, ":")
          testColon1: if(pos > 0) then
             testContib: if(value(:pos-1) == "contribute to forces") then

             ! Determine whether or not this family contributes to the
             ! forces and moments.

                value = value(pos+1:)
                value  = adjustl(value)

                select case (value)
                case ("yes")
                   do nn=1,nFamLocated
                      mm = indFamLoc(nn)
                      cgnsDoms(nZone)%bocoInfo(mm)%contributeToForce = .true.
                   end do
                case ("no")
                   do nn=1,nFamLocated
                      mm = indFamLoc(nn)
                      cgnsDoms(nZone)%bocoInfo(mm)%contributeToForce  = .false.
                   end do
                case default
                   write(string,"(3a)") "Boundary family, ", &
                        trim(famName),       &
                        ": Contribute to forces should &
                        &be either yes or no in the &
                        &parameter file."
                   if(myID == 0) &
                        call terminate("overwriteBoundaryInfo", string)
                   call mpi_barrier(SUmb_comm_world, ierr)
                end select

             ! Information has been extracted

             endif testContib
          endif testColon1
          
       endif wallBCfound

       deallocate(indFamLoc)
       return
       
   end subroutine overwriteWallBCInfo
