!
!      ******************************************************************
!      *                                                                *
!      * File:          readLevel0CoolingParameters.f90                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-28-2005                                      *
!      * Last modified: 03-29-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readLevel0CoolingParameters
!
!      ******************************************************************
!      *                                                                *
!      * readLevel0CoolingParameters reads the parameters for the level *
!      * 0 turbine cooling model from the parameter file and stores the *
!      * info in the derived data type for this model.                  *
!      *                                                                *
!      * This model has been developed by Pratt and Whitney and should  *
!      * not be given to third parties.                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use cgnsGrid
       use communication
       use coolingModelLevel0
       use inputIO
       use inputIteration
       implicit none
!
!      Local parameter.
!
       integer, parameter :: readUnit = 32
!
!      Local variables.
!
       integer :: ios, pos, ierr

       integer(kind=intType) :: nn, mm, kk, ll
       integer(kind=intType) :: nLevels, nZones, nSubFaces, zoneID
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
       integer(kind=intType) :: dirDownStream

       integer(kind=intType), dimension(0:nDom) :: nBlocksPerZone
       integer(kind=intType), dimension(nDom)   :: blocksPerZone, tmp
       integer(kind=intType), dimension(nDom)   :: zoneNumbers

       character(len=maxCGNSNameLen) :: zoneName
       character(len=2*maxStringLen) :: errorMessage
       character(len=maxStringLen)   :: keyword, value
       character(len=512)            :: string

       character(len=512), dimension(:), allocatable :: strings

       character(len=maxCGNSNameLen), dimension(nDom) :: zoneNames
!
!      Function definitions
!
       integer(kind=intType) :: bsearchStrings
       integer(kind=intType) :: bsearchIntegers
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine for global blockID's/names for the locally stored    *
!      * blocks and the corresponding ID's. As it is possible that      *
!      * blocks are split also the multiplicity must be taken into      *
!      * account.                                                       *
!      *                                                                *
!      ******************************************************************
!
       select case(fileFormatRead)
         case (cgnsFormat)

           ! CGNS format used. Blocks are given by names.
           ! Copy the global block names into zoneNames, sort them and
           ! remove the possible ties. 

           do nn=1,nDom
             mm = flowDoms(nn,1,1)%cgnsBlockID
             zoneNames(nn) = cgnsDoms(mm)%zoneName
           enddo

           nZones = nDom
           call qsortStrings(zoneNames, nZones)
           call removeTiesStrings(zoneNames, nBlocksPerZone, nZones)

           ! Determine the blocks per zone.

           do nn=1,nZones
             tmp(nn) = nBlocksPerZone(nn-1)
           enddo

           do nn=1,nDom
             mm = flowDoms(nn,1,1)%cgnsBlockID
             mm = bsearchStrings(cgnsDoms(mm)%zoneName, zoneNames, &
                                 nZones)
             tmp(mm) = tmp(mm) + 1
             blocksPerZone(tmp(mm)) = nn
           enddo

         !===============================================================

         case (plot3DFormat)

           ! Plot 3D format. Blocks are given by numbers.
           ! Copy the global block numbers into zoneNumbers, sort them
           ! and remove the possible ties. 

           do nn=1,nDom
             zoneNumbers(nn) = flowDoms(nn,1,1)%cgnsBlockID
           enddo

           nZones = nDom
           call qsortIntegers(zoneNumbers, nZones)
           call removeTiesIntegers(zoneNumbers, nBlocksPerZone, nZones)

           ! Determine the blocks per zone.

           do nn=1,nZones
             tmp(nn) = nBlocksPerZone(nn-1)
           enddo

           do nn=1,nDom
             mm = bsearchIntegers(flowDoms(nn,1,1)%cgnsBlockID, &
                                  zoneNumbers, nZones)
             tmp(mm) = tmp(mm) + 1
             blocksPerZone(tmp(mm)) = nn
           enddo

       end select
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of cooling planes.                        *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the number of cooling planes to 0.

       nPlanesLevel0CoolingModel = 0

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
           call terminate("readLevel0CoolingParameters", errorMessage)
         endif

         call mpi_barrier(SUmb_comm_world, ierr)

       endif

       ! Loop to read the data

       dataLoop: do

         ! Read a string from the file. In case the end of the file
         ! has been reached, exit the loop.

         call getNextInfoLineCooling(string)
         if(ios /= 0) exit dataLoop

         ! Search for the : in the string. If not present, continue
         ! with the next line.

         pos = index(string, ":")
         if(pos == 0) cycle

         ! Create the string keyword get rid of the leading and trailing
         ! spaces.

         keyword = string(:pos-1)
         keyword = trim(keyword)

         ! Check if the line corresponds to the keyword for the
         ! number of cooling planes.

         call convertToLowerCase(keyword)

         if(keyword == "number of level 0 cooling planes") then

           ! Create the string value and read the number it contains.

           value = string(pos+1:)
           read(value,*) nPlanesLevel0CoolingModel

           ! The rest of the information should be present after
           ! this line. Exit the data loop to read this info.

           exit dataLoop

         endif

       enddo dataLoop

       ! Allocate the memory for level0Cooling, the array of the
       ! derived datatype to store the information.

       nLevels = max(nMGLevels, mgStartlevel)
       allocate(level0Cooling(nPlanesLevel0CoolingModel, nLevels), &
                stat=ierr)
       if(ierr /= 0)                                   &
         call terminate("readLevel0CoolingParameters", &
                        "Memory allocation failure for level0Cooling")
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local subfaces for each of the cooling planes.   *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of cooling planes specified.

       coolingLoop: do nn=1,nPlanesLevel0CoolingModel

         ! Find the line with the number of subfaces for this plane.
         ! If not successful, terminate.

         call getNextInfoLineCooling(string)
         if(ios /= 0) then
           if(myID == 0) &
             call terminate("readLevel0CoolingParameters", &
                            "Number of subfaces not found for &
                            &cooling model")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Search for the : in the string and create the string keyword.
         ! If the search was not successful an error message will be
         ! printed later.

         pos = max(index(string, ":"),1)
         keyword = string(:pos-1)
         keyword = trim(keyword)

         ! Check if the line corresponds to the keyword for the
         ! number of subfaces for the cooling plane. If not terminate.

         call convertToLowerCase(keyword)
         if(keyword /= "number of subfaces cooling plane") then
           if(myID == 0) &
             call terminate("readLevel0CoolingParameters", &
                            "Number of subfaces not found for &
                            &cooling model")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Create the string value and read the number it contains.
         ! Allocate the memory for strings accordingly.

         value = string(pos+1:)
         read(value,*) nSubfaces

         allocate(strings(nSubfaces), stat=ierr)
         if(ierr /= 0)                                   &
           call terminate("readLevel0CoolingParameters", &
                          "Memory allocation failure for strings")

         ! Find the next line with info, which should contain the
         ! mass flow, pressure and temperature loss parameters.

         call getNextInfoLineCooling(string)
         if(ios /= 0) then
           if(myID == 0) &
             call terminate("readLevel0CoolingParameters", &
                            "Parameters not found for cooling model")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Read the parameters.

         read(string,*) level0Cooling(nn,1)%mDotRatio, &
                        level0Cooling(nn,1)%dpLog,     &
                        level0Cooling(nn,1)%dTLog

         ! Loop over the number of subfaces and store the parameters
         ! in strings. They will be analyzed later. 

         do mm=1,nSubfaces

           ! Find the next line with info and check if it was found.

           call getNextInfoLineCooling(strings(mm))
           if(ios /= 0) then
             if(myID == 0) &
               call terminate("readLevel0CoolingParameters", &
                              "Parameters not found for cooling model")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

         enddo

         ! Loop over the strings just read and determine the local number
         ! of subfaces belonging to this cooling plane. Make a
         ! distinction between the two grid formats.

         kk = 0
         select case(fileFormatRead)
           case (cgnsFormat)

             ! CGNS format. The block should be identified by a string.

             subfaceLoop1CGNS: do mm=1,nSubfaces

               read(strings(mm),*,iostat=ios) zoneName, iBeg, iEnd,   &
                                              jBeg, jEnd, kBeg, kEnd, &
                                              dirDownStream
               if(ios /= 0) then
                 if(myID == 0)                                   &
                   call terminate("readLevel0CoolingParameters", &
                                  "Something wrong when reading &
                                  &subface parameters.")
                 call mpi_barrier(SUmb_comm_world, ierr)
               endif

               ! Search for the name in the sorted global names
               ! corresponding to the blocks stored on this processor.
               ! Update the number of local subfaces.

               ll = bsearchStrings(zoneName, zoneNames, nZones)
               kk = kk + nLocalSubfacesCoolingSubface(ll)

             enddo subfaceLoop1CGNS

           !=============================================================

           case (plot3DFormat)

             ! Plot 3D format.
             ! The block should be identified by a number.

             subfaceLoop1Plot3D: do mm=1,nSubfaces

               read(strings(mm),*,iostat=ios) zoneID, iBeg, iEnd,     &
                                              jBeg, jEnd, kBeg, kEnd, &
                                              dirDownStream
               if(ios /= 0) then
                 if(myID == 0)                                   &
                   call terminate("readLevel0CoolingParameters", &
                                  "Something wrong when reading &
                                  &subface parameters.")
                 call mpi_barrier(SUmb_comm_world, ierr)
               endif

               ! Search for the zone ID in the sorted zone numbers
               ! corresponding to the blocks stored on this processor.
               ! Update the number of local subfaces.

               ll = bsearchIntegers(zoneID, zoneNumbers, nZones)
               kk = kk + nLocalSubfacesCoolingSubface(ll)

             enddo subfaceLoop1Plot3D

         end select

         ! Allocate the memory to store the local subfaces of this
         ! cooling plane.

         level0Cooling(nn,1)%nSubfaces = kk

         allocate(level0Cooling(nn,1)%blockID(kk),  &
                  level0Cooling(nn,1)%indexDir(kk), &
                  level0Cooling(nn,1)%indSol(kk),   &
                  level0Cooling(nn,1)%indNorm(kk),  &
                  level0Cooling(nn,1)%indX1(kk),    &
                  level0Cooling(nn,1)%indX2(kk),    &
                  level0Cooling(nn,1)%jcBeg(kk),    &
                  level0Cooling(nn,1)%jcEnd(kk),    &
                  level0Cooling(nn,1)%icBeg(kk),    &
                  level0Cooling(nn,1)%icEnd(kk),    &
                  stat=ierr)
         if(ierr /= 0)                                   &
           call terminate("readLevel0CoolingParameters", &
                          "Memory allocation failure for subface data")

         ! Repeat the loop over the number of global subfaces for the
         ! cooling plane, but now store the local subface information.

         kk = 0
         select case(fileFormatRead)
           case (cgnsFormat)

             ! CGNS format. The block is identified by a string. This
             ! has already been checked, so no need to do again.

             subfaceLoop2CGNS: do mm=1,nSubfaces

               read(strings(mm),*) zoneName, iBeg, iEnd, jBeg, jEnd, &
                                             kBeg, kEnd, dirDownStream

               ! Search for the name in the sorted global names
               ! corresponding to the blocks stored on this processor.
               ! Store the local subfaces.

               ll = bsearchStrings(zoneName, zoneNames, nZones)
               call localSubfacesCoolingSubface(ll, kk, &
                                                level0Cooling(nn,1))

             enddo subfaceLoop2CGNS

           !=============================================================

           case (plot3DFormat)

             ! Plot 3D format. The block is identified by a number. This
             ! has already been checked, so no need to do again.

             subfaceLoop2Plot3D: do mm=1,nSubfaces

               read(strings(mm),*) zoneID, iBeg, iEnd, jBeg, jEnd, &
                                           kBeg, kEnd, dirDownStream

               ! Search for the zone ID in the sorted zone numbers
               ! corresponding to the blocks stored on this processor.
               ! Store the local subfaces.

               ll = bsearchIntegers(zoneID, zoneNumbers, nZones)
               call localSubfacesCoolingSubface(ll, kk, &
                                                level0Cooling(nn,1))

             enddo subfaceLoop2Plot3D

         end select

         ! Release the memory of strings.

         deallocate(strings, stat=ierr)
         if(ierr /= 0)                                   &
           call terminate("readLevel0CoolingParameters", &
                          "Deallocation failure for strings")

       enddo coolingLoop

       ! Close the parameter file.

       close(unit=readUnit)

       !=================================================================

       contains

       !=================================================================

         subroutine getNextInfoLineCooling(stringBuf)
!
!        ****************************************************************
!        *                                                              *
!        * getNextInfoLineCooling determines the next line with         *
!        * with information.                                            *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         character(len=512), intent(out) :: stringBuf
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
  99     continue

         ! Try to read the next line. Return if it fails.

         read(unit=readUnit, fmt="(a512)", iostat=ios) stringBuf
         if(ios /= 0) return

         ! Replace all the tab and return characters by spaces and
         ! get rid of the leading and trailing spaces in stringBuf.

         call replaceTabsAndReturns(stringBuf)
         stringBuf = adjustl(stringBuf)
         stringBuf = trim(stringBuf)

         ! In case this is an empty stringBuf or if the first character
         ! is a comment sign, try the next line

         if(len_trim(stringBuf) == 0) goto 99
         if(stringBuf(:1) == "#")     goto 99

         ! Find a possible comment sign somewhere in the stringBuf.
         ! If present the info following the comment sign is ignored.

         pos = index(stringBuf, "#")
         if(pos > 0) then
           stringBuf = stringBuf(:pos-1)
           stringBuf = trim(stringBuf)
         endif

         end subroutine getNextInfoLineCooling

         !===============================================================

         function nLocalSubfacesCoolingSubface(indZone)
!
!        ****************************************************************
!        *                                                              *
!        * nLocalSubfacesCoolingSubface determines the number of local  *
!        * subfaces for the currently global subface. This subface      *
!        * belongs to a global block ID, which corresponds to indZone.  *
!        * If this index is 0 it means that none of the local blocks is *
!        * a subblock of the global block ID and consequently the       *
!        * cooling subface will have 0 local subfaces.                  *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function type.
!
         integer(kind=intType) :: nLocalSubfacesCoolingSubface
!
!        Function arguments.
!
         integer(kind=intType), intent(in) :: indZone
!
!        Local variables.
!
         integer(kind=intType) :: i, j
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Initialize the number of local subfaces to 0.

         nLocalSubfacesCoolingSubface = 0

         ! Check if local blocks are present, which are a subblock of
         ! the global block currently investiged. This is the case if
         ! indZone > 0.

         blocksPresent: if(indZone > 0) then

           ! Subblocks are present. Correct the downstream direction,
           ! if needed.

           if(dirDownStream > 0) then
             dirDownStream  = 1
           else
             dirDownStream = -1
           endif

           ! Loop over the local subblocks to check how many have a
           ! local subface, which is part of the cooling face.

           do i=(nBlocksPerZone(indZone-1)+1), nBlocksPerZone(indZone)

             ! Store the local block ID a bit easier and check if a
             ! local subface is present.
             ! If so, update nLocalSubfacesCoolingSubface.

             j = blocksPerZone(i)

             if(iBeg == iEnd) then  ! I-face

               if(flowDoms(j,1,1)%iBegor <= (iEnd+dirDownStream) .and. &
                  flowDoms(j,1,1)%iEndor >= (iBeg+dirDownStream) .and. &
                  flowDoms(j,1,1)%jBegor < jEnd                  .and. &
                  flowDoms(j,1,1)%jEndor > jBeg                  .and. &
                  flowDoms(j,1,1)%kBegor < kEnd                  .and. &
                  flowDoms(j,1,1)%kEndor > kBeg)                       &
                 nLocalSubfacesCoolingSubface =                       &
                   nLocalSubfacesCoolingSubface + 1

             else if(jBeg == jEnd) then ! J-face

               if(flowDoms(j,1,1)%jBegor <= (jEnd+dirDownStream) .and. &
                  flowDoms(j,1,1)%jEndor >= (jBeg+dirDownStream) .and. &
                  flowDoms(j,1,1)%iBegor < iEnd                  .and. &
                  flowDoms(j,1,1)%iEndor > iBeg                  .and. &
                  flowDoms(j,1,1)%kBegor < kEnd                  .and. &
                  flowDoms(j,1,1)%kEndor > kBeg)                       &
                 nLocalSubfacesCoolingSubface =                       &
                   nLocalSubfacesCoolingSubface + 1

             else if(kBeg == kEnd) then ! K-face

               if(flowDoms(j,1,1)%kBegor <= (kEnd+dirDownStream) .and. &
                  flowDoms(j,1,1)%kEndor >= (kBeg+dirDownStream) .and. &
                  flowDoms(j,1,1)%iBegor < iEnd                  .and. &
                  flowDoms(j,1,1)%iEndor > iBeg                  .and. &
                  flowDoms(j,1,1)%jBegor < jEnd                  .and. &
                  flowDoms(j,1,1)%jEndor > jBeg)                       &
                 nLocalSubfacesCoolingSubface =                       &
                   nLocalSubfacesCoolingSubface + 1

             else

               ! No constant index found for subface. Processor 0
               ! prints an error message and the code terminates.

               if(myID == 0)                                    &
                 call terminate("nLocalSubfacesCoolingSubface", &
                                "No constant index found for cooling &
                                &subface.")
                 call mpi_barrier(SUmb_comm_world, ierr)

             endif

           enddo

         endif blocksPresent

         end function nLocalSubfacesCoolingSubface

         !===============================================================

         subroutine localSubfacesCoolingSubface(indZone, nSub, &
                                                coolingPlane)
!
!        ****************************************************************
!        *                                                              *
!        * localSubfacesCoolingSubface stores the local of subfaces,    *
!        * which are part of the global subface of the currently active *
!        * cooling plane.                                               *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) ::    indZone
         integer(kind=intType), intent(inout) :: nSub

         type(level0CoolingType), intent(inout) :: coolingPlane
!
!        Local variables.
!
         integer(kind=intType) :: i, j, jj, cellOffset
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Check if local blocks are present, which are a subblock of
         ! the global block currently investiged. This is the case if
         ! indZone > 0.

         blocksPresent: if(indZone > 0) then

           ! Subblocks are present. Correct the downstream direction,
           ! if needed and set the cellOffset.

           if(dirDownStream > 0) then
             dirDownStream  = 1
             cellOffset     = 1
           else
             dirDownStream = -1
             cellOffset     = 0
           endif

           ! Loop over the local subblocks to check how many have a
           ! local subface, which is part of the cooling face.

           do i=(nBlocksPerZone(indZone-1)+1), nBlocksPerZone(indZone)

             ! Store the local block ID a bit easier and check if a
             ! local subface is present. If so, store the data.

             j = blocksPerZone(i)

             if(iBeg == iEnd) then  ! I-face

               if(flowDoms(j,1,1)%iBegor <= (iEnd+dirDownStream) .and. &
                  flowDoms(j,1,1)%iEndor >= (iBeg+dirDownStream) .and. &
                  flowDoms(j,1,1)%jBegor < jEnd                  .and. &
                  flowDoms(j,1,1)%jEndor > jBeg                  .and. &
                  flowDoms(j,1,1)%kBegor < kEnd                  .and. &
                  flowDoms(j,1,1)%kEndor > kBeg) then

                 ! Overlapping I-face. Store all the info.

                 nSub = nSub + 1

                 coolingPlane%blockID(nSub)  = j
                 coolingPlane%indexDir(nSub) = iMin
                 coolingPlane%indNorm(nSub)  = iBeg + 1                   & 
                                             - flowDoms(j,1,1)%iBegor
                 coolingPlane%indSol(nSub)   = coolingPlane%indNorm(nSub) &
                                             + cellOffset
                 coolingPlane%indX1(nSub)    = coolingPlane%indNorm(nSub)
                 coolingPlane%indX2(nSub)    = coolingPlane%indX1(nSub)   &
                                             + dirDownStream

                 ! The face range on the plane.

                 jj = jBeg - flowDoms(j,1,1)%jBegor + 1
                 coolingPlane%icBeg(nSub) = max(jj, 1_intType) + 1

                 jj = jEnd - flowDoms(j,1,1)%jBegor + 1
                 coolingPlane%icEnd(nSub) = min(jj, flowDoms(j,1,1)%jl)

                 jj = kBeg - flowDoms(j,1,1)%kBegor + 1
                 coolingPlane%jcBeg(nSub) = max(jj, 1_intType) + 1

                 jj = kEnd - flowDoms(j,1,1)%kBegor + 1
                 coolingPlane%jcEnd(nSub) = min(jj, flowDoms(j,1,1)%kl)

               endif

             else if(jBeg == jEnd) then ! J-face

               if(flowDoms(j,1,1)%jBegor <= (jEnd+dirDownStream) .and. &
                  flowDoms(j,1,1)%jEndor >= (jBeg+dirDownStream) .and. &
                  flowDoms(j,1,1)%iBegor < iEnd                  .and. &
                  flowDoms(j,1,1)%iEndor > iBeg                  .and. &
                  flowDoms(j,1,1)%kBegor < kEnd                  .and. &
                  flowDoms(j,1,1)%kEndor > kBeg) then

                 ! Overlapping J-face. Store all the info.

                 nSub = nSub + 1

                 coolingPlane%blockID(nSub)  = j
                 coolingPlane%indexDir(nSub) = jMin
                 coolingPlane%indNorm(nSub)  = jBeg + 1                   & 
                                             - flowDoms(j,1,1)%jBegor
                 coolingPlane%indSol(nSub)   = coolingPlane%indNorm(nSub) &
                                             + cellOffset
                 coolingPlane%indX1(nSub)    = coolingPlane%indNorm(nSub)
                 coolingPlane%indX2(nSub)    = coolingPlane%indX1(nSub)   &
                                             + dirDownStream

                 ! The face range on the plane.

                 jj = iBeg - flowDoms(j,1,1)%iBegor + 1
                 coolingPlane%icBeg(nSub) = max(jj, 1_intType) + 1

                 jj = iEnd - flowDoms(j,1,1)%iBegor + 1
                 coolingPlane%icEnd(nSub) = min(jj, flowDoms(j,1,1)%il)

                 jj = kBeg - flowDoms(j,1,1)%kBegor + 1
                 coolingPlane%jcBeg(nSub) = max(jj, 1_intType) + 1

                 jj = kEnd - flowDoms(j,1,1)%kBegor + 1
                 coolingPlane%jcEnd(nSub) = min(jj, flowDoms(j,1,1)%kl)

               endif

             else   ! K-face; no need to test, because this has
                    ! already been done.

               if(flowDoms(j,1,1)%kBegor <= (kEnd+dirDownStream) .and. &
                  flowDoms(j,1,1)%kEndor >= (kBeg+dirDownStream) .and. &
                  flowDoms(j,1,1)%iBegor < iEnd                  .and. &
                  flowDoms(j,1,1)%iEndor > iBeg                  .and. &
                  flowDoms(j,1,1)%jBegor < jEnd                  .and. &
                  flowDoms(j,1,1)%jEndor > jBeg) then

                 ! Overlapping K-face. Store all the info.

                 nSub = nSub + 1

                 coolingPlane%blockID(nSub)  = j
                 coolingPlane%indexDir(nSub) = kMin
                 coolingPlane%indNorm(nSub)  = kBeg + 1                   & 
                                             - flowDoms(j,1,1)%kBegor
                 coolingPlane%indSol(nSub)   = coolingPlane%indNorm(nSub) &
                                             + cellOffset
                 coolingPlane%indX1(nSub)    = coolingPlane%indNorm(nSub)
                 coolingPlane%indX2(nSub)    = coolingPlane%indX1(nSub)   &
                                             + dirDownStream

                 ! The face range on the plane.

                 jj = iBeg - flowDoms(j,1,1)%iBegor + 1
                 coolingPlane%icBeg(nSub) = max(jj, 1_intType) + 1

                 jj = iEnd - flowDoms(j,1,1)%iBegor + 1
                 coolingPlane%icEnd(nSub) = min(jj, flowDoms(j,1,1)%il)

                 jj = jBeg - flowDoms(j,1,1)%jBegor + 1
                 coolingPlane%jcBeg(nSub) = max(jj, 1_intType) + 1

                 jj = jEnd - flowDoms(j,1,1)%jBegor + 1
                 coolingPlane%jcEnd(nSub) = min(jj, flowDoms(j,1,1)%jl)

               endif

             endif

           enddo

         endif blocksPresent

         end subroutine localSubfacesCoolingSubface

       end subroutine readLevel0CoolingParameters
