!
!      ******************************************************************
!      *                                                                *
!      * File:          determineInterfaceIDs.f90                       *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 02-18-2005                                      *
!      * Last modified: 09-19-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineInterfaceIDs
!
!      ******************************************************************
!      *                                                                *
!      * DetermineInterfaceIDs determines more information for both the *
!      * sliding mesh and domain interfaces with other codes, which are *
!      * both specified as user defined boundary conditions in CGNS. In *
!      * particular the number of sliding mesh interfaces and their     *
!      * pairings, and number of interfaces with other codes are        *
!      * determined.                                                    *
!      *                                                                *
!      * There are some parts in the coupler API routines where         *
!      * the family-specified domain interfaces are implicitly assumed. *
!      * Therefore, it is recommended for the time being that all the   *
!      * domain interfaces should be family-specified.                  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use iteration
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, mm, nn
       integer(kind=intType) :: nSlidingFam, nSlidingBC, nSliding
       integer(kind=intType) :: nDomainFam,  nDomainBC

       integer(kind=intType), dimension(:),   allocatable :: famSlidingID
       integer(kind=intType), dimension(:,:), allocatable :: bcSlidingID
       integer(kind=intType), dimension(:),   allocatable :: orID

       character(len=maxStringLen) :: errorMessage

       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                            namesSliding
       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                            namesSorted

       logical :: validInterface
!
!      Function definition.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Count the total number of each type of user-defined BC that
       ! needs to be treated here. Note that if a BC is specified by the
       ! the family it is not counted, such that BCs making up a user-
       ! defined family are only counted once.

       nSlidingFam = 0
       nSlidingBC  = 0
       nDomainFam  = 0
       nDomainBC   = 0

       do nn=1,cgnsNFamilies

         select case (cgnsFamilies(nn)%BCType)
           case (SlidingInterface)
             nSlidingFam = nSlidingFam + 1
           case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                 DomainInterfaceP,   DomainInterfaceRho,    &
                 DomainInterfaceTotal)
             nDomainFam = nDomainFam + 1
         end select

       enddo

       do nn=1,cgnsNDom
         do mm=1,cgnsDoms(nn)%nBocos
           if (cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
               cgnsDoms(nn)%bocoInfo(mm)%familyID == 0) then

             select case (cgnsDoms(nn)%bocoInfo(mm)%BCType)
               case (SlidingInterface)
                 nSlidingBC = nSlidingBC + 1
               case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                     DomainInterfaceP  , DomainInterfaceRho,    &
                     DomainInterfaceTotal)
                 nDomainBC = nDomainBC + 1
             end select

           end if
         end do
       end do

       ! Domain interfaces are only allowed when the code is run in a
       ! multi-disciplinary environment. Check this.

       cgnsNDomainInterfaces = nDomainFam + nDomainBC

       if(standAloneMode .and. cgnsNDomainInterfaces > 0) then
         if(myID == 0) &
           call terminate("determineInterfaceIDs", &
                          "Domain interfaces are not allowed in &
                          &stand alone mode.")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! The number of sliding mesh interfaces must be even, because each
       ! sliding interface should have two sides. Check this.

       nSliding = nSlidingFam + nSlidingBC
       if(mod(nSliding,2) == 1) then
         if(myID == 0)                              &
           call terminate("determineInterfaceIDs", &
                          "Odd number of sliding mesh families found")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Allocate memory to store the names and IDs for the interfaces.

       allocate(famIDsDomainInterfaces(nDomainFam),            &
                bcIDsDomainInterfaces(2,nDomainBC),            &
                namesSliding(nSliding), namesSorted(nSliding), &
                orID(nSliding), famSlidingID(nSlidingFam),     &
                bcSlidingID(2,nSlidingBC), stat=ierr)
       if(ierr /= 0) &
         call terminate("determineInterfaceIDs", &
                        "Memory allocation failure for names, IDs")

       ! Loop back over the families again and store the names and
       ! IDs this time around. Note the ID is just the index of the
       ! corresponding family.

       ii = 0
       jj = 0

       do nn=1,cgnsNFamilies

         select case (cgnsFamilies(nn)%BCType)
           case (SlidingInterface)
             ii = ii + 1
             namesSliding(ii) = cgnsFamilies(nn)%familyName
             famSlidingID(ii) = nn
           case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                 DomainInterfaceP,   DomainInterfaceRho,    &
                 DomainInterfaceTotal)
             jj = jj + 1
             famIDsDomainInterfaces(jj) = nn
         end select

       enddo

       ! Loop back over the boundary conditions again and store the
       ! names and IDs this time around. Note the ID has two parts:
       ! the domain and the index of the BC info in that domain.

       jj = 0

       do nn=1,cgnsNDom
         do mm=1,cgnsDoms(nn)%nBocos
           if (cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
               cgnsDoms(nn)%bocoInfo(mm)%familyID == 0) then

             select case (cgnsDoms(nn)%bocoInfo(mm)%BCType)
               case (SlidingInterface)
                 ii = ii + 1
                 namesSliding(ii) = cgnsDoms(nn)%bocoInfo(mm)%bocoName
                 bcSlidingID(:,ii-nSlidingFam) = (/ nn, mm /)
               case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                     DomainInterfaceP,   DomainInterfaceRho,    &
                     DomainInterfaceTotal)
                 jj = jj + 1
                 bcIDsDomainInterfaces(:,jj) = (/ nn, mm /)
             end select

           end if
         end do
       end do

       ! Initialize orID to -1, which serves as a check later on, and
       ! copy the names of the sliding mesh families in namesSorted.

       do ii=1,nSliding
         orID(ii)        = -1
         namesSorted(ii) = namesSliding(ii)
       enddo

       ! Sort the names of the sliding mesh families in increasing
       ! order and find the corresponding entry in namesSliding.

       call qsortStrings(namesSorted, nSliding)

       do ii=1,nSliding

         ! Search the sorted strings.

         mm = bsearchStrings(namesSliding(ii), namesSorted, nSliding)

         if(orID(mm) /= -1) then

           ! Family name occurs more than once. This is not allowed.

           write(errorMessage,101) trim(namesSliding(ii))
           if(myID == 0) &
             call terminate("determineInterfaceIDs", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)

         endif

         ! Store the entry in orID.

         orID(mm) = ii

       enddo

       ! Set the number of sliding mesh interfaces and allocate the
       ! memory for famIDsSliding.

       cgnsNSliding = nSliding/2
       allocate(famIDsSliding(cgnsNSliding,2), stat=ierr)
       if(ierr /= 0) &
         call terminate("determineInterfaceIDs", &
                        "Memory allocation failure for famIDsSliding")

       ! Check if the sorted family names indeed form a
       ! sliding mesh interface.

       do ii=1,cgnsNSliding

         ! Store the entries in namesSorted, which should form a sliding
         ! interface, a bit easier.

         nn = 2*ii
         mm = nn - 1

         ! Store the original values in famIDsSliding.

         jj = orID(nn)/2

         famIDsSliding(jj,1) = famSlidingID(orID(mm))
         famIDsSliding(jj,2) = famSlidingID(orID(nn))

         ! Check if the names form a valid sliding interface.

         validInterface = .true.
         if(len_trim(namesSorted(nn)) /= len_trim(namesSorted(mm))) &
           validInterface = .false.

         jj = len_trim(namesSorted(nn)) - 2
         if(jj <= 0) then
           validInterface = .false.
         else if(namesSorted(nn)(:jj) /= namesSorted(mm)(:jj)) then
           validInterface = .false.
         endif

         ! Print an error message and exit if the two families do not
         ! form a valid sliding interface.

         if(.not. validInterface) then
           write(errorMessage,102) trim(namesSliding(nn)), &
                                   trim(namesSliding(mm))
           if(myID == 0) &
             call terminate("determineInterfaceIDs", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! The two names form a valid sliding interface. Store the
         ! sliding interface ID for the original family or BC. Note that
         ! one gets a positive ID and the other one a negative ID.
         ! In this way a distinction is made between the two sides.

         if (orID(mm) > nSlidingFam) then
           jj = bcSlidingID(2,orID(mm))
           mm = bcSlidingID(1,orID(mm))
           cgnsDoms(mm)%bocoInfo(jj)%slidingID = -ii
         else
           mm = famSlidingID(orID(mm))
           cgnsFamilies(mm)%slidingID = -ii
         end if

         if (orID(nn) > nSlidingFam) then
           jj = bcSlidingID(2,orID(nn))
           nn = bcSlidingID(1,orID(nn))
           cgnsDoms(nn)%bocoInfo(jj)%slidingID = ii
         else
           nn = famSlidingID(orID(nn))
           cgnsFamilies(nn)%slidingID = ii
         end if

       enddo

       ! Loop over all the boundary conditions again, this time
       ! copying the sliding interface IDs for those BCs that
       ! were part of a sliding family.

       do nn=1,cgnsNDom
         do mm=1,cgnsDoms(nn)%nBocos
           if (cgnsDoms(nn)%bocoInfo(mm)%actualFace .and.   &
               cgnsDoms(nn)%bocoInfo(mm)%familyID > 0 .and. &
               cgnsDoms(nn)%bocoInfo(mm)%BCType == SlidingInterface) then

             cgnsDoms(nn)%bocoInfo(mm)%slidingID = &
             cgnsFamilies(cgnsDoms(nn)%bocoInfo(mm)%familyID)%slidingID

           end if
         end do
       end do

       ! Deallocate the memory used to determine the sliding IDs.

       deallocate(namesSliding, namesSorted, orID, &
                  famSlidingID, bcSlidingID, stat=ierr)
       if(ierr /= 0) &
         call terminate("determineInterfaceIDs", &
                        "Deallocation failure for names, IDs")

       ! Format statements.

 101   format("Family name", 1X,A,1X,"occurs more than once.")
 102   format("Family names", 1X,A,1X,"and",1X,A,1X, &
              "do not form a valid sliding mesh interface")

       end subroutine determineInterfaceIDs
