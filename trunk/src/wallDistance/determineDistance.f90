!
!      ******************************************************************
!      *                                                                *
!      * File:          determineDistance.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-03-2003                                      *
!      * Last modified: 02-10-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineDistance(level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * determineDistance determines the distance from the center      *
!      * of the cell to the nearest viscous wall for owned cells.       *
!      *                                                                *
!      ******************************************************************
!
       use adtAPI
       use blockPointers
       use communication
       use inputPhysics
       use section
       use viscSurface
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sps
!
!      Local parameter, which defines the name of the adt to be create
!      and used here. Only needed because of the api of the adt library
!
       character(len=10), parameter :: viscAdt = "ViscousADT"
!
!      Local variables.
!
       integer :: ierr

       integer, dimension(:), allocatable :: procID

       integer(kind=adtIntType) :: nCell, nCellPer, nTria

       integer(kind=adtIntType), dimension(1,1) :: connTria
       real(kind=adtRealType),   dimension(3,2) :: dummy

       integer(kind=adtIntType), dimension(:), allocatable :: elementID

       real(kind=adtRealType), dimension(:),   allocatable :: dist2
       real(kind=adtRealType), dimension(:),   allocatable :: dist2per
       real(kind=adtRealType), dimension(:,:), allocatable :: coor, uvw
       real(kind=adtRealType), dimension(:,:), allocatable :: coorPer

       integer(kind=adtElementType), dimension(:), allocatable :: &
                                                            elementType

       integer(kind=intType) :: nn, mm, ll, ii, jj, i, j, k

       real(kind=realType), dimension(3) :: xc
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Build the adt of the surface grid. As the api requires to
       ! specify both the quadrilateral and the triangular connectivity,
       ! some dummy variables must be passed.

       nTria    = 0
       connTria = 0
       dummy    = zero

       call adtBuildSurfaceADT(nTria,    nquadVisc, nNodeVisc,       &
                               coorVisc, connTria,  connVisc,        &
                               dummy,    .false.,   SUmb_comm_world, &
                               viscAdt)

       ! Determine the number of cell centers for which the distance
       ! must be computed. Also determine how many of them are periodic.

       nCell    = 0
       nCellPer = 0

       do nn=1,nDom
         ll    = flowDoms(nn,level,sps)%nx * flowDoms(nn,level,sps)%ny &
               * flowDoms(nn,level,sps)%nz
         nCell = nCell + ll

         mm = flowDoms(nn,level,sps)%sectionID
         if( sections(mm)%periodic ) nCellPer = nCellPer + ll
       enddo

       ! Allocate the memory for the arrays needed by the ADT.

       allocate(coor(3,nCell), procID(nCell), elementType(nCell), &
                elementID(nCell), uvw(3,nCell), dist2(nCell),      &
                coorPer(3,nCellPer), dist2Per(nCellPer), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("determineDistance", &
                        "Memory allocation failure for the variables &
                        &needed by the adt.")
!
!      ******************************************************************
!      *                                                                *
!      * Step 1: The search of the original coordinates; possibly a     *
!      *         rotational periodic transformation is applied to align *
!      *         the sections.                                          *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the domains to store the coordinates of the cell
       ! centers of the owned cells. Apply the transformation such that
       ! the sections are aligned.

       mm = 0
       domains: do nn=1,nDom

         ! Set the pointers for this block, store the section id a bit
         ! easier in ll and loop over the cell centers.

         call setPointers(nn, level, sps)
         ll = sectionID

         do k=2,kl
           do j=2,jl
             do i=2,il

               ! Compute the coordinates of the cell center relative
               ! to the rotation center of this section.

               xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                     +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                     +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                     +         x(i-1,j,  k,  1) + x(i,j,  k,  1)) &
                     - sections(ll)%rotCenter(1)

               xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                     +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                     +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                     +         x(i-1,j,  k,  2) + x(i,j,  k,  2)) &
                     - sections(ll)%rotCenter(2)

               xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                     +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                     +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                     +         x(i-1,j,  k,  3) + x(i,j,  k,  3)) &
                     - sections(ll)%rotCenter(3)

               ! Apply the periodic transformation for this section to
               ! align it with other sections and store this coordinate
               ! in the appropriate place in coor.

               mm = mm + 1
               coor(1,mm) = rotMatrixSections(ll,1,1)*xc(1) &
                          + rotMatrixSections(ll,1,2)*xc(2) &
                          + rotMatrixSections(ll,1,3)*xc(3) &
                          + sections(ll)%rotCenter(1)

               coor(2,mm) = rotMatrixSections(ll,2,1)*xc(1) &
                          + rotMatrixSections(ll,2,2)*xc(2) &
                          + rotMatrixSections(ll,2,3)*xc(3) &
                          + sections(ll)%rotCenter(2)

               coor(3,mm) = rotMatrixSections(ll,3,1)*xc(1) &
                          + rotMatrixSections(ll,3,2)*xc(2) &
                          + rotMatrixSections(ll,3,3)*xc(3) &
                          + sections(ll)%rotCenter(3)

               ! Initialize the distance squared, because this is an
               ! inout argument in the call to adtMinDistanceSearch.

               dist2(mm) = d2Wall(i,j,k)

             enddo
           enddo
         enddo
       enddo domains

       ! Perform the search. As no no interpolations are required,
       ! some dummies are passed.

       call adtMinDistanceSearch(nCell,  coor,        viscAdt,      &
                                 procID, elementType, elementID,    &
                                 uvw,    dist2,       0_adtIntType, &
                                 dummy,  dummy)
!
!      ******************************************************************
!      *                                                                *
!      * Step 2: For periodic sections the nearest wall may be in the   *
!      *         periodic part of the grid that is not stored.          *
!      *         Therefore apply the periodic transformation to the     *
!      *         node and compute the minimum distance for this         *
!      *         coordinate.                                            *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the counters mm and ii. Mm is the counter for coor;
       ! ii is the counter for coorPer.

       mm = 0
       ii = 0

       ! Loop over the domains and find the periodic ones.

       domainsPer1: do nn=1,nDom
         jj    = flowDoms(nn,level,sps)%nx * flowDoms(nn,level,sps)%ny &
               * flowDoms(nn,level,sps)%nz

         ll = flowDoms(nn,level,sps)%sectionID

         ! Check if the section is periodic.

         if( sections(ll)%periodic ) then

           ! Loop over the corresponding entries in coor of this block
           ! and apply the periodic transformation. The transformed
           ! coordinates are stored in coorPer. Also initialize
           ! the wall distance squared to the value just computed.

           do i=1,jj
             mm = mm + 1
             ii = ii + 1

             xc(1) = coor(1,mm) - sections(ll)%rotCenter(1)
             xc(2) = coor(2,mm) - sections(ll)%rotCenter(2)
             xc(3) = coor(3,mm) - sections(ll)%rotCenter(3)

             coorPer(1,ii) = sections(ll)%rotMatrix(1,1)*xc(1) &
                           + sections(ll)%rotMatrix(1,2)*xc(2) &
                           + sections(ll)%rotMatrix(1,3)*xc(3) &
                           + sections(ll)%rotCenter(1)         &
                           + sections(ll)%translation(1)

             coorPer(2,ii) = sections(ll)%rotMatrix(2,1)*xc(1) &
                           + sections(ll)%rotMatrix(2,2)*xc(2) &
                           + sections(ll)%rotMatrix(2,3)*xc(3) &
                           + sections(ll)%rotCenter(2)         &
                           + sections(ll)%translation(2)

             coorPer(3,ii) = sections(ll)%rotMatrix(3,1)*xc(1) &
                           + sections(ll)%rotMatrix(3,2)*xc(2) &
                           + sections(ll)%rotMatrix(3,3)*xc(3) &
                           + sections(ll)%rotCenter(3)         &
                           + sections(ll)%translation(3)

             dist2Per(ii) = dist2(mm)
           enddo

         else
           ! Section is not periodic. Update the counter mm.

           mm = mm + jj
         endif

       enddo domainsPer1

       ! Perform the adt search of this set of periodic coordinates.
       ! As no no interpolations are required, some dummies are passed.

       call adtMinDistanceSearch(nCellPer, coorPer,     viscAdt,      &
                                 procID,   elementType, elementID,    &
                                 uvw,      dist2Per,    0_adtIntType, &
                                 dummy,    dummy)
!
!      ******************************************************************
!      *                                                                *
!      * Step 3: Also apply the inverse periodic transformation.        *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the counters mm and ii. Mm is the counter for coor;
       ! ii is the counter for coorPer.

       mm = 0
       ii = 0

       ! Loop over the domains and find the periodic ones.

       domainsPer2: do nn=1,nDom
         jj    = flowDoms(nn,level,sps)%nx * flowDoms(nn,level,sps)%ny &
               * flowDoms(nn,level,sps)%nz

         ll = flowDoms(nn,level,sps)%sectionID

         ! Check if the section is periodic.

         if( sections(ll)%periodic ) then

           ! Loop over the corresponding entries in coor of this block
           ! and apply the inverse periodic transformation. Again the
           ! transformed coordinates are stored in coorPer. Note that
           ! the inverse of the rotation matrix is the transpose and
           ! that the translation vector should be multiplied by the
           ! inverse of the rotation matrix. Note that the wall distance
           ! has already been initialized in the previous periodic
           ! search.

           do i=1,jj
             mm = mm + 1
             ii = ii + 1

             xc(1) = coor(1,mm) - sections(ll)%rotCenter(1) &
                   -              sections(ll)%translation(1)
             xc(2) = coor(2,mm) - sections(ll)%rotCenter(2) &
                   -              sections(ll)%translation(2)
             xc(3) = coor(3,mm) - sections(ll)%rotCenter(3) &
                   -              sections(ll)%translation(3)

             coorPer(1,ii) = sections(ll)%rotMatrix(1,1)*xc(1) &
                           + sections(ll)%rotMatrix(2,1)*xc(2) &
                           + sections(ll)%rotMatrix(3,1)*xc(3) &
                           + sections(ll)%rotCenter(1)

             coorPer(2,ii) = sections(ll)%rotMatrix(1,2)*xc(1) &
                           + sections(ll)%rotMatrix(2,2)*xc(2) &
                           + sections(ll)%rotMatrix(3,2)*xc(3) &
                           + sections(ll)%rotCenter(2)

             coorPer(3,ii) = sections(ll)%rotMatrix(1,3)*xc(1) &
                           + sections(ll)%rotMatrix(2,3)*xc(2) &
                           + sections(ll)%rotMatrix(3,3)*xc(3) &
                           + sections(ll)%rotCenter(3)
           enddo

         else
           ! Section is not periodic. Update the counter mm.

           mm = mm + jj
         endif

       enddo domainsPer2

       ! Perform the adt search of this set of periodic coordinates.
       ! As no no interpolations are required, some dummies are passed.

       call adtMinDistanceSearch(nCellPer, coorPer,     viscAdt,      &
                                 procID,   elementType, elementID,    &
                                 uvw,      dist2Per,    0_adtIntType, &
                                 dummy,    dummy)
!
!      ******************************************************************
!      *                                                                *
!      * Step 4: Store the minimum distance in the block type.          *
!      *                                                                *
!      ******************************************************************
!
       mm = 0
       ii = 0

       domainsStore: do nn=1,nDom

         ! Set the pointers for this block and store the section id a
         ! bit easier in ll.

         call setPointers(nn, level, sps)
         ll = sectionID

         ! Check if the section is periodic. If so the distance is set
         ! to the minimum of the value computed in step 1 and the
         ! periodic values. Note that mm should not be updated in this
         ! loop, because it is updated in the loop over the cell centers
         ! of this block. Instead the counter j is used.

         if( sections(ll)%periodic ) then

           jj = nx*ny*nz

           j = mm
           do i=1,jj
             j  = j  + 1
             ii = ii + 1

             dist2(j) = min(dist2(j), dist2Per(ii))
           enddo
         endif

         ! Loop over the cell centers of the block to store the wall
         ! distance. Note that dist2 stores the distance squared and
         ! thus a square root must be taken.
         ! Add a possible offset for the debugging of wall functions.

         do k=2,kl
           do j=2,jl
             do i=2,il
               mm = mm + 1
               d2Wall(i,j,k) = sqrt(dist2(mm)) + wallOffset
             enddo
           enddo
         enddo

       enddo domainsStore

       ! Release the memory of the ADT and the arrays of the module
       ! viscSurface.

       call adtDeallocateADTs(viscAdt)

       deallocate(connVisc, coorVisc, rotMatrixSections, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("determineDistance", &
                         "Deallocation error for the arrays &
                         &of viscSurface")

       ! Release the variables needed by the ADT.

       deallocate(coor, procID, elementType, elementID, uvw, dist2, &
                  coorPer, dist2Per, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("determineDistance", &
                        "Deallocation failure for the variables &
                        &needed by the adt.")

       end subroutine determineDistance
