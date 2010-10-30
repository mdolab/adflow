!
!      ******************************************************************
!      *                                                                *
!      * File:          interpolateSlide.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-04-2003                                      *
!      * Last modified: 02-10-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       module storeEntity
!
!      ******************************************************************
!      *                                                                *
!      * Local module to store an entity to be interpolated such that   *
!      * the same source code can be used for both the nodes and the    *
!      * faces to be interpolated.                                      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the derived data type to store the data which    *
!      * determines the arrays sent to the adt routines.                *
!      *                                                                *
!      ******************************************************************
!
       type storeEntityType

         ! nEntities:    Number of entities on the subface.
         ! searchEntity: Whether or not the entity must be interpolated.
         ! coor:         The coordinates to be interpolated.

         integer(kind=intType)                        :: nEntities
         logical,             dimension(:),   pointer :: searchEntity
         real(kind=realType), dimension(:,:), pointer :: coor

       end type storeEntityType

       end module storeEntity

       !=================================================================

       subroutine interpolateSlide(nMySubfaces, mySubfaces, nNode,     &
                                   nQuad,       thetapMin,  thetapMax, &
                                   thetanMin,   thetanMax,  conn,      &
                                   coor,        coorInt,    level,     &
                                   sps,         color,      interType, &
                                   nSlices)
!
!      ******************************************************************
!      *                                                                *
!      * interpolateSlide interpolates either the coordinates,          *
!      * interType == 1, or the face centers, interType /= 1, of        *
!      * mySubfaces in the surface mesh given by conn and coor.         *
!      * Possible periodicity is taken into account via the four angles *
!      * thetapMin, thetapMax, thetanMin and thetanMax.                 *
!      *                                                                *
!      ******************************************************************
!
       use adtAPI
       use block
       use interfaceGroups
       use localSubfacesMod
       use storeEntity
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=adtIntType), intent(in) :: nNode, nQuad

       integer(kind=intType),   intent(in) :: nMySubfaces, color
       integer(kind=intType),   intent(in) :: level, sps, nSlices
       integer(kind=intType),   intent(in) :: interType

       real(kind=realType),     intent(in) :: thetapMin, thetapMax
       real(kind=realType),     intent(in) :: thetanMin, thetanMax

       integer(kind=adtIntType), dimension(4,nQuad), intent(in) :: conn
       real(kind=adtRealType),   dimension(3,nNode), intent(in) :: coor

       real(kind=adtRealType), dimension(3,nNode), intent(in) :: coorInt

       type(localSubfaceType), dimension(nMySubfaces), &
                                           intent(inout) :: mySubfaces
!
!      Local parameter, which defines the name of the ADT to be created
!      and used here. Only needed because of the API of the ADT library.
!
       character(len=10), parameter :: slideADT = "SlidingADT"
!
!      Local variables.
!
       integer :: ierr, comm, myID, nProc
       integer, dimension(:), allocatable :: procID

       integer(kind=adtIntType) :: nInter, nTria

       integer(kind=adtIntType), dimension(1,1) :: connTria
       real(kind=adtRealType),   dimension(3,2) :: dummy

       integer(kind=adtIntType), dimension(:), allocatable :: elementID
       real(kind=adtRealType), dimension(:,:), allocatable :: uvw
       real(kind=adtRealType), dimension(:,:), allocatable :: bufInterpol

       integer(kind=intType) :: bb, nn, mm, ii, jj, i, j, k

       integer(kind=intType), dimension(:), allocatable :: subfaceID
       integer(kind=intType), dimension(:), allocatable :: entityID
       integer(kind=intType), dimension(:), allocatable :: rotIndex

       real(kind=realType) :: theta, perAngle, perAngleInv
       real(kind=realType) :: thetaMin, thetaMax
       real(kind=realType) :: scaleInv
       real(kind=realType) :: ax, r, r1, r2

       real(kind=realType), dimension(3) :: xx

       real(kind=realType), dimension(:),   allocatable :: rotAngles

       logical :: lineThetaPiCrossed

       type(storeEntityType), dimension(nMySubfaces) :: tmpSubfaces
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the communicator, myID and the number of processors
       ! of this sliding mesh a bit easier.

       comm  = myInterfaces(color)%commSlide
       myID  = myInterfaces(color)%myIDSlide
       nProc = myInterfaces(color)%nProcSlide

       ! Initialize some variables to do the interpolation, especially
       ! the continuous range of the polar angle.

       call initInterpol

       ! Build the adt of the surface grid. As the API requires to
       ! specify both the quadrilateral and the triangular connectivity,
       ! some dummy variables must be passed.

       ntria = 0
       connTria = 0

       call adtBuildSurfaceADT(nTria, nQuad, nNode,   coor, connTria, &
                               conn,  dummy, .false., comm, slideADT)

       ! Determine whether the coordinates or the face centers must be
       ! interpolated.

       testInterType: select case (interType)

         case (1_intType)
!
!          **************************************************************
!          *                                                            *
!          * Interpolation of the halo coordinates, i.e. the            *
!          * coordinates one layer outside the original block. This is  *
!          * accomplished by interpolating the coordinates at the       *
!          * interface and assuming that the grids are normal to the    *
!          * interface, i.e. the interpolation weights of the surface   *
!          * nodes can be used for the halo nodes.                      *
!          *                                                            *
!          **************************************************************
!
           ! Set the pointers in tmpSubfaces, such that they correspond
           ! to nodal data.

           do nn=1,nMySubfaces
             tmpSubfaces(nn)%nEntities     = mySubfaces(nn)%nNode
             tmpSubfaces(nn)%searchEntity => mySubfaces(nn)%searchNode
             tmpSubfaces(nn)%coor         => mySubfaces(nn)%coorN
           enddo

           ! Determine the interpolation data for the nodes.
           ! The coordinates of the halo nodes are interpolated.

           call getInterpolationData(3_adtIntType, coorInt)

           ! Release the memory of the ADT, because all the searches
           ! have been done.

           call adtDeallocateADTs(slideADT)

           ! Copy the scaled cylindrical coordinates of the halo nodes
           ! into the block derived datatype.

           haloNodeLoop: do nn=1,nInter

             ! Abbreviate the subface and entity ID in jj and ii
             ! respectively.

             jj = subfaceID(nn)
             ii = entityID(nn)

             ! Determine the unscaled cylindrical coordinates.
             ! Note that the rotation angle must now be substracted to
             ! obtain the correct polar angle.

             ax    = myInterfaces(color)%axmin &
                   + bufInterpol(1,nn)*scaleInv
             r     = myInterfaces(color)%rmin  &
                   + bufInterpol(2,nn)*scaleInv
             theta = bufInterpol(3,nn) - rotAngles(nn)

             ! Determine the Cartesian coordinates. As the transformation
             ! from the Cartesian frame to the local cylindrical frame is
             ! orthonormal, the inverse is the transpose. After that the
             ! coordinates of the center of rotation must be added to
             ! obtain the correct coordinates.

             r1 = r*cos(theta)
             r2 = r*sin(theta)

             xx(1) = ax*myInterfaces(color)%rotAxis(1) &
                   + r1*myInterfaces(color)%radVec1(1) &
                   + r2*myInterfaces(color)%radVec2(1) &
                   + myInterfaces(color)%rotCenter(1)
             xx(2) = ax*myInterfaces(color)%rotAxis(2) &
                   + r1*myInterfaces(color)%radVec1(2) &
                   + r2*myInterfaces(color)%radVec2(2) &
                   + myInterfaces(color)%rotCenter(2)
             xx(3) = ax*myInterfaces(color)%rotAxis(3) &
                   + r1*myInterfaces(color)%radVec1(3) &
                   + r2*myInterfaces(color)%radVec2(3) &
                   + myInterfaces(color)%rotCenter(3)

             ! Determine the block ID and the i, j and k indices of the
             ! halo node and store the coordinates in the block derived
             ! datatype.

             bb = mySubfaces(jj)%blockID
             i  = mySubfaces(jj)%indHaloN(ii,1)
             j  = mySubfaces(jj)%indHaloN(ii,2)
             k  = mySubfaces(jj)%indHaloN(ii,3)

             flowDoms(bb,level,sps)%x(i,j,k,1) = xx(1)
             flowDoms(bb,level,sps)%x(i,j,k,2) = xx(2)
             flowDoms(bb,level,sps)%x(i,j,k,3) = xx(3)

           enddo haloNodeLoop                            

         !===============================================================

         case default
!
!          **************************************************************
!          *                                                            *
!          * Interpolation of the face centers.                         *
!          *                                                            *
!          **************************************************************
!
           ! Set the pointers in tmpSubfaces, such that they correspond
           ! to face data.

           do nn=1,nMySubfaces
             tmpSubfaces(nn)%nEntities     = mySubfaces(nn)%nQuad
             tmpSubfaces(nn)%searchEntity => mySubfaces(nn)%searchQuad
             tmpSubfaces(nn)%coor         => mySubfaces(nn)%coorQuad
           enddo

           ! Determine the interpolation data for the faces.
           ! No data needs to be interpolated.

           call getInterpolationData(0_adtIntType, dummy)

           ! Release the memory of the ADT, because all the searches
           ! have been done.

           call adtDeallocateADTs(slideADT)

           ! Loop over the subfaces and allocate the memory to store
           ! the interpolation data.

           do nn=1,nMySubfaces
             mm = mySubfaces(nn)%nQuad

             allocate(mySubfaces(nn)%u(mm), &
                      mySubfaces(nn)%v(mm), &
                      mySubfaces(nn)%nRotations(mm), &
                      mySubfaces(nn)%donorProc(mm), stat=ierr)
             if(ierr /= 0)                        &
               call terminate("interpolateSlide", &
                              "Memory allocation failure for &
                              &face stuff.")
           enddo

           ! Loop over the number of faces interpolated and store the
           ! data in mySubfaces. Note that the absolute value of the
           ! elementID must be taken. The reason is that a minimum
           ! distance search has been performed and it is possible
           ! that, due to round off or a small physical gap, the
           ! distance is not exactly zero. In that case the ADT search
           ! routine returns a negative value to indicate this.

           do ii=1,nInter
             nn = subfaceID(ii)
             mm = entityID(ii)

             mySubfaces(nn)%donorDualQ(mm) = abs(elementID(ii))
             mySubfaces(nn)%donorProc(mm)  = procID(ii)
             mySubfaces(nn)%u(mm)          = uvw(1,ii)
             mySubfaces(nn)%v(mm)          = uvw(2,ii)
             mySubfaces(nn)%nRotations(mm) = rotIndex(ii)
           enddo

       end select testInterType

       ! Deallocate the memory of the arrays allocated in
       ! getInterpolationData.

       deallocate(subfaceID, entityID, rotIndex,  rotAngles,   &
                  uvw,       procID,   elementID, bufInterpol, &
                  stat=ierr)
       if(ierr /= 0)                        &
         call terminate("interpolateSlide", &
                        "Deallocation failure for arrays allocated &
                        &in getInterpolationData")

       ! Synchronize the processors. Just to be sure, because wild
       ! cards have been used in the receives.

       call mpi_barrier(comm, ierr)

       !=================================================================

       contains

         !===============================================================

         subroutine initInterpol
!
!        ****************************************************************
!        *                                                              *
!        * InitInterpol performs a couple of initializations needed     *
!        * for the interpolation.                                       *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        ****************************************************************
!        *                                                              *
!        * Begin executation.                                           *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the periodic angle and its inverse and the scale
         ! factor to transform to the original cylindrical coordinates.

         perAngle    = two*pi/nSlices
         perAngleInv = one/perAngle
         scaleInv    = one/myInterfaces(color)%scale

         ! Determine the range of angles of the surface grid in the
         ! continuous interval [-pi..2*pi]. The range [pi..2*pi] is only
         ! used if the line theta == pi is crossed. The addition and
         ! substraction of 1.e-8 is there to avoid possible numerical
         ! problems.

         lineThetaPiCrossed = .false.

         if(thetanMin > thetanMax) then

           ! All the angles are positive.

           thetaMin = thetapMin - 1.e-8_realType
           thetaMax = thetapMax + 1.e-8_realType

         else if(thetapMin > thetapMax) then

           ! All the angles are negative.

           thetaMin = thetanMin - 1.e-8_realType
           thetaMax = thetanMax + 1.e-8_realType

         else if(thetanMin == -pi) then

           ! The line theta == pi is crossed. Add two*pi to the maximum
           ! of the negative range to obtain a continuous angle interval
           ! between [0..2*pi].

           thetaMin  = thetapMin - 1.e-8_realType
           thetaMax  = thetanMax + two*pi + 1.e-8_realType

           lineThetaPiCrossed = .true.

         else

           ! The line theta == 0 is crossed. The minimum is the negative
           ! minimum and the maximum the positive maximum.

           thetaMin = thetanMin - 1.e-8_realType
           thetaMax = thetapMax + 1.e-8_realType

         endif

         end subroutine initInterpol

         !===============================================================

         subroutine getInterpolationData(nInterpol, donorData)
!
!        ****************************************************************
!        *                                                              *
!        * getInterpolationData determines the interpolation data for   *
!        * the coordinates stored in tmpSubfaces. The data in this      *
!        * derived data type either points to nodal data or face data.  *
!        * Furthermore it interpolates nInterpol, possibly 0, values    *
!        * from the given donorData.                                    *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=adtIntType), intent(in) :: nInterpol

         real(kind=adtRealType), dimension(:,:), intent(in) :: donorData
!
!        Local variables.
!
         integer(kind=intType) :: nn, mm, ii, jj, nLow, nUpp, nPosRot

         integer(kind=intType), dimension(:), allocatable :: nRotations

         integer(kind=adtIntType) :: nSearch

         integer(kind=adtElementType), dimension(:), allocatable :: &
                                                            elementType

         real(kind=realType) :: theta, thetaNew, angleRot, d2

         real(kind=adtRealType), dimension(:),   allocatable :: dist2
         real(kind=adtRealType), dimension(:,:), allocatable :: coorEnt

         logical :: thetaCorrected
!
!        ****************************************************************
!        *                                                              *
!        * Begin executation.                                           *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the number of local nodes to be interpolated.

         nInter = 0
         do nn=1,nMySubfaces
           do mm=1,tmpSubfaces(nn)%nEntities
             if( tmpSubfaces(nn)%searchEntity(mm) ) nInter = nInter + 1
           enddo
         enddo

         ! Allocate the memory for the subface id, entity id and number
         ! of possible rotations. The last array is in cumulative storage
         ! format, such that the cylindrical coordinates can be stored in
         ! a 1D array later on.

         allocate(subfaceID(nInter),    entityID(nInter), &
                  nRotations(0:nInter), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("getInterpolationData", &
                          "Memory allocation failure for subfaceID, &
                          &entityID and nRotations.")

         ! Repeat the loop over the subfaces, but now store the info in
         ! the arrays just allocated. Make a distinction between 1 slice
         ! present, i.e. the entire wheel, or the case when periodicity
         ! is used, i.e. multiple slices.

         nRotations(0) = 0
         nInter = 0

         testOneSlice1: if(nSlices == 1) then

           ! Only 1 slice present, which is the entire wheel. Therefore
           ! the number of possible rotations is only 1 as well.

           do nn=1,nMySubfaces
             do mm=1,tmpSubfaces(nn)%nEntities
               if( tmpSubfaces(nn)%searchEntity(mm) ) then
                 nInter = nInter + 1
                 subfaceID(nInter)  = nn
                 entityID(nInter)   = mm
                 nRotations(nInter) = nInter
               endif
             enddo
           enddo

         else testOneSlice1

           ! Multiple slices to form the entire wheel. It is possible
           ! that multiple rotations over the periodic angle cover the
           ! range of the surface grid. The usual cause is curved
           ! periodic boundaries.

           do nn=1,nMySubfaces
             do mm=1,tmpSubfaces(nn)%nEntities
               if( tmpSubfaces(nn)%searchEntity(mm) ) then
                 nInter = nInter + 1
                 subfaceID(nInter) = nn
                 entityID(nInter)  = mm

                 ! Determine the possible number of rotations, such that
                 ! the angle theta (plus a certain number times the
                 ! periodic angle) is covered by the range of the surface
                 ! grid stored in the adt.

                 theta = tmpSubfaces(nn)%coor(mm,3)

                 ! Add 2*pi to theta if the negative range is shifted
                 ! and if theta was in this range.

                 if(lineThetaPiCrossed .and. theta <= thetanMax) &
                   theta = theta + two*pi

                 ! Determine the number of possibilities.

                 nLow = (thetaMin - theta)*perAngleInv
                 nUpp = (thetaMax - theta)*perAngleInv

                 nPosRot = 0
                 do ii=nLow,nUpp

                   ! Determine the angle obtained by applying ii
                   ! rotations of the periodic angle to theta. If this
                   ! angle is in the range of the surface mesh, update
                   ! the number of possibilities.

                   thetaNew = theta + ii*perAngle
                   if(thetaNew >= thetaMin .and. &
                      thetaNew <= thetaMax) nPosRot = nPosRot + 1
                 enddo

                 ! Store the number of rotations for this node in
                 ! nRotations, which is in cumulative storage format.

                 nRotations(nInter) = nRotations(nInter-1) + nPosRot

               endif
             enddo
           enddo

         endif testOneSlice1

         ! Check in debug mode if at least one possibility was found
         ! for every node.

         if( debug ) then

           do nn=1,nInter
             if(nRotations(nn) == nRotations(nn-1)) exit
           enddo

           if(nn <= nInter)                         &
             call terminate("getInterpolationData", &
                            "No rotation angle found.")
         endif

         ! Allocate the memory for the arrays needed for the
         ! interpolation.

         nn = nRotations(nInter)
         mm = nInterpol

         allocate(coorEnt(3,nn), procID(nn),    elementType(nn),    &
                  elementID(nn), uvw(3,nn),     dist2(nn),          &
                  rotIndex(nn),  rotAngles(nn), bufInterpol(mm,nn), &
                  stat=ierr)
         if(ierr /= 0)                            &
           call terminate("getInterpolationData", &
                          "Memory allocation failure for &
                          &interpolation arrays.")

         ! Again loop over the nodes to be interpolated. Now store the
         ! scaled cylindrical coordinates and the rotation angles for
         ! all the possibilities. And again make a distinction between
         ! one and multiple slices present.

         nInter = 0

         testOneSlice2: if(nSlices == 1) then

           ! One slice present. No rotations need to be applied.
           ! Simply store the scaled cylindrical coordinates and
           ! initialize the wall distance squared to a large value.
           ! The latter is an inout argument to the ADT search routine
           ! and therefore has to be initialized.

           do nn=1,nMySubfaces
             do mm=1,tmpSubfaces(nn)%nEntities
               if( tmpSubfaces(nn)%searchEntity(mm) ) then
                 nInter = nInter + 1
                 coorEnt(1,nInter) = tmpSubfaces(nn)%coor(mm,1)
                 coorEnt(2,nInter) = tmpSubfaces(nn)%coor(mm,2)
                 coorEnt(3,nInter) = tmpSubfaces(nn)%coor(mm,3)

                 rotIndex(nInter)  = 0
                 rotAngles(nInter) = zero
                 dist2(nInter)     = large
               endif
             enddo
           enddo

         else testOneSlice2

           ! Multiple slices present. Determine the possible rotation
           ! angles and corresponding indices.

           do nn=1,nMySubfaces
             do mm=1,tmpSubfaces(nn)%nEntities
               if( tmpSubfaces(nn)%searchEntity(mm) ) then

                 ! Determine the possible rotation indices and angles
                 ! for this point.

                 theta = tmpSubfaces(nn)%coor(mm,3)
                 ii    = nRotations(nInter) + 1

                 ! Add 2*pi to theta if the negative range is shifted
                 ! and if theta was in this range. Set thetaCorrected
                 ! to .True. In this case.

                 thetaCorrected = .false.
                 if(lineThetaPiCrossed .and. theta <= thetanMax) then
                   theta = theta + two*pi
                   thetaCorrected = .true.
                 endif

                 ! Determine the lower and upper end of the rotation
                 ! index.

                 nLow = (thetaMin - theta)*perAngleInv
                 nUpp = (thetaMax - theta)*perAngleInv

                 ! Loop over this range of rotation indices and store
                 ! the ones that are real possibilities.

                 nPosRot = nRotations(nInter)
                 do jj=nLow,nUpp

                   angleRot = jj*perAngle
                   thetaNew = theta + angleRot
                   if(thetaNew >= thetaMin .and. &
                      thetaNew <= thetaMax) then

                     ! Angle is in the range. Update nPosRot and
                     ! store the rotation index jj.

                     nPosRot = nPosRot + 1
                     rotIndex(nPosRot) = jj

                     ! Adapt angleRot if thetaNew is outside the range
                     ! (-pi..Pi). Note that thetaNew >= -pi by
                     ! definition, such that it is only needed to test
                     ! the upper boundary. If theta was corrected (by
                     ! adding 2pi) this must also be done to angleRot
                     ! as the new angle for the search is computed
                     ! relative to the original angle and not to theta.

                     if(thetaNew > pi)  angleRot = angleRot - two*pi
                     if(thetaCorrected) angleRot = angleRot + two*pi

                     rotAngles(nPosRot) = angleRot

                     ! Make sure that zero rotation angle is always the
                     ! first possibility to be checked. In this way no
                     ! rotation is applied in case of ties, which saves
                     ! CPU time during the cell centered exchange.

                     if(jj == 0) then
                       rotIndex(nPosRot) = rotIndex(ii)
                       rotIndex(ii)      = jj

                       rotAngles(nPosRot) = rotAngles(ii)
                       rotAngles(ii)      = angleRot
                     endif
                   endif
                 enddo

                 ! Loop over the possible number of rotations for this
                 ! point and store the cylindrical coordinates. Note the
                 ! offset added to the angle to make sure that the point
                 ! coincides with the surface mesh stored in the adt.
                 ! Also initialize the distance squared to a large value,
                 ! because this variable is an inout argument in the adt
                 ! search routine.

                 do jj=ii,nRotations(nInter+1)
                   coorEnt(1,jj) = tmpSubfaces(nn)%coor(mm,1)
                   coorEnt(2,jj) = tmpSubfaces(nn)%coor(mm,2)
                   coorEnt(3,jj) = tmpSubfaces(nn)%coor(mm,3) &
                                 + rotAngles(jj)

                   dist2(jj) = large
                 enddo

                 ! Update the counter nInter.

                 nInter = nInter + 1

               endif
             enddo
           enddo

         endif testOneSlice2

         ! Determine the elements which minimize the distance
         ! to the points.

         nSearch = nRotations(nInter)
         call adtMinDistanceSearch(nSearch,   coorEnt,     slideADT,  &
                                   procID,    elementType, elementID, &
                                   uvw,       dist2,       nInterpol, &
                                   donorData, bufInterpol)

         ! Determine for the interpolated points which of the possible
         ! rotations should be taken.

         do nn=1,nInter

           ! Determine the index which has the minimal distance.
           ! This value is stored in mm.

           mm = nRotations(nn-1) + 1
           d2 = dist2(mm)

           do ii=(nRotations(nn-1)+2), nRotations(nn)
             if(dist2(ii) < d2) then
               d2 = dist2(ii)
               mm = ii
             endif
           enddo

           ! Copy the relevant info into the position nn.

           procID(nn)    = procID(mm)
           elementID(nn) = elementID(mm)
           uvw(1,nn)     = uvw(1,mm)
           uvw(2,nn)     = uvw(2,mm)
           uvw(3,nn)     = uvw(3,mm)
           rotIndex(nn)  = rotIndex(mm)
           rotAngles(nn) = rotAngles(mm)

           do ii=1,nInterpol
             bufInterpol(ii,nn) = bufInterpol(ii,mm)
           enddo

         enddo

         ! Release the local memory of this internal subroutine.

         deallocate(nRotations, coorEnt, elementType, dist2, stat=ierr)
         if(ierr /= 0)                            &
           call terminate("getInterpolationData", &
                          "Deallocation failure for nRotations, &
                          &coorEnt, elementType and dist2.")

         end subroutine getInterpolationData

       end subroutine interpolateSlide
