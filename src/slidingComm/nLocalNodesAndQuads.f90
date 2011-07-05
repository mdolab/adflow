!
!      ******************************************************************
!      *                                                                *
!      * File:          nLocalNodesAndQuads.F90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-21-2005                                      *
!      * Last modified: 02-10-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine nLocalNodesAndQuads(nMySubfaces, mySubfaces, &
                                      gridType, nQuad, nNode)
!
!      ******************************************************************
!      *                                                                *
!      * nLocalNodesAndQuads determines the number of local nodes       *
!      * and quadrilaters present in the surface grid to be given to    *
!      * the ADT. The information stored in mySubfaces needs to be      *
!      * modified, because in mySubfaces the halo faces are taken into  *
!      * account and it is possible that some faces and thus nodes must *
!      * be duplicated, because of crossing the line theta = pi and     *
!      * radial singularities.                                          *
!      * For gridType == 1 the connectivity of the primary grid must be *
!      * constructed, while otherwise the dual grid connectivity has to *
!      * be created.                                                    *
!      *                                                                *
!      ******************************************************************
!
       use adtAPI
       use localSubfacesMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nMySubfaces, gridType
       type(localSubfaceType), dimension(nMySubfaces), &
                                              intent(in) :: mySubfaces

       integer(kind=adtIntType), intent(out) :: nQuad, nNode
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, ii, jj
       integer(kind=intType) :: nFace, nPoint

       integer(kind=intType), dimension(:), allocatable :: nodeFlag

       integer(kind=intType), dimension(:,:), pointer :: connFace
       integer(kind=porType), dimension(:),   pointer :: statusFace

       logical, dimension(:), pointer :: storeFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of subfaces of the interface stored on
       ! this processor to determine the local number of nodes and
       ! quads of this interface.

       nQuad = 0
       nNode = 0

       loopSubfaces: do nn=1,nMySubfaces

         ! Set the variables depending on the type of grid that has to
         ! be constructed.

         select case (gridType)

           case (1_intType)

             ! Primary grid must be constructed.

             nFace  = mySubfaces(nn)%nQuad
             nPoint = mySubfaces(nn)%nNode

             storeFace  => mySubfaces(nn)%storeQuad
             statusFace => mySubfaces(nn)%statusQuad
             connFace   => mySubfaces(nn)%connQuad

           !=============================================================

           case default

             ! Dual grid must be constructed.

             nFace  = mySubfaces(nn)%nDual
             nPoint = mySubfaces(nn)%nQuad

             storeFace  => mySubfaces(nn)%storeDual
             statusFace => mySubfaces(nn)%statusDual
             connFace   => mySubfaces(nn)%connDual

         end select

         ! Allocate the memory for nodeFlag, which indicates whether or
         ! not a node belongs to an owned dual face or not. Initialize
         ! nodeFlag to 0 to indicate that it does not belong to an
         ! owned face; 1 means that it does belong to such a face and a
         ! value greater than 1 means it should be stored several times,
         ! because of crossing the line theta == pi or because of a
         ! polar singularity.

         allocate(nodeFlag(nPoint), stat=ierr)
         if(ierr /= 0)                           &
           call terminate("nLocalNodesAndQuads", &
                          "Memory allocation failure for nodeFlag.")

         nodeFlag = 0

         ! Loop over the number quadrilaterals of this face.

         quadLoop: do mm=1,nFace

           ! Check if the face must be stored in the surface mesh.

           if( storeFace(mm) ) then

             ! Face must be stored. Check what kind of face we are
             ! dealing with here.

             select case (statusFace(mm))

               case (normalQuad)

                 ! Normal quadrilateral. It must be stored only once.

                 ii = 1
                 jj = 1

               case (piCrossed)

                 ! Quadrilateral intersects the line theta == pi.
                 ! It must be stored twice.

                 ii = 2
                 jj = 2

               case (polarSingular)

                 ! Quadrilateral contains a polar singularity. The face
                 ! is stored 5 times and the nodes 3 times.

                 ii = 5
                 jj = 3

             end select

             ! Update nQuad and set the nodeFlag of its 4 nodes
             ! to the maximum of jj and the currently stored value.

             nQuad = nQuad + ii

             ii = connFace(mm,1); nodeFlag(ii) = max(nodeFlag(ii),jj)
             ii = connFace(mm,2); nodeFlag(ii) = max(nodeFlag(ii),jj)
             ii = connFace(mm,3); nodeFlag(ii) = max(nodeFlag(ii),jj)
             ii = connFace(mm,4); nodeFlag(ii) = max(nodeFlag(ii),jj)

           endif
         enddo quadLoop

         ! Loop over the number of points of this subface
         ! and update nNode.

         do mm=1,nPoint
           nNode = nNode + nodeFlag(mm)
         enddo

         ! Release the memory of nodeFlag.

         deallocate(nodeFlag, stat=ierr)
         if(ierr /= 0)                           &
           call terminate("nLocalNodesAndQuads", &
                          "Deallocation error for nodeFlag.")

       enddo loopSubfaces

       end subroutine nLocalNodesAndQuads
