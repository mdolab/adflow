!
!      ******************************************************************
!      *                                                                *
!      * File:          localNodesAndQuads.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-03-2003                                      *
!      * Last modified: 02-10-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine localNodesAndQuads(nMySubfaces, mySubfaces, &
                                     nQuad,        nNode,     &
                                     gridType,     conn,      &
                                     subface,      quadID,    &
                                     coor,         coorInt,   &
                                     thetapMin,    thetanMin, &
                                     thetapMax,    thetanMax)
!
!      ******************************************************************
!      *                                                                *
!      * localNodesAndQuads stores the owned part of mySubfaces in      *
!      * the given arrays. Crossing of the line theta == pi is taken    *
!      * into account by storing the corresponding face twice. Also     *
!      * polar singularities for faces are taken into account.          *
!      * For gridType == 1 the connectivity of the primary grid is      *
!      * constructed, while otherwise the dual grid connectivity is     *
!      * be created.                                                    *
!      *                                                                *
!      ******************************************************************
!
       use adtAPI
       use constants
       use localSubfacesMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nMySubfaces, gridType

       type(localSubfaceType), dimension(*), intent(inout) :: mySubfaces

       real(kind=realType), intent(out) :: thetapMin, thetanMin
       real(kind=realType), intent(out) :: thetapMax, thetanMax

       integer(kind=adtIntType), intent(in) :: nQuad, nNode

       integer(kind=intType), dimension(nQuad), intent(out) :: subface
       integer(kind=intType), dimension(nQuad), intent(out) :: quadID

       integer(kind=adtIntType), dimension(4,nQuad), intent(out) :: conn
       real(kind=adtRealType),   dimension(3,nNode), intent(out) :: coor

       real(kind=realType), dimension(3,nNode), intent(out) :: coorInt
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, ii, jj, kk, ll
       integer(kind=intType) :: iNode, iQuad
       integer(kind=intType) :: n1, n2, n3, n4, m1, m2, m3, m4
       integer(kind=intType) :: nFace, nPoint

       integer(kind=intType), dimension(:), allocatable :: nodeFlag
       integer(kind=intType), dimension(:), allocatable :: nodeDupl

       integer(kind=intType), dimension(:,:), pointer :: connFace
       integer(kind=porType), dimension(:),   pointer :: statusFace

       real(kind=realType) :: val

       real(kind=realType), dimension(:,:), pointer :: coorPoint
       real(kind=realType), dimension(:,:), pointer :: coorPointInt

       logical :: thetaPiCrossed, thetaZeroCrossed, storeIntCoor

       logical, dimension(:), pointer :: storeFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Some initializations.

       iNode = 0
       iQuad = 0

       thetapMin =  pi
       thetapMax =  zero
       thetanMin =  zero
       thetanMax = -pi

       ! Loop over the subfaces of the sliding interrface to store
       ! surface grid info in the arrays subface, quadID, conn and
       ! coor. Also the angles theta are updated.

       loopSubfaces: do nn=1,nMySubfaces

         ! Set the variables depending on the type of grid that has to
         ! be constructed.

         select case (gridType)

           case (1_intType)

             ! Primary grid must be constructed. The coordinates
             ! one layer into the block must be stored.

             nFace  = mySubfaces(nn)%nQuad
             nPoint = mySubfaces(nn)%nNode

             thetaPiCrossed   = mySubfaces(nn)%piCrossedQuad
             thetaZeroCrossed = mySubfaces(nn)%zeroCrossedQuad

             storeFace    => mySubfaces(nn)%storeQuad
             statusFace   => mySubfaces(nn)%statusQuad
             connFace     => mySubfaces(nn)%connQuad
             coorPoint    => mySubfaces(nn)%coorN
             coorPointInt => mySubfaces(nn)%coorNInt

             storeIntCoor = .true.

           !=============================================================

           case default

             ! Dual grid must be constructed. The coordinates one layer
             ! into the block are not present and therefore they cannot
             ! be stored.

             nFace  = mySubfaces(nn)%nDual
             nPoint = mySubfaces(nn)%nQuad

             thetaPiCrossed   = mySubfaces(nn)%piCrossedDual
             thetaZeroCrossed = mySubfaces(nn)%zeroCrossedDual

             storeFace  => mySubfaces(nn)%storeDual
             statusFace => mySubfaces(nn)%statusDual
             connFace   => mySubfaces(nn)%connDual
             coorPoint  => mySubfaces(nn)%coorQuad

             nullify(coorPointInt)
             storeIntCoor = .false.

         end select

         ! Allocate the memory for nodeFlag and nodeDupl;
         ! initialize both arrays to 0.

         allocate(nodeFlag(nPoint), nodeDupl(nPoint), stat=ierr)
         if(ierr /= 0)                          &
           call terminate("localNodesAndQuads", &
                          "Memory allocation failure for nodeFlag &
                          &and nodeDupl.")
         nodeFlag = 0
         nodeDupl = 0

         ! Loop over the faces to flag the nodes, which must be stored.

         connLoop1: do mm=1,nFace

           ! Check if the face must be stored.

           if( storeFace(mm) ) then

             ! Set nodeFlag of all these nodes to 1, indicating that
             ! it must be stored.

             nodeFlag(connFace(mm,1)) = 1
             nodeFlag(connFace(mm,2)) = 1
             nodeFlag(connFace(mm,3)) = 1
             nodeFlag(connFace(mm,4)) = 1

             ! For a face that intersects the line theta == pi the
             ! nodes must be duplicated.

             if(statusFace(mm) == piCrossed) then
               nodeDupl(connFace(mm,1)) = 1
               nodeDupl(connFace(mm,2)) = 1
               nodeDupl(connFace(mm,3)) = 1
               nodeDupl(connFace(mm,4)) = 1
             endif

             ! For a polar singularity two nodes must be duplicated.

             if(statusFace(mm) == polarSingular) then

               ! Find out the nodes in the second and third quadrant.
               ! These must be duplicated.

               do ii=1,4
                 jj = connFace(mm,ii)

                 if(coorPoint(jj,3) >   half*pi) nodeDupl(jj) = 1
                 if(coorPoint(jj,3) <= -half*pi) nodeDupl(jj) = 1
               enddo

             endif
           endif
         enddo connLoop1

         ! Loop over the points of this subface and store the new
         ! connectivity in nodeFlag. Also copy the coordinates.

         pointLoop1: do mm=1,nPoint
           if(nodeFlag(mm) == 1) then

             ! Node must be stored. Update the counter and store the
             ! node number for the connectivity in nodeFlag.

             iNode = iNode + 1
             nodeFlag(mm) = iNode

             ! Copy the coordinates of the node on the face and the
             ! node 1 layer into the block. The latter is only needed
             ! for the primary grid.

             coor(1,iNode) = coorPoint(mm,1)
             coor(2,iNode) = coorPoint(mm,2)
             coor(3,iNode) = coorPoint(mm,3)

             if( storeIntCoor ) then
               coorInt(1,iNode) = coorPointInt(mm,1)
               coorInt(2,iNode) = coorPointInt(mm,2)
               coorInt(3,iNode) = coorPointInt(mm,3)
             endif

             ! Check the minimum and maximum angle.

             if(coor(3,iNode) >= zero) then
               thetapMax = max(thetapMax, &
                               real(coor(3,iNode),kind=realType))
               thetapMin = min(thetapMin, &
                               real(coor(3,iNode),kind=realType))
             else
               thetanMax = max(thetanMax, &
                               real(coor(3,iNode),kind=realType))
               thetanMin = min(thetanMin, &
                               real(coor(3,iNode),kind=realType))
             endif

           endif
         enddo pointLoop1

         ! Loop again over the points to store the points that must
         ! be duplicated.

         pointLoop2: do mm=1,nPoint
           if(nodeDupl(mm) == 1) then

             ! Node must be duplicated. Update the counter and store
             ! the new node number in nodeDupl(mm). Needed for the
             ! connectivity later on.

             iNode = iNode + 1
             nodeDupl(mm) = iNode

             ! Store the axial and radial coordinate of the node.

             coor(1,iNode) = coorPoint(mm,1)
             coor(2,iNode) = coorPoint(mm,2)

             ! Add or substract 2*pi from the angle, depending on the
             ! situation.

             val = two*pi
             if(coorPoint(mm,3) > zero) val = -two*pi

             coor(3,iNode)    = coorPoint(mm,3)    + val

             ! Store the coordinates one layer into the block, if needed.

             if( storeIntCoor ) then
               coorInt(1,iNode) = coorPointInt(mm,1)
               coorInt(2,iNode) = coorPointInt(mm,2)
               coorInt(3,iNode) = coorPointInt(mm,3) + val
             endif

           endif
         enddo pointLoop2

         ! Correct the angles if the lines theta == pi and/or theta == 0
         ! intersect with the current subface.

         if( thetaPiCrossed ) then
           thetapMax =  pi
           thetanMin = -pi
         endif

         if( thetaZeroCrossed ) then
           thetapMin = zero
           thetanMax = zero
         endif

         ! Loop over the faces to store the info of the owned faces in
         ! conn. Also update the coordinates for a polar singular quad.

         connLoop2: do mm=1,nFace

           ! Test if the face must be stored.

           testStoreFace: if( storeFace(mm) ) then

             ! Face must be stored. Check what kind of face we are
             ! dealing with here.

             select case (statusFace(mm))

               case (normalQuad)

                 ! Normal quad. It should be stored once. Set the counter
                 ! jj to 1 and store the connectivity.

                 jj = 1

                 ii = iQuad + 1
                 conn(1,ii) = nodeFlag(connFace(mm,1))
                 conn(2,ii) = nodeFlag(connFace(mm,2))
                 conn(3,ii) = nodeFlag(connFace(mm,3))
                 conn(4,ii) = nodeFlag(connFace(mm,4))
 
               case (piCrossed)

                 ! Quadrilateral intersects the line theta == pi.
                 ! It must be stored twice. Store the connectivity.

                 ii = iQuad + 1
                 jj = iQuad + 2

                 ! Loop over the 4 nodes of the connectivity.

                 do kk=1,4

                   ! In place ii the quad is stored with the positive
                   ! angles; in jj the one with the negative angles.

                   ll = connFace(mm,kk)
                   if(coorPoint(ll,3) > zero) then
                     conn(kk,ii) = nodeFlag(ll)
                     conn(kk,jj) = nodeDupl(ll)
                   else
                     conn(kk,ii) = nodeDupl(ll)
                     conn(kk,jj) = nodeFlag(ll)
                   endif

                 enddo

                 ! And set the counter jj to 2 to indicate that 2
                 ! faces have been stored.

                 jj = 2

               case (polarSingular)

                 ! Quadrilateral contains a polar singularity. Without
                 ! doing something special interpolation will go wrong
                 ! here. Therefore the quadrilateral in cylindrical
                 ! coordinates is replaced by 5 quadrilaterals. For
                 ! every quadrilateral two nodes are real nodes and
                 ! two other nodes are ghost nodes with the same angle
                 ! but with a negative radius. In this way a kind of
                 ! a cartesian interpolation is recovered on this face.

                 ! First determine the node in the second quadrant;
                 ! store the index in n2.

                 do n2=1,4
                   ii = connFace(mm,n2)
                   if(coorPoint(ii,3) > half*pi) exit
                 enddo

                 ! Idem for n3 for the third quadrant.

                 do n3=1,4
                   ii = connFace(mm,n3)
                   if(coorPoint(ii,3) <= -half*pi) exit
                 enddo

                 ! Check if the search was successful. If not terminate.

                 if(n2 > 4 .or. n3 > 4)                 &
                   call terminate("localNodesAndQuads", &
                                  "No points in second or third &
                                  &quadrant for polar singular quad")

                 ! Determine the indices n1 and n4.

                 n1 = 2*n2 - n3
                 if(n1 > 4) n1 = 1
                 if(n1 < 1) n1 = 4

                 n4 = 2*n3 - n2
                 if(n4 > 4) n4 = 1
                 if(n4 < 1) n4 = 4

                 ! Store the local node numbers in m1, m2, m3 and m4.

                 m1 = connFace(mm,n1)
                 m2 = connFace(mm,n2)
                 m3 = connFace(mm,n3)
                 m4 = connFace(mm,n4)

                 ! Store 5 quads to represent this "singular" quad. This
                 ! is just a trick to avoid problems. The sequence of the
                 ! connectivities is such that a Cartesian interpolation
                 ! is achieved as closely as possible.

                 ! Quad 1. This is the left most in terms of theta. Two
                 ! additional nodes need to be introduced.

                 ii = iQuad + 1
                 iNode = iNode + 2

                 conn(n1,ii) = iNode - 1
                 conn(n2,ii) = nodeDupl(m2)
                 conn(n3,ii) = nodeFlag(m3)
                 conn(n4,ii) = iNode

                 coor(1,iNode-1) =  coor(1,nodeDupl(m2))
                 coor(2,iNode-1) = -coor(2,nodeDupl(m2))
                 coor(3,iNode-1) =  coor(3,nodeDupl(m2))

                 coorInt(1,iNode-1) =  coorInt(1,nodeDupl(m2))
                 coorInt(2,iNode-1) = -coorInt(2,nodeDupl(m2))
                 coorInt(3,iNode-1) =  coorInt(3,nodeDupl(m2))

                 coor(1,iNode) =  coor(1,nodeFlag(m3))
                 coor(2,iNode) = -coor(2,nodeFlag(m3))
                 coor(3,iNode) =  coor(3,nodeFlag(m3))

                 coorInt(1,iNode) =  coorInt(1,nodeFlag(m3))
                 coorInt(2,iNode) = -coorInt(2,nodeFlag(m3))
                 coorInt(3,iNode) =  coorInt(3,nodeFlag(m3))

                 ! Quad 2. Adjacent to quad 1. One additional node is
                 ! introduced.
 
                 ii = iQuad + 2
                 iNode = iNode + 1

                 conn(n1,ii) = iNode
                 conn(n2,ii) = iNode - 1
                 conn(n3,ii) = nodeFlag(m3)
                 conn(n4,ii) = nodeFlag(m4)

                 coor(1,iNode) =  coor(1,nodeFlag(m4))
                 coor(2,iNode) = -coor(2,nodeFlag(m4))
                 coor(3,iNode) =  coor(3,nodeFlag(m4))

                 coorInt(1,iNode) =  coorInt(1,nodeFlag(m4))
                 coorInt(2,iNode) = -coorInt(2,nodeFlag(m4))
                 coorInt(3,iNode) =  coorInt(3,nodeFlag(m4))

                 ! Quad 3. Adjacent to quad 2. One additional node is
                 ! introduced.
 
                 ii = iQuad + 3
                 iNode = iNode + 1

                 conn(n1,ii) = nodeFlag(m1)
                 conn(n2,ii) = iNode
                 conn(n3,ii) = iNode - 1
                 conn(n4,ii) = nodeFlag(m4)

                 coor(1,iNode) =  coor(1,nodeFlag(m1))
                 coor(2,iNode) = -coor(2,nodeFlag(m1))
                 coor(3,iNode) =  coor(3,nodeFlag(m1))

                 coorInt(1,iNode) =  coorInt(1,nodeFlag(m1))
                 coorInt(2,iNode) = -coorInt(2,nodeFlag(m1))
                 coorInt(3,iNode) =  coorInt(3,nodeFlag(m1))

                 ! Quad 4. Adjacent to quad 3. One additional node is
                 ! introduced.
 
                 ii = iQuad + 4
                 iNode = iNode + 1

                 conn(n1,ii) = nodeFlag(m1)
                 conn(n2,ii) = nodeFlag(m2)
                 conn(n3,ii) = iNode
                 conn(n4,ii) = iNode - 1

                 coor(1,iNode) =  coor(1,nodeFlag(m2))
                 coor(2,iNode) = -coor(2,nodeFlag(m2))
                 coor(3,iNode) =  coor(3,nodeFlag(m2))

                 coorInt(1,iNode) =  coorInt(1,nodeFlag(m2))
                 coorInt(2,iNode) = -coorInt(2,nodeFlag(m2))
                 coorInt(3,iNode) =  coorInt(3,nodeFlag(m2))

                 ! Quad 5. Adjacent to quad 3 and the most right one in
                 ! theta. One additional node is introduced.
 
                 ii = iQuad + 5
                 iNode = iNode + 1

                 conn(n1,ii) = iNode - 1
                 conn(n2,ii) = nodeFlag(m2)
                 conn(n3,ii) = nodeDupl(m3)
                 conn(n4,ii) = iNode

                 coor(1,iNode) =  coor(1,nodeDupl(m3))
                 coor(2,iNode) = -coor(2,nodeDupl(m3))
                 coor(3,iNode) =  coor(3,nodeDupl(m3))

                 coorInt(1,iNode) =  coorInt(1,nodeDupl(m3))
                 coorInt(2,iNode) = -coorInt(2,nodeDupl(m3))
                 coorInt(3,iNode) =  coorInt(3,nodeDupl(m3))

                 ! And set the counter.

                 jj = 5

             end select

             ! Loop over the number of times this face is stored.

             do ii=1,jj

               ! Update the counter iQuad and store the subface id and
               ! quad id.

               iQuad = iQuad + 1
               subface(iQuad) = nn
               quadID(iQuad)  = mm

             enddo

           endif testStoreFace
         enddo connLoop2

         ! Deallocate the memory of nodeFlag and nodeDupl.

         deallocate(nodeFlag, nodeDupl, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("localNodesAndQuads", &
                          "Deallocation error for nodeFlag &
                          &and nodeDupl.")

       enddo loopSubfaces

       ! Check if the numbers of nodes and cells are okay in debug mode.

       if( debug ) then

         if(nQuad /= iQuad) &
           call terminate("localNodesAndQuads", &
                          "nQuad and iQuad differ")

         if(nNode /= iNode) &
           call terminate("localNodesAndQuads", &
                          "nNode and iNode differ")

       endif

       end subroutine localNodesAndQuads
