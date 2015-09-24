!
!      ******************************************************************
!      *                                                                *
!      * File:          storeSurfsolInBuffer.f90                        *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 05-19-2003                                      *
!      * Last modified: 07-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine storeSurfsolInBuffer2(sps, buffer, nn, blockID,   &
     faceID, nodeRange, cellRange, solName, &
     viscousSubface)
!
!      ******************************************************************
!      *                                                                *
!      * storeSurfsolInBuffer stores the variable indicated by          *
!      * solName of the given block ID in the buffer. The solution is   *
!      * stored at the vertices of the boundary. To do so, a weighting  * 
!      * is determined for each node by performing a Newton search for  *
!      * the parametric coordinates of the node in the cell formed by   *
!      * cell centers. This scheme exactly reproducs a linear function. *
!      * Since the value is actually required on the surface, the first *
!      * halo and the first real cell are averaged to compute the value *
!      * on the surface. Note that this is different that the original  *
!      * cell centered implementation. The reason for the modification  *
!      * is that it allows for more accurate contour lines in post      *
!      * processing which is necessary for overset solutions.           *
!      *                                                                *
!      ******************************************************************
!

  use blockPointers
  use BCTypes
  use cgnsNames
  use constants
  use flowVarRefState
  use inputPhysics
  use inputIO
  use communication 
  use costFunctions
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in)    :: sps, blockID, faceID
  integer(kind=intType), intent(inout) :: nn
  integer(kind=intType), dimension(3,2), intent(in) :: nodeRange
  integer(kind=intType), dimension(3,2), intent(in) :: cellRange
  real(kind=realType), dimension(*), intent(out) :: buffer
  character(len=*), intent(in) :: solName
  logical, intent(in) :: viscousSubface
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, sigma, ll
  integer(kind=intType) :: ii, jj, mm, iiMax, jjMax

  integer(kind=intType), dimension(2,2) :: rangeFace, rangeFaceCell
  integer(kind=intType), dimension(3,2) :: rangeNode
  integer(kind=intType), dimension(3,2) :: rangeCell

  integer(kind=intType), dimension(:,:), pointer :: viscPointer
  integer(kind=intType), dimension(:,:), pointer :: iblank2

  real(kind=realType) :: fact, gm1, ptotInf, ptot, psurf, rsurf
  real(kind=realType) :: usurf, vsurf, wsurf, m2surf, musurf
  real(kind=realType) :: fx, fy, fz, fn, a2Tot, a2, qw
  real(kind=realType) :: tauxx, tauyy, tauzz
  real(kind=realType) :: tauxy, tauxz, tauyz, top , bottom
  real(kind=realType) :: pm1, scaleDim, sensor, plocal, sensor1
  real(kind=realType), dimension(3) :: norm

  real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
  real(kind=realType), dimension(:,:,:), pointer :: ss1, ss2, ss
  real(kind=realType), dimension(:,:),   pointer :: pp1, pp2

  real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
  real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
  real(kind=realType), dimension(:,:),   pointer :: dd2Wall
  real(kind=realType), dimension(:,:,:),   pointer :: xx
  real(kind=realType), dimension(3) :: n1, n2, n3, n4, x1, x21, x41, x3142, xf
  real(kind=realType), dimension(3) :: an, bn, vt, b, a, vf
  real(kind=realType) :: u, v, du, dv, uv, uold, vvold, val, vn, invLen

  integer(kind=intType), parameter :: iterMax   = 15
  real(kind=realType),   parameter :: adteps    = 1.e-25_realType
  real(kind=realType),   parameter :: thresConv = 1.e-10_realType
  real(kind=realType), dimension(:, :, :), allocatable :: nodeWeights
  integer(kind=intType), dimension(:, :), allocatable :: iblankTmp

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Set the pointers to this block.

  call setPointers(blockID, 1_intType, sps)

  ! Set the offset for the viscous data, such that the range is
  ! limited to the actual physical face. Viscous data, like skin
  ! friction, need gradIent information, which is not available
  ! in the halo's.

  ! nodeRange contains the range of the current block in the
  ! original cgns block. Substract the offset and store the local
  ! range in rangeNode.
  
  rangeNode(1,1) = nodeRange(1,1) - iBegor + 1
  rangeNode(1,2) = nodeRange(1,2) - iBegor + 1

  rangeNode(2,1) = nodeRange(2,1) - jBegor + 1
  rangeNode(2,2) = nodeRange(2,2) - jBegor + 1

  rangeNode(3,1) = nodeRange(3,1) - kBegor + 1
  rangeNode(3,2) = nodeRange(3,2) - kBegor + 1

  ! CellRange contains the range of the current block in the
  ! original cgns block. Substract the offset and store the local
  ! range in rangeCell.
  
  rangeCell(1,1) = cellRange(1,1) - iBegor + 1 !+ 1
  rangeCell(1,2) = cellRange(1,2) - iBegor + 1 !- 1
  
  rangeCell(2,1) = cellRange(2,1) - jBegor + 1 !+ 1
  rangeCell(2,2) = cellRange(2,2) - jBegor + 1 !- 1
  
  rangeCell(3,1) = cellRange(3,1) - kBegor + 1 !+ 1
  rangeCell(3,2) = cellRange(3,2) - kBegor + 1 !- 1


  !
  !      ******************************************************************
  !      *                                                                *
  !      * Determine the face on which the subface is located and set     *
  !      * a couple of variables accordingly. In this way a generic       *
  !      * treatment is possible and there is no need to repeat the code  *
  !      * for each of the six block faces.                               *
  !      * Note that for dd2Wall a slightly different notation must be    *
  !      * used. Reason is that d2Wall starts at index 2, rather than 0.  *
  !      *                                                                *
  !      ******************************************************************
  !

  select case (faceID)
     
  case (iMin)
     rangeFace(1,1:2) = rangeNode(2,1:2)
     rangeFace(2,1:2) = rangeNode(3,1:2)
     rangeFaceCell(1,1:2) = rangeCell(2,1:2)
     rangeFaceCell(2,1:2) = rangeCell(3,1:2)

     iiMax = jl; jjMax = kl

     ww1    => w(1,1:,1:,:);   ww2    => w(2,1:,1:,:)
     pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)
     ss => si(1,:,:,:) ; fact = -one

     pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)           
     gamma1 => gamma(1,1:,1:); gamma2 => gamma(2,1:,1:)
     xx => x(1, :, :, :)
     if( blockIsMoving)then
        ss1    => s(1,1:,1:,:);   ss2    => s(2,1:,1:,:)
     end if

     iblank2 => iblank(2,1:,1:)
     viscPointer => viscIminPointer

     if( viscous ) then
        rlv1 => rlv(1,1:,1:); rlv2 => rlv(2,1:,1:)
     endif

     if(equations == RANSEquations) dd2Wall => d2Wall(2,:,:)

     !===============================================================

  case (iMax)
     rangeFace(1,1:2) = rangeNode(2,1:2)
     rangeFace(2,1:2) = rangeNode(3,1:2)
     rangeFaceCell(1,1:2) = rangeCell(2,1:2)
     rangeFaceCell(2,1:2) = rangeCell(3,1:2)

     iiMax = jl; jjMax = kl

     ww1    => w(ie,1:,1:,:);   ww2    => w(il,1:,1:,:)
     ss => si(il,:,:,:) ; fact = one
     xx => x(il, :, :, :)

     pp1    => p(ie,1:,1:);     pp2    => p(il,1:,1:)
     gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)

     if( blockIsMoving)then
        ss1    => s(ie-1,1:,1:,:);   ss2    => s(ie,1:,1:,:)
     end if

     iblank2 => iblank(il,1:,1:)
     viscPointer => viscImaxPointer

     if( viscous ) then
        rlv1 => rlv(ie,1:,1:); rlv2 => rlv(il,1:,1:)
     endif

     if(equations == RANSEquations) dd2Wall => d2Wall(il,:,:)

     !===============================================================

  case (jMin)
     rangeFace(1,1:2) = rangeNode(1,1:2)
     rangeFace(2,1:2) = rangeNode(3,1:2)
     rangeFaceCell(1,1:2) = rangeCell(1,1:2)
     rangeFaceCell(2,1:2) = rangeCell(3,1:2)
     iiMax = il; jjMax = kl

     ww1    => w(1:,1,1:,:);   ww2    => w(1:,2,1:,:)
     ss => sj(:,1,:,:) ; fact = -one
     xx => x(:, 1, :, :)

     pp1    => p(1:,1,1:);     pp2    => p(1:,2,1:)
     gamma1 => gamma(1:,1,1:); gamma2 => gamma(1:,2,1:)

     if( blockIsMoving)then
        ss1    => s(1:,1,1:,:);   ss2    => s(1:,2,1:,:)
     end if

     iblank2 => iblank(1:,2,1:)
     viscPointer => viscJminPointer

     if( viscous ) then
        rlv1 => rlv(1:,1,1:); rlv2 => rlv(1:,2,1:)
     endif

     if(equations == RANSEquations) dd2Wall => d2Wall(:,2,:)

     !===============================================================

  case (jMax)
     rangeFace(1,1:2) = rangeNode(1,1:2)
     rangeFace(2,1:2) = rangeNode(3,1:2)
     rangeFaceCell(1,1:2) = rangeCell(1,1:2)
     rangeFaceCell(2,1:2) = rangeCell(3,1:2)
     iiMax = il; jjMax = kl

     ww1    => w(1:,je,1:,:);   ww2    => w(1:,jl,1:,:)
     ss => sj(:,jl,:,:); fact = one
     xx => x(:, jl, :, :)

     pp1    => p(1:,je,1:);     pp2    => p(1:,jl,1:)
     gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)

     if( blockIsMoving)then
        ss1    => s(1:,je-1,1:,:);   ss2    => s(1:,je,1:,:)
     end if

     iblank2 => iblank(1:,jl,1:)
     viscPointer => viscJmaxPointer

     if( viscous ) then
        rlv1 => rlv(1:,je,1:); rlv2 => rlv(1:,jl,1:)
     endif

     if(equations == RANSEquations) dd2Wall => d2Wall(:,jl,:)

     !===============================================================

  case (kMin)
     rangeFace(1,1:2) = rangeNode(1,1:2)
     rangeFace(2,1:2) = rangeNode(2,1:2)
     rangeFaceCell(1,1:2) = rangeCell(1,1:2)
     rangeFaceCell(2,1:2) = rangeCell(2,1:2)
     iiMax = il; jjMax = jl

     ww1    => w(1:,1:,1,:);   ww2    => w(1:,1:,2,:)
     ss => sk(:,:,1,:);  fact = -one
     xx => x(:, :, 1, :)

     pp1    => p(1:,1:,1);     pp2    => p(1:,1:,2)
     gamma1 => gamma(1:,1:,1); gamma2 => gamma(1:,1:,2)

     if( blockIsMoving)then
        ss1    => s(1:,1:,1,:);   ss2    => s(1:,1:,2,:)
     end if

     iblank2 => iblank(1:,1:,2)
     viscPointer => viscKminPointer

     if( viscous ) then
        rlv1 => rlv(1:,1:,1); rlv2 => rlv(1:,1:,2)
     endif

     if(equations == RANSEquations) dd2Wall => d2Wall(:,:,2)

     !===============================================================

  case (kMax)
     rangeFace(1,1:2) = rangeNode(1,1:2)
     rangeFace(2,1:2) = rangeNode(2,1:2)
     rangeFaceCell(1,1:2) = rangeCell(1,1:2)
     rangeFaceCell(2,1:2) = rangeCell(2,1:2)
     iiMax = il; jjMax = jl

     ww1    => w(1:,1:,ke,:);   ww2    => w(1:,1:,kl,:)
     ss => sk(:,:,kl,:);  fact = one
     xx => x(:, :, kl, :)

     pp1    => p(1:,1:,ke);     pp2    => p(1:,1:,kl)
     gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)

     if( blockIsMoving)then
        ss1    => s(1:,1:,ke-1,:);   ss2    => s(1:,1:,ke,:)
     end if

     iblank2 => iblank(1:,1:,kl)
     viscPointer => viscKmaxPointer

     if( viscous ) then
        rlv1 => rlv(1:,1:,ke); rlv2 => rlv(1:,1:,kl)
     endif

     if(equations == RANSEquations) dd2Wall => d2Wall(:,:,kl)

  end select

  !allocate space for the weights:
  allocate(nodeWeights(4,  rangeFace(1, 1):rangeFace(1,2), rangeFace(2, 1):rangeFace(2, 2)))
  call computeNodeWeighting

  !
  !      ******************************************************************
  !      *                                                                *
  !      * The actual part for storing the data. Determine the variable   *
  !      * to be written and loop over the boundary faces of the subface. *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the variable to be written.

  varName: select case (solName)

  case (cgnsDensity)
     
     ! Loop over the given range of nodes
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           nn = nn + 1

           top = nodeWeights(1, i, j) * ww2(i  , j  , irho) + & 
                 nodeWeights(2, i, j) * ww2(i+1, j  , irho) + & 
                 nodeWeights(3, i, j) * ww2(i  , j+1, irho) + & 
                 nodeWeights(4, i, j) * ww2(i+1, j+1, irho) 

           bottom = nodeWeights(1, i, j) * ww1(i  , j  , irho) + & 
                    nodeWeights(2, i, j) * ww1(i+1, j  , irho) + & 
                    nodeWeights(3, i, j) * ww1(i  , j+1, irho) + & 
                    nodeWeights(4, i, j) * ww1(i+1, j+1, irho) 
           buffer(nn) = half*(top + bottom)
        end do
     end do


  case (cgnsPressure)
     
     ! Loop over the given range of nodes
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           nn = nn + 1

           top = nodeWeights(1, i, j) * pp2(i  , j) + & 
                 nodeWeights(2, i, j) * pp2(i+1, j) + & 
                 nodeWeights(3, i, j) * pp2(i  , j+1) + & 
                 nodeWeights(4, i, j) * pp2(i+1, j+1) 

           bottom = nodeWeights(1, i, j) * pp1(i  , j) + & 
                    nodeWeights(2, i, j) * pp1(i+1, j) + & 
                    nodeWeights(3, i, j) * pp1(i  , j+1) + & 
                    nodeWeights(4, i, j) * pp1(i+1, j+1) 
           buffer(nn) = half*(top + bottom)
        end do
     end do

  case (cgnsVelx)
     
     ! Loop over the given range of nodes
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           nn = nn + 1

           top = nodeWeights(1, i, j) * ww2(i  , j  , iVx) + & 
                 nodeWeights(2, i, j) * ww2(i+1, j  , iVx) + & 
                 nodeWeights(3, i, j) * ww2(i  , j+1, iVx) + & 
                 nodeWeights(4, i, j) * ww2(i+1, j+1, iVx) 

           bottom = nodeWeights(1, i, j) * ww1(i  , j  , iVx) + & 
                    nodeWeights(2, i, j) * ww1(i+1, j  , iVx) + & 
                    nodeWeights(3, i, j) * ww1(i  , j+1, iVx) + & 
                    nodeWeights(4, i, j) * ww1(i+1, j+1, iVx) 
           buffer(nn) = half*(top + bottom)
        end do
     end do

  case (cgnsVely)
     
     ! Loop over the given range of nodes
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           nn = nn + 1

           top = nodeWeights(1, i, j) * ww2(i  , j  , iVy) + & 
                 nodeWeights(2, i, j) * ww2(i+1, j  , iVy) + & 
                 nodeWeights(3, i, j) * ww2(i  , j+1, iVy) + & 
                 nodeWeights(4, i, j) * ww2(i+1, j+1, iVy) 

           bottom = nodeWeights(1, i, j) * ww1(i  , j  , iVy) + & 
                    nodeWeights(2, i, j) * ww1(i+1, j  , iVy) + & 
                    nodeWeights(3, i, j) * ww1(i  , j+1, iVy) + & 
                    nodeWeights(4, i, j) * ww1(i+1, j+1, iVy) 
           buffer(nn) = half*(top + bottom)
        end do
     end do


  case (cgnsVelz)
     
     ! Loop over the given range of nodes
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           nn = nn + 1

           top = nodeWeights(1, i, j) * ww2(i  , j  , iVz) + & 
                 nodeWeights(2, i, j) * ww2(i+1, j  , iVz) + & 
                 nodeWeights(3, i, j) * ww2(i  , j+1, iVz) + & 
                 nodeWeights(4, i, j) * ww2(i+1, j+1, iVz) 

           bottom = nodeWeights(1, i, j) * ww1(i  , j  , iVz) + & 
                    nodeWeights(2, i, j) * ww1(i+1, j  , iVz) + & 
                    nodeWeights(3, i, j) * ww1(i  , j+1, iVz) + & 
                    nodeWeights(4, i, j) * ww1(i+1, j+1, iVz) 
           buffer(nn) = half*(top + bottom)
        end do
     end do

  case (cgnsBlank)

     ! This one is a litle tricker. The procedure is:
     ! 1. Set all iblanks to 1
     ! 2. Loop over cells, check if iblank = 0, set all nodes around that cell to 0
     ! 3. Loop over cells, check if iblank = -1, set all nodes around that cell to -1
     
     allocate(iblankTmp(rangeFace(1,1):rangeFace(1,2), rangeFace(2,1):rangeFace(2,2)))
     iBlankTmp = 1

     ! Loop over the cells checking for 0:
     do j=rangeFaceCell(2,1), rangeFaceCell(2,2)
        do i=rangeFaceCell(1,1), rangeFaceCell(1,2)
           if (iBlank2(i, j) == 0) then 
              iBlankTmp(i-1, j-1) = 0
              iBlankTmp(i  , j-1) = 0
              iBlankTmp(i-1, j  ) = 0
              iBlankTmp(i  , j  ) = 0
           end if
        end do
     end do

     ! Loop over the cells checking for -1
     do j=rangeFaceCell(2,1), rangeFaceCell(2,2)
        do i=rangeFaceCell(1,1), rangeFaceCell(1,2)
           if (iBlank2(i, j) == -1) then 
              iBlankTmp(i-1, j-1) = -1
              iBlankTmp(i  , j-1) = -1
              iBlankTmp(i-1, j  ) = -1
              iBlankTmp(i  , j  ) = -1
           end if
        end do
     end do

     ! We can now just copy in
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           nn = nn + 1
           buffer(nn) = real(iBlankTmp(i, j), realType)
        enddo
     enddo

     deallocate(iblankTmp)

  case (cgnsTemp, cgnsRelVelX, cgnsRelVely, cgnsRelVelz, &
       cgnsPtotLoss, cgnsCp, cgnsMach, cgnsRelMach, &
       cgnsLift, cgnsDrag, cgnsSepSensor, cgnsCavitation)

     print *,' Not implemented yet'
     stop
     
  case (cgnsSkinFmag, cgnsYplus, &
       cgnsSkinFx, cgnsSkinFy, cgnsSkinFz)

     print *,'not done for node based output'
     stop
  case (cgnsStanton)
     print *,'not done for node based output'
     stop

  end select varName
  
  ! Clean up memory
  deallocate(nodeWeights)

contains
  
  subroutine computeNodeWeighting


    ! Loop over the given range of nodes
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)

           ! Do this correctly. Extract out the dual mesh surrounding
           ! this node. Everything is shifted by 1 due to pointer offset

           n1 = fourth*(xx(i  , j  , :) + xx(i+1, j, :) + xx(i+1, j+1, :) + xx(i  , j+1, :))
           n2 = fourth*(xx(i+1, j  , :) + xx(i+2, j, :) + xx(i+2, j+1, :) + xx(i+1, j+1, :))
           n3 = fourth*(xx(i  , j+1, :) + xx(i+1, j+1, :) + xx(i+1, j+2, :) + xx(i  , j+2, :))
           n4 = fourth*(xx(i+1, j+1, :) + xx(i+2, j+1, :) + xx(i+2, j+2, :) + xx(i+1, j+2, :))

           ! We now need to find the weights for point xx(i, j, :) in
           ! the quad defined n1, through n4 (cyclic). This will tell
           ! us the 4 weights each of the 4 cells. 

           x1(1) = n1(1)
           x1(2) = n1(2)
           x1(3) = n1(3)

           x21(1) = n2(1) - x1(1)
           x21(2) = n2(2) - x1(2)
           x21(3) = n2(3) - x1(3)

           x41(1) = n3(1) - x1(1)
           x41(2) = n3(2) - x1(2)
           x41(3) = n3(3) - x1(3)

           x3142(1) = n4(1) - x1(1) - x21(1) - x41(1)
           x3142(2) = n4(2) - x1(2) - x21(2) - x41(2)
           x3142(3) = n4(3) - x1(3) - x21(3) - x41(3)

           ! Initialize u and v to 0.5 and determine the
           ! corresponding coordinates on the face, which is the
           ! centroid.

           u  = half
           v  = half
           uv = u*v

           xf(1) = x1(1) + u*x21(1) + v*x41(1) + uv*x3142(1)
           xf(2) = x1(2) + u*x21(2) + v*x41(2) + uv*x3142(2)
           xf(3) = x1(3) + u*x21(3) + v*x41(3) + uv*x3142(3)
           
           ! Newton loop to determine the point on the surface,
           ! which minimizes the distance to the given coordinate.
           
           NewtonQuads: do ll=1,iterMax

              ! Store the current values of u and v for a stop
              ! criterion later on.
              
              uold = u
              vvold = v
              
              ! Determine the vector vf from xf to given coordinate.

              vf(1) = xx(i+1, j+1, 1) - xf(1)
              vf(2) = xx(i+1, j+1, 2) - xf(2)
              vf(3) = xx(i+1, j+1, 3) - xf(3)
              
              ! Determine the tangent vectors in u- and v-direction.
              ! Store these in a and b respectively.
              
              a(1) = x21(1) + v*x3142(1)
              a(2) = x21(2) + v*x3142(2)
              a(3) = x21(3) + v*x3142(3)
              
              b(1) = x41(1) + u*x3142(1)
              b(2) = x41(2) + u*x3142(2)
              b(3) = x41(3) + u*x3142(3)
              
              ! Determine the normal vector of the face by taking the
              ! cross product of a and b. Afterwards this vector will
              ! be scaled to a unit vector.
              
              norm(1) = a(2)*b(3) - a(3)*b(2)
              norm(2) = a(3)*b(1) - a(1)*b(3)
              norm(3) = a(1)*b(2) - a(2)*b(1)
              
              invLen = one/max(adtEps, sqrt(norm(1)*norm(1) &
                   +                        norm(2)*norm(2) &
                   +                        norm(3)*norm(3)))
              
              norm(1) = norm(1)*invLen
              norm(2) = norm(2)*invLen
              norm(3) = norm(3)*invLen

              ! Determine the projection of the vector vf onto
              ! the face.
              
              vn = vf(1)*norm(1) + vf(2)*norm(2) + vf(3)*norm(3)
              vt(1) = vf(1) - vn*norm(1)
              vt(2) = vf(2) - vn*norm(2)
              vt(3) = vf(3) - vn*norm(3)
              
              ! The vector vt points from the current point on the
              ! face to the new point. However this new point lies on
              ! the plane determined by the vectors a and b, but not
              ! necessarily on the face itself. The new point on the
              ! face is obtained by projecting the point in the a-b
              ! plane onto the face. this can be done by determining
              ! the coefficients du and dv, such that vt = du*a + dv*b.
              ! To solve du and dv the vectors normal to a and b
              ! inside the plane ab are needed.
              
              an(1) = a(2)*norm(3) - a(3)*norm(2)
              an(2) = a(3)*norm(1) - a(1)*norm(3)
              an(3) = a(1)*norm(2) - a(2)*norm(1)
              
              bn(1) = b(2)*norm(3) - b(3)*norm(2)
              bn(2) = b(3)*norm(1) - b(1)*norm(3)
              bn(3) = b(1)*norm(2) - b(2)*norm(1)
              
              ! Solve du and dv. the clipping of vn should not be
              ! active, as this would mean that the vectors a and b
              ! are parallel. This corresponds to a quad degenerated
              ! to a line, which should not occur in the surface mesh.
              
              vn = a(1)*bn(1) + a(2)*bn(2) + a(3)*bn(3)
              vn = sign(max(eps,abs(vn)),vn)
              du = (vt(1)*bn(1) + vt(2)*bn(2) + vt(3)*bn(3))/vn
              
              vn = b(1)*an(1) + b(2)*an(2) + b(3)*an(3)
              vn = sign(max(eps,abs(vn)),vn)
              dv = (vt(1)*an(1) + vt(2)*an(2) + vt(3)*an(3))/vn
              
              ! Determine the new parameter values uu and vv. These
              ! are limited to 0 <= (uu,vv) <= 1.
              
              u = u + du; u = min(one,max(zero,u))
              v = v + dv; v = min(one,max(zero,v))
              
              ! Determine the final values of the corrections.

              du = abs(u-uold)
              dv = abs(v-vvold)
              
              ! Determine the new coordinates of the point xf.
              
              uv  = u*v
              xf(1) = x1(1) + u*x21(1) + v*x41(1) + uv*x3142(1)
              xf(2) = x1(2) + u*x21(2) + v*x41(2) + uv*x3142(2)
              xf(3) = x1(3) + u*x21(3) + v*x41(3) + uv*x3142(3)
              
              ! Exit the loop if the update of the parametric
              ! weights is below the threshold
              
              val = sqrt(du*du + dv*dv)
              if(val <= thresConv) exit NewtonQuads
              
           enddo NewtonQuads
         
           ! Finally store the weights
           nodeWeights(1, i, j) = (one - u)*(one -v)
           nodeWeights(2, i, j) = (      u)*(one -v)
           nodeWeights(3, i, j) = (one - u)*(     v)
           nodeWeights(4, i, j) = (      u)*(     v)

        end do
     end do
  end subroutine computeNodeWeighting

 
end subroutine storeSurfsolInBuffer2
