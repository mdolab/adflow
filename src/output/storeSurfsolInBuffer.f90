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
subroutine storeSurfsolInBuffer(sps, buffer, nn, blockID,   &
     faceID, nodeRange, cellRange, &
     solName, viscousSubface)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * storeSurfsolInBuffer stores the variable indicated by          *
  !      * solName of the given block ID in the buffer. As the solution   *
  !      * must be stored in the center of the boundary face the average  *
  !      * value of the first internal cell and its corresponding halo is *
  !      * computed. The counter nn is updated in this routine. However   *
  !      * it is not initialized, because multiple contributions may be   *
  !      * stored in buffer.                                              *
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
  integer(kind=intType) :: i, j, k, bufSize
  integer(kind=intType) :: ii, jj, mm, ll, iiMax, jjMax, offVis

  integer(kind=intType), dimension(2,2) :: rangeFace, rangeFaceNode
  integer(kind=intType), dimension(3,2) :: rangeCell, rangeNode

  integer(kind=intType), dimension(:,:), pointer :: viscPointer
  integer(kind=intType), dimension(:,:), pointer :: iblank2

  real(kind=realType) :: fact, gm1, ptotInf, ptot, psurf, rsurf, sigma
  real(kind=realType) :: usurf, vsurf, wsurf, m2surf, musurf
  real(kind=realType) :: fx, fy, fz, fn, a2Tot, a2, qw
  real(kind=realType) :: tauxx, tauyy, tauzz
  real(kind=realType) :: tauxy, tauxz, tauyz
  real(kind=realType) :: pm1, scaleDim, a, sensor, plocal, sensor1
  real(kind=realType), dimension(3) :: norm, V

  real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
  real(kind=realType), dimension(:,:,:), pointer :: ss1, ss2, ss
  real(kind=realType), dimension(:,:),   pointer :: pp1, pp2

  real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
  real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
  real(kind=realType), dimension(:,:),   pointer :: dd2Wall
  real(kind=realType), dimension(:, :), allocatable :: tmpBuffer
  real(kind=realType), dimension(:, :, :), pointer :: nodeWeights
  integer(kind=intType), dimension(:, :), allocatable :: iblankTmp

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Starting counter for the temporary buffer
  ll = 0

  ! Set the pointers to this block.

  call setPointers(blockID, 1_intType, sps)

  ! Set the offset for the viscous data, such that the range is
  ! limited to the actual physical face. Viscous data, like skin
  ! friction, need gradIent information, which is not available
  ! in the halo's.

  offVis = 0
  if( storeRindLayer .or. nodalOutput ) offVis = 1

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

  rangeCell(1,1) = cellRange(1,1) - iBegor + 1
  rangeCell(1,2) = cellRange(1,2) - iBegor + 1

  rangeCell(2,1) = cellRange(2,1) - jBegor + 1
  rangeCell(2,2) = cellRange(2,2) - jBegor + 1

  rangeCell(3,1) = cellRange(3,1) - kBegor + 1
  rangeCell(3,2) = cellRange(3,2) - kBegor + 1

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

  nodeWeights => flowDoms(blockID, 1_intType, sps)%nodalWeights(faceID)%weight

  select case (faceID)

  case (iMin)
     rangeFace(1,1:2) = rangeCell(2,1:2)
     rangeFace(2,1:2) = rangeCell(3,1:2)
     rangeFaceNode(1,1:2) = rangeNode(2,1:2)
     rangeFaceNode(2,1:2) = rangeNode(3,1:2)

     iiMax = jl; jjMax = kl

     ww1    => w(1,1:,1:,:);   ww2    => w(2,1:,1:,:)
     pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)
     ss => si(1,:,:,:) ; fact = -one

     pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)           
     gamma1 => gamma(1,1:,1:); gamma2 => gamma(2,1:,1:)

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
     rangeFace(1,1:2) = rangeCell(2,1:2)
     rangeFace(2,1:2) = rangeCell(3,1:2)
     rangeFaceNode(1,1:2) = rangeNode(2,1:2)
     rangeFaceNode(2,1:2) = rangeNode(3,1:2)
     iiMax = jl; jjMax = kl

     ww1    => w(ie,1:,1:,:);   ww2    => w(il,1:,1:,:)
     ss => si(il,:,:,:) ; fact = one

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
     rangeFace(1,1:2) = rangeCell(1,1:2)
     rangeFace(2,1:2) = rangeCell(3,1:2)
     rangeFaceNode(1,1:2) = rangeNode(1,1:2)
     rangeFaceNode(2,1:2) = rangeNode(3,1:2)
     iiMax = il; jjMax = kl
   
     ww1    => w(1:,1,1:,:);   ww2    => w(1:,2,1:,:)
     ss => sj(:,1,:,:) ; fact = -one

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
     rangeFace(1,1:2) = rangeCell(1,1:2)
     rangeFace(2,1:2) = rangeCell(3,1:2)
     rangeFaceNode(1,1:2) = rangeNode(1,1:2)
     rangeFaceNode(2,1:2) = rangeNode(3,1:2)
     iiMax = il; jjMax = kl

     ww1    => w(1:,je,1:,:);   ww2    => w(1:,jl,1:,:)
     ss => sj(:,jl,:,:); fact = one

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
     rangeFace(1,1:2) = rangeCell(1,1:2)
     rangeFace(2,1:2) = rangeCell(2,1:2)
     rangeFaceNode(1,1:2) = rangeNode(1,1:2)
     rangeFaceNode(2,1:2) = rangeNode(2,1:2)
     iiMax = il; jjMax = jl

     ww1    => w(1:,1:,1,:);   ww2    => w(1:,1:,2,:)
     ss => sk(:,:,1,:);  fact = -one

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
     rangeFace(1,1:2) = rangeCell(1,1:2)
     rangeFace(2,1:2) = rangeCell(2,1:2)
     rangeFaceNode(1,1:2) = rangeNode(1,1:2)
     rangeFaceNode(2,1:2) = rangeNode(2,1:2)
     iiMax = il; jjMax = jl

     ww1    => w(1:,1:,ke,:);   ww2    => w(1:,1:,kl,:)
     ss => sk(:,:,kl,:);  fact = one

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

  ! Allocate space for the temporary buffer
  allocate(tmpBuffer(&
       rangeFace(1,1):rangeFace(1,2), &
       rangeFace(2,1):rangeFace(2,2)))

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

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           tmpBuffer(i, j) = half*(ww1(i,j,irho) + ww2(i,j,irho))
        enddo
     enddo

     !===============================================================

  case (cgnsPressure)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           tmpBuffer(i, j) = half*(pp1(i,j) + pp2(i,j))
        enddo
     enddo

     !===============================================================

  case (cgnsTemp)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           tmpBuffer(i, j) = (pp1(i,j) + pp2(i,j)) &
                / (RGas*(ww1(i,j,irho) + ww2(i,j,irho)))
        enddo
     enddo

     !===============================================================

  case (cgnsVelx)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           if (viscousSurfaceVelocities .and. viscous) then
              tmpBuffer(i, j) = ww2(i,j,ivx)
           else
              tmpBuffer(i, j) = half*(ww1(i,j,ivx) + ww2(i,j,ivx))
           end if
        enddo
     enddo

     !===============================================================

  case (cgnsVely)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           if (viscousSurfaceVelocities .and. viscous) then
              tmpBuffer(i, j) = ww2(i,j,ivy)
           else
              tmpBuffer(i, j) = half*(ww1(i,j,ivy) + ww2(i,j,ivy))
           end if
        enddo
     enddo

     !===============================================================

  case (cgnsVelz)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           if (viscousSurfaceVelocities .and. viscous) then
              tmpBuffer(i, j) = ww2(i,j,ivz)
           else
              tmpBuffer(i, j) = half*(ww1(i,j,ivz) + ww2(i,j,ivz))
           end if

        enddo
     enddo

     !===============================================================

  case (cgnsRelVelx)
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           if (viscousSurfaceVelocities .and. viscous) then
              tmpBuffer(i, j) = ww2(i,j,ivx) - ss2(i,j,1)
           else
              tmpBuffer(i, j) = half*(ww1(i,j,ivx) + ww2(i,j,ivx))-half*(ss1(i,j,1) + ss2(i,j,1))
           end if
        enddo
     enddo

     !===============================================================

  case (cgnsRelVely)
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           if (viscousSurfaceVelocities .and. viscous) then
              tmpBuffer(i, j) = ww2(i,j,ivy) - ss2(i,j,2)
           else
              tmpBuffer(i, j) = half*(ww1(i,j,ivy) + ww2(i,j,ivy))-half*(ss1(i,j,2) + ss2(i,j,2))
           end if
        enddo
     enddo

     !===============================================================

  case (cgnsRelVelz)
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           if (viscousSurfaceVelocities .and. viscous) then
              tmpBuffer(i, j) = ww2(i,j,ivz) - ss2(i,j,3)
           else
              tmpBuffer(i, j) = half*(ww1(i,j,ivz) + ww2(i,j,ivz))-half*(ss1(i,j,3) + ss2(i,j,3))
           end if
        enddo
     enddo


     !================================================================

  case (cgnsCp)

     ! Factor multiplying p-pInf

     fact = two/(gammaInf*pInf*MachCoef*MachCoef)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           tmpBuffer(i, j) = fact*(half*(pp1(i,j) + pp2(i,j)) - pInf)
        enddo
     enddo

     !===============================================================

  case (cgnsPtotloss)

     ! First compute the total pressure of the free stream.

     call computePtot(rhoInf, uInf, zero, zero, &
          pInf, ptotInf, 1_intType)
     ptotInf = one/ptotInf

     ! Loop over the faces and compute the total pressure loss.

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)

           psurf = half*(pp1(i,j) + pp2(i,j))
           rsurf = half*(ww1(i,j,irho) + ww2(i,j,irho))
           usurf = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
           vsurf = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
           wsurf = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))

           call computePtot(rsurf, usurf, vsurf, wsurf, &
                psurf, ptot, 1_intType)

           tmpBuffer(i, j) = one - ptot*ptotInf
        enddo
     enddo

     !===============================================================

  case (cgnsMach)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)

           psurf  = half*(pp1(i,j) + pp2(i,j))
           rsurf  = half*(ww1(i,j,irho) + ww2(i,j,irho))
           usurf  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
           vsurf  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
           wsurf  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))
           m2surf = rsurf*(usurf**2 + vsurf**2 + wsurf**2) &
                / (half*(gamma1(i,j) + gamma2(i,j))*psurf)
           tmpBuffer(i, j) = sqrt(m2surf)
        enddo
     enddo


     !===============================================================

  case (cgnsRelMach)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)

           psurf  = half*(pp1(i,j) + pp2(i,j))
           rsurf  = half*(ww1(i,j,irho) + ww2(i,j,irho))
           usurf  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))-half*(ss1(i,j,1) + ss2(i,j,1))
           vsurf  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))-half*(ss1(i,j,2) + ss2(i,j,2))
           wsurf  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))-half*(ss1(i,j,3) + ss2(i,j,3))
           m2surf = rsurf*(usurf**2 + vsurf**2 + wsurf**2) &
                / (half*(gamma1(i,j) + gamma2(i,j))*psurf)
           tmpBuffer(i, j) = sqrt(m2surf)
        enddo
     enddo

     !        ================================================================

  case (cgnsSkinFmag, cgnsYplus, &
       cgnsSkinFx, cgnsSkinFy, cgnsSkinFz)

     viscSubface1: if(.not. viscousSubface) then
        do j=rangeFace(2,1), rangeFace(2,2)
           do i=rangeFace(1,1), rangeFace(1,2)
              tmpBuffer(i, j) = zero
           end do
        end do
     else

        ! To avoid a lot of code duplication these 5 variables are
        ! treated together.

        ! Multiplication factor to obtain the skin friction from
        ! the wall shear stress.

        fact = two/(gammaInf*pInf*MachCoef*MachCoef)

        ! Loop over the given range of faces. As the viscous data is
        ! only present in the owned faces, the values of the halo's
        ! are set equal to the nearest physical face. Therefore the
        ! working indices are ii and jj.

        do j=rangeFace(2,1), rangeFace(2,2)
           if(j == rangeFace(2,1)) then
              jj = j + offVis
           else if(j == rangeFace(2,2)) then
              jj = j - offVis
           else
              jj = j
           endif

           do i=rangeFace(1,1), rangeFace(1,2)
              if(i == rangeFace(1,1)) then
                 ii = i + offVis
              else if(i == rangeFace(1,2)) then
                 ii = i - offVis
              else
                 ii = i
              endif

              ! Determine the viscous subface on which this
              ! face is located.

              mm = viscPointer(ii,jj)

              ! Store the 6 components of the viscous stress tensor
              ! a bit easier.

              tauxx = viscSubface(mm)%tau(ii,jj,1)
              tauyy = viscSubface(mm)%tau(ii,jj,2)
              tauzz = viscSubface(mm)%tau(ii,jj,3)
              tauxy = viscSubface(mm)%tau(ii,jj,4)
              tauxz = viscSubface(mm)%tau(ii,jj,5)
              tauyz = viscSubface(mm)%tau(ii,jj,6)

              ! Compute the "unit" force on this face. The unit normal
              ! is outward pointing per definition. A minus sign is
              ! present, because of the definition of the viscous
              ! stress tensor. Note that in the normal the indices i
              ! and j could be used. However this is not done.

              norm(1) = BCData(mm)%norm(ii,jj,1)
              norm(2) = BCData(mm)%norm(ii,jj,2)
              norm(3) = BCData(mm)%norm(ii,jj,3)

              fx = -(tauxx*norm(1) + tauxy*norm(2) + tauxz*norm(3))
              fy = -(tauxy*norm(1) + tauyy*norm(2) + tauyz*norm(3))
              fz = -(tauxz*norm(1) + tauyz*norm(2) + tauzz*norm(3))

              fn = fx*norm(1) + fy*norm(2) + fz*norm(3)

              fx = fx - fn*norm(1)
              fy = fy - fn*norm(2)
              fz = fz - fn*norm(3)

              ! Determine the variable to be stored and compute it.
              ! Note that an offset of -1 must be used in dd2Wall,
              ! because the original array, d2Wall, starts at 2.
              ! First update the counter ll.

              select case (solName)
              case (cgnsSkinFmag)
                 tmpBuffer(i, j) = fact*sqrt(fx*fx + fy*fy + fz*fz)

              case (cgnsSkinFx)
                 tmpBuffer(i, j) = fact*fx

              case (cgnsSkinFy)
                 tmpBuffer(i, j) = fact*fy

              case (cgnsSkinFz)
                 tmpBuffer(i, j) = fact*fz

              case (cgnsYplus)
                 rsurf      = half*(ww1(ii,jj,irho) + ww2(ii,jj,irho))
                 musurf     = half*(rlv1(ii,jj)     + rlv2(ii,jj))
                 tmpBuffer(i, j) = sqrt(rsurf*sqrt(fx*fx + fy*fy + fz*fz)) &
                      * dd2Wall(ii-1,jj-1)/musurf
              end select

           enddo
        enddo
     end if viscSubface1
     !        ================================================================

  case (cgnsStanton)

     viscSubface2: if(.not. viscousSubface) then
        do j=rangeFace(2,1), rangeFace(2,2)
           do i=rangeFace(1,1), rangeFace(1,2)
              tmpBuffer(i, j) = zero
           end do
        end do
     else
        ! Some constants needed to compute the stanton number.

        gm1   = gammaInf - one
        a2Tot = gammaInf*pInf*(one + half*gm1*MachCoef*MachCoef) &
             / rhoInf
        fact   = MachCoef*sqrt(gammaInf*pInf*rhoInf)/gm1

        ! Loop over the given range of faces. As the viscous data is
        ! only present in the owned faces, the values of the halo's
        ! are set equal to the nearest physical face. Therefore the
        ! working indices are ii and jj.

        do j=rangeFace(2,1), rangeFace(2,2)
           if(j == rangeFace(2,1)) then
              jj = j + offVis
           else if(j == rangeFace(2,2)) then
              jj = j - offVis
           else
              jj = j
           endif

           do i=rangeFace(1,1), rangeFace(1,2)
              if(i == rangeFace(1,1)) then
                 ii = i + offVis
              else if(i == rangeFace(1,2)) then
                 ii = i - offVis
              else
                 ii = i
              endif

              ! Determine the viscous subface on which this
              ! face is located.

              mm = viscPointer(ii,jj)

              ! Compute the heat flux. Multipy with the sign of the
              ! normal to obtain the correct value.

              qw = viscSubface(mm)%q(ii,jj,1)*BCData(mm)%norm(ii,jj,1) &
                   + viscSubface(mm)%q(ii,jj,2)*BCData(mm)%norm(ii,jj,2) &
                   + viscSubface(mm)%q(ii,jj,3)*BCData(mm)%norm(ii,jj,3)

              ! Compute the speed of sound squared at the wall and
              ! the stanton number, which is stored in tmpBuffer.

              a2 = half*(gamma1(ii,jj)   + gamma2(ii,jj)) &
                   *      (pp1(ii,jj)      + pp2(ii,jj))    &
                   /      (ww1(ii,jj,irho) + ww2(ii,jj,irho))

              tmpBuffer(i, j) = qw/(fact*(a2Tot-a2))

           enddo
        enddo
     end if viscSubface2
     !        ================================================================

  case (cgnsBlank)

     ! Loop over the given range of faces. 

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           tmpBuffer(i, j) = real(iblank2(i,j), realType)
        enddo
     enddo

  case (cgnsLift,cgnsDrag) 

     fact = fact*two/(gammaInf*pInf*MachCoef*MachCoef)
     scaleDim = pRef/pInf
     do j=rangeFace(2,1), rangeFace(2,2)
        if(j == rangeFace(2,1)) then
           jj = j + offVis
        else if(j == rangeFace(2,2)) then
           jj = j - offVis
        else
           jj = j
        endif

        do i=rangeFace(1,1), rangeFace(1,2)
           if(i == rangeFace(1,1)) then
              ii = i + offVis
           else if(i == rangeFace(1,2)) then
              ii = i - offVis
           else
              ii = i
           endif

           norm(1) = ss(ii,jj,1)
           norm(2) = ss(ii,jj,2)
           norm(3) = ss(ii,jj,3)

           ! Compute inviscid force
           pm1 = fact*(half*(pp2(i,j) + pp1(i,j)) - pInf)*scaleDim
           fx = pm1*norm(1)
           fy = pm1*norm(2)
           fz = pm1*norm(3)

           ! Compute possible viscous force
           if (viscousSubface) then

              ! Determine the viscous subface on which this
              ! face is located.

              mm = viscPointer(ii,jj)

              ! Store the 6 components of the viscous stress tensor
              ! a bit easier.

              tauxx = viscSubface(mm)%tau(ii,jj,1)
              tauyy = viscSubface(mm)%tau(ii,jj,2)
              tauzz = viscSubface(mm)%tau(ii,jj,3)
              tauxy = viscSubface(mm)%tau(ii,jj,4)
              tauxz = viscSubface(mm)%tau(ii,jj,5)
              tauyz = viscSubface(mm)%tau(ii,jj,6)

              fx = fx -(tauxx*norm(1) + tauxy*norm(2) + tauxz*norm(3))*fact*scaleDim
              fy = fy -(tauxy*norm(1) + tauyy*norm(2) + tauyz*norm(3))*fact*scaleDim
              fz = fz -(tauxz*norm(1) + tauyz*norm(2) + tauzz*norm(3))*fact*scaleDim

           end if

           ! Next we get the traction by dividing by the area
           a =  sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)
           fx = fx / a
           fy = fy / a
           fz = fz / a
           select case(solName)
           case (cgnsLift)
              ! Take dot-product with lift vector
              tmpBuffer(i, j) = fx*liftDirection(1) +fy*liftDirection(2) + fz*liftDirection(3)
           case(cgnsDrag)
              tmpBuffer(i, j) = fx*dragDirection(1) +fy*dragDirection(2) + fz*dragDirection(3)
           end select
        end do
     end do

  case (cgnsSepSensor)

     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)

           ! Get normalized surface velocity:
           v(1) = ww2(i, j, ivx)
           v(2) = ww2(i, j, ivy)
           v(3) = ww2(i, j, ivz)

           ! Normalize
           v = v / (sqrt(v(1)**2 + v(2)**2 + v(3)**2) + 1e-16)

           ! Dot product with free stream
           sensor = -dot_product(v, velDirFreeStream)

           !Now run through a smooth heaviside function:
           sensor = one/(one + exp(-2*sepSensorSharpness*(sensor - sepSensorOffset)))
           tmpBuffer(i, j) = sensor
        enddo
     enddo

  case (cgnsCavitation)
     fact = two/(gammaInf*pInf*MachCoef*MachCoef)
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)

           ! Get local pressure
           plocal = half*(pp1(i,j) + pp2(i,j))

           sigma = 1.4 
           sensor1 = (-(fact)*(plocal-pInf))- sigma
           sensor1 = one/(one + exp(-2*10*sensor1))
           tmpBuffer(i, j) = sensor1

        enddo
     enddo
  end select varName


  ! We have now computed the cell centered values in tmpBuf. If this
  ! is what we want, we can simply copy the values into the actual
  ! buffer:

  if (.not. nodalOutput) then 
     do j=rangeFace(2,1), rangeFace(2,2)
        do i=rangeFace(1,1), rangeFace(1,2)
           nn = nn + 1
           buffer(nn) = tmpBuffer(i, j)
        end do
     end do
  else

     ! For the nodal output we have to use the precomputed nodal
     ! weights to average the values back to the nodes. We have to do
     ! the iblank variables slightly differently since averaging with
     ! the weights makes no sense for that variable. 

     varName2: select case (solName)

        case (cgnsBlank) 

           ! This one is a litle tricker. The procedure is:
           ! 1. Set all iblanks to 1
           ! 2. Loop over REAL cells, check if cell iblank = 0, set all nodes around that cell to 0
           ! 3. Loop over REAL cells, check if cell iblank = -1, set all nodes around that cell to -1

           ! Note that in the loops below we have +1 on the low end
           ! and -1 on the high end. This is because the cell indices
           ! include the rind cells (due to requesting nodal output),
           ! but here we want to loop only over the real cells. 
  
           allocate(iblankTmp(rangeFaceNode(1,1):rangeFaceNode(1,2), &
                              rangeFaceNode(2,1):rangeFaceNode(2,2)))
           iBlankTmp = 1

           ! Loop over the cells checking for 0:
           do j=rangeFace(2,1)+1, rangeFace(2,2)-1
              do i=rangeFace(1,1)+1, rangeFace(1,2)-1
                 if (iBlank2(i, j) == 0) then 
                    iBlankTmp(i-1, j-1) = 0
                    iBlankTmp(i  , j-1) = 0
                    iBlankTmp(i-1, j  ) = 0
                    iBlankTmp(i  , j  ) = 0
                 end if
              end do
           end do
           
           ! Loop over the cells checking for -1
           do j=rangeFace(2,1)+1, rangeFace(2,2)-1
              do i=rangeFace(1,1)+1, rangeFace(1,2)-1
                 if (iBlank2(i, j) == -1) then 
                    iBlankTmp(i-1, j-1) = -1
                    iBlankTmp(i  , j-1) = -1
                    iBlankTmp(i-1, j  ) = -1
                    iBlankTmp(i  , j  ) = -1
                 end if
              end do
           end do
           
           ! We can now just copy in
           do j=rangeFaceNode(2,1), rangeFaceNode(2,2)
              do i=rangeFaceNode(1,1), rangeFaceNode(1,2)
                 nn = nn + 1
                 buffer(nn) = real(iBlankTmp(i, j), realType)
              enddo
           enddo

           deallocate(iblankTmp)

        case default

           do j=rangeFaceNode(2,1), rangeFaceNode(2,2)
              do i=rangeFaceNode(1,1), rangeFaceNode(1,2)
                 nn = nn + 1
                 buffer(nn) = & 
                      nodeWeights(1, i, j) * tmpBuffer(i  , j  ) + &
                      nodeWeights(2, i, j) * tmpBuffer(i+1, j  ) + &
                      nodeWeights(3, i, j) * tmpBuffer(i  , j+1) + &
                      nodeWeights(4, i, j) * tmpBuffer(i+1, j+1) 
              end do
           end do
        end select varName2
     end if

     ! Cleanup the temporary buffer
  deallocate(tmpBuffer)
end subroutine storeSurfsolInBuffer
