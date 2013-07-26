subroutine getNPatches(nPatches)

  ! Return the number of wall patches we have on this proc

  use BCTypes
  use blockPointers
 
  implicit none

  ! Output Variables
  integer(kind=intType),intent(out) :: nPatches

  ! Working Variables
  integer(kind=intType) :: nn, mm

  nPatches = 0_intType

  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           nPatches = nPatches + 1
        end if
     end do bocos
  end do domains
end subroutine getNPatches

subroutine getPatchName(iPatch, patchName) 
  ! Get names one at a time since f2py can't do arrays of strings (nicely)
  use BCTypes
  use blockPointers
  use cgnsGrid

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  character*(*), intent(inout) :: patchName

  ! Working
  integer(kind=intType) :: nn, mm, patchCount, cgb

  patchCount = 0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           cgb = flowDoms(nn,1,1)%cgnsBlockID
           if (patchCount + 1 == iPatch) then
              patchName = cgnsDoms(cgb)%bocoInfo(cgnsSubface(mm))%wallBCName
           end if
           patchCount = patchCount + 1
        end if
     end do bocos
  end do domains
end subroutine getPatchName

subroutine getPatchSize(iPatch, patchSize)
  ! Get size of one patach
  use BCTypes
  use blockPointers
  use cgnsGrid

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  integer(kind=intType), intent(out) :: patchSize(2)

  ! Working
  integer(kind=intType) :: nn, mm, patchCount, iBeg, iEnd, jBeg, jEnd

  patchCount = 0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           if (patchCount + 1 == iPatch) then
              jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
              patchSize(1) = (iEnd - iBeg + 1)
              patchSize(2) = (jEnd - jBeg + 1)
           end if
           patchCount = patchCount + 1
        end if
     end do bocos
  end do domains
end subroutine getPatchSize

subroutine getForceSize(size, sizeCell)
  ! Compute the number of points that will be returned from getForces
  ! or getForcePoints
  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none

  integer(kind=intType),intent(out) :: size, sizeCell
  integer(kind=intType) :: nn,mm
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd

  size = 0_intType
  sizeCell = 0_intType
 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           size = size + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           sizeCell = sizeCell + (iEnd - iBeg)*(jEnd - jBeg)
        end if
     end do bocos
  end do domains
end subroutine getForceSize

subroutine getForceConnectivity(conn, ncell)
  ! Return the connectivity list for the each of the patches

  use BCTypes
  use blockPointers
  use inputPhysics
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: ncell
  integer(kind=intType), intent(inout) :: conn(4*ncell)

  ! Working
  integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  
  cellCount = 0
  nodeCount = 0

  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           
           ! Loop over generic face size...Note we are doing zero
           ! based ordering!
           do j=0,nj-2
              do i=0,ni-2
                 conn(4*cellCount+1) = nodeCount + (j  )*ni + i
                 conn(4*cellCount+2) = nodeCount + (j  )*ni + i + 1
                 conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1
                 conn(4*cellCount+4) = nodeCount + (j+1)*ni + i  
                 cellCount = cellCount + 1
              end do
           end do
           nodeCount = nodeCount + ni*nj
        end if
     end do bocos
  end do domains
end subroutine getForceConnectivity

subroutine getForcePoints(points, npts, sps_in)

  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts,sps_in
  real(kind=realType), intent(inout) :: points(3,npts)
  integer :: ierr

  integer(kind=intType) :: mm, nn, i, j, ii,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  sps = sps_in

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           
           ! NODE Based
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           
           do j=jBeg, jEnd ! This is a node loop
              do i=iBeg, iEnd ! This is a node loop
                 ii = ii +1
                 select case(BCFaceID(mm))
                    
                 case(imin)
                    points(:,ii) = x(1,i,j,:)
                 case(imax)
                    points(:,ii) = x(il,i,j,:)
                 case(jmin) 
                    points(:,ii) = x(i,1,j,:)
                 case(jmax) 
                    points(:,ii) = x(i,jl,j,:)
                 case(kmin) 
                    points(:,ii) = x(i,j,1,:)
                 case(kmax) 
                    points(:,ii) = x(i,j,kl,:)
                 end select    
              end do
           end do
        end if
     end do bocos
  end do domains
  
end subroutine getForcePoints

subroutine getForces(forcesP, forcesV, npts, sps_in)

  use BCTypes
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts, sps_in
  real(kind=realType), intent(out) :: forcesP(3,npts), forcesV(3, nPts)

  integer :: ierr
  real(kind=realType) :: area(npts) ! Dual area's
  integer(kind=intType) :: mm, nn, i, j, ipt, ii, jj,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: sss(3),v2(3),v1(3), qa
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yplusmax, qf(3)

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  forcesp = zero
  forcesv = zero
  area = zero
  sps = sps_in

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax)
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd
           
           ! Compute the inviscid force on each of the faces and
           ! scatter it to the 4 nodes, whose forces are updated.
           
           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 
                 
                 lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                 lower_right = lower_left + 1
                 upper_left  = lower_right + iend - ibeg + 1
                 upper_right = upper_left + 1

                 ! Assign the quarter of the pressure forces to each node
                 qf = bcData(mm)%Fp(i, j, :) * fourth

                 forcesp(:, lower_left)  = forcesp(:, lower_left)  + qf
                 forcesp(:, lower_right) = forcesp(:, lower_right) + qf
                 forcesp(:, upper_left)  = forcesp(:, upper_left)  + qf
                 forcesp(:, upper_right) = forcesp(:, upper_right) + qf

                 ! Assign the quarter of the viscous forces to each node
                 qf = bcData(mm)%Fv(i, j, :) * fourth

                 forcesv(:, lower_left)  = forcesv(:, lower_left)  + qf
                 forcesv(:, lower_right) = forcesv(:, lower_right) + qf
                 forcesv(:, upper_left)  = forcesv(:, upper_left)  + qf
                 forcesv(:, upper_right) = forcesv(:, upper_right) + qf

              end do
           end do
           
           if (forcesAsTractions) then
              jj = 1
              do j=jBeg-1, jEnd ! This is a NODE loop
                 do i=iBeg-1, iEnd ! This is a NODE loop 
                    forcesp(:, ii + jj) = forcesp(:, ii + jj) * bcData(mm)%oArea(i, j)
                    forcesv(:, ii + jj) = forcesv(:, ii + jj) * bcData(mm)%oArea(i, j)
                    jj = jj + 1
                 end do
              end do
           end if

           ! Note how iBeg,iBeg is defined above... it is one MORE
           ! then the starting node (used for looping over faces, not
           ! nodes)
           ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)
           
        end if
     end do bocos
  end do domains
end subroutine getForces
