subroutine getNPatches(nPatches)

  ! Return the number of wall patches we have on this proc

  use BCTypes
  use blockPointers
 
  implicit none

  ! Output Variables
  integer(kind=intType),intent(out) :: nPatches

  ! Working Variables
  integer(kind=intType) :: nn, mm, cgb

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
  ! Get names one at a time since f2py can't do arrays of strings
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
           if (patchCount == iPatch) then
              patchName = cgnsDoms(cgb)%bocoInfo(cgnsSubface(mm))%wallBCName
           end if
           patchCount = patchCount + 1
        end if
     end do bocos
  end do domains
end subroutine getPatchName

subroutine getPatchSize(iPatch, patchSize)
  ! Get names one at a time since f2py can't do arrays of strings
  use BCTypes
  use blockPointers
  use cgnsGrid

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  integer(kind=intType), intent(out) :: patchSize(2)

  ! Working
  integer(kind=intType) :: nn, mm, patchCount, cgb, iBeg, iEnd, jBeg, jEnd

  patchCount = 0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           cgb = flowDoms(nn,1,1)%cgnsBlockID
           if (patchCount == iPatch) then
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

subroutine getForceSize(size, sizeCell, nTS)
  ! Compute the number of points that will be returned from getForces
  ! or getForcePoints
  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none

  integer(kind=intType),intent(out) :: size, sizeCell, nTS
  integer(kind=intType) :: nn,mm
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd

  size = 0_intType
  sizeCell = 0_intType
  nTS = nTimeIntervalsSpectral
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

subroutine getForceConnectivitySize(size)
  ! Compute the number of cellss that will be returned from
  ! getForceConnectivity
  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none
  
  ! Input/Output
  integer(kind=intType), intent(out) :: size

  ! Working
  integer(kind=intType) :: nn, mm
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  size = 0_intType
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           size = size + (iEnd - iBeg)*(jEnd - jBeg)
        end if
     end do bocos
  end do domains
end subroutine getForceConnectivitySize

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

subroutine getForces(forces, pts, npts, sps_in)

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
  real(kind=realType), intent(in)  :: pts(3,npts)
  real(kind=realType), intent(out) :: forces(3,npts)

  integer :: ierr
  real(kind=realType) :: area(npts) ! Dual area's

  integer(kind=intType) :: mm, nn, i, j, ipt, ii, jj,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: scaleDim, fact, fact2, pp, fx,fy,fz

  real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
  real(kind=realType) :: sss(3),v2(3),v1(3), qa
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  logical :: viscousSubFace
  real(kind=realType) :: tauXx, tauYy, tauZz
  real(kind=realType) :: tauXy, tauXz, tauYz
  real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yplusmax, qf(3)
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  ! Compute the scaling factor to create the correct dimensional
  ! force in newton. As the coordinates are already in meters,
  ! this scaling factor is pRef.

  scaleDim = pRef/pInf
  forces = zero
  area = zero

  ! Convert to fortran numbering
  sps = sps_in+ 1

  ! Compute the local forces (or tractions). Take the scaling
  ! factor into account to obtain the forces in SI-units,
  ! i.e. Newton.
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

                 ! Assign the quarter of the forces to each node
                 qf = bcData(mm)%F(i, j, :)*fourth

                 forces(:, lower_left)  = forces(:, lower_left)  + qf
                 forces(:, lower_right) = forces(:, lower_right) + qf
                 forces(:, upper_left)  = forces(:, upper_left)  + qf
                 forces(:, upper_right) = forces(:, upper_right) + qf
              end do
           end do
           
           if (forcesAsTractions) then
              jj = 1
              do j=jBeg-1, jEnd ! This is a NODE loop
                 do i=iBeg-1, iEnd ! This is a NODE loop 
                    forces(:, ii + jj) = forces(:, ii + jj) * bcData(mm)%oArea(i, j)
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

subroutine getForcePoints(points,npts,nTS)

  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts,nTS
  real(kind=realType), intent(inout) :: points(3,npts,nTS)
  integer :: ierr

  integer(kind=intType) :: mm, nn, i, j, ii,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  do sps = 1,nTimeIntervalsSpectral
     ii = 0 
     domains: do nn=1,nDom
        call setPointers(nn,1_intType,sps)

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
                       points(:,ii,sps) = x(1,i,j,:)
                    case(imax)
                       points(:,ii,sps) = x(il,i,j,:)
                    case(jmin) 
                       points(:,ii,sps) = x(i,1,j,:)
                    case(jmax) 
                       points(:,ii,sps) = x(i,jl,j,:)
                    case(kmin) 
                       points(:,ii,sps) = x(i,j,1,:)
                    case(kmax) 
                       points(:,ii,sps) = x(i,j,kl,:)
                    end select
                 end do
              end do
           end if
        end do bocos
     end do domains
  enddo
end subroutine getForcePoints


