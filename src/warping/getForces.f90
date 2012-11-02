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

subroutine getForceSize(size, nTS)
  ! Compute the number of points that will be returned from getForces
  ! or getForcePoints
  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none

  integer(kind=intType),intent(out) :: size,nTS
  integer(kind=intType) :: nn,mm
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd

  size = 0_intType
  nTS = nTimeIntervalsSpectral
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           size = size + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
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
  integer(kind=intType) :: nn,mm
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd

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
           
           ! Loop over generic face size...Note we are doing zero based ordering!
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
     if (flowDoms(nn,1_intType,sps)%rightHanded) then
        fact2 = half
     else
        fact2 = -half
     end if
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           
           viscousSubface = .true.
           if(BCType(mm) == EulerWall) viscousSubface = .false.
           
           select case (BCFaceID(mm))
              
              
              ! NOTE: The 'fact' here are NOT the same as you will
              ! find in ForcesAndMoment.f90. The reason is that, we
              ! are not using points to si, sj, sk. Those have teh
              ! normals pointing in the direction of increasing
              ! {i,j,k}. Here we are evaluating the normal from
              ! directly from the coordinates on the faces. As it
              ! happens, the normals for the jMin and jMax faces are
              ! flipped. 
           case (iMin)
              pp2 => p( 2,1:,1:); pp1 => p( 1,1:,1:)
              fact = -one 
              
           case (iMax)
              pp2 => p(il,1:,1:); pp1 => p(ie,1:,1:)
              fact = one
              
           case (jMin)
              pp2 => p(1:, 2,1:); pp1 => p(1:, 1,1:)
              fact = one
              
           case (jMax)
              pp2 => p(1:,jl,1:); pp1 => p(1:,je,1:)
              fact = -one
              
           case (kMin)
              pp2 => p(1:,1:, 2); pp1 => p(1:,1:, 1)
              fact = -one
              
           case (kMax)
              pp2 => p(1:,1:,kl); pp1 => p(1:,1:,ke)
              fact = one
              
           end select
           
           ! Store the cell range of the subfaces a bit easier.
           ! As only owned faces must be considered the nodal range
           ! in BCData must be used to obtain this data.
           
           jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd
           
           ! Compute the inviscid force on each of the faces and
           ! scatter it to the 4 nodes, whose forces are updated.
           
           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 
                 
                 ! Compute the pressure in the center of the boundary
                 ! face, which is an average between pp2 and pp1. The
                 ! value of pp is multiplied by 1/4 (the factor to
                 ! scatter to its 4 nodes, by scaleDim (to obtain the
                 ! correct dimensional value) and by fact (which takes
                 ! the possibility of inward or outward pointing normal
                 ! into account).
                 
                 pp = half*(pp2(i,j) + pp1(i,j))-Pinf
                 pp = fourth*fact*scaleDim*pp
                 
                 ! Compute Normal
                 
                 lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                 lower_right = lower_left + 1
                 upper_left  = lower_right + iend - ibeg + 1
                 upper_right = upper_left + 1

                 v1(:) = pts(:,upper_right)-pts(:,lower_left)
                 v2(:) = pts(:,upper_left)-pts(:,lower_right)
                 ! The face normal, which is the cross product of the two
                 ! diagonal vectors times fact; remember that fact2 is
                 ! either -0.5 or 0.5.
                 
                 sss(1) = fact2*(v1(2)*v2(3) - v1(3)*v2(2))
                 sss(2) = fact2*(v1(3)*v2(1) - v1(1)*v2(3))
                 sss(3) = fact2*(v1(1)*v2(2) - v1(2)*v2(1))
                 
                 ! Compute 1/4 of the area of the cell:
                 qa = fourth*sqrt(sss(1)*sss(1) + sss(2)*sss(2) + sss(3)*sss(3))
                 
                 fx = pp*sss(1)
                 fy = pp*sss(2)
                 fz = pp*sss(3)
                 
                 ! If we have viscous forces, add these:
                 if (viscousSubface) then
                    
                    ! Store the viscous stress tensor a bit easier.
                    tauXx = viscSubface(mm)%tau(i,j,1)
                    tauYy = viscSubface(mm)%tau(i,j,2)
                    tauZz = viscSubface(mm)%tau(i,j,3)
                    tauXy = viscSubface(mm)%tau(i,j,4)
                    tauXz = viscSubface(mm)%tau(i,j,5)
                    tauYz = viscSubface(mm)%tau(i,j,6)
                    
                    ! Compute the viscous force on the face. A minus sign
                    ! is now present, due to the definition of this force.
                    ! Also must multiply by scattering factor of 1/4
                    fx = fx -fact*(tauXx*sss(1) + tauXy*sss(2) &
                         +        tauXz*sss(3))*scaleDim*fourth
                    
                    fy = fy -fact*(tauXy*sss(1) + tauYy*sss(2) &
                         +        tauYz*sss(3))*scaleDim*fourth
                    
                    fz = fz -fact*(tauXz*sss(1) + tauYz*sss(2) &
                         +        tauZz*sss(3))*scaleDim*fourth
                    
                 end if
                 
                 ! Assign the quarter of the forces to each node
                 forces(1,lower_left)  = forces(1,lower_left)  + fx
                 forces(2,lower_left)  = forces(2,lower_left)  + fy
                 forces(3,lower_left)  = forces(3,lower_left)  + fz
                 
                 forces(1,lower_right) = forces(1,lower_right) + fx
                 forces(2,lower_right) = forces(2,lower_right) + fy
                 forces(3,lower_right) = forces(3,lower_right) + fz
                 
                 forces(1,upper_left)  = forces(1,upper_left)  + fx
                 forces(2,upper_left)  = forces(2,upper_left)  + fy
                 forces(3,upper_left)  = forces(3,upper_left)  + fz
                 
                 forces(1,upper_right) = forces(1,upper_right) + fx
                 forces(2,upper_right) = forces(2,upper_right) + fy
                 forces(3,upper_right) = forces(3,upper_right) + fz
                    
                 area(lower_left ) = area(lower_left ) + qa
                 area(lower_right) = area(lower_right) + qa
                 area(upper_left ) = area(upper_left ) + qa
                 area(upper_right) = area(upper_right) + qa

              end do
           end do
           
           ! Note how iBeg,iBeg is defined above... it is one MORE
           ! then the starting node (used for looping over faces, not
           ! nodes)
           ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)
           
        end if
     end do bocos
  end do domains
  
  ! If we want tractions...simply divide the forces by the dual areas
  if (forcesAsTractions) then
     do i=1,npts
        forces(:, i) = forces(:, i) / area(i)
     end do
  end if

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


