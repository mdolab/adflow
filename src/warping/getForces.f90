subroutine getForceSize(size,nTS)

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

subroutine getForces(forces,pts,npts,nTS)

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
  integer(kind=intType), intent(in) :: npts,nTS
  real(kind=realType), intent(in)  :: pts(3,npts,nTS)
  real(kind=realType), intent(out) :: forces(3,npts,nTS)

  integer :: ierr

  integer(kind=intType) :: mm, nn, i, j, ipt, ii, jj,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  real(kind=realType) :: scaleDim, fact, fact2, pp, fx,fy,fz

  real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
  real(kind=realType) :: sss(3),v2(3),v1(3)
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  logical :: viscousSubFace
  real(kind=realType) :: tauXx, tauYy, tauZz
  real(kind=realType) :: tauXy, tauXz, tauYz
  real(kind=realType) :: tmp, sss_mag
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  ! Compute the scaling factor to create the correct dimensional
  ! force in newton. As the coordinates are already in meters,
  ! this scaling factor is pRef.

  scaleDim = pRef/pInf
  forces = 0.0
  do sps = 1,nTimeIntervalsSpectral

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
                    pp = fact*scaleDim*pp

                    ! Compute Normal

                    lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                    lower_right = lower_left + 1
                    upper_left  = lower_right + iend - ibeg + 1
                    upper_right = upper_left + 1

                    v1(:) = pts(:,upper_right,sps)-pts(:,lower_left,sps)
                    v2(:) = pts(:,upper_left,sps)-pts(:,lower_right,sps)
                    ! The face normal, which is the cross product of the two
                    ! diagonal vectors times fact; remember that fact2 is
                    ! either -0.5 or 0.5.

                    sss(1) = fact2*(v1(2)*v2(3) - v1(3)*v2(2))
                    sss(2) = fact2*(v1(3)*v2(1) - v1(1)*v2(3))
                    sss(3) = fact2*(v1(1)*v2(2) - v1(2)*v2(1))

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
                            +        tauXz*sss(3))*scaleDim

                       fy = fy -fact*(tauXy*sss(1) + tauYy*sss(2) &
                            +        tauYz*sss(3))*scaleDim

                       fz = fz -fact*(tauXz*sss(1) + tauYz*sss(2) &
                            +        tauZz*sss(3))*scaleDim

                    end if
                       
                    if (.not. forcesAsTractions) then

                       ! This is the normal computation. Divide the
                       ! forces by 4 and scatter accordingly
                       fx = fx * fourth
                       fy = fy * fourth
                       fz = fz * fourth
                    else
                       
                       ! We want to get tractions (force per unit
                       ! area) and not the forces. We just have to
                       ! take fx and divide by the magnitude of ss
                       ! Note the ACTUAL tractions taken from this
                       ! routine are USELESS!!!! The nodal tractions
                       ! MUST be divided by their nodal multiplicity
                       ! (valance) in order to get the correct
                       ! values. These should ONLY be used in
                       ! conjunction with pyWarp that performs the
                       ! required multiplications. 

                       sss_mag = sqrt(sss(1)*sss(1) + &
                                      sss(2)*sss(2) + &
                                      sss(3)*sss(3))

                       fx = fx / sss_mag
                       fy = fy / sss_mag
                       fz = fz / sss_mag
                    end if

                    ! Assigning the forces is the same in either case
                    
                    forces(1,lower_left,sps)  = forces(1,lower_left,sps)  + fx
                    forces(2,lower_left,sps)  = forces(2,lower_left,sps)  + fy
                    forces(3,lower_left,sps)  = forces(3,lower_left,sps)  + fz
                    
                    forces(1,lower_right,sps) = forces(1,lower_right,sps) + fx
                    forces(2,lower_right,sps) = forces(2,lower_right,sps) + fy
                    forces(3,lower_right,sps) = forces(3,lower_right,sps) + fz
                    
                    forces(1,upper_left,sps)  = forces(1,upper_left,sps)  + fx
                    forces(2,upper_left,sps)  = forces(2,upper_left,sps)  + fy
                    forces(3,upper_left,sps)  = forces(3,upper_left,sps)  + fz
                    
                    forces(1,upper_right,sps) = forces(1,upper_right,sps) + fx
                    forces(2,upper_right,sps) = forces(2,upper_right,sps) + fy
                    forces(3,upper_right,sps) = forces(3,upper_right,sps) + fz
                       
                 end do
              end do

              ! Note how iBeg,iBeg is defined above... it is one MORE
              ! then the starting node (used for looping over faces, not
              ! nodes)
              ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)

           end if
        end do bocos
     end do domains
  end do
 
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
  real(kind=realType), intent(out) :: points(3,npts,nTS)
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

