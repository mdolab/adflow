subroutine getAreas(areas, pts, npts, sps_in, axis)

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
  real(kind=realType), intent(in)  :: pts(3,npts), axis(3)
  real(kind=realType), intent(out) :: areas(3,npts)

  integer :: ierr
  integer(kind=intType) :: mm, nn, i, j, ipt, ii, jj,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType) :: da, fact, fact2

  areas = zero

  ! Convert to fortran numbering
  sps = sps_in+ 1

  ! Compute the local forces (or tractions). Take the scaling
  ! factor into account to obtain the forces in SI-units,
  ! i.e. Newton.
  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     if (flowDoms(nn,1_intType,sps)%rightHanded) then
        fact2 = one
     else
        fact2 = -one
     end if
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           
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
              fact = -one 
           case (iMax)
              fact = one
           case (jMin)
              fact = one
           case (jMax)
              fact = -one
           case (kMin)
              fact = -one
           case (kMax)
              fact = one
           end select
           
           ! Store the cell range of the subfaces a bit easier.
           ! As only owned faces must be considered the nodal range
           ! in BCData must be used to obtain this data.
           
           jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd
           
           ! Compute the dual area at each node. Just store in first dof
           
           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 
                 
                 ! Compute Normal
                 
                 lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                 lower_right = lower_left + 1
                 upper_left  = lower_right + iend - ibeg + 1
                 upper_right = upper_left + 1

                 call quad_area(&
                      pts(:, lower_left), pts(:, lower_right), &
                      pts(:, upper_left), pts(:, upper_right), &
                      axis, da)
                 da = fourth *da * fact * fact2

                 if (da > zero) then
                    ! Scatter to nodes
                    areas(1,lower_left)  = areas(1,lower_left)  + da
                    areas(1,lower_right) = areas(1,lower_right) + da
                    areas(1,upper_left)  = areas(1,upper_left)  + da
                    areas(1,upper_right) = areas(1,upper_right) + da
                 end if
              end do
           end do
           
           ! Note how iBeg,iBeg is defined above... it is one MORE
           ! then the starting node (used for looping over faces, not
           ! nodes)
           ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)
           
        end if
     end do bocos
  end do domains
end subroutine getAreas

subroutine getAreaSensitivity(darea, pts, npts, sps_in, axis)

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
  real(kind=realType), intent(in)  :: pts(3,npts), axis(3)
  real(kind=realType), intent(out) :: darea(3,npts)

  integer :: ierr
  integer(kind=intType) :: mm, nn, i, j, ipt, ii, jj,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType) :: area, areab, pt1b(3), pt2b(3), pt3b(3), pt4b(3)
  real(kind=realType) :: fact, fact2, da

  darea = zero

  ! Convert to fortran numbering
  sps = sps_in+ 1

  ! Compute the local forces (or tractions). Take the scaling
  ! factor into account to obtain the forces in SI-units,
  ! i.e. Newton.
  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     if (flowDoms(nn,1_intType,sps)%rightHanded) then
        fact2 = one
     else
        fact2 = -one
     end if
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           
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
              fact = -one 
           case (iMax)
              fact = one
           case (jMin)
              fact = one
           case (jMax)
              fact = -one
           case (kMin)
              fact = -one
           case (kMax)
              fact = one
           end select
           
           ! Store the cell range of the subfaces a bit easier.
           ! As only owned faces must be considered the nodal range
           ! in BCData must be used to obtain this data.
           
           jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd
           
           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 
                 
                 ! Extract 4 corner points
                 lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                 lower_right = lower_left + 1
                 upper_left  = lower_right + iend - ibeg + 1
                 upper_right = upper_left + 1

                 ! Compute actual area since we need to know to
                 ! include or not. The reverse mode calc does NOT
                 ! compute area
                 call quad_area(&
                      pts(:, lower_left), pts(:, lower_right), &
                      pts(:, upper_left), pts(:, upper_right), &
                      axis, da)

                 da = fourth *da * fact * fact2
                 if (da > zero) then

                    areab = one
                    call quad_area_b(&
                         pts(:, lower_left) , pt1b, &
                         pts(:, lower_right), pt2b, &
                         pts(:, upper_left) , pt3b, &
                         pts(:, upper_right), pt4b, &
                         axis, area, areab)

                    darea(:,lower_left)  = darea(:, lower_left)  + pt1b*fact*fact2
                    darea(:,lower_right) = darea(:, lower_right) + pt2b*fact*fact2
                    darea(:,upper_left)  = darea(:, upper_left)  + pt3b*fact*fact2
                    darea(:,upper_right) = darea(:, upper_right) + pt4b*fact*fact2
                 end if
              end do
           end do
           
           ! Note how iBeg,iBeg is defined above... it is one MORE
           ! then the starting node (used for looping over faces, not
           ! nodes)
           ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)
           
        end if
     end do bocos
  end do domains
end subroutine getAreaSensitivity

subroutine quad_area(p1, p2, p3, p4, axis, area)
  ! Kernel-level function to get area of quad defined by 4 points
  ! projected onto plane defined by axis. Only +ve areas are computed.
  use constants
  implicit none

  ! I/O
  real(kind=realType), intent(in) :: p1(3), p2(3), p3(3), p4(3), axis(3)
  real(kind=realType), intent(out) :: area

  ! Working
  real(kind=realType) :: v1(3), v2(3), sss(3)
  ! Vectors for Cross Product

  v1(:) = p4 - p1
  v2(:) = p3 - p2

  ! Cross Product
  sss(1) = half*(v1(2)*v2(3) - v1(3)*v2(2))
  sss(2) = half*(v1(3)*v2(1) - v1(1)*v2(3))
  sss(3) = half*(v1(1)*v2(2) - v1(2)*v2(1))

  area = sss(1)*axis(1)+ sss(2)*axis(2)+ sss(3)*axis(3)

end subroutine quad_area

!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!  Differentiation of quad_area in reverse (adjoint) mode:
!   gradient     of useful results: area
!   with respect to varying inputs: area p1 p2 p3 p4
!   RW status of diff variables: area:in-zero p1:out p2:out p3:out
!                p4:out
SUBROUTINE QUAD_AREA_B(p1, p1b, p2, p2b, p3, p3b, p4, p4b, axis, area, &
&  areab)
  use constants
  IMPLICIT NONE
! Kernel-level function to get area of quad defined by 4 points
! projected onto plane defined by axis. Only +ve areas are computed.
! I/O
  REAL(kind=realtype), INTENT(IN) :: p1(3), p2(3), p3(3), p4(3), axis(3)
  REAL(kind=realtype) :: p1b(3), p2b(3), p3b(3), p4b(3)
  REAL(kind=realtype) :: area
  REAL(kind=realtype) :: areab
! Working
  REAL(kind=realtype) :: v1(3), v2(3), sss(3)
  REAL(kind=realtype) :: v1b(3), v2b(3), sssb(3)
  REAL(kind=realtype) :: tempb1
  REAL(kind=realtype) :: tempb0
  INTRINSIC ABS
  REAL(kind=realtype) :: tempb
! Vectors for Cross Product
  v1(:) = p4 - p1
  v2(:) = p3 - p2
! Cross Product
  sss(1) = half*(v1(2)*v2(3)-v1(3)*v2(2))
  sss(2) = half*(v1(3)*v2(1)-v1(1)*v2(3))
  sss(3) = half*(v1(1)*v2(2)-v1(2)*v2(1))
  IF (sss(1)*axis(1) + sss(2)*axis(2) + sss(3)*axis(3) .GE. 0.) THEN
    sssb = 0.0
    sssb(1) = axis(1)*areab
    sssb(2) = axis(2)*areab
    sssb(3) = axis(3)*areab
  ELSE
    sssb = 0.0
    sssb(1) = -(axis(1)*areab)
    sssb(2) = -(axis(2)*areab)
    sssb(3) = -(axis(3)*areab)
  END IF
  v1b = 0.0
  v2b = 0.0
  tempb = half*sssb(3)
  v1b(1) = v2(2)*tempb
  v2b(2) = v1(1)*tempb
  v1b(2) = -(v2(1)*tempb)
  sssb(3) = 0.0
  tempb0 = half*sssb(2)
  v2b(1) = v1(3)*tempb0 - v1(2)*tempb
  v1b(3) = v1b(3) + v2(1)*tempb0
  v1b(1) = v1b(1) - v2(3)*tempb0
  sssb(2) = 0.0
  tempb1 = half*sssb(1)
  v2b(3) = v2b(3) + v1(2)*tempb1 - v1(1)*tempb0
  v1b(2) = v1b(2) + v2(3)*tempb1
  v1b(3) = v1b(3) - v2(2)*tempb1
  v2b(2) = v2b(2) - v1(3)*tempb1
  p2b = 0.0
  p3b = 0.0
  p3b = v2b(:)
  p2b = -v2b(:)
  p1b = 0.0
  p4b = 0.0
  p4b = v1b(:)
  p1b = -v1b(:)
  areab = 0.0
END SUBROUTINE QUAD_AREA_B
