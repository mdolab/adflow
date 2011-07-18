subroutine metric_block(nn,level,sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * metric computes the face normals and the volume for the given  *
  !      * grid level for all spectral solutions. First the volumes are   *
  !      * computed assuming that the block is right handed. Then the     *
  !      * number of positive and negative volumes are determined. If all *
  !      * volumes are positive the block is indeed right handed; if all  *
  !      * volumes are negative the block is left handed and both the     *
  !      * volumes and the normals must be negated (for the normals this  *
  !      * is done by the introduction of fact, which is either -0.5 or   *
  !      * 0.5); if there are both positive and negative volumes the mesh *
  !      * is not valid.                                                  *
  !      *                                                                *
  !      ******************************************************************
  !
  use BCTypes
  use blockPointers
  use cgnsGrid
  use communication
  use inputTimeSpectral

  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level
  !
  !      Local parameter.
  !
  real(kind=realType), parameter :: thresVolume = 1.e-2_realType
  !
  !      Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: i, j, k, n, m, l
  integer(kind=intType) :: nn, mm, sps
  integer(kind=intType) :: nVolNeg,   nVolPos
  integer(kind=intType) :: nVolBad,   nVolBadGlobal
  integer(kind=intType) :: nBlockBad, nBlockBadGlobal

  real(kind=realType) :: fact, mult
  real(kind=realType) :: xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6

  real(kind=realType), dimension(3) :: v1, v2

  real(kind=realType), dimension(:,:,:), pointer :: ss

  character(len=10) :: integerString

  logical :: checkK, checkJ, checkI, checkAll
  logical :: badVolume

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
 
  ! Compute the volumes. The hexahedron is split into 6 pyramids
  ! whose volumes are computed. The volume is positive for a
  ! right handed block.
  ! Initialize the volumes to zero. The reasons is that the second
  ! level halo's must be initialized to zero and for convenience
  ! all the volumes are set to zero.

  vol = zero

  do k=1,ke
     n = k -1

     checkK = .true.
     if(k == 1 .or. k == ke) checkK = .false.

     do j=1,je
        m = j -1

        checkJ = .true.
        if(j == 1 .or. j == je) checkJ = .false.

        do i=1,ie
           l = i -1

           checkI = .true.
           if(i == 1 .or. i == ie) checkI = .false.

           ! Determine whether or not the voluem must be checked for
           ! quality. Only owned volumes are checked, not halo's.

           checkAll = .false.
           if(checkK .and. checkJ .and. checkI) checkAll = .true.

           ! Compute the coordinates of the center of gravity.

!!$           xp = eighth*(x(i,j,k,1) + x(i,m,k,1) &
!!$                +         x(i,m,n,1) + x(i,j,n,1) &
!!$                +         x(l,j,k,1) + x(l,m,k,1) &
!!$                +         x(l,m,n,1) + x(l,j,n,1))
!!$           yp = eighth*(x(i,j,k,2) + x(i,m,k,2) &
!!$                +         x(i,m,n,2) + x(i,j,n,2) &
!!$                +         x(l,j,k,2) + x(l,m,k,2) &
!!$                +         x(l,m,n,2) + x(l,j,n,2))
!!$           zp = eighth*(x(i,j,k,3) + x(i,m,k,3) &
!!$                +         x(i,m,n,3) + x(i,j,n,3) &
!!$                +         x(l,j,k,3) + x(l,m,k,3) &
!!$                +         x(l,m,n,3) + x(l,j,n,3))


           ! Compute the volumes of the 6 sub pyramids. The
           ! arguments of volpym must be such that for a (regular)
           ! right handed hexahedron all volumes are positive.

           vp1 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                x(i,m,k,1), x(i,m,k,2), x(i,m,k,3))

           vp2 = volpym(x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                x(l,j,n,1), x(l,j,n,2), x(l,j,n,3))

           vp3 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                x(i,j,n,1), x(i,j,n,2), x(i,j,n,3))

           vp4 = volpym(x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                x(l,m,k,1), x(l,m,k,2), x(l,m,k,3))

           vp5 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                x(l,j,k,1), x(l,j,k,2), x(l,j,k,3))

           vp6 = volpym(x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                x(i,m,n,1), x(i,m,n,2), x(i,m,n,3))

           ! Set the volume to 1/6 of the sum of the volumes of the
           ! pyramid. Remember that volpym computes 6 times the
           ! volume.

           vol(i,j,k) = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)

           ! Check the volume and update the number of positive
           ! and negative volumes if needed.

           ! Set the volume to the absolute value.

           vol(i,j,k) = abs(vol(i,j,k))

        enddo
     enddo
  enddo

!!$  ! Some additional safety stuff for halo volumes.
!!$
!!$  do k=2,kl
!!$     do j=2,jl
!!$        if(vol(1, j,k) <= eps) vol(1, j,k) = vol(2, j,k)
!!$        if(vol(ie,j,k) <= eps) vol(ie,j,k) = vol(il,j,k)
!!$     enddo
!!$  enddo
!!$
!!$  do k=2,kl
!!$     do i=1,ie
!!$        if(vol(i,1, k) <= eps) vol(i,1, k) = vol(i,2, k)
!!$        if(vol(i,je,k) <= eps) vol(i,je,k) = vol(i,jl,k)
!!$     enddo
!!$  enddo
!!$
!!$  do j=1,je
!!$     do i=1,ie
!!$        if(vol(i,j,1)  <= eps) vol(i,j,1)  = vol(i,j,2)
!!$        if(vol(i,j,ke) <= eps) vol(i,j,ke) = vol(i,j,kl)
!!$     enddo
!!$  enddo
!!$
!!$  ! Set the factor in the surface normals computation. For a
!!$  ! left handed block this factor is negative, such that the
!!$  ! normals still point in the direction of increasing index.
!!$  ! The formulae used later on assume a right handed block
!!$  ! and fact is used to correct this for a left handed block,
!!$  ! as well as the scaling factor of 0.5
!!$
!!$  if (rightHanded) then!( flowDoms(nn,level,sps)%rightHanded ) then
!!$     fact =  half
!!$  else
!!$     fact = -half
!!$  endif
!!$
!!$  ! Check if both positive and negative volumes occur. If so,
!!$  ! the block is bad and the counter nBlockBad is updated.
!!$
!!$  !
!!$  !          **************************************************************
!!$  !          *                                                            *
!!$  !          * Computation of the face normals in i-, j- and k-direction. *
!!$  !          * Formula's are valid for a right handed block; for a left   *
!!$  !          * handed block the correct orientation is obtained via fact. *
!!$  !          * The normals point in the direction of increasing index.    *
!!$  !          * The absolute value of fact is 0.5, because the cross       *
!!$  !          * product of the two diagonals is twice the normal vector.   *
!!$  !          *                                                            *
!!$  !          * Note that also the normals of the first level halo cells   *
!!$  !          * are computed. These are needed for the viscous fluxes.     *
!!$  !          *                                                            *
!!$  !          **************************************************************
!!$  !
!!$  ! Projected areas of cell faces in the i direction.
!!$
!!$  do k=1,ke
!!$     n = k -1
!!$     do j=1,je
!!$        m = j -1
!!$        do i=0,ie
!!$
!!$           ! Determine the two diagonal vectors of the face.
!!$
!!$           v1(1) = x(i,j,n,1) - x(i,m,k,1)
!!$           v1(2) = x(i,j,n,2) - x(i,m,k,2)
!!$           v1(3) = x(i,j,n,3) - x(i,m,k,3)
!!$
!!$           v2(1) = x(i,j,k,1) - x(i,m,n,1)
!!$           v2(2) = x(i,j,k,2) - x(i,m,n,2)
!!$           v2(3) = x(i,j,k,3) - x(i,m,n,3)
!!$
!!$           ! The face normal, which is the cross product of the two
!!$           ! diagonal vectors times fact; remember that fact is
!!$           ! either -0.5 or 0.5.
!!$        
!!$           si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
!!$           si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
!!$           si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
!!$         
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  ! Projected areas of cell faces in the j direction.
!!$
!!$  do k=1,ke
!!$     n = k -1
!!$     do j=0,je
!!$        do i=1,ie
!!$           l = i -1
!!$
!!$           ! Determine the two diagonal vectors of the face.
!!$
!!$           v1(1) = x(i,j,n,1) - x(l,j,k,1)
!!$           v1(2) = x(i,j,n,2) - x(l,j,k,2)
!!$           v1(3) = x(i,j,n,3) - x(l,j,k,3)
!!$
!!$           v2(1) = x(l,j,n,1) - x(i,j,k,1)
!!$           v2(2) = x(l,j,n,2) - x(i,j,k,2)
!!$           v2(3) = x(l,j,n,3) - x(i,j,k,3)
!!$
!!$           ! The face normal, which is the cross product of the two
!!$           ! diagonal vectors times fact; remember that fact is
!!$           ! either -0.5 or 0.5.
!!$
!!$           sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
!!$           sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
!!$           sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
!!$
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  ! Projected areas of cell faces in the k direction.
!!$
!!$  do k=0,ke
!!$     do j=1,je
!!$        m = j -1
!!$        do i=1,ie
!!$           l = i -1
!!$
!!$           ! Determine the two diagonal vectors of the face.
!!$
!!$           v1(1) = x(i,j,k,1) - x(l,m,k,1)
!!$           v1(2) = x(i,j,k,2) - x(l,m,k,2)
!!$           v1(3) = x(i,j,k,3) - x(l,m,k,3)
!!$
!!$           v2(1) = x(l,j,k,1) - x(i,m,k,1)
!!$           v2(2) = x(l,j,k,2) - x(i,m,k,2)
!!$           v2(3) = x(l,j,k,3) - x(i,m,k,3)
!!$
!!$           ! The face normal, which is the cross product of the two
!!$           ! diagonal vectors times fact; remember that fact is
!!$           ! either -0.5 or 0.5.
!!$
!!$           sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
!!$           sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
!!$           sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
!!$
!!$        enddo
!!$     enddo
!!$  enddo
!!$  !
!!$  !          **************************************************************
!!$  !          *                                                            *
!!$  !          * The unit normals on the boundary faces. These always point *
!!$  !          * out of the domain, so a multiplication by -1 is needed for *
!!$  !          * the iMin, jMin and kMin boundaries.                        *
!!$  !          *                                                            *
!!$  !          **************************************************************
!!$  !
!!$  ! Loop over the boundary subfaces of this block.
!!$
!!$  bocoLoop: do mm=1,nBocos
!!$
!!$     ! Determine the block face on which this subface is located
!!$     ! and set ss and mult accordingly.
!!$
!!$     select case (BCFaceID(mm))
!!$
!!$     case (iMin)
!!$        mult = -one; ss => si(1,:,:,:)
!!$
!!$     case (iMax)
!!$        mult = one;  ss => si(il,:,:,:)
!!$
!!$     case (jMin)
!!$        mult = -one; ss => sj(:,1,:,:)
!!$
!!$     case (jMax)
!!$        mult = one;  ss => sj(:,jl,:,:)
!!$
!!$     case (kMin)
!!$        mult = -one; ss => sk(:,:,1,:)
!!$
!!$     case (kMax)
!!$        mult = one;  ss => sk(:,:,kl,:)
!!$
!!$     end select
!!$
!!$     ! Loop over the boundary faces of the subface.
!!$
!!$     do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
!!$        do i=BCData(mm)%icBeg, BCData(mm)%icEnd
!!$
!!$           ! Compute the inverse of the length of the normal vector
!!$           ! and possibly correct for inward pointing.
!!$
!!$           xp = ss(i,j,1);  yp = ss(i,j,2);  zp = ss(i,j,3)
!!$           fact = sqrt(xp*xp + yp*yp + zp*zp)
!!$           if(fact > zero) fact = mult/fact
!!$
!!$           ! Compute the unit normal.
!!$
!!$           BCData(mm)%norm(i,j,1) = fact*xp
!!$           BCData(mm)%norm(i,j,2) = fact*yp
!!$           BCData(mm)%norm(i,j,3) = fact*zp
!!$
!!$        enddo
!!$     enddo
!!$
!!$  enddo bocoLoop


contains

  !        ================================================================

  function volpym(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd)
    !
    !        ****************************************************************
    !        *                                                              *
    !        * volpym computes 6 times the volume of a pyramid. Node p,     *
    !        * whose coordinates are set in the subroutine metric itself,   *
    !        * is the top node and a-b-c-d is the quadrilateral surface.    *
    !        * It is assumed that the cross product vCa * vDb points in     *
    !        * the direction of the top node. Here vCa is the diagonal      *
    !        * running from node c to node a and vDb the diagonal from      *
    !        * node d to node b.                                            *
    !        *                                                              *
    !        ****************************************************************
    !
    use precision
    implicit none
    !
    !        Function type.
    !
    real(kind=realType) :: volpym
    !
    !        Function arguments.
    !
    real(kind=realType), intent(in) :: xa, ya, za, xb, yb, zb
    real(kind=realType), intent(in) :: xc, yc, zc, xd, yd, zd
    !
    !        ****************************************************************
    !        *                                                              *
    !        * Begin execution                                              *
    !        *                                                              *
    !        ****************************************************************
    !
    volpym = (xp - fourth*(xa + xb  + xc + xd))              &
         * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
         (yp - fourth*(ya + yb  + yc + yd))              &
         * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
         (zp - fourth*(za + zb  + zc + zd))              &
         * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

  end function volpym

end subroutine metric_block
