subroutine volume_block
  
  ! This is COPY of metric.f90. It was necessary to copy this file
  ! since there is debugging stuff in the original that is not
  ! necessary for AD.

  use BCTypes
  use blockPointers
  use cgnsGrid
  use communication
  use inputTimeSpectral

  implicit none
  !
  !      Local parameter.
  !
  real(kind=realType), parameter :: thresVolume = 1.e-2_realType
  real(kind=realType), parameter :: haloCellRatio = 1e-10_realType
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, n, m, l, ii
  integer(kind=intType) :: mm
  real(kind=realType) :: fact, mult
  real(kind=realType) :: xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6
  real(kind=realType) :: xxp, yyp, zzp
  real(kind=realType), dimension(3) :: v1, v2

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

  do k=1, ke
     n = k - 1
     do j=1, je
        m = j - 1
        do i=1, ie

              l = i - 1
              
              ! Compute the coordinates of the center of gravity.

              xp = eighth*(x(i,j,k,1) + x(i,m,k,1) &
                   +         x(i,m,n,1) + x(i,j,n,1) &
                   +         x(l,j,k,1) + x(l,m,k,1) &
                   +         x(l,m,n,1) + x(l,j,n,1))
              yp = eighth*(x(i,j,k,2) + x(i,m,k,2) &
                   +         x(i,m,n,2) + x(i,j,n,2) &
                   +         x(l,j,k,2) + x(l,m,k,2) &
                +         x(l,m,n,2) + x(l,j,n,2))
              zp = eighth*(x(i,j,k,3) + x(i,m,k,3) &
                   +         x(i,m,n,3) + x(i,j,n,3) &
                   +         x(l,j,k,3) + x(l,m,k,3) &
                   +         x(l,m,n,3) + x(l,j,n,3))
              
              
              ! Compute the volumes of the 6 sub pyramids. The
              ! arguments of volpym must be such that for a (regular)
              ! right handed hexahedron all volumes are positive.
              
              call volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                   x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                   x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                   x(i,m,k,1), x(i,m,k,2), x(i,m,k,3),vp1)
              
              call volpym(x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                   x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                   x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                   x(l,j,n,1), x(l,j,n,2), x(l,j,n,3),vp2)
              
              call volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                   x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                   x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                   x(i,j,n,1), x(i,j,n,2), x(i,j,n,3),vp3)
              
              call volpym(x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                   x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                   x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                   x(l,m,k,1), x(l,m,k,2), x(l,m,k,3),vp4)
              
              call volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                   x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                   x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                   x(l,j,k,1), x(l,j,k,2), x(l,j,k,3),vp5)
              
              call volpym(x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                   x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                   x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                   x(i,m,n,1), x(i,m,n,2), x(i,m,n,3),vp6)
              
              ! Set the volume to 1/6 of the sum of the volumes of the
              ! pyramid. Remember that volpym computes 6 times the
              ! volume.
              
              vol(i,j,k) = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)
              
              ! Set the volume to the absolute value.
              vol(i, j, k) = abs(vol(i, j, k))
        enddo
     enddo
  enddo

  ! Some additional safety stuff for halo volumes.
  
  do k=2,kl
     do j=2,jl
        if(vol(1, j,k)/vol(2, j, k) < haloCellRatio) then 
           vol(1, j,k) = vol(2, j,k)
        end if
        if(vol(ie,j,k)/vol(il,j,k)  < haloCellRatio) then 
           vol(ie,j,k) = vol(il,j,k)
        end if
     enddo
  enddo

  do k=2,kl
     do i=1,ie
        if(vol(i,1, k)/vol(i,2,k) < haloCellRatio) then 
           vol(i,1, k) = vol(i,2, k)
        end if
        if(vol(i,je,k)/voL(i,jl,k) < haloCellRatio) then 
           vol(i,je,k) = vol(i,jl,k)
        end if
     enddo
  enddo
  
  do j=1,je
     do i=1,ie
        if(vol(i,j,1)/vol(i,j,2)  < haloCellRatio) then 
           vol(i,j,1)  = vol(i,j,2)
        end if
        if(vol(i,j,ke)/vol(i,j,kl) < haloCellRatio) then 
           vol(i,j,ke) = vol(i,j,kl)
        end if
     enddo
  enddo


contains

  subroutine volpym(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,volume)
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
    real(kind=realType) :: volume
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
    volume = (xp - fourth*(xa + xb  + xc + xd))              &
         * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
         (yp - fourth*(ya + yb  + yc + yd))              &
         * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
         (zp - fourth*(za + zb  + zc + zd))              &
         * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

  end subroutine volpym
end subroutine volume_block

subroutine metric_block
  
  use blockPointers
  implicit none

  ! Local variables.
  integer(kind=intType) :: i, j, k, n, m, l, ii
  real(kind=realType) :: fact
  real(kind=realType) :: xxp, yyp, zzp
  real(kind=realType), dimension(3) :: v1, v2

  ! Set the factor in the surface normals computation. For a
  ! left handed block this factor is negative, such that the
  ! normals still point in the direction of increasing index.
  ! The formulae used later on assume a right handed block
  ! and fact is used to correct this for a left handed block,
  ! as well as the scaling factor of 0.5

  if (rightHanded) then
     fact =  half
  else
     fact = -half
  endif

  !
  ! **************************************************************
  ! *                                                            *
  ! * Computation of the face normals in i-, j- and k-direction. *
  ! * Formula's are valid for a right handed block; for a left   *
  ! * handed block the correct orientation is obtained via fact. *
  ! * The normals point in the direction of increasing index.    *
  ! * The absolute value of fact is 0.5, because the cross       *
  ! * product of the two diagonals is twice the normal vector.   *
  ! *                                                            *
  ! * Note that also the normals of the first level halo cells   *
  ! * are computed. These are needed for the viscous fluxes.     *
  ! *                                                            *
  ! **************************************************************
  !
  ! Projected areas of cell faces in the i direction.
  !$AD II-LOOP
  do ii=0,ke*je*(ie+1)-1
     i = mod(ii, ie+1) + 0 ! 0:ie
     j = mod(ii/(ie+1), je) + 1 !1:je
     k = ii/((ie+1)*je) + 1 !1:ke

     n = k -1
     m = j -1

     ! Determine the two diagonal vectors of the face.

     v1(1) = x(i,j,n,1) - x(i,m,k,1)
     v1(2) = x(i,j,n,2) - x(i,m,k,2)
     v1(3) = x(i,j,n,3) - x(i,m,k,3)
     
     v2(1) = x(i,j,k,1) - x(i,m,n,1)
     v2(2) = x(i,j,k,2) - x(i,m,n,2)
     v2(3) = x(i,j,k,3) - x(i,m,n,3)

     ! The face normal, which is the cross product of the two
     ! diagonal vectors times fact; remember that fact is
     ! either -0.5 or 0.5.
     
     si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
     si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
     si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
  enddo
 

  ! Projected areas of cell faces in the j direction
  !$AD II-LOOP
  do ii=0,ke*(je+1)*ie-1
     i = mod(ii, ie) + 1 ! 1:ie
     j = mod(ii/ie, je+1) + 0 !0:je
     k = ii/(ie*(je+1)) + 1 !1:ke
     n = k -1
     l = i -1

     ! Determine the two diagonal vectors of the face.

     v1(1) = x(i,j,n,1) - x(l,j,k,1)
     v1(2) = x(i,j,n,2) - x(l,j,k,2)
     v1(3) = x(i,j,n,3) - x(l,j,k,3)
     
     v2(1) = x(l,j,n,1) - x(i,j,k,1)
     v2(2) = x(l,j,n,2) - x(i,j,k,2)
     v2(3) = x(l,j,n,3) - x(i,j,k,3)

     ! The face normal, which is the cross product of the two
     ! diagonal vectors times fact; remember that fact is
     ! either -0.5 or 0.5.
     
     sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
     sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
     sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
     
  enddo

  ! Projected areas of cell faces in the k direction.
  !$AD II-LOOP
  do ii=0,(ke+1)*je*ie-1
     i = mod(ii, ie) + 1 ! 1:ie
     j = mod(ii/ie, je) + 1 !1:je
     k = ii/(ie*je) + 0 !0:ke
     m = j -1
     l = i -1

     ! Determine the two diagonal vectors of the face.

     v1(1) = x(i,j,k,1) - x(l,m,k,1)
     v1(2) = x(i,j,k,2) - x(l,m,k,2)
     v1(3) = x(i,j,k,3) - x(l,m,k,3)
     
     v2(1) = x(l,j,k,1) - x(i,m,k,1)
     v2(2) = x(l,j,k,2) - x(i,m,k,2)
     v2(3) = x(l,j,k,3) - x(i,m,k,3)
     
     ! The face normal, which is the cross product of the two
     ! diagonal vectors times fact; remember that fact is
     ! either -0.5 or 0.5.

     sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
     sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
     sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))
  enddo 
end subroutine metric_block

subroutine boundaryNormals

  ! **************************************************************
  ! *                                                            *
  ! * The unit normals on the boundary faces. These always point *
  ! * out of the domain, so a multiplication by -1 is needed for *
  ! * the iMin, jMin and kMin boundaries.                        *
  ! *                                                            *
  ! **************************************************************
  !
  use BCTypes
  use blockPointers
  use cgnsGrid
  use communication
  use inputTimeSpectral
  implicit none

  ! Local variables.
  integer(kind=intType) :: i, j, ii
  integer(kind=intType) :: mm
  real(kind=realType) :: fact, mult
  real(kind=realType) :: xxp, yyp, zzp

  !Loop over the boundary subfaces of this block.
  !$AD II-LOOP
  bocoLoop: do mm=1,nBocos
    
     ! Loop over the boundary faces of the subface.
     !$AD II-LOOP
     do ii=0,(BCData(mm)%jcEnd - bcData(mm)%jcBeg + 1)*(BCData(mm)%icEnd - BCData(mm)%icBeg + 1) - 1 
        i = mod(ii, (BCData(mm)%icEnd - BCData(mm)%icBeg + 1)) + BCData(mm)%icBeg
        j = ii/(BCData(mm)%icEnd - BCData(mm)%icBeg + 1) + BCData(mm)%jcBeg
        
        select case (BCFaceID(mm))
        case (iMin)
           mult = -one
           xxp = si(1,i,j,1);  yyp = si(1,i,j,2);  zzp = si(1,i,j,3)
        case (iMax)
           mult = one
           xxp = si(il,i,j,1);  yyp = si(il,i,j,2);  zzp = si(il,i,j,3)
        case (jMin)
           mult = -one
           xxp = sj(i,1,j,1);  yyp = sj(i,1,j,2);  zzp = sj(i,1,j,3)
        case (jMax)
           mult = one
           xxp = sj(i,jl,j,1);  yyp = sj(i,jl,j,2);  zzp = sj(i,jl,j,3)
        case (kMin)
           mult = -one
           xxp = sk(i,j,1,1);  yyp = sk(i,j,1,2);  zzp = sk(i,j,1,3)
        case (kMax)
           mult = one
           xxp = sk(i,j,kl,1);  yyp = sk(i,j,kl,2);  zzp = sk(i,j,kl,3)
        end select
        
        ! Compute the inverse of the length of the normal vector
        ! and possibly correct for inward pointing.
        
        fact = sqrt(xxp*xxp + yyp*yyp + zzp*zzp)
        if(fact > zero) fact = mult/fact
        
        ! Compute the unit normal.
        
        BCData(mm)%norm(i,j,1) = fact*xxp
        BCData(mm)%norm(i,j,2) = fact*yyp
        BCData(mm)%norm(i,j,3) = fact*zzp
     end do
  enddo bocoLoop
end subroutine boundaryNormals
