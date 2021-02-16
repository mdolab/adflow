module adjointExtra

contains

  subroutine volume_block

    ! This is COPY of metric.f90. It was necessary to copy this file
    ! since there is debugging stuff in the original that is not
    ! necessary for AD.
    use constants
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
      !         volpym computes 6 times the volume of a pyramid. Node p,
      !         whose coordinates are set in the subroutine metric itself,
      !         is the top node and a-b-c-d is the quadrilateral surface.
      !         It is assumed that the cross product vCa * vDb points in
      !         the direction of the top node. Here vCa is the diagonal
      !         running from node c to node a and vDb the diagonal from
      !         node d to node b.
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

      volume = (xp - fourth*(xa + xb  + xc + xd))              &
           * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
           (yp - fourth*(ya + yb  + yc + yd))              &
           * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
           (zp - fourth*(za + zb  + zc + zd))              &
           * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

    end subroutine volpym
  end subroutine volume_block

  subroutine metric_block
    use constants
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
    !  Computation of the face normals in i-, j- and k-direction.
    !  Formula's are valid for a right handed block; for a left
    !  handed block the correct orientation is obtained via fact.
    !  The normals point in the direction of increasing index.
    !  The absolute value of fact is 0.5, because the cross
    !  product of the two diagonals is twice the normal vector.
    !  Note that also the normals of the first level halo cells
    !  are computed. These are needed for the viscous fluxes.
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

    !  The unit normals on the boundary faces. These always point
    !  out of the domain, so a multiplication by -1 is needed for
    !  the iMin, jMin and kMin boundaries.
    !
    use constants
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

  subroutine xhalo_block
    !
    !       xhalo determines the coordinates of the nodal halo's.
    !       First it sets all halo coordinates by simple extrapolation,
    !       then the symmetry planes are treated (also the unit normal of
    !       symmetry planes are determined) and finally an exchange is
    !       made for the internal halo's.
    !
    use constants
    use blockPointers
    use communication
    use inputTimeSpectral
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: mm, i, j, k
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iiMax, jjMax
    logical err
    real(kind=realType) :: length, dot
    real(kind=realType), dimension(3) :: v1, v2, norm

    ! Extrapolation in i-direction.

    do k=1,kl
       do j=1,jl
          x(0,j,k,1) = two*x(1,j,k,1) - x(2,j,k,1)
          x(0,j,k,2) = two*x(1,j,k,2) - x(2,j,k,2)
          x(0,j,k,3) = two*x(1,j,k,3) - x(2,j,k,3)

          x(ie,j,k,1) = two*x(il,j,k,1) - x(nx,j,k,1)
          x(ie,j,k,2) = two*x(il,j,k,2) - x(nx,j,k,2)
          x(ie,j,k,3) = two*x(il,j,k,3) - x(nx,j,k,3)
       enddo
    enddo

    ! Extrapolation in j-direction.

    do k=1,kl
       do i=0,ie
          x(i,0,k,1) = two*x(i,1,k,1) - x(i,2,k,1)
          x(i,0,k,2) = two*x(i,1,k,2) - x(i,2,k,2)
          x(i,0,k,3) = two*x(i,1,k,3) - x(i,2,k,3)

          x(i,je,k,1) = two*x(i,jl,k,1) - x(i,ny,k,1)
          x(i,je,k,2) = two*x(i,jl,k,2) - x(i,ny,k,2)
          x(i,je,k,3) = two*x(i,jl,k,3) - x(i,ny,k,3)
       enddo
    enddo

    ! Extrapolation in k-direction.

    do j=0,je
       do i=0,ie
          x(i,j,0,1) = two*x(i,j,1,1) - x(i,j,2,1)
          x(i,j,0,2) = two*x(i,j,1,2) - x(i,j,2,2)
          x(i,j,0,3) = two*x(i,j,1,3) - x(i,j,2,3)

          x(i,j,ke,1) = two*x(i,j,kl,1) - x(i,j,nz,1)
          x(i,j,ke,2) = two*x(i,j,kl,2) - x(i,j,nz,2)
          x(i,j,ke,3) = two*x(i,j,kl,3) - x(i,j,nz,3)
       enddo
    enddo
    !
    !           Mirror the halo coordinates adjacent to the symmetry
    !           planes
    !
    ! Loop over boundary subfaces.

    loopBocos: do mm=1,nBocos

       ! The actual correction of the coordinates only takes
       ! place for symmetry planes.

       testSymmetry: if(BCType(mm) == Symm) then

          ! Set some variables, depending on the block face on
          ! which the subface is located.
          norm(1) = bcData(mm)%symNorm(1)
          norm(2) = bcData(mm)%symNorm(2)
          norm(3) = bcData(mm)%symNorm(3)

          length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)

          ! Compute the unit normal of the subface.

          norm(1) = norm(1)/length
          norm(2) = norm(2)/length
          norm(3) = norm(3)/length
          ! See xhalo_block for comments for below:
          testSingular: if(length > eps) then

             select case (BCFaceID(mm))
             case (iMin)
                iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl

                if(iBeg == 1)     iBeg = 0
                if(iEnd == iiMax) iEnd = iiMax + 1

                if(jBeg == 1)     jBeg = 0
                if(jEnd == jjMax) jEnd = jjMax + 1

                do j=jBeg,jEnd
                   do i=iBeg,iEnd
                      v1(1) = x(1, i,j,1) - x(2, i,j,1)
                      v1(2) = x(1, i,j,2) - x(2, i,j,2)
                      v1(3) = x(1, i,j,3) - x(2, i,j,3)
                      dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                           +      v1(3)*norm(3))
                      x(0,i,j,1) = x(2,i,j,1) + dot*norm(1)
                      x(0,i,j,2) = x(2,i,j,2) + dot*norm(2)
                      x(0,i,j,3) = x(2,i,j,3) + dot*norm(3)
                   enddo
                enddo

             case (iMax)
                iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl

                if(iBeg == 1)     iBeg = 0
                if(iEnd == iiMax) iEnd = iiMax + 1

                if(jBeg == 1)     jBeg = 0
                if(jEnd == jjMax) jEnd = jjMax + 1

                do j=jBeg,jEnd
                   do i=iBeg,iEnd
                      v1(1) = x(il, i,j,1) - x(nx, i,j,1)
                      v1(2) = x(il, i,j,2) - x(nx, i,j,2)
                      v1(3) = x(il, i,j,3) - x(nx, i,j,3)
                      dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                           +      v1(3)*norm(3))
                      x(ie,i,j,1) = x(nx,i,j,1) + dot*norm(1)
                      x(ie,i,j,2) = x(nx,i,j,2) + dot*norm(2)
                      x(ie,i,j,3) = x(nx,i,j,3) + dot*norm(3)
                   enddo
                enddo

             case (jMin)
                iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl

                if(iBeg == 1)     iBeg = 0
                if(iEnd == iiMax) iEnd = iiMax + 1

                if(jBeg == 1)     jBeg = 0
                if(jEnd == jjMax) jEnd = jjMax + 1

                do j=jBeg,jEnd
                   do i=iBeg,iEnd
                      v1(1) = x(i,1,j,1) - x(i,2,j,1)
                      v1(2) = x(i,1,j,2) - x(i,2,j,2)
                      v1(3) = x(i,1,j,3) - x(i,2,j,3)
                      dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                           +      v1(3)*norm(3))
                      x(i,0,j,1) = x(i,2,j,1) + dot*norm(1)
                      x(i,0,j,2) = x(i,2,j,2) + dot*norm(2)
                      x(i,0,j,3) = x(i,2,j,3) + dot*norm(3)
                   enddo
                enddo

             case (jMax)
                iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl

                if(iBeg == 1)     iBeg = 0
                if(iEnd == iiMax) iEnd = iiMax + 1

                if(jBeg == 1)     jBeg = 0
                if(jEnd == jjMax) jEnd = jjMax + 1

                do j=jBeg,jEnd
                   do i=iBeg,iEnd
                      v1(1) = x(i,jl,j,1) - x(i,ny,j,1)
                      v1(2) = x(i,jl,j,2) - x(i,ny,j,2)
                      v1(3) = x(i,jl,j,3) - x(i,ny,j,3)
                      dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                           +      v1(3)*norm(3))
                      x(i,je,j,1) = x(i,ny,j,1) + dot*norm(1)
                      x(i,je,j,2) = x(i,ny,j,2) + dot*norm(2)
                      x(i,je,j,3) = x(i,ny,j,3) + dot*norm(3)
                   enddo
                enddo

             case (kMin)
                iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl

                if(iBeg == 1)     iBeg = 0
                if(iEnd == iiMax) iEnd = iiMax + 1

                if(jBeg == 1)     jBeg = 0
                if(jEnd == jjMax) jEnd = jjMax + 1

                do j=jBeg,jEnd
                   do i=iBeg,iEnd
                      v1(1) = x(i,j,1,1) - x(i,j,2,1)
                      v1(2) = x(i,j,1,2) - x(i,j,2,2)
                      v1(3) = x(i,j,1,3) - x(i,j,2,3)
                      dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                           +      v1(3)*norm(3))
                      x(i,j,0,1) = x(i,j,2,1) + dot*norm(1)
                      x(i,j,0,2) = x(i,j,2,2) + dot*norm(2)
                      x(i,j,0,3) = x(i,j,2,3) + dot*norm(3)
                   enddo
                enddo

             case (kMax)
                iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl

                if(iBeg == 1)     iBeg = 0
                if(iEnd == iiMax) iEnd = iiMax + 1

                if(jBeg == 1)     jBeg = 0
                if(jEnd == jjMax) jEnd = jjMax + 1

                do j=jBeg,jEnd
                   do i=iBeg,iEnd
                      v1(1) = x(i,j,kl,1) - x(i,j,nz,1)
                      v1(2) = x(i,j,kl,2) - x(i,j,nz,2)
                      v1(3) = x(i,j,kl,3) - x(i,j,nz,3)
                      dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                           +      v1(3)*norm(3))
                      x(i,j,ke,1) = x(i,j,nz,1) + dot*norm(1)
                      x(i,j,ke,2) = x(i,j,nz,2) + dot*norm(2)
                      x(i,j,ke,3) = x(i,j,nz,3) + dot*norm(3)
                   enddo
                enddo
             end select
          endif testSingular
       end if testSymmetry
    enddo loopBocos
  end subroutine xhalo_block

  subroutine resScale

    use constants
    use blockPointers, only : il, jl, kl, nx, ny, nz, volRef, dw
    use flowVarRefState, only : nwf, nt1, nt2
    use inputIteration, only : turbResScale
    implicit none

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii, nTurb
    real(kind=realType) :: ovol

    ! Divide through by the reference volume
    nTurb = nt2-nt1+1
#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do ii=0,nx*ny*nz-1
       i = mod(ii, nx) + 2
       j = mod(ii/nx, ny) + 2
       k = ii/(nx*ny) + 2
#else
       do k=2,kl
          do j=2,jl
             do i=2,il
#endif
                oVol = one/volRef(i,j,k)
                dw(i, j, k, 1:nwf) = dw(i,j, k, 1:nwf)* ovol
                dw(i, j, k, nt1:nt2) = dw(i, j, k, nt1:nt2) * ovol * turbResScale(1:nTurb)
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine resScale

  subroutine sumDwAndFw

    use constants
    use blockPointers, only :il, jl, kl, dw, fw, iBlank
    use flowVarRefState, only : nwf

    implicit none

    ! Local Variables
    integer(kind=intType) :: i, j, k, l

    do l=1, nwf
       do k=2, kl
          do j=2, jl
             do i=2, il
                dw(i,j,k,l) = (dw(i,j,k,l) + fw(i,j,k,l)) &
                     * max(real(iblank(i,j,k), realType), zero)
             end do
          end do
       end do
    end do
  end subroutine sumDwAndFw

  ! EXPL: The three following functions; gridVelocitiesFineLevel_block(),
  !       slipVelocitiesFineLevel_block() and normalVelocities_block()
  !       all had to bee rewritten due the use of pointers in the three
  !       functions. Now, having removed the pointers we can feed these
  !       functions to Tapenade and get differentiated code that can be
  !       readily compiled.
  !            An instructive example already present in the code for
  !       this re-write is to compare the subroutine metric() from
  !       prepocessingAPI.F90 and the 'improved' subroutine called
  !       'boundaryNormals()' which can be found in the present file.
  !       One can basically compare the two 'bocoLoop'-loops.
  !
  subroutine gridVelocitiesFineLevel_block(useOldCoor, t, sps)
    !
    !       gridVelocitiesFineLevel computes the grid velocities for
    !       the cell centers and the normal grid velocities for the faces
    !       of moving blocks for the currently finest grid, i.e.
    !       groundLevel. The velocities are computed at time t for
    !       spectral mode sps. If useOldCoor is .true. the velocities
    !       are determined using the unsteady time integrator in
    !       combination with the old coordinates; otherwise the analytic
    !       form is used.
    !
    use blockPointers
    use cgnsGrid
    use flowVarRefState
    use inputMotion
    use inputUnsteady
    use iteration
    use inputPhysics
    use inputTSStabDeriv
    use monitor
    use communication ! gives myID and adflow_comm_world
    use flowUtils, only :  derivativeRotMatrixRigid, getDirVector
    use utils, only : setCoefTimeIntegrator,tsAlpha, tsBeta, tsMach, terminate, &
         rotMatrixRigidBody, getDirAngle
    

    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor

    real(kind=realType), dimension(*), intent(in) :: t
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm
    integer(kind=intType) :: i, j, k, ii, iie, jje, kke, ii_

    real(kind=realType) :: oneOver4dt, oneOver8dt
    real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
    real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

    real(kind=realType), dimension(3) :: sc, xc, xxc
    real(kind=realType), dimension(3) :: rotCenter, rotRate

    real(kind=realType), dimension(3)   :: rotationPoint
    real(kind=realType), dimension(3,3) :: rotationMatrix,&
         derivRotationMatrix

    real(kind=realType) :: tNew, tOld
    ! mham: removing pointers for Tapenade
    ! this line was commented because the version of tapenade used could not handle 
    ! pointers used like this. 
    ! As a workaround, we use sFace_jk variable, which is a regular fortran variable. 
    ! The math is then reformulated to use this variable, 
    ! and the version of tapenade used works well with this. 
    ! real(kind=realType), dimension(:,:), pointer :: sFace
    ! sFace points to e.g. sFace => sFaceI(i,:,:). sFace has the
    ! shape: sFaceI(0:ie,je,ke) (see l. 617 in initializeFlow.F90)
    real(kind=realType) :: sFace_jk ! fake pointer

    ! mham: removing pointers for Tapenade
    ! real(kind=realType), dimension(:,:,:),   pointer :: xx, ss
    real(kind=realType) :: xx_x, xx_y, xx_z ! fake pointer
    real(kind=realType) :: xx_x_jk_10, xx_y_jk_10, xx_z_jk_10, &
         xx_x_jk_01, xx_y_jk_01, xx_z_jk_01, &
         xx_x_jk_11, xx_y_jk_11, xx_z_jk_11
    real(kind=realType) :: ss_x, ss_y, ss_z ! fake pointer

    ! only used when useOldCoor==.TRUE.
    ! real(kind=realType), dimension(:,:,:,:), pointer :: xxOld

    real(kind=realType) :: intervalMach,alphaTS,alphaIncrement,&
         betaTS,betaIncrement
    real(kind=realType), dimension(3) ::velDir
    real(kind=realType), dimension(3) :: refDirection


    ! Compute the mesh velocity from the given mesh Mach number.

    ! vel{x,y,z}Grid0 is the ACTUAL velocity you want at the
    ! geometry.
    aInf = sqrt(gammaInf*pInf/rhoInf)
    velxGrid0 = (aInf*machgrid)*(-velDirFreestream(1))
    velyGrid0 = (aInf*machgrid)*(-velDirFreestream(2))
    velzGrid0 = (aInf*machgrid)*(-velDirFreestream(3))

    ! Compute the derivative of the rotation matrix and the rotation
    ! point; needed for velocity due to the rigid body rotation of
    ! the entire grid. It is assumed that the rigid body motion of
    ! the grid is only specified if there is only 1 section present.

    call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, t(1))

    !compute the rotation matrix to update the velocities for the time
    !spectral stability derivative case...

    if(TSStability)then
       ! Determine the time values of the old and new time level.
       ! It is assumed that the rigid body rotation of the mesh is only
       ! used when only 1 section is present.

       tNew = timeUnsteady + timeUnsteadyRestart
       tOld = tNew - t(1)

       if(TSpMode.or. TSqMode .or.TSrMode) then
          ! Compute the rotation matrix of the rigid body rotation as
          ! well as the rotation point; the latter may vary in time due
          ! to rigid body translation.

          call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

          if(TSAlphaFollowing) then

             velxgrid0 = rotationMatrix(1,1)*velxgrid0 &
                  + rotationMatrix(1,2)*velygrid0 &
                  + rotationMatrix(1,3)*velzgrid0
             velygrid0 = rotationMatrix(2,1)*velxgrid0 &
                  + rotationMatrix(2,2)*velygrid0 &
                  + rotationMatrix(2,3)*velzgrid0
             velzgrid0 = rotationMatrix(3,1)*velxgrid0 &
                  + rotationMatrix(3,2)*velygrid0 &
                  + rotationMatrix(3,3)*velzgrid0

          end if

       elseif(tsAlphaMode)then

          !Determine the alpha for this time instance
          alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
               degreeFourAlpha,  omegaFourAlpha,     &
               cosCoefFourAlpha, sinCoefFourAlpha, t(1))

          alphaTS = alpha+alphaIncrement
          !Determine the grid velocity for this alpha
          refDirection(:) = zero
          refDirection(1) = one
          call getDirVector(refDirection, alphaTS, beta, velDir, liftIndex)

          !do I need to update the lift direction and drag direction as well?
          !set the effictive grid velocity for this time interval
          velxGrid0 = (aInf*machgrid)*(-velDir(1))
          velyGrid0 = (aInf*machgrid)*(-velDir(2))
          velzGrid0 = (aInf*machgrid)*(-velDir(3))

       elseif(tsBetaMode)then

          !Determine the alpha for this time instance
          betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
               degreeFourBeta,  omegaFourBeta,     &
               cosCoefFourBeta, sinCoefFourBeta, t(1))

          betaTS = beta+betaIncrement
          !Determine the grid velocity for this alpha
          refDirection(:) = zero
          refDirection(1) = one
          call getDirVector(refDirection, alpha, betaTS, velDir, liftIndex)

          !do I need to update the lift direction and drag direction as well?
          !set the effictive grid velocity for this time interval
          velxGrid0 = (aInf*machgrid)*(-velDir(1))
          velyGrid0 = (aInf*machgrid)*(-velDir(2))
          velzGrid0 = (aInf*machgrid)*(-velDir(3))
       elseif(TSMachMode)then
          !determine the mach number at this time interval
          IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
               degreeFourMach,  omegaFourMach,     &
               cosCoefFourMach, sinCoefFourMach, t(1))
          !set the effective grid velocity
          velxGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(1))
          velyGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(2))
          velzGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(3))

       elseif(TSAltitudeMode)then
          call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
       else
          call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
       end if
    endif

    testMoving: if( blockIsMoving ) then
       ! Determine the situation we are having here.

       testUseOldCoor: if( useOldCoor ) then
          ! nothing here... we call the rotational functions with the
          ! flag send to useOldCoor=.FALSE.
          !
          ! OBS! this should of course be improved!
          !
       else testUseOldCoor
          !
          !             The velocities must be determined analytically.
          !
          ! Store the rotation center and determine the
          ! nonDimensional rotation rate of this block. As the
          ! reference length is 1 timeRef == 1/uRef and at the end
          ! the nonDimensional velocity is computed.

          j = nbkGlobal

          rotCenter = cgnsDoms(j)%rotCenter
          rotRate   = timeRef*cgnsDoms(j)%rotRate

          velXgrid = velXGrid0
          velYgrid = velYGrid0
          velZgrid = velZGrid0
          !
          !             Grid velocities of the cell centers, including the
          !             1st level halo cells.
          !
          ! Loop over the cells, including the 1st level halo's.

          do k=1,ke
             do j=1,je
                do i=1,ie

                   ! Determine the coordinates of the cell center,
                   ! which are stored in xc.

                   xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1) &
                        +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1) &
                        +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1) &
                        +         x(i-1,j,  k,  1) + x(i,j,  k,  1))
                   xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2) &
                        +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2) &
                        +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2) &
                        +         x(i-1,j,  k,  2) + x(i,j,  k,  2))
                   xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3) &
                        +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3) &
                        +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3) &
                        +         x(i-1,j,  k,  3) + x(i,j,  k,  3))

                   ! Determine the coordinates relative to the
                   ! center of rotation.

                   xxc(1) = xc(1) - rotCenter(1)
                   xxc(2) = xc(2) - rotCenter(2)
                   xxc(3) = xc(3) - rotCenter(3)

                   ! Determine the rotation speed of the cell center,
                   ! which is omega*r.

                   sc(1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
                   sc(2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
                   sc(3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

                   ! Determine the coordinates relative to the
                   ! rigid body rotation point.

                   xxc(1) = xc(1) - rotationPoint(1)
                   xxc(2) = xc(2) - rotationPoint(2)
                   xxc(3) = xc(3) - rotationPoint(3)

                   ! Determine the total velocity of the cell center.
                   ! This is a combination of rotation speed of this
                   ! block and the entire rigid body rotation.

                   s(i,j,k,1) = sc(1) + velxGrid           &
                        + derivRotationMatrix(1,1)*xxc(1) &
                        + derivRotationMatrix(1,2)*xxc(2) &
                        + derivRotationMatrix(1,3)*xxc(3)
                   s(i,j,k,2) = sc(2) + velyGrid           &
                        + derivRotationMatrix(2,1)*xxc(1) &
                        + derivRotationMatrix(2,2)*xxc(2) &
                        + derivRotationMatrix(2,3)*xxc(3)
                   s(i,j,k,3) = sc(3) + velzGrid           &
                        + derivRotationMatrix(3,1)*xxc(1) &
                        + derivRotationMatrix(3,2)*xxc(2) &
                        + derivRotationMatrix(3,3)*xxc(3)
                enddo
             enddo
          enddo
          !
          !             Normal grid velocities of the faces.
          !
          ! Loop over the three directions.

          loopDirection: do mm=1,3

             ! Set the upper boundaries depending on the direction.

             select case (mm)
             case (1_intType)       ! Normals in i-direction
                iie = ie; jje = je; kke = ke

             case (2_intType)       ! Normals in j-direction
                iie = je; jje = ie; kke = ke

             case (3_intType)       ! Normals in k-direction
                iie = ke; jje = ie; kke = je
             end select
             !
             !               Normal grid velocities in generalized i-direction.
             !               mm == 1: i-direction
             !               mm == 2: j-direction
             !               mm == 3: k-direction
             !
             
             ! mham: here we insert the the fake double loop
             ! to avoid use of pointers
             ! mham: removing old code ->
             ! do i=0,iie
             ! Set the pointers for the coordinates, normals and
             ! normal velocities for this generalized i-plane.
             ! This depends on the value of mm.
             ! mham: inserting new code ->
             FakePointerLoop: do ii=0,( (iie+1)*(kke)*(jje) )-1
                ! outer loop                  
                i = ii/(jje*(kke)) + 0   ! 0:iie ! 
                ! middle loop                      
                k = mod(ii/jje, kke) + 1 ! 1:kke ! 
                ! inner loop                       
                j = mod(ii, jje) + 1     ! 1:jje ! 

                select case (mm)
                   ! mham: Take a look at the indices below. All 'i' and 'j'
                   ! indices for x(:,:,:,:) have been moved one back, i.e. we
                   ! have added a minus 1. WHY? because before (in the original
                   ! code) the pointers DID NOT INHERET THE BOUNDS of the target
                   ! array. This means that for x(:,:,:,:) which has lower
                   ! bound of 0 the pointer to x, i.e. xx would have lower
                   ! bound of 1. Now, in the present code we DO NOT use pointers
                   ! so we still have lower bound of 0. We therefore shift the
                   ! 'i' and 'j' (which we designed for pointers) one down to
                   ! acces the right data.
                   !
                   ! We do NOT just shift the entire 'i' and 'j' indices since
                   ! they are also used for (fake) si, sj and sk pointers and
                   ! here we know from preprocessing/preprocessingAPI.F90:
                   !    allocate(flowDoms(nn,level,sps)%si(0:ie,1:je,1:ke,3),
                   ! that si(:,:,:,:) for example does not start from 0 for 'i'
                   ! and 'j' but it starts from 1...
                case (1_intType)       ! normals in i-direction
                   ! xx =>  x(i,:,:,:)
                   ! mham: for xx we note that it is allocated in
                   !       partitioning.F90 around l. 1758:
                   ! allocate(flowDoms(nn,1,mm)%x(0:ie,0:je,0:ke,3), stat=ierr)
                   !       clearly, it is the last index that refers to dim.
                   xx_x=x(i,j-1,k-1,1); xx_y=x(i,j-1,k-1,2); xx_z=x(i,j-1,k-1,3);
                   ! j+=1,k+=0
                   xx_x_jk_10=x(i,j+1-1,k-1,1)
                   xx_y_jk_10=x(i,j+1-1,k-1,2)
                   xx_z_jk_10=x(i,j+1-1,k-1,3)
                   ! j+=0,k+=1
                   xx_x_jk_01=x(i,j-1,k+1-1,1)
                   xx_y_jk_01=x(i,j-1,k+1-1,2)
                   xx_z_jk_01=x(i,j-1,k+1-1,3)
                   ! j+=1,k+=1
                   xx_x_jk_11=x(i,j+1-1,k+1-1,1)
                   xx_y_jk_11=x(i,j+1-1,k+1-1,2)
                   xx_z_jk_11=x(i,j+1-1,k+1-1,3)
                   ! ss => si(i,:,:,:);
                   ss_x=si(i,j,k,1); ss_y=si(i,j,k,2); ss_z=si(i,j,k,3);
                   ! mham: we do NOT make fake pointers to the sFaceI/J/K
                   !       since they are not used in the computation.
                   !       Instead, we have a simply container that receives
                   !       the values from the computation. We then must
                   !       remember to insert these values in sFaceI/J/K
                   !       that we ideally should have simply pointed to...
                   ! sFace => sFaceI(i,:,:)
                case (2_intType)       ! normals in j-direction
                   ! xx =>  x(:,i,:,:)
                   xx_x=x(j-1,i,k-1,1); xx_y=x(j-1,i,k-1,2); xx_z=x(j-1,i,k-1,3); 
                   ! j+=1,k+=0
                   xx_x_jk_10=x(j+1-1,i,k-1,1)
                   xx_y_jk_10=x(j+1-1,i,k-1,2)
                   xx_z_jk_10=x(j+1-1,i,k-1,3)
                   ! j+=0,k+=1
                   xx_x_jk_01=x(j-1,i,k+1-1,1)
                   xx_y_jk_01=x(j-1,i,k+1-1,2)
                   xx_z_jk_01=x(j-1,i,k+1-1,3)
                   ! j+=1,k+=1
                   xx_x_jk_11=x(j+1-1,i,k+1-1,1)
                   xx_y_jk_11=x(j+1-1,i,k+1-1,2)
                   xx_z_jk_11=x(j+1-1,i,k+1-1,3)
                   ! ss => sj(:,i,:,:)
                   ss_x=sj(j,i,k,1); ss_y=sj(j,i,k,2); ss_z=sj(j,i,k,3);
                   ! sFace => sFaceJ(:,i,:)
                case (3_intType)       ! normals in k-direction
                   ! xx =>  x(:,:,i,:)
                   xx_x=x(j-1,k-1,i,1); xx_y=x(j-1,k-1,i,2); xx_z=x(j-1,k-1,i,3); 
                   ! j+=1,k+=0
                   xx_x_jk_10=x(j+1-1,k-1,i,1)
                   xx_y_jk_10=x(j+1-1,k-1,i,2)
                   xx_z_jk_10=x(j+1-1,k-1,i,3)
                   ! j+=0,k+=1
                   xx_x_jk_01=x(j-1,k+1-1,i,1)
                   xx_y_jk_01=x(j-1,k+1-1,i,2)
                   xx_z_jk_01=x(j-1,k+1-1,i,3)
                   ! j+=1,k+=1
                   xx_x_jk_11=x(j+1-1,k+1-1,i,1)
                   xx_y_jk_11=x(j+1-1,k+1-1,i,2)
                   xx_z_jk_11=x(j+1-1,k+1-1,i,3)
                   ! ss => sk(:,:,i,:);
                   ss_x=sk(j,k,i,1); ss_y=sk(j,k,i,2); ss_z=sk(j,k,i,3);
                   ! sFace => sFaceK(:,:,i)
                end select

                ! Loop over the k and j-direction of this generalized
                ! i-face. Note that due to the usage of the pointer
                ! xx an offset of +1 must be used in the coordinate
                ! array, because x originally starts at 0 for the
                ! i, j and k indices.

                ! mham: commnet out use of pointer nested loop
                ! do k=1,kke
                !    do j=1,jje

                ! Determine the coordinates of the face center,
                ! which are stored in xc.

                ! xc(1) = fourth*(xx(j+1,k+1,1) + xx(j,k+1,1) &
                !      +         xx(j+1,k,  1) + xx(j,k,  1))
                ! xc(2) = fourth*(xx(j+1,k+1,2) + xx(j,k+1,2) &
                !      +         xx(j+1,k,  2) + xx(j,k,  2))
                ! xc(3) = fourth*(xx(j+1,k+1,3) + xx(j,k+1,3) &
                !      +         xx(j+1,k,  3) + xx(j,k,  3))
                xc(1) = fourth*(xx_x_jk_11 + xx_x_jk_01 &
                     +         xx_x_jk_10 + xx_x)
                xc(2) = fourth*(xx_y_jk_11 + xx_y_jk_01 &
                     +         xx_y_jk_10 + xx_y)
                xc(3) = fourth*(xx_z_jk_11 + xx_z_jk_01 &
                     +         xx_z_jk_10 + xx_z)

                ! Determine the coordinates relative to the
                ! center of rotation.

                xxc(1) = xc(1) - rotCenter(1)
                xxc(2) = xc(2) - rotCenter(2)
                xxc(3) = xc(3) - rotCenter(3)

                ! Determine the rotation speed of the face center,
                ! which is omega*r.

                sc(1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
                sc(2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
                sc(3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

                ! Determine the coordinates relative to the
                ! rigid body rotation point.

                xxc(1) = xc(1) - rotationPoint(1)
                xxc(2) = xc(2) - rotationPoint(2)
                xxc(3) = xc(3) - rotationPoint(3)

                ! Determine the total velocity of the cell face.
                ! This is a combination of rotation speed of this
                ! block and the entire rigid body rotation.

                sc(1) = sc(1) + velxGrid           &
                     + derivRotationMatrix(1,1)*xxc(1) &
                     + derivRotationMatrix(1,2)*xxc(2) &
                     + derivRotationMatrix(1,3)*xxc(3)
                sc(2) = sc(2) + velyGrid           &
                     + derivRotationMatrix(2,1)*xxc(1) &
                     + derivRotationMatrix(2,2)*xxc(2) &
                     + derivRotationMatrix(2,3)*xxc(3)
                sc(3) = sc(3) + velzGrid           &
                     + derivRotationMatrix(3,1)*xxc(1) &
                     + derivRotationMatrix(3,2)*xxc(2) &
                     + derivRotationMatrix(3,3)*xxc(3)

                ! Store the dot product of grid velocity sc and
                ! the normal ss in sFace.

                ! sFace(j,k) = sc(1)*ss(j,k,1) + sc(2)*ss(j,k,2) &
                !      + sc(3)*ss(j,k,3)
                sFace_jk = sc(1)*ss_x + sc(2)*ss_y &
                     + sc(3)*ss_z

                ! mham: we have to remember to manually insert 
                !       the sFace_jk in the correct container since
                !       we removed the pointers...
                select case (mm)
                case (1_intType)       ! normals in i-direction
                   sFaceI(i,j,k) = sFace_jk 
                case (2_intType)       ! normals in j-direction
                   sFaceJ(j,i,k) = sFace_jk
                case (3_intType)       ! normals in k-direction
                   sFaceK(j,k,i) = sFace_jk
                end select

                ! mham: remove the nested pointer loop 
                !    enddo
                ! enddo
                ! mham: remember to close the fake double loop
                ! enddo FakeNestedLoop
             enddo FakePointerLoop

          enddo loopDirection
       endif testUseOldCoor
    endif testMoving

  end subroutine gridVelocitiesFineLevel_block
  ! mham addition:
  ! this routine was added using the previous adjoint implementation in SUMB.
  ! most of the comments and variables are to go around the tapenade issues with pointers.
  subroutine slipVelocitiesFineLevel_block(useOldCoor, t, sps)
    !
    !       slipVelocitiesFineLevel computes the slip velocities for
    !       viscous subfaces on all viscous boundaries on groundLevel for
    !       the given spectral solution. If useOldCoor is .true. the
    !       velocities are determined using the unsteady time integrator;
    !       otherwise the analytic form is used.
    !
    use constants
    use inputTimeSpectral
    use blockPointers
    use cgnsGrid
    use flowVarRefState
    use inputMotion
    use inputUnsteady
    use iteration
    use inputPhysics
    use inputTSStabDeriv
    use monitor
    use communication
    use flowUtils, only :  derivativeRotMatrixRigid, getDirVector
    use utils, only : tsAlpha, tsBeta, tsMach, terminate, rotMatrixRigidBody, &
         setCoefTimeIntegrator, getDirAngle
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    logical,               intent(in) :: useOldCoor

    real(kind=realType), dimension(*), intent(in) :: t
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, i, j, level

    real(kind=realType) :: oneOver4dt
    real(kind=realType) :: velxGrid, velyGrid, velzGrid,ainf
    real(kind=realType) :: velxGrid0, velyGrid0, velzGrid0

    real(kind=realType), dimension(3) :: xc, xxc
    real(kind=realType), dimension(3) :: rotCenter, rotRate

    real(kind=realType), dimension(3)   :: rotationPoint
    real(kind=realType), dimension(3,3) :: rotationMatrix,&
         derivRotationMatrix

    real(kind=realType) :: tNew, tOld

    ! MHAM corrections for Tapenade use
    ! real(kind=realType), dimension(:,:,:),   pointer :: uSlip
    ! real(kind=realType), dimension(:,:,:),   pointer :: xFace
    ! mham: we do not need xFaceOld...
    ! real(kind=realType), dimension(:,:,:,:), pointer :: xFaceOld
    !
    ! l. 2241 in BCData.F90
    ! BCData(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3)
    ! This one we simply insert below. It is one-to-one. Code will be longer
    ! / less readable but effectively unaltered
    !
    real(kind=realType)   :: xFace_x,xFace_y,xFace_z
    real(kind=realType)   :: xFace_x_ij_10,xFace_y_ij_10,xFace_z_ij_10 ,&
         xFace_x_ij_01,xFace_y_ij_01,xFace_z_ij_01 ,&
         xFace_x_ij_11,xFace_y_ij_11,xFace_z_ij_11 
    ! real(kind=realType)   :: xFaceOld
    ! mham: now, we must remember to declare new counter used in the fake-loop
    integer :: ii,i_,j_
    !        
    real(kind=realType) :: intervalMach,alphaTS,alphaIncrement,&
         betaTS,betaIncrement
    real(kind=realType), dimension(3) ::velDir
    real(kind=realType), dimension(3) :: refDirection

    ! Determine the situation we are having here.

    testUseOldCoor: if( useOldCoor ) then
       ! mham: Nothing here, since we set useOldCoor=.FALSE.
       !       Everything has been cut out.
    else testUseOldCoor

       ! The velocities must be determined analytically.

       ! Compute the mesh velocity from the given mesh Mach number.

       !  aInf = sqrt(gammaInf*pInf/rhoInf)
       !  velxGrid = aInf*MachGrid(1)
       !  velyGrid = aInf*MachGrid(2)
       !  velzGrid = aInf*MachGrid(3)

       aInf = sqrt(gammaInf*pInf/rhoInf)
       velxGrid0 = (aInf*machgrid)*(-velDirFreestream(1))
       velyGrid0 = (aInf*machgrid)*(-velDirFreestream(2))
       velzGrid0 = (aInf*machgrid)*(-velDirFreestream(3))

       ! Compute the derivative of the rotation matrix and the rotation
       ! point; needed for velocity due to the rigid body rotation of
       ! the entire grid. It is assumed that the rigid body motion of
       ! the grid is only specified if there is only 1 section present.

       call derivativeRotMatrixRigid(derivRotationMatrix, rotationPoint, &
            t(1))

       !compute the rotation matrix to update the velocities for the time
       !spectral stability derivative case...

       if(TSStability)then
          ! Determine the time values of the old and new time level.
          ! It is assumed that the rigid body rotation of the mesh is only
          ! used when only 1 section is present.

          tNew = timeUnsteady + timeUnsteadyRestart
          tOld = tNew - t(1)

          if(TSpMode.or. TSqMode .or.TSrMode) then
             ! Compute the rotation matrix of the rigid body rotation as
             ! well as the rotation point; the latter may vary in time due
             ! to rigid body translation.

             call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

             if(TSAlphaFollowing) then

                velxgrid0 = rotationMatrix(1,1)*velxgrid0 &
                     + rotationMatrix(1,2)*velygrid0 &
                     + rotationMatrix(1,3)*velzgrid0
                velygrid0 = rotationMatrix(2,1)*velxgrid0 &
                     + rotationMatrix(2,2)*velygrid0 &
                     + rotationMatrix(2,3)*velzgrid0
                velzgrid0 = rotationMatrix(3,1)*velxgrid0 &
                     + rotationMatrix(3,2)*velygrid0 &
                     + rotationMatrix(3,3)*velzgrid0

             endif
          elseif(tsAlphaMode)then
             !Determine the alpha for this time instance
             alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
                  degreeFourAlpha,  omegaFourAlpha,     &
                  cosCoefFourAlpha, sinCoefFourAlpha, t(1))

             alphaTS = alpha+alphaIncrement
             !Determine the grid velocity for this alpha
             refDirection(:) = zero
             refDirection(1) = one
             call getDirVector(refDirection, alphaTS, beta, velDir, liftIndex)

             !do I need to update the lift direction and drag direction as well?
             !set the effictive grid velocity for this time interval
             velxGrid0 = (aInf*machgrid)*(-velDir(1))
             velyGrid0 = (aInf*machgrid)*(-velDir(2))
             velzGrid0 = (aInf*machgrid)*(-velDir(3))

          elseif(tsBetaMode)then

             !Determine the alpha for this time instance
             betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
                  degreeFourBeta,  omegaFourBeta,     &
                  cosCoefFourBeta, sinCoefFourBeta, t(1))

             betaTS = beta+betaIncrement
             !Determine the grid velocity for this alpha
             refDirection(:) = zero
             refDirection(1) = one
             call getDirVector(refDirection, alpha, betaTS, velDir, liftIndex)

             !do I need to update the lift direction and drag direction as well?
             !set the effictive grid velocity for this time interval
             velxGrid0 = (aInf*machgrid)*(-velDir(1))
             velyGrid0 = (aInf*machgrid)*(-velDir(2))
             velzGrid0 = (aInf*machgrid)*(-velDir(3))
          elseif(TSMachMode)then
             !determine the mach number at this time interval
             IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
                  degreeFourMach,  omegaFourMach,     &
                  cosCoefFourMach, sinCoefFourMach, t(1))
             !set the effective grid velocity
             velxGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(1))
             velyGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(2))
             velzGrid0 = (aInf*(IntervalMach+machgrid))*(-velDirFreestream(3))

          elseif(TSAltitudeMode)then
             call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
          else
             call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
          end if
       endif

       ! Loop over the number of viscous subfaces.

       bocoLoop2: do mm=1,nViscBocos
          do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
             do i=BCData(mm)%icBeg, BCData(mm)%icEnd



             ! Determine the grid face on which the subface is located
             ! and set some variables accordingly.

             select case (BCFaceID(mm))

             case (iMin)
                ! xFace => x(1,:,:,:)
                xFace_x =  x(1,i-1,j-1,1);
                xFace_y =  x(1,i-1,j-1,2);
                xFace_z =  x(1,i-1,j-1,3);             
                xFace_x_ij_10 =  x(1,i+1-1,j+0-1,1);
                xFace_y_ij_10 =  x(1,i+1-1,j+0-1,2);
                xFace_z_ij_10 =  x(1,i+1-1,j+0-1,3);             
                xFace_x_ij_01 =  x(1,i+0-1,j+1-1,1);
                xFace_y_ij_01 =  x(1,i+0-1,j+1-1,2);
                xFace_z_ij_01 =  x(1,i+0-1,j+1-1,3);             
                xFace_x_ij_11 =  x(1,i+1-1,j+1-1,1);
                xFace_y_ij_11 =  x(1,i+1-1,j+1-1,2);
                xFace_z_ij_11 =  x(1,i+1-1,j+1-1,3);             

             case (iMax)
                ! xFace => x(il,:,:,:)
                xFace_x =  x(il,i-1,j-1,1);
                xFace_y =  x(il,i-1,j-1,2);
                xFace_z =  x(il,i-1,j-1,3);             
                xFace_x_ij_10 =  x(il,i+1-1,j+0-1,1);
                xFace_y_ij_10 =  x(il,i+1-1,j+0-1,2);
                xFace_z_ij_10 =  x(il,i+1-1,j+0-1,3);             
                xFace_x_ij_01 =  x(il,i+0-1,j+1-1,1);
                xFace_y_ij_01 =  x(il,i+0-1,j+1-1,2);
                xFace_z_ij_01 =  x(il,i+0-1,j+1-1,3);             
                xFace_x_ij_11 =  x(il,i+1-1,j+1-1,1);
                xFace_y_ij_11 =  x(il,i+1-1,j+1-1,2);
                xFace_z_ij_11 =  x(il,i+1-1,j+1-1,3);             

             case (jMin)
                ! xFace => x(:,1,:,:)
                xFace_x =  x(i-1,1,j-1,1);
                xFace_y =  x(i-1,1,j-1,2);
                xFace_z =  x(i-1,1,j-1,3);             
                xFace_x_ij_10 =  x(i+1-1,1,j+0-1,1);
                xFace_y_ij_10 =  x(i+1-1,1,j+0-1,2);
                xFace_z_ij_10 =  x(i+1-1,1,j+0-1,3);             
                xFace_x_ij_01 =  x(i+0-1,1,j+1-1,1);
                xFace_y_ij_01 =  x(i+0-1,1,j+1-1,2);
                xFace_z_ij_01 =  x(i+0-1,1,j+1-1,3);             
                xFace_x_ij_11 =  x(i+1-1,1,j+1-1,1);
                xFace_y_ij_11 =  x(i+1-1,1,j+1-1,2);
                xFace_z_ij_11 =  x(i+1-1,1,j+1-1,3);             

             case (jMax)
                ! xFace => x(:,jl,:,:)
                xFace_x =  x(i-1,jl,j-1,1);
                xFace_y =  x(i-1,jl,j-1,2);
                xFace_z =  x(i-1,jl,j-1,3);             
                xFace_x_ij_10 =  x(i+1-1,jl,j+0-1,1);
                xFace_y_ij_10 =  x(i+1-1,jl,j+0-1,2);
                xFace_z_ij_10 =  x(i+1-1,jl,j+0-1,3);             
                xFace_x_ij_01 =  x(i+0-1,jl,j+1-1,1);
                xFace_y_ij_01 =  x(i+0-1,jl,j+1-1,2);
                xFace_z_ij_01 =  x(i+0-1,jl,j+1-1,3);             
                xFace_x_ij_11 =  x(i+1-1,jl,j+1-1,1);
                xFace_y_ij_11 =  x(i+1-1,jl,j+1-1,2);
                xFace_z_ij_11 =  x(i+1-1,jl,j+1-1,3);             

             case (kMin)
                ! xFace => x(:,:,1,:)
                xFace_x =  x(i-1,j-1,1,1);
                xFace_y =  x(i-1,j-1,1,2);
                xFace_z =  x(i-1,j-1,1,3);             
                xFace_x_ij_10 =  x(i+1-1,j+0-1,1,1);
                xFace_y_ij_10 =  x(i+1-1,j+0-1,1,2);
                xFace_z_ij_10 =  x(i+1-1,j+0-1,1,3);             
                xFace_x_ij_01 =  x(i+0-1,j+1-1,1,1);
                xFace_y_ij_01 =  x(i+0-1,j+1-1,1,2);
                xFace_z_ij_01 =  x(i+0-1,j+1-1,1,3);             
                xFace_x_ij_11 =  x(i+1-1,j+1-1,1,1);
                xFace_y_ij_11 =  x(i+1-1,j+1-1,1,2);
                xFace_z_ij_11 =  x(i+1-1,j+1-1,1,3);             

             case (kMax)
                ! xFace => x(:,:,kl,:)
                xFace_x =  x(i-1,j-1,kl,1);
                xFace_y =  x(i-1,j-1,kl,2);
                xFace_z =  x(i-1,j-1,kl,3);             
                xFace_x_ij_10 =  x(i+1-1,j+0-1,kl,1);
                xFace_y_ij_10 =  x(i+1-1,j+0-1,kl,2);
                xFace_z_ij_10 =  x(i+1-1,j+0-1,kl,3);             
                xFace_x_ij_01 =  x(i+0-1,j+1-1,kl,1);
                xFace_y_ij_01 =  x(i+0-1,j+1-1,kl,2);
                xFace_z_ij_01 =  x(i+0-1,j+1-1,kl,3);             
                xFace_x_ij_11 =  x(i+1-1,j+1-1,kl,1);
                xFace_y_ij_11 =  x(i+1-1,j+1-1,kl,2);
                xFace_z_ij_11 =  x(i+1-1,j+1-1,kl,3);             

             end select

             ! Store the rotation center and the rotation rate
             ! for this subface.

             ! mham: Notice we had to rename these two placeholders, (i,j)->
             ! (i_,j_) since i and j are now inside the i,j-double loop
             j_ = nbkGlobal
             i_ = cgnsSubface(mm)

             rotCenter = cgnsDoms(j_)%bocoInfo(i_)%rotCenter
             rotRate   = timeRef*cgnsDoms(j_)%bocoInfo(i_)%rotRate

             ! useWindAxis should go back here!
             velXgrid = velXGrid0
             velYgrid = velYGrid0
             velZgrid = velZGrid0

             ! Loop over the quadrilateral faces of the viscous
             ! subface.

             ! mham: this nested loop has been replaced by FakeNestedLoop2
             ! do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
             !    do i=BCData(mm)%icBeg, BCData(mm)%icEnd

             ! Compute the coordinates of the centroid of the face.
             ! Normally this would be an average of i-1 and i, but
             ! due to the usage of the pointer xFace and the fact
             ! that x starts at index 0 this is shifted 1 index.

             ! mham: this region is deprecated and has been replaced to
             ! allow for a more smooth Tapenade differentiation w.o. any pointers
             ! xc(1) = fourth*(xFace(i+1,j+1,1) + xFace(i+1,j,1) &
             !      +         xFace(i,  j+1,1) + xFace(i,  j,1))
             ! xc(2) = fourth*(xFace(i+1,j+1,2) + xFace(i+1,j,2) &
             !      +         xFace(i,  j+1,2) + xFace(i,  j,2))
             ! xc(3) = fourth*(xFace(i+1,j+1,3) + xFace(i+1,j,3) &
             !      +         xFace(i,  j+1,3) + xFace(i,  j,3))
             xc(1) = fourth*(xFace_x_ij_11 + xFace_x_ij_10 &
                  +         xFace_x_ij_01 + xFace_x)
             xc(2) = fourth*(xFace_y_ij_11 + xFace_y_ij_10 &
                  +         xFace_y_ij_01 + xFace_y)
             xc(3) = fourth*(xFace_z_ij_11 + xFace_z_ij_10 &
                  +         xFace_z_ij_01 + xFace_z)

             ! Determine the coordinates relative to the center
             ! of rotation.

             xxc(1) = xc(1) - rotCenter(1)
             xxc(2) = xc(2) - rotCenter(2)
             xxc(3) = xc(3) - rotCenter(3)

             ! Compute the velocity, which is the cross product
             ! of rotRate and xc.

             BCData(mm)%uSlip(i,j,1) = rotRate(2)*xxc(3) - rotRate(3)*xxc(2)
             BCData(mm)%uSlip(i,j,2) = rotRate(3)*xxc(1) - rotRate(1)*xxc(3)
             BCData(mm)%uSlip(i,j,3) = rotRate(1)*xxc(2) - rotRate(2)*xxc(1)

             ! Determine the coordinates relative to the
             ! rigid body rotation point.

             xxc(1) = xc(1) - rotationPoint(1)
             xxc(2) = xc(2) - rotationPoint(2)
             xxc(3) = xc(3) - rotationPoint(3)

             ! Determine the total velocity of the cell center.
             ! This is a combination of rotation speed of this
             ! block and the entire rigid body rotation.

             BCData(mm)%uSlip(i,j,1) = BCData(mm)%uSlip(i,j,1) + velxGrid    &
                  + derivRotationMatrix(1,1)*xxc(1) &
                  + derivRotationMatrix(1,2)*xxc(2) &
                  + derivRotationMatrix(1,3)*xxc(3)
             BCData(mm)%uSlip(i,j,2) = BCData(mm)%uSlip(i,j,2) + velyGrid    &
                  + derivRotationMatrix(2,1)*xxc(1) &
                  + derivRotationMatrix(2,2)*xxc(2) &
                  + derivRotationMatrix(2,3)*xxc(3)
             BCData(mm)%uSlip(i,j,3) = BCData(mm)%uSlip(i,j,3) + velzGrid    &
                  + derivRotationMatrix(3,1)*xxc(1) &
                  + derivRotationMatrix(3,2)*xxc(2) &
                  + derivRotationMatrix(3,3)*xxc(3)
             ! mham: these nested loops have been replaced by FakeNestedLoop2
             !    enddo
             ! enddo
          enddo
       end do
    enddo bocoLoop2

       
    endif testUseOldCoor

  end subroutine slipVelocitiesFineLevel_block

    ! mham addition
    ! The newly added variables are needed to avoid tapenade issues with pointers.
  subroutine normalVelocities_block(sps)
    !
    !       normalVelocitiesAllLevels computes the normal grid
    !       velocities of some boundary faces of the moving blocks for
    !       spectral mode sps. All grid levels from ground level to the
    !       coarsest level are considered.
    !
    use constants
    use blockPointers, only : il, jl, kl, addGridVelocities, nBocos, BCData, &
         sfaceI, sfaceJ, sfaceK, bcFaceID, si, sj, sk
    !use iteration
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps
    !
    !      Local variables.
    !
    integer(kind=intType) :: mm
    integer(kind=intType) :: i, j

    real(kind=realType) :: weight, mult

    ! mham: comment out deprecated pointers
    ! real(kind=realType), dimension(:,:),   pointer :: sFace  
    ! real(kind=realType), dimension(:,:,:), pointer :: ss
    real(kind=realType) :: sFace_ij
    real(kind=realType) :: ss_x,ss_y,ss_z
    ! mham: we also need a pointer for the fake double loop:
    integer(kind=intType) :: ii_

    ! Check for a moving block. As it is possible that in a
    ! multidisicplinary environment additional grid velocities
    ! are set, the test should be done on addGridVelocities
    ! and not on blockIsMoving.

    testMoving: if( addGridVelocities ) then
       !
       !             Determine the normal grid velocities of the boundaries.
       !             As these values are based on the unit normal. A division
       !             by the length of the normal is needed.
       !             Furthermore the boundary unit normals are per definition
       !             outward pointing, while on the iMin, jMin and kMin
       !             boundaries the face normals are inward pointing. This
       !             is taken into account by the factor mult.
       !
       ! Loop over the boundary subfaces.

       bocoLoop: do mm=1,nBocos

          ! Check whether rFace is allocated.

          testAssoc: if( associated(BCData(mm)%rFace) ) then

             ! Determine the block face on which the subface is
             ! located and set some variables accordingly.


            FakeDoubleLoop: do ii_=0,(BCData(mm)%jcEnd - bcData(mm)%jcBeg + 1)*(BCData(mm)%icEnd - BCData(mm)%icBeg + 1) - 1
               j = ii_/(BCData(mm)%icEnd - BCData(mm)%icBeg + 1) + BCData(mm)%jcBeg
               i = mod(ii_, (BCData(mm)%icEnd - BCData(mm)%icBeg + 1)) + BCData(mm)%icBeg

             select case (BCFaceID(mm))
             case (iMin)
                mult = -one
                ! ss => si(1,:,:,:)
                ss_x=si(1,i,j,1); ss_y=si(1,i,j,2); ss_z=si(1,i,j,3)
                ! sFace => sFaceI(1,:,:)
                sFace_ij = sFaceI(1,i,j)
             case (iMax)
                mult = one
                ! ss => si(il,:,:,:)
                ss_x=si(il,i,j,1); ss_y=si(il,i,j,2); ss_z=si(il,i,j,3)
                ! sFace => sFaceI(il,:,:)
                sFace_ij = sFaceI(il,i,j)
             case (jMin)
                mult = -one
                ! ss => sj(:,1,:,:) 
                ss_x=sj(i,1,j,1); ss_y=sj(i,1,j,2); ss_z=sj(i,1,j,3)
                ! sFace => sFaceJ(:,1,:)
                sFace_ij = sFaceJ(i,1,j)
             case (jMax)
                mult = one
                ! ss => sj(:,jl,:,:)
                ss_x=sj(i,jl,j,1); ss_y=sj(i,jl,j,2); ss_z=sj(i,jl,j,3)
                ! sFace => sFaceJ(:,jl,:)
                sFace_ij = sFaceJ(i,jl,j)
             case (kMin)
                mult = -one
                ! ss => sk(:,:,1,:)
                ss_x=sk(i,j,1,1); ss_y=sk(i,j,1,2); ss_z=sk(i,j,1,3)
                ! sFace => sFaceK(:,:,1)
                sFace_ij = sFaceK(i,j,1)
             case (kMax)
                mult = one
                ! ss => sk(:,:,kl,:)
                ss_x=sk(i,j,kl,1); ss_y=sk(i,j,kl,2); ss_z=sk(i,j,kl,3)
                ! sFace => sFaceK(:,:,kl)
                sFace_ij = sFaceK(i,j,kl)
             end select

             ! Loop over the faces of the subface.
             
             ! mham: remove nested pointer loop
             ! do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
             !    do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                   ! Compute the inverse of the length of the normal
                   ! vector and possibly correct for inward pointing.

             ! mham: remove old pointer version
             ! weight = sqrt(ss(i,j,1)**2 + ss(i,j,2)**2 &
             !      +      ss(i,j,3)**2)
             ! mham: insert new fake pointers
                   weight = sqrt(ss_x**2 + ss_y**2 &
                        +      ss_z**2)

                   if(weight > zero) weight = mult/weight

                   ! Compute the normal velocity based on the outward
                   ! pointing unit normal.

             ! mham: remove old pointer version
                   ! BCData(mm)%rFace(i,j) = weight*sFace(i,j)
             ! mham: insert new fake pointers
                   BCData(mm)%rFace(i,j) = weight*sFace_ij

                   ! mham: remember to remove old enddo's from
                   !       nested pointer loop
             !    enddo
             ! enddo
                   ! mham: remember to close the fake loop
                enddo FakeDoubleLoop
          endif testAssoc
       enddo bocoLoop

    else testMoving

       ! Block is not moving. Loop over the boundary faces and set
       ! the normal grid velocity to zero if allocated.

       do mm=1,nBocos
          if( associated(BCData(mm)%rFace) ) &
               BCData(mm)%rFace = zero
       enddo

    endif testMoving

  end subroutine normalVelocities_block

end module adjointExtra
