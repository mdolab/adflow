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
                     * real(iblank(i,j,k), realType)
             end do
          end do
       end do
    end do
  end subroutine sumDwAndFw

  subroutine getCostFunctions(globalVals)

    use constants
    use costFunctions
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : pRef, rhoRef, tRef, LRef, gammaInf
    use inputPhysics, only : liftDirection, dragDirection, surfaceRef, machCoef, lengthRef
    use inputTSStabDeriv, only : TSstability
    implicit none

    ! Input 
    real(kind=realType), intent(in), dimension(nLocalValues, nTimeIntervalsSpectral) :: globalVals

    ! Working
    real(kind=realType) :: fact, factMoment, ovrNTS
    real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment, cForce, cMoment
    real(kind=realType) ::  mAvgPtot, mAvgTtot, mAvgPs, mFlow
    integer(kind=intType) :: sps

    ! Factor used for time-averaged quantities.
    ovrNTS = one/nTimeIntervalsSpectral

    ! Sum pressure and viscous contributions
    Force = globalvals(iFp:iFp+2, :) + globalvals(iFv:iFv+2, :)
    Moment = globalvals(iMp:iMp+2, :) + globalvals(iMv:iMv+2, :)

    fact = two/(gammaInf*MachCoef*MachCoef &
         *surfaceRef*LRef*LRef*pRef)
    cForce = fact*force

    ! Moment factor has an extra lengthRef
    fact = fact/(lengthRef*LRef)
    cMoment = fact*Moment

    ! Zero values since we are summing.
    funcValues = zero

    ! Here we finally assign the final function values
    !$AD II-LOOP
    do sps=1, nTimeIntervalsSpectral
       funcValues(costFuncForceX) = funcValues(costFuncForceX) + ovrNTS*force(1, sps)
       funcValues(costFuncForceY) = funcValues(costFuncForceY) + ovrNTS*force(2, sps)
       funcValues(costFuncForceZ) = funcValues(costFuncForceZ) + ovrNTS*force(3, sps)

       funcValues(costFuncForceXCoef) = funcValues(costFuncForceXCoef) + ovrNTS*cForce(1, sps)
       funcValues(costFuncForceYCoef) = funcValues(costFuncForceYCoef) + ovrNTS*cForce(2, sps)
       funcValues(costFuncForceZCoef) = funcValues(costFuncForceZCoef) + ovrNTS*cForce(3, sps)

       funcValues(costFuncMomX) = funcValues(costFuncMomX) + ovrNTS*moment(1, sps)
       funcValues(costFuncMomY) = funcValues(costFuncMomY) + ovrNTS*moment(2, sps)
       funcValues(costFuncMomZ) = funcValues(costFuncMomZ) + ovrNTS*moment(3, sps)

       funcValues(costFuncMomXCoef) = funcValues(costFuncMomXCoef) + ovrNTS*cMoment(1, sps)
       funcValues(costFuncMomYCoef) = funcValues(costFuncMomYCoef) + ovrNTS*cMoment(2, sps)
       funcValues(costFuncMomZCoef) = funcValues(costFuncMomZCoef) + ovrNTS*cMoment(3, sps)

       funcValues(costFuncSepSensor) = funcValues(costFuncSepSensor) + ovrNTS*globalVals(iSepSensor, sps)
       funcValues(costFuncCavitation) = funcValues(costFuncCavitation) + ovrNTS*globalVals(iCavitation, sps)
       funcValues(costFuncSepSensorAvgX) = funcValues(costFuncSepSensorAvgX) + ovrNTS*globalVals(iSepAvg  , sps)
       funcValues(costFuncSepSensorAvgY) = funcValues(costFuncSepSensorAvgY) + ovrNTS*globalVals(iSepAvg+1, sps)
       funcValues(costFuncSepSensorAvgZ) = funcValues(costFuncSepSensorAvgZ) + ovrNTS*globalVals(iSepAvg+2, sps)

       ! Mass flow like objective
       mFlow = globalVals(iMassFlow, sps)
       if (mFlow /= zero) then 
          mAvgPtot = globalVals(iMassPtot, sps)/mFlow*pRef
          mAvgTtot = globalVals(iMassTtot, sps)/mFlow*tRef
          mAvgPs   = globalVals(iMassPs, sps)/mFlow*pRef
          mFlow = globalVals(iMassFlow, sps)*sqrt(Pref/rhoRef)
       else
          mAvgPtot = zero
          mAvgTtot = zero
          mAvgPs = zero
       end if

       funcValues(costFuncMdot)      = funcValues(costFuncMdot) + ovrNTS*mFlow
       funcValues(costFuncMavgPtot ) = funcValues(costFuncMavgPtot) + ovrNTS*mAvgPtot
       funcValues(costFuncMavgPtot)  = funcValues(costFuncMavgTtot) + ovrNTS*mAvgTtot
       funcValues(costFuncMavgPs)    = funcValues(costFuncMAvgPs) + ovrNTS*mAvgPs

       ! Bending moment calc - also broken. 
       ! call computeRootBendingMoment(cForce, cMoment, liftIndex, bendingMoment)
       ! funcValues(costFuncBendingCoef) = funcValues(costFuncBendingCoef) + ovrNTS*bendingMoment

    end do

    ! Lift and Drag (coefficients): Dot product with the lift/drag direction.
    funcValues(costFuncLift) = &
         funcValues(costFuncForceX)*liftDirection(1) + &
         funcValues(costFuncForceY)*liftDirection(2) + &
         funcValues(costFuncForceZ)*liftDirection(3)

    funcValues(costFuncDrag) = &
         funcValues(costFuncForceX)*dragDirection(1) + &
         funcValues(costFuncForceY)*dragDirection(2) + &
         funcValues(costFuncForceZ)*dragDirection(3)

    funcValues(costFuncLiftCoef) = &
         funcValues(costFuncForceXCoef)*liftDirection(1) + &
         funcValues(costFuncForceYCoef)*liftDirection(2) + &
         funcValues(costFuncForceZCoef)*liftDirection(3)

    funcValues(costFuncDragCoef) = &
         funcValues(costFuncForceXCoef)*dragDirection(1) + &
         funcValues(costFuncForceYCoef)*dragDirection(2) + &
         funcValues(costFuncForceZCoef)*dragDirection(3)

    ! -------------------- Time Spectral Objectives ------------------

    if (TSSTability) then 
       print *,'Error: TSStabilityDerivatives are *BROKEN*. They need to be '&
            &'completely verifed from scratch'
       stop
    end if

  end subroutine getCostFunctions
  
end module adjointExtra
