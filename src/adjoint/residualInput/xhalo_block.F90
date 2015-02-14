!
!      ******************************************************************
!      *                                                                *
!      * File:          xhalo.f90                                       *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader            *
!      * Starting date: 02-23-2003                                      *
!      * Last modified: 08-12-2009                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine xhalo_block
  !
  !      ******************************************************************
  !      *                                                                *
  !      * xhalo determines the coordinates of the nodal halo's.          *
  !      * First it sets all halo coordinates by simple extrapolation,    *
  !      * then the symmetry planes are treated (also the unit normal of  *
  !      * symmetry planes are determined) and finally an exchange is     *
  !      * made for the internal halo's.                                  *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use BCTypes
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
  !          **************************************************************
  !          *                                                            *
  !          * Mirror the halo coordinates adjacent to the symmetry       *
  !          * planes                                                     *
  !          *                                                            *
  !          **************************************************************
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

