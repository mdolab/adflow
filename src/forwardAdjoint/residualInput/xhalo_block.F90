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
  real(kind=realType), dimension(:,:,:), pointer :: x0, x1, x2
  logical err
  real(kind=realType) :: length, dot
  logical :: iMinInternal, jMinInternal, kMinInternal 
  logical :: iMaxInternal, jMaxInternal, kMaxInternal
  real(kind=realType), dimension(3) :: v1, v2, norm
 

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  !
  !          **************************************************************
  !          *                                                            *
  !          * Extrapolation of the coordinates. First extrapolation in   *
  !          * i-direction, without halo's, followed by extrapolation in  *
  !          * j-direction, with i-halo's and finally extrapolation in    *
  !          * k-direction, with both i- and j-halo's. In this way also   *
  !          * the indirect halo's get a value, albeit a bit arbitrary.   *
  !          *                                                            *
  !          **************************************************************
  !


  iMinInternal=.false.
  jMinInternal=.false.
  kMinInternal=.false.
  iMaxInternal=.false.
  jMaxInternal=.false.
  kMaxInternal=.false.

  ! Loop over all the subfaces to determine which ones do NOT need to be extrapolated.
  loopSubface: do mm=1,nSubface
     testInternal: if(BCType(mm) ==B2BMatch)then
        select case (BCFaceID(mm))
        case (iMin)
           iMinInternal =.true.
        case (iMax)
           iMaxInternal =.true.
        case (jMin)
           jMinInternal =.true.
        case (jMax)
           jMaxInternal =.true. 
        case (kMin)
           kMinInternal =.true.
        case (kMax)
           kMaxInternal =.true.
        end select
     end if testInternal
  end do loopSubface

  ! Re-loop back over and see if any subface that is NOT B2BMatch is
  ! on the same logical face as a Block2Block. We cannot deal with
  ! properly so will print an error and quit
  err = .False.
  loopSubface2: do mm=1,nSubface
     testNotInternal: if(BCType(mm) .ne. B2BMatch)then
        select case (BCFaceID(mm))
        case (iMin)
           if (iMinInternal) then
              err = .True.
           end if
        case (iMax)
           if (iMaxInternal) then
              err = .True.
           end if
        case (jMin)
           if (jMinInternal) then
              err = .True.
           end if
        case (jMax)
           if (jMaxInternal) then
              err = .True.
           end if
        case (kMin)
           if (kMinInternal) then
              err = .True.
           end if
        case (kMax)
           if (kMaxInternal) then
              err = .True.
           end if
        end select
     end if testNotInternal
  end do loopSubface2

  if (err) then
     print *,'Detected a block-to-block boundary condition on the same face'
     print *,'as another boundary condition. This is not supported.'
     stop
  end if

  ! Extrapolation in i-direction.

  do k=1,kl
     do j=1,jl
        if (.not. iMinInternal) then
           x(0,j,k,1) = two*x(1,j,k,1) - x(2,j,k,1)
           x(0,j,k,2) = two*x(1,j,k,2) - x(2,j,k,2)
           x(0,j,k,3) = two*x(1,j,k,3) - x(2,j,k,3)
        end if

        if (.not. iMaxInternal) then
           x(ie,j,k,1) = two*x(il,j,k,1) - x(nx,j,k,1)
           x(ie,j,k,2) = two*x(il,j,k,2) - x(nx,j,k,2)
           x(ie,j,k,3) = two*x(il,j,k,3) - x(nx,j,k,3)
        end if
     enddo
  enddo

  ! Extrapolation in j-direction.

  do k=1,kl
     do i=0,ie  
        if (.not. jMinInternal) then
           x(i,0,k,1) = two*x(i,1,k,1) - x(i,2,k,1)
           x(i,0,k,2) = two*x(i,1,k,2) - x(i,2,k,2)
           x(i,0,k,3) = two*x(i,1,k,3) - x(i,2,k,3)
        end if

        if (.not. jMaxInternal) then
           x(i,je,k,1) = two*x(i,jl,k,1) - x(i,ny,k,1)
           x(i,je,k,2) = two*x(i,jl,k,2) - x(i,ny,k,2)
           x(i,je,k,3) = two*x(i,jl,k,3) - x(i,ny,k,3)
        end if
     enddo
  enddo

  ! Extrapolation in k-direction.

  do j=0,je
     do i=0,ie
        if (.not. kMinInternal) then
           x(i,j,0,1) = two*x(i,j,1,1) - x(i,j,2,1)
           x(i,j,0,2) = two*x(i,j,1,2) - x(i,j,2,2)
           x(i,j,0,3) = two*x(i,j,1,3) - x(i,j,2,3)
        end if

        if (.not. kMaxInternal) then
           x(i,j,ke,1) = two*x(i,j,kl,1) - x(i,j,nz,1)
           x(i,j,ke,2) = two*x(i,j,kl,2) - x(i,j,nz,2)
           x(i,j,ke,3) = two*x(i,j,kl,3) - x(i,j,nz,3)
        end if
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

        select case (BCFaceID(mm))
        case (iMin)
           iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
           jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
           x0 => x(0,:,:,:); x1 => x(1,:,:,:); x2 => x(2,:,:,:)

        case (iMax)
           iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
           jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
           x0 => x(ie,:,:,:); x1 => x(il,:,:,:); x2 => x(nx,:,:,:)

        case (jMin)
           iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
           jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
           x0 => x(:,0,:,:); x1 => x(:,1,:,:); x2 => x(:,2,:,:)

        case (jMax)
           iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
           jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
           x0 => x(:,je,:,:); x1 => x(:,jl,:,:); x2 => x(:,ny,:,:)

        case (kMin)
           iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
           jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
           x0 => x(:,:,0,:); x1 => x(:,:,1,:); x2 => x(:,:,2,:)

        case (kMax)
           iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
           jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
           x0 => x(:,:,ke,:); x1 => x(:,:,kl,:); x2 => x(:,:,nz,:)
        end select


        if (.not. bcData(mm)%symNormSet) then
           ! This code technically should not run. symNormSet should
           ! already be set from the regular Xhao on the
           ! first call.

           ! Determine the vector from the lower left corner to
           ! the upper right corner. Due to the usage of pointers an
           ! offset of +1 must be used, because the original array x
           ! start at 0.

           v1(1) = x1(iimax+1,jjmax+1,1) - x1(1+1,1+1,1)
           v1(2) = x1(iimax+1,jjmax+1,2) - x1(1+1,1+1,2)
           v1(3) = x1(iimax+1,jjmax+1,3) - x1(1+1,1+1,3)

           ! And the vector from the upper left corner to the
           ! lower right corner.
           
           v2(1) = x1(iimax+1,1+1,1) - x1(1+1,jjmax+1,1)
           v2(2) = x1(iimax+1,1+1,2) - x1(1+1,jjmax+1,2)
           v2(3) = x1(iimax+1,1+1,3) - x1(1+1,jjmax+1,3)
           
           ! Determine the normal of the face by taking the cross
           ! product of v1 and v2 and add it to norm.
           
           norm(1) = v1(2)*v2(3) - v1(3)*v2(2)
           norm(2) = v1(3)*v2(1) - v1(1)*v2(3)
           norm(3) = v1(1)*v2(2) - v1(2)*v2(1)
           
           bcData(mm)%symNorm(1) = norm(1)
           bcData(mm)%symNorm(2) = norm(2)
           bcData(mm)%symNorm(3) = norm(3)
        else
           ! Copy out the saved symNorm
           norm(1) = bcData(mm)%symNorm(1)
           norm(2) = bcData(mm)%symNorm(2)
           norm(3) = bcData(mm)%symNorm(3)
        end if
           
        ! Compute the length of the normal and test if this is
        ! larger than eps. If this is the case this means that
        ! it is a nonsingular subface and the coordinates are
        ! corrected.

        length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)

        testSingular: if(length > eps) then

           ! Compute the unit normal of the subface.

           norm(1) = norm(1)/length
           norm(2) = norm(2)/length
           norm(3) = norm(3)/length

           ! Add an overlap to the symmetry subface if the
           ! boundaries coincide with the block boundaries.
           ! This way the indirect halo's are treated properly.

           if(iBeg == 1)     iBeg = 0
           if(iEnd == iiMax) iEnd = iiMax + 1

           if(jBeg == 1)     jBeg = 0
           if(jEnd == jjMax) jEnd = jjMax + 1

           ! Loop over the nodes of the subface and set the
           ! corresponding halo coordinates.

           do j=jBeg,jEnd
              do i=iBeg,iEnd

                 ! Determine the vector from the internal node to the
                 ! node on the face. Again an offset of +1 must be
                 ! used, due to the usage of pointers.

                 v1(1) = x1(i+1,j+1,1) - x2(i+1,j+1,1)
                 v1(2) = x1(i+1,j+1,2) - x2(i+1,j+1,2)
                 v1(3) = x1(i+1,j+1,3) - x2(i+1,j+1,3)

                 ! Determine two times the normal component of this
                 ! vector; this vector must be added to the
                 ! coordinates of the internal node to obtain the
                 ! halo coordinates. Again the offset of +1.

                 dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                      +      v1(3)*norm(3))

                 x0(i+1,j+1,1) = x2(i+1,j+1,1) + dot*norm(1)
                 x0(i+1,j+1,2) = x2(i+1,j+1,2) + dot*norm(2)
                 x0(i+1,j+1,3) = x2(i+1,j+1,3) + dot*norm(3)


              enddo
           enddo
        endif testSingular
     endif testSymmetry
  enddo loopBocos

end subroutine xhalo_block






