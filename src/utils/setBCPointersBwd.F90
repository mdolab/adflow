!
!      ******************************************************************
!      *                                                                *
!      * File:          setBcPointers.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-17-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setBCPointersBwd(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                rev1, rev2, offset)
!
!      ******************************************************************
!      *                                                                *
!      * setBCPointers sets the pointers needed for the boundary        *
!      * condition treatment on a general face, such that the boundary  *
!      * routines are only implemented once instead of 6 times.         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, offset

       real(kind=realType), dimension(imaxDim,jmaxDim,nw) :: ww1, ww2
       real(kind=realType), dimension(imaxDim,jmaxDim) :: pp1, pp2
       real(kind=realType), dimension(imaxDim,jmaxDim) :: rlv1, rlv2
       real(kind=realType), dimension(imaxDim,jmaxDim) :: rev1, rev2
!
!      Local variables
!
       integer(kind=intType) :: id, ih, ierr, i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the face id on which the subface is located and set
       ! the pointers accordinly.

       select case (BCFaceID(nn))

         case (iMin)
          
           id = 2 + offset;   ih = 1 - offset
           do k = 1,kl+1
              do j = 1,jl+1
                 ww1(j,k,:) = w(ih,j,k,:)
                 ww2(j,k,:) = w(id,j,k,:)
                 pp1(j,k)   = p(ih,j,k)
                 pp2(j,k)   = p(id,j,k)
                 
                 if( viscous ) then
                    rlv1(j,k) = rlv(ih,j,k)
                    rlv2(j,k) = rlv(id,j,k)
                 endif

                 if( eddyModel ) then
                    rev1(j,k) = rev(ih,j,k)
                    rev2(j,k) = rev(id,j,k)
                 endif
              end do
           end do
         !===============================================================

         case (iMax)
           
           id = il - offset;  ih = ie + offset
           do k = 1,kl+1
              do j = 1,jl+1
                 ww1(j,k,:) = w(ih,j,k,:)
                 ww2(j,k,:) = w(id,j,k,:)
                 pp1(j,k)   = p(ih,j,k)
                 pp2(j,k)   = p(id,j,k)
                 
                 if( viscous ) then
                    rlv1(j,k) = rlv(ih,j,k)
                    rlv2(j,k) = rlv(id,j,k)
                 endif

                 if( eddyModel ) then
                    rev1(j,k) = rev(ih,j,k)
                    rev2(j,k) = rev(id,j,k)
                 endif
              end do
           end do
         !===============================================================

         case (jMin)
        
           id = 2 + offset;   ih = 1 - offset
           do k = 1,kl+1
              do i = 1,il+1
                 ww1(i,k,:) = w(i,ih,k,:)
                 ww2(i,k,:) = w(i,id,k,:)
                 pp1(i,k)   = p(i,ih,k)
                 pp2(i,k)   = p(i,id,k)
                 
                 if( viscous ) then
                    rlv1(i,k) = rlv(i,ih,k)
                    rlv2(i,k) = rlv(i,id,k)
                 endif

                 if( eddyModel ) then
                    rev1(i,k) = rev(i,ih,k)
                    rev2(i,k) = rev(i,id,k)
                 endif
              end do
           end do
         !===============================================================

         case (jMax)
          
           id = jl - offset;  ih = je + offset
           do k = 1,kl+1
              do i = 1,il+1
                 ww1(i,k,:) = w(i,ih,k,:)
                 ww2(i,k,:) = w(i,id,k,:)
                 pp1(i,k)   = p(i,ih,k)
                 pp2(i,k)   = p(i,id,k)
                 
                 if( viscous ) then
                    rlv1(i,k) = rlv(i,ih,k)
                    rlv2(i,k) = rlv(i,id,k)
                 endif

                 if( eddyModel ) then
                    rev1(i,k) = rev(i,ih,k)
                    rev2(i,k) = rev(i,id,k)
                 endif
              end do
           end do
         !===============================================================

         case (kMin)
           
           id = 2 + offset;   ih = 1 - offset
           do j = 1,jl+1
              do i = 1,il+1
                 ww1(i,j,:) = w(i,j,ih,:)
                 ww2(i,j,:) = w(i,j,id,:)
                 pp1(i,j)   = p(i,j,ih)
                 pp2(i,j)   = p(i,j,id)
                 
                 if( viscous ) then
                    rlv1(i,j) = rlv(i,j,ih)
                    rlv2(i,j) = rlv(i,j,id)
                 endif

                 if( eddyModel ) then
                    rev1(i,j) = rev(i,j,ih)
                    rev2(i,j) = rev(i,j,id)
                 endif
              end do
           end do
         !===============================================================

         case (kMax)
        
           id = kl - offset;  ih = ke + offset
           do j = 1,jl+1
              do i = 1,il+1
                 ww1(i,j,:) = w(i,j,ih,:)
                 ww2(i,j,:) = w(i,j,id,:)
                 pp1(i,j)   = p(i,j,ih)
                 pp2(i,j)   = p(i,j,id)
                 
                 if( viscous ) then
                    rlv1(i,j) = rlv(i,j,ih)
                    rlv2(i,j) = rlv(i,j,id)
                 endif

                 if( eddyModel ) then
                    rev1(i,j) = rev(i,j,ih)
                    rev2(i,j) = rev(i,j,id)
                 endif
              end do
           end do
       end select

       end subroutine setBCPointersBwd
