!
!      ******************************************************************
!      *                                                                *
!      * File:          resetBcPointersBwd.f90                          *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-02-2014                                      *
!      * Last modified: 10-02-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetBCPointersBwd(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
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
       ! Copy back the values from the pointers

       select case (BCFaceID(nn))

         case (iMin)
          
           id = 2 + offset;   ih = 1 - offset
           do k = 1,kl+1
              do j = 1,jl+1
                 w(ih,j,k,:) = ww1(j,k,:)
                 w(id,j,k,:) = ww2(j,k,:) 
                 p(ih,j,k)  = pp1(j,k)   
                 p(id,j,k)  = pp2(j,k)   
                 
                 if( viscous ) then
                    rlv(ih,j,k) = rlv1(j,k)
                    rlv(id,j,k) = rlv2(j,k)
                 endif

                 if( eddyModel ) then
                     rev(ih,j,k) = rev1(j,k)
                     rev(id,j,k) = rev2(j,k)
                 endif
              end do
           end do

         !===============================================================

         case (iMax)
           
           id = il - offset;  ih = ie + offset
           do k = 1,kl+1
              do j = 1,jl+1
                 w(ih,j,k,:) = ww1(j,k,:)
                 w(id,j,k,:) = ww2(j,k,:) 
                 p(ih,j,k)  = pp1(j,k)   
                 p(id,j,k)  = pp2(j,k)   
                 
                 if( viscous ) then
                    rlv(ih,j,k) = rlv1(j,k)
                    rlv(id,j,k) = rlv2(j,k)
                 endif

                 if( eddyModel ) then
                     rev(ih,j,k) = rev1(j,k)
                     rev(id,j,k) = rev2(j,k)
                 endif
              end do
           end do
         !===============================================================

         case (jMin)
        
           id = 2 + offset;   ih = 1 - offset
           do k = 1,kl+1
              do i = 1,il+1
                 w(i,ih,k,:) = ww1(i,k,:)
                 w(i,id,k,:) = ww2(i,k,:)
                 p(i,ih,k)   = pp1(i,k)   
                 p(i,id,k)   = pp2(i,k)
                 
                 if( viscous ) then
                     rlv(i,ih,k) = rlv1(i,k)
                     rlv(i,id,k) = rlv2(i,k)
                 endif

                 if( eddyModel ) then
                     rev(i,ih,k) = rev1(i,k)
                     rev(i,id,k) = rev2(i,k)
                 endif
              end do
           end do
         !===============================================================

         case (jMax)
          
           id = jl - offset;  ih = je + offset
           do k = 1,kl+1
              do i = 1,il+1
                 w(i,ih,k,:) = ww1(i,k,:)
                 w(i,id,k,:) = ww2(i,k,:)
                 p(i,ih,k)   = pp1(i,k)   
                 p(i,id,k)   = pp2(i,k)
                 
                 if( viscous ) then
                     rlv(i,ih,k) = rlv1(i,k)
                     rlv(i,id,k) = rlv2(i,k)
                 endif

                 if( eddyModel ) then
                     rev(i,ih,k) = rev1(i,k)
                     rev(i,id,k) = rev2(i,k)
                 endif
              end do
           end do
         !===============================================================

         case (kMin)
           
           id = 2 + offset;   ih = 1 - offset
           do j = 1,jl+1
              do i = 1,il+1
                 w(i,j,ih,:) = ww1(i,j,:)
                 w(i,j,id,:) = ww2(i,j,:)
                 p(i,j,ih)   = pp1(i,j)
                 p(i,j,id)   = pp2(i,j)
                 
                 if( viscous ) then
                    rlv(i,j,ih) = rlv1(i,j)
                    rlv(i,j,id) = rlv2(i,j)
                 endif

                 if( eddyModel ) then
                    rev(i,j,ih) = rev1(i,j)
                    rev(i,j,id) = rev2(i,j)
                 endif
              end do
           end do
         !===============================================================

         case (kMax)
        
           id = kl - offset;  ih = ke + offset
           do j = 1,jl+1
              do i = 1,il+1
                 w(i,j,ih,:) = ww1(i,j,:)
                 w(i,j,id,:) = ww2(i,j,:)
                 p(i,j,ih)   = pp1(i,j)
                 p(i,j,id)   = pp2(i,j)
                 
                 if( viscous ) then
                    rlv(i,j,ih) = rlv1(i,j)
                    rlv(i,j,id) = rlv2(i,j)
                 endif

                 if( eddyModel ) then
                    rev(i,j,ih) = rev1(i,j)
                    rev(i,j,id) = rev2(i,j)
                 endif
              end do
           end do
       end select

       end subroutine resetBCPointersBwd
