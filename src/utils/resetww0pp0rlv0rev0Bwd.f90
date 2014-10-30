!
!      ******************************************************************
!      *                                                                *
!      * File:          resetww0pp0rlv0rev0Bwd.f90                      *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-21-2014                                      *
!      * Last modified: 10-21-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetww0pp0rlv0rev0Bwd(nn, idim, ddim, ww0, pp0, rlv0, rev0)
       
       use BCTypes
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       integer(kind=intType) :: idim, ddim

       real(kind=realType), dimension(imaxDim,jmaxDim,nw) :: ww0
       real(kind=realType), dimension(imaxDim,jmaxDim) :: pp0
       real(kind=realType), dimension(imaxDim,jmaxDim) :: rlv0
       real(kind=realType), dimension(imaxDim,jmaxDim) :: rev0

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
           w(0,1:je,1:ke,:) = ww0(1:je,1:ke,:)
           p(0,1:je,1:ke) = pp0(1:je,1:ke)
           if( viscous )   rlv(0,1:je,1:ke) = rlv0(1:je,1:ke)
           if( eddyModel ) rev(0,1:je,1:ke) = rev0(1:je,1:ke)
           idim = 1; ddim = 0

         case (iMax)
           w(ib,1:je,1:ke,:) = ww0(1:je,1:ke,:)
           p(ib,1:je,1:ke) = pp0(1:je,1:ke)
           if( viscous )   rlv(ib,1:je,1:ke) = rlv0(1:je,1:ke)
           if( eddyModel ) rev(ib,1:je,1:ke) = rev0(1:je,1:ke)
           idim = 1; ddim = ib

         case (jMin)
           w(1:ie,0,1:ke,:) = ww0(1:ie,1:ke,:)
           p(1:ie,0,1:ke) =  pp0(1:ie,1:ke)
           if( viscous )   rlv(1:ie,0,1:ke) = rlv0(1:ie,1:ke)
           if( eddyModel ) rev(1:ie,0,1:ke) = rev0(1:ie,1:ke)
           idim = 2; ddim = 0

         case (jMax)
           w(1:ie,jb,1:ke,:) = ww0(1:ie,1:ke,:)
           p(1:ie,jb,1:ke) = pp0(1:ie,1:ke)
           if( viscous )   rlv(1:ie,jb,1:ke) = rlv0(1:ie,1:ke)
           if( eddyModel ) rev(1:ie,jb,1:ke) = rev0(1:ie,1:ke)
           idim = 2; ddim = jb

         case (kMin)
           w(1:ie,1:je,0,:) = ww0(1:ie,1:je,:)
           p(1:ie,1:je,0) = pp0(1:ie,1:je)
           if( viscous )   rlv(1:ie,1:je,0) = rlv0(1:ie,1:je)
           if( eddyModel ) rev(1:ie,1:je,0) = rev0(1:ie,1:je)
           idim = 3; ddim = 0

         case (kMax)
           w(1:ie,1:je,kb,:) = ww0(1:ie,1:je,:)
           p(1:ie,1:je,kb) = pp0(1:ie,1:je)
           if( viscous )   rlv(1:ie,1:je,kb) = rlv0(1:ie,1:je)
           if( eddyModel ) rev(1:ie,1:je,kb) = rev0(1:ie,1:je)
           idim = 3; ddim = kb

       end select
       end subroutine resetww0pp0rlv0rev0Bwd   
       
       



