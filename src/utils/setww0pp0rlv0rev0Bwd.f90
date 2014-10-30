!
!      ******************************************************************
!      *                                                                *
!      * File:          setww0pp0rlv0rev0Bwd.f90                        *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-21-2014                                      *
!      * Last modified: 10-21-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setww0pp0rlv0rev0Bwd(nn, idim, ddim, ww0, pp0, rlv0, rev0)
       
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
           ww0(1:je,1:ke,:) = w(0,1:je,1:ke,:)
           pp0(1:je,1:ke) = p(0,1:je,1:ke)
           if( viscous )   rlv0(1:je,1:ke) = rlv(0,1:je,1:ke)
           if( eddyModel ) rev0(1:je,1:ke) = rev(0,1:je,1:ke)
           idim = 1; ddim = 0

         case (iMax)
           ww0(1:je,1:ke,:) = w(ib,1:je,1:ke,:)
           pp0(1:je,1:ke) = p(ib,1:je,1:ke)
           if( viscous )   rlv0(1:je,1:ke) = rlv(ib,1:je,1:ke)
           if( eddyModel ) rev0(1:je,1:ke) = rev(ib,1:je,1:ke)
           idim = 1; ddim = ib

         case (jMin)
           ww0(1:ie,1:ke,:) = w(1:ie,0,1:ke,:)
           pp0(1:ie,1:ke) = p(1:ie,0,1:ke)
           if( viscous )   rlv0(1:ie,1:ke) = rlv(1:ie,0,1:ke)
           if( eddyModel ) rev0(1:ie,1:ke) = rev(1:ie,0,1:ke)
           idim = 2; ddim = 0

         case (jMax)
           ww0(1:ie,1:ke,:) = w(1:ie,jb,1:ke,:)
           pp0(1:ie,1:ke) = p(1:ie,jb,1:ke)
           if( viscous )   rlv0(1:ie,1:ke) = rlv(1:ie,jb,1:ke)
           if( eddyModel ) rev0(1:ie,1:ke) = rev(1:ie,jb,1:ke)
           idim = 2; ddim = jb

         case (kMin)
           ww0(1:ie,1:je,:) = w(1:ie,1:je,0,:)
           pp0(1:ie,1:je) = p(1:ie,1:je,0)
           if( viscous )   rlv0(1:ie,1:je) = rlv(1:ie,1:je,0)
           if( eddyModel ) rev0(1:ie,1:je) = rev(1:ie,1:je,0)
           idim = 3; ddim = 0

         case (kMax)
           ww0(1:ie,1:je,:) = w(1:ie,1:je,kb,:)
           pp0(1:ie,1:je) = p(1:ie,1:je,kb)
           if( viscous )   rlv0(1:ie,1:je) = rlv(1:ie,1:je,kb)
           if( eddyModel ) rev0(1:ie,1:je) = rev(1:ie,1:je,kb)
           idim = 3; ddim = kb

       end select


       end subroutine setww0pp0rlv0rev0Bwd   
       
       



