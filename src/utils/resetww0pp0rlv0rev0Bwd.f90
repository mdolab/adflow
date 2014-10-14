!
!      ******************************************************************
!      *                                                                *
!      * File:          resetww0pp0rlv0rev0Bwd.f90                           *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
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
           w(0,1:,1:,:) = ww0
           p(0,1:,1:) = pp0
           if( viscous )   rlv(0,1:,1:) = rlv0
           if( eddyModel ) rev(0,1:,1:) = rev0
           idim = 1; ddim = 0

         case (iMax)
           w(ib,1:,1:,:) = ww0
           p(ib,1:,1:) = pp0
           if( viscous )   rlv(ib,1:,1:) = rlv0
           if( eddyModel ) rev(ib,1:,1:) = rev0
           idim = 1; ddim = ib

         case (jMin)
           w(1:,0,1:,:) = ww0
           p(1:,0,1:) =  pp0
           if( viscous )   rlv(1:,0,1:) = rlv0
           if( eddyModel ) rev(1:,0,1:) = rev0
           idim = 2; ddim = 0

         case (jMax)
           w(1:,jb,1:,:) = ww0
           p(1:,jb,1:) = pp0
           if( viscous )   rlv(1:,jb,1:) = rlv0
           if( eddyModel ) rev(1:,jb,1:) = rev0
           idim = 2; ddim = jb

         case (kMin)
           w(1:,1:,0,:) = ww0 
           p(1:,1:,0) = pp0
           if( viscous )   rlv(1:,1:,0) = rlv0
           if( eddyModel ) rev(1:,1:,0) = rev0
           idim = 3; ddim = 0

         case (kMax)
           w(1:,1:,kb,:) = ww0
           p(1:,1:,kb) = pp0
           if( viscous )   rlv(1:,1:,kb) = rlv0
           if( eddyModel ) rev(1:,1:,kb) = rev0
           idim = 3; ddim = kb

       end select


       end subroutine resetww0pp0rlv0rev0Bwd   
       
       



