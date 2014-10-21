!
!      ******************************************************************
!      *                                                                *
!      * File:          setww0pp0rlv0rev0Bwd.f90                        *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
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
           ww0 = w(0,1:,1:,:); pp0 = p(0,1:,1:)
           if( viscous )   rlv0 = rlv(0,1:,1:)
           if( eddyModel ) rev0 = rev(0,1:,1:)
           idim = 1; ddim = 0

         case (iMax)
           ww0 = w(ib,1:,1:,:); pp0 = p(ib,1:,1:)
           if( viscous )   rlv0 = rlv(ib,1:,1:)
           if( eddyModel ) rev0 = rev(ib,1:,1:)
           idim = 1; ddim = ib

         case (jMin)
           ww0 = w(1:,0,1:,:); pp0 = p(1:,0,1:)
           if( viscous )   rlv0 = rlv(1:,0,1:)
           if( eddyModel ) rev0 = rev(1:,0,1:)
           idim = 2; ddim = 0

         case (jMax)
           ww0 = w(1:,jb,1:,:); pp0 = p(1:,jb,1:)
           if( viscous )   rlv0 = rlv(1:,jb,1:)
           if( eddyModel ) rev0 = rev(1:,jb,1:)
           idim = 2; ddim = jb

         case (kMin)
           ww0 = w(1:,1:,0,:); pp0 = p(1:,1:,0)
           if( viscous )   rlv0 = rlv(1:,1:,0)
           if( eddyModel ) rev0 = rev(1:,1:,0)
           idim = 3; ddim = 0

         case (kMax)
           ww0 = w(1:,1:,kb,:); pp0 = p(1:,1:,kb)
           if( viscous )   rlv0 = rlv(1:,1:,kb)
           if( eddyModel ) rev0 = rev(1:,1:,kb)
           idim = 3; ddim = kb

       end select


       end subroutine setww0pp0rlv0rev0Bwd   
       
       



