!
!      ******************************************************************
!      *                                                                *
!      * File:          turb2ndHalo.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-16-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine turb2ndHalo(nn)
!
!      ******************************************************************
!      *                                                                *
!      * turb2ndHalo sets the turbulent variables in the second halo    *
!      * cell for the given subface. Simple constant extrapolation is   *
!      * used to avoid problems.                                        *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l

       real(kind=realType), dimension(:,:,:), pointer :: ww0,   ww1
       real(kind=realType), dimension(:,:),   pointer :: rrev0, rrev1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the face on which this subface is located and set
       ! some pointers accordingly.

       select case (BCFaceID(nn))
         case (iMin)
           ww0 => w(0,1:,1:,:); ww1 => w(1,1:,1:,:)

           if( eddyModel ) then
             rrev0 => rev(0,1:,1:); rrev1 => rev(1,1:,1:)
           endif

         !===============================================================

         case (iMax)
           ww0 => w(ib,1:,1:,:); ww1 => w(ie,1:,1:,:)

           if( eddyModel ) then
             rrev0 => rev(ib,1:,1:); rrev1 => rev(ie,1:,1:)
           endif

         !===============================================================

         case (jMin)
           ww0 => w(1:,0,1:,:); ww1 => w(1:,1,1:,:)

           if( eddyModel ) then
             rrev0 => rev(1:,0,1:); rrev1 => rev(1:,1,1:)
           endif

         !===============================================================

         case (jMax)
           ww0 => w(1:,jb,1:,:); ww1 => w(1:,je,1:,:)

           if( eddyModel ) then
             rrev0 => rev(1:,jb,1:); rrev1 => rev(1:,je,1:)
           endif

         !===============================================================

         case (kMin)
           ww0 => w(1:,1:,0,:); ww1 => w(1:,1:,1,:)

           if( eddyModel ) then
             rrev0 => rev(1:,1:,0); rrev1 => rev(1:,1:,1)
           endif

         !===============================================================

         case (kMax)
           ww0 => w(1:,1:,kb,:); ww1 => w(1:,1:,ke,:)

           if( eddyModel ) then
             rrev0 => rev(1:,1:,kb); rrev1 => rev(1:,1:,ke)
           endif
       end select

       ! Loop over the faces of the subface.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
         do i=BCData(nn)%icBeg, BCData(nn)%icEnd

           ! Loop over the turbulent variables and set the second halo
           ! value. If this is an eddy model, also set the eddy viscosity.

           do l=nt1,nt2
             ww0(i,j,l) = ww1(i,j,l)
           enddo

           if( eddyModel ) rrev0(i,j) = rrev1(i,j)

         enddo
       enddo

       end subroutine turb2ndHalo
