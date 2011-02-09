!
!      ******************************************************************
!      *                                                                *
!      * File:          turbBCNSWall.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-30-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine turbBCNSWall(secondHalo)
!
!      ******************************************************************
!      *                                                                *
!      * turbBCNSWall applies the viscous wall boundary conditions      *
!      * of the turbulent transport equations to a block. It is assumed *
!      * that the pointers in blockPointers are already set to the      *
!      * correct block on the correct grid level.                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Subroutine argument.
!
       logical, intent(in) :: secondHalo
!
!      Local variables.
!
       integer(kind=intType) :: nn, i, j, l, m

       real(kind=realType), dimension(:,:,:,:), pointer :: bmt
       real(kind=realType), dimension(:,:,:), pointer :: bvt
       real(kind=realType), dimension(:,:,:), pointer :: ww0, ww1, ww2
       real(kind=realType), dimension(:,:), pointer :: rev0, rev1, rev2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the viscous subfaces of this block.

       bocos: do nn=1,nViscBocos

         ! Set the corresponding arrays.

         call BCTurbWall(nn)

         ! Determine the block face on which this subface is located
         ! and set some pointers accordingly.

         select case (BCFaceID(nn))
           case (iMin)
             bmt => bmti1; bvt => bvti1
             ww0 => w(0 ,1:,1:,:); ww1 => w(1 ,1:,1:,:)
             ww2 => w(2 ,1:,1:,:)

             if( eddyModel ) then
               rev0 => rev(0 ,1:,1:); rev1 => rev(1 ,1:,1:)
               rev2 => rev(2 ,1:,1:)
             endif

           case (iMax)
             bmt => bmti2; bvt => bvti2
             ww0 => w(ib,1:,1:,:); ww1 => w(ie,1:,1:,:)
             ww2 => w(il,1:,1:,:)

             if( eddyModel ) then
               rev0 => rev(ib,1:,1:); rev1 => rev(ie,1:,1:)
               rev2 => rev(il,1:,1:)
             endif

           case (jMin)
             bmt => bmtj1; bvt => bvtj1
             ww0 => w(1:,0 ,1:,:); ww1 => w(1:,1 ,1:,:)
             ww2 => w(1:,2 ,1:,:)

             if( eddyModel ) then
               rev0 => rev(1:,0 ,1:); rev1 => rev(1:,1 ,1:)
               rev2 => rev(1:,2 ,1:)
             endif

           case (jMax)
             bmt => bmtj2; bvt => bvtj2
             ww0 => w(1:,jb,1:,:); ww1 => w(1:,je,1:,:)
             ww2 => w(1:,jl,1:,:)

             if( eddyModel ) then
               rev0 => rev(1:,jb,1:); rev1 => rev(1:,je,1:)
               rev2 => rev(1:,jl,1:)
             endif

           case (kMin)
             bmt => bmtk1; bvt => bvtk1
             ww0 => w(1:,1:,0 ,:); ww1 => w(1:,1:,1 ,:)
             ww2 => w(1:,1:,2 ,:)

             if( eddyModel ) then
               rev0 => rev(1:,1:, 0); rev1 => rev(1:,1:, 1)
               rev2 => rev(1:,1:, 2)
             endif

           case (kMax)
             bmt => bmtk2; bvt => bvtk2
             ww0 => w(1:,1:,kb,:); ww1 => w(1:,1:,ke,:)
             ww2 => w(1:,1:,kl,:)

             if( eddyModel ) then
               rev0 => rev(1:,1:,kb); rev1 => rev(1:,1:,ke)
               rev2 => rev(1:,1:,kl)
             endif
         end select

         ! Loop over the faces and set the state in
         ! the turbulent halo cells.

         do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
             do l=nt1,nt2
               ww1(i,j,l) = bvt(i,j,l)
               do m=nt1,nt2
                 ww1(i,j,l) = ww1(i,j,l) - bmt(i,j,l,m)*ww2(i,j,m)
               enddo
             enddo
           enddo
         enddo

         ! Use constant extrapolation if the state in the second halo
         ! must be computed.

         if( secondHalo ) then
           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
               do l=nt1,nt2
                 ww0(i,j,l) = ww1(i,j,l)
               enddo
             enddo
           enddo
         endif

         ! Set the eddy viscosity for an eddy model.

         if( eddyModel ) then
           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
               rev1(i,j) = -rev2(i,j)
             enddo
           enddo

           if( secondHalo ) then
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                 rev0(i,j) = rev1(i,j)
               enddo
             enddo
           endif
         endif

       enddo bocos

       end subroutine turbBCNSWall
