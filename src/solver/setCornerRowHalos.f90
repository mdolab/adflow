!
!      ******************************************************************
!      *                                                                *
!      * File:          setCornerRowHalos.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-20-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setCornerRowHalos(nVar)
!
!      ******************************************************************
!      *                                                                *
!      * setCornerRowHalos initializes the halo's next to corner row    *
!      * halo's, such that it contains some values. Otherwise it may    *
!      * be uninitialized or cause a floating point exception, as this  *
!      * memory is also used to compute the mg corrections.             *
!      * It is assumed that the pointers in blockPointers already       *
!      * point to the correct block.                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: nVar
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l, mm, ll
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Halo's on the i=iMin and i=iMax plane.                         *
!      *                                                                *
!      ******************************************************************
!
       ! K-rows.

       mm = min(3_intType,jl)
       ll = max(2_intType,ny)

       do k=2,kl
         do l=1,nVar
           w(1, 2, k,l) = w(2, 2, k,l)
           w(1, mm,k,l) = w(2, mm,k,l)
           w(1, jl,k,l) = w(2, jl,k,l)
           w(1, ll,k,l) = w(2, ll,k,l)
           w(ie,2, k,l) = w(il,2, k,l)
           w(ie,mm,k,l) = w(il,mm,k,l)
           w(ie,jl,k,l) = w(il,jl,k,l)
           w(ie,ll,k,l) = w(il,ll,k,l)
         enddo

         p(1, 2, k) = p(2, 2, k)
         p(1, mm,k) = p(2, mm,k)
         p(1, jl,k) = p(2, jl,k)
         p(1, ll,k) = p(2, ll,k)
         p(ie,2, k) = p(il,2, k)
         p(ie,mm,k) = p(il,mm,k)
         p(ie,jl,k) = p(il,jl,k)
         p(ie,ll,k) = p(il,ll,k)

         if( viscous) then
           rlv(1, 2, k) = rlv(2, 2, k)
           rlv(1, mm,k) = rlv(2, mm,k)
           rlv(1, jl,k) = rlv(2, jl,k)
           rlv(1, ll,k) = rlv(2, ll,k)
           rlv(ie,2, k) = rlv(il,2, k)
           rlv(ie,mm,k) = rlv(il,mm,k)
           rlv(ie,jl,k) = rlv(il,jl,k)
           rlv(ie,ll,k) = rlv(il,ll,k)
         endif

         if( eddyModel ) then
           rev(1, 2, k) = rev(2, 2, k)
           rev(1, mm,k) = rev(2, mm,k)
           rev(1, jl,k) = rev(2, jl,k)
           rev(1, ll,k) = rev(2, ll,k)
           rev(ie,2, k) = rev(il,2, k)
           rev(ie,mm,k) = rev(il,mm,k)
           rev(ie,jl,k) = rev(il,jl,k)
           rev(ie,ll,k) = rev(il,ll,k)
         endif
       enddo

       ! J-rows; no need to include the corners. These have been set in
       ! the previous k-loop.

       mm = min(3_intType,kl)
       ll = max(2_intType,nz)

       do j=3,ny
         do l=1,nVar
           w(1, j,2, l) = w(2, j,2, l)
           w(1, j,mm,l) = w(2, j,mm,l)
           w(1, j,kl,l) = w(2, j,kl,l)
           w(1, j,ll,l) = w(2, j,ll,l)
           w(ie,j,2, l) = w(il,j,2, l)
           w(ie,j,mm,l) = w(il,j,mm,l)
           w(ie,j,kl,l) = w(il,j,kl,l)
           w(ie,j,ll,l) = w(il,j,ll,l)
         enddo

         p(1, j, 2) = p(2, j, 2)
         p(1, j,mm) = p(2, j,mm)
         p(1, j,kl) = p(2, j,kl)
         p(1, j,ll) = p(2, j,ll)
         p(ie,j, 2) = p(il,j, 2)
         p(ie,j,mm) = p(il,j,mm)
         p(ie,j,kl) = p(il,j,kl)
         p(ie,j,ll) = p(il,j,ll)

         if( viscous) then
           rlv(1, j, 2) = rlv(2, j, 2)
           rlv(1, j,mm) = rlv(2, j,mm)
           rlv(1, j,kl) = rlv(2, j,kl)
           rlv(1, j,ll) = rlv(2, j,ll)
           rlv(ie,j, 2) = rlv(il,j, 2)
           rlv(ie,j,mm) = rlv(il,j,mm)
           rlv(ie,j,kl) = rlv(il,j,kl)
           rlv(ie,j,ll) = rlv(il,j,ll)
         endif

         if( eddyModel ) then
           rev(1, j, 2) = rev(2, j, 2)
           rev(1, j,mm) = rev(2, j,mm)
           rev(1, j,kl) = rev(2, j,kl)
           rev(1, j,ll) = rev(2, j,ll)
           rev(ie,j, 2) = rev(il,j, 2)
           rev(ie,j,mm) = rev(il,j,mm)
           rev(ie,j,kl) = rev(il,j,kl)
           rev(ie,j,ll) = rev(il,j,ll)
         endif
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Halo's on the j=jMin and j=jMax plane.                         *
!      *                                                                *
!      ******************************************************************
!
       ! K-rows; no need to include the corners; this is done in the
       ! next i-loop.

       mm = min(3_intType,il)
       ll = max(2_intType,nx)

       do k=3,nz
         do l=1,nVar
           w(2, 1, k,l) = w(2, 2, k,l)
           w(mm,1, k,l) = w(mm,2, k,l)
           w(il,1, k,l) = w(il,2, k,l)
           w(ll,1, k,l) = w(ll,2, k,l)
           w(2, je,k,l) = w(2, jl,k,l)
           w(mm,je,k,l) = w(mm,jl,k,l)
           w(il,je,k,l) = w(il,jl,k,l)
           w(ll,je,k,l) = w(ll,jl,k,l)
         enddo

         p(2, 1, k) = p(2, 2, k)
         p(mm,1, k) = p(mm,2, k)
         p(il,1, k) = p(il,2, k)
         p(ll,1, k) = p(ll,2, k)
         p(2, je,k) = p(2, jl,k)
         p(mm,je,k) = p(mm,jl,k)
         p(il,je,k) = p(il,jl,k)
         p(ll,je,k) = p(ll,jl,k)

         if( viscous) then
           rlv(2, 1, k) = rlv(2, 2, k)
           rlv(mm,1, k) = rlv(mm,2, k)
           rlv(il,1, k) = rlv(il,2, k)
           rlv(ll,1, k) = rlv(ll,2, k)
           rlv(2, je,k) = rlv(2, jl,k)
           rlv(mm,je,k) = rlv(mm,jl,k)
           rlv(il,je,k) = rlv(il,jl,k)
           rlv(ll,je,k) = rlv(ll,jl,k)
         endif

         if( eddyModel ) then
           rev(2, 1, k) = rev(2, 2, k)
           rev(mm,1, k) = rev(mm,2, k)
           rev(il,1, k) = rev(il,2, k)
           rev(ll,1, k) = rev(ll,2, k)
           rev(2, je,k) = rev(2, jl,k)
           rev(mm,je,k) = rev(mm,jl,k)
           rev(il,je,k) = rev(il,jl,k)
           rev(ll,je,k) = rev(ll,jl,k)
         endif
       enddo

       ! I-rows, including halo's set on the iMin and iMax plane.

       mm = min(3_intType,kl)
       ll = max(2_intType,nz)

       do i=1,ie
         do l=1,nVar
           w(i,1, 2, l) = w(i,2, 2, l)
           w(i,1, mm,l) = w(i,2, mm,l)
           w(i,1, kl,l) = w(i,2, kl,l)
           w(i,1, ll,l) = w(i,2, ll,l)
           w(i,je,2, l) = w(i,jl,2, l)
           w(i,je,mm,l) = w(i,jl,mm,l)
           w(i,je,kl,l) = w(i,jl,kl,l)
           w(i,je,ll,l) = w(i,jl,ll,l)
         enddo

         p(i, 1, 2) = p(i, 2, 2)
         p(i, 1,mm) = p(i, 2,mm)
         p(i, 1,kl) = p(i, 2,kl)
         p(i, 1,ll) = p(i, 2,ll)
         p(i,je, 2) = p(i,jl, 2)
         p(i,je,mm) = p(i,jl,mm)
         p(i,je,kl) = p(i,jl,kl)
         p(i,je,ll) = p(i,jl,ll)

         if( viscous) then
           rlv(i, 1, 2) = rlv(i, 2, 2)
           rlv(i, 1,mm) = rlv(i, 2,mm)
           rlv(i, 1,kl) = rlv(i, 2,kl)
           rlv(i, 1,ll) = rlv(i, 2,ll)
           rlv(i,je, 2) = rlv(i,jl, 2)
           rlv(i,je,mm) = rlv(i,jl,mm)
           rlv(i,je,kl) = rlv(i,jl,kl)
           rlv(i,je,ll) = rlv(i,jl,ll)
         endif

         if( eddyModel ) then
           rev(i, 1, 2) = rev(i, 2, 2)
           rev(i, 1,mm) = rev(i, 2,mm)
           rev(i, 1,kl) = rev(i, 2,kl)
           rev(i, 1,ll) = rev(i, 2,ll)
           rev(i,je, 2) = rev(i,jl, 2)
           rev(i,je,mm) = rev(i,jl,mm)
           rev(i,je,kl) = rev(i,jl,kl)
           rev(i,je,ll) = rev(i,jl,ll)
         endif
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Halo's on the k=kMin and k=kMax plane.                         *
!      *                                                                *
!      ******************************************************************
!
       ! J-rows, including halo's set on the jMin and jMax plane.

       mm = min(3_intType,il)
       ll = max(2_intType,nx)

       do j=1,je
         do l=1,nVar
           w(2, j,1, l) = w(2, j,2, l)
           w(mm,j,1, l) = w(mm,j,2, l)
           w(il,j,1, l) = w(il,j,2, l)
           w(ll,j,1, l) = w(ll,j,2, l)
           w(2, j,ke,l) = w(2, j,kl,l)
           w(mm,j,ke,l) = w(mm,j,kl,l)
           w(il,j,ke,l) = w(il,j,kl,l)
           w(ll,j,ke,l) = w(ll,j,kl,l)
         enddo

         p( 2,j, 1) = p( 2,j, 2)
         p(mm,j, 1) = p(mm,j, 2)
         p(il,j, 1) = p(il,j, 2)
         p(ll,j, 1) = p(ll,j, 2)
         p( 2,j,ke) = p( 2,j,kl)
         p(mm,j,ke) = p(mm,j,kl)
         p(il,j,ke) = p(il,j,kl)
         p(ll,j,ke) = p(ll,j,kl)

         if( viscous) then
           rlv( 2,j, 1) = rlv( 2,j, 2)
           rlv(mm,j, 1) = rlv(mm,j, 2)
           rlv(il,j, 1) = rlv(il,j, 2)
           rlv(ll,j, 1) = rlv(ll,j, 2)
           rlv( 2,j,ke) = rlv( 2,j,kl)
           rlv(mm,j,ke) = rlv(mm,j,kl)
           rlv(il,j,ke) = rlv(il,j,kl)
           rlv(ll,j,ke) = rlv(ll,j,kl)
         endif

         if( eddyModel ) then
           rev( 2,j, 1) = rev( 2,j, 2)
           rev(mm,j, 1) = rev(mm,j, 2)
           rev(il,j, 1) = rev(il,j, 2)
           rev(ll,j, 1) = rev(ll,j, 2)
           rev( 2,j,ke) = rev( 2,j,kl)
           rev(mm,j,ke) = rev(mm,j,kl)
           rev(il,j,ke) = rev(il,j,kl)
           rev(ll,j,ke) = rev(ll,j,kl)
         endif
       enddo

       ! I-rows, including halo's set on the iMin and iMax plane.

       mm = min(3_intType,jl)
       ll = max(2_intType,ny)

       do i=1,ie
         do l=1,nVar
           w(i, 2, 1,l) = w(i, 2, 2,l)
           w(i,mm, 1,l) = w(i,mm, 2,l)
           w(i,jl, 1,l) = w(i,jl, 2,l)
           w(i,ll, 1,l) = w(i,ll, 2,l)
           w(i, 2,ke,l) = w(i, 2,kl,l)
           w(i,mm,ke,l) = w(i,mm,kl,l)
           w(i,jl,ke,l) = w(i,jl,kl,l)
           w(i,ll,ke,l) = w(i,ll,kl,l)
         enddo

         p(i, 2, 1) = p(i, 2, 2)
         p(i,mm, 1) = p(i,mm, 2)
         p(i,jl, 1) = p(i,jl, 2)
         p(i,ll, 1) = p(i,ll, 2)
         p(i, 2,ke) = p(i, 2,kl)
         p(i,mm,ke) = p(i,mm,kl)
         p(i,jl,ke) = p(i,jl,kl)
         p(i,ll,ke) = p(i,ll,kl)

         if( viscous) then
           rlv(i, 2, 1) = rlv(i, 2, 2)
           rlv(i,mm, 1) = rlv(i,mm, 2)
           rlv(i,jl, 1) = rlv(i,jl, 2)
           rlv(i,ll, 1) = rlv(i,ll, 2)
           rlv(i, 2,ke) = rlv(i, 2,kl)
           rlv(i,mm,ke) = rlv(i,mm,kl)
           rlv(i,jl,ke) = rlv(i,jl,kl)
           rlv(i,ll,ke) = rlv(i,ll,kl)
         endif

         if( eddyModel ) then
           rev(i, 2, 1) = rev(i, 2, 2)
           rev(i,mm, 1) = rev(i,mm, 2)
           rev(i,jl, 1) = rev(i,jl, 2)
           rev(i,ll, 1) = rev(i,ll, 2)
           rev(i, 2,ke) = rev(i, 2,kl)
           rev(i,mm,ke) = rev(i,mm,kl)
           rev(i,jl,ke) = rev(i,jl,kl)
           rev(i,ll,ke) = rev(i,ll,kl)
         endif
       enddo

       end subroutine setCornerRowHalos
