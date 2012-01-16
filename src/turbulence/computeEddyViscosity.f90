!
!      ******************************************************************
!      *                                                                *
!      * File:          computeEddyViscosity.f90                        *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeEddyViscosity
!
!      ******************************************************************
!      *                                                                *
!      * computeEddyViscosity computes the eddy viscosity in the        *
!      * owned cell centers of the given block. It is assumed that the  *
!      * pointes already point to the correct block before entering     *
!      * this subroutine.                                               *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputPhysics
       use iteration
       use blockPointers
       implicit none
!
!      Local variables.
!
       logical :: returnImmediately
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if an immediate return can be made.

       if( eddyModel ) then
         if((currentLevel <= groundLevel) .or. turbCoupled) then
           returnImmediately = .false.
         else
           returnImmediately = .true.
         endif
       else
         returnImmediately = .true.
       endif

       if( returnImmediately ) return

       ! Determine the turbulence model and call the appropriate
       ! routine to compute the eddy viscosity.

       select case (turbModel)

         case (baldwinLomax)
           call blEddyViscosity

         case (spalartAllmaras, spalartAllmarasEdwards)
           call saEddyViscosity

         case (komegaWilcox, komegaModified)
           call kwEddyViscosity

         case (menterSST)
           call SSTEddyViscosity

         case (ktau)
           call ktEddyViscosity

         case (v2f)
           call vfEddyViscosity

         case default
           call terminate("computeEddyViscosity", &
                          "Turbulence model not implemented yet")

       end select

       end subroutine computeEddyViscosity

!      ==================================================================

       subroutine vfEddyViscosity
!
!      ******************************************************************
!      *                                                                *
!      * vfEddyViscosity computes the eddy-viscosity according to the   *
!      * v2f turbulence model for the block given in blockPointers.     *
!      * This routine is for both the n=1 and n=6 version.              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use paramTurb
       use turbMod
       use bcTypes
       use inputPhysics
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, nn
       real(kind=realType)   :: tke, tep, tkea, tepa, tepl, tv2, tv2a
       real(kind=realType)   :: yp, utau

       real(kind=realType), dimension(itu1:itu5) :: tup
       real(kind=realType), dimension(:,:,:), pointer :: ww
       real(kind=realType), dimension(:,:),   pointer :: rrlv, rrev
       real(kind=realType), dimension(:,:),   pointer :: dd2Wall

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute time and length scale
 
       call vfScale
 
       ! Loop over the cells of this block and compute the eddy viscosity.
       ! Do not include halo's.

       do k=2,kl
         do j=2,jl
           do i=2,il
             tke   = w(i,j,k,itu1)
             tep   = w(i,j,k,itu2)
             tv2   = w(i,j,k,itu3)
             tkea  = abs(tke)
             tepa  = abs(tep)
             tv2a  = abs(tv2)
             tepl  = max(tepa,rvfLimitE)

             rev(i,j,k) = rvfCmu*w(i,j,k,irho)*tv2a/tepl*sct(i,j,k)

           enddo
         enddo
       enddo

       ! Modify the rhs of the 1st internal cell, if wall functions
       ! are used; their value is determined by the table.

       testWallFunctions: if( wallFunctions ) then

         bocos: do nn=1,nViscBocos

           ! Determine the block face on which the subface is located
           ! and set some variables. As flag points to the entire array
           ! flagI2, etc., its starting indices are the starting indices
           ! of its target and not 1.

           select case (bcFaceid(nn))
             case (iMin)
               ww      => w(2,1:,1:,1:);   rrlv => rlv(2,1:,1:)
               dd2Wall => d2Wall(2,:,:);   rrev => rev(2,1:,1:)

             case (iMax)
               ww      => w(il,1:,1:,1:);   rrlv => rlv(il,1:,1:)
               dd2Wall => d2Wall(il,:,:);   rrev => rev(il,1:,1:)

             case (jMin)
               ww      => w(1:,2,1:,1:);   rrlv => rlv(1:,2,1:)
               dd2Wall => d2Wall(:,2,:);   rrev => rev(1:,2,1:)

             case (jMax)
               ww      => w(1:,jl,1:,1:);   rrlv => rlv(1:,jl,1:)
               dd2Wall => d2Wall(:,jl,:);   rrev => rev(1:,jl,1:)

             case (kMin)
               ww      => w(1:,1:,2,1:);   rrlv => rlv(1:,1:,2)
               dd2Wall => d2Wall(:,:,2);   rrev => rev(1:,1:,2)

             case (kMax)
               ww      => w(1:,1:,kl,1:);   rrlv => rlv(1:,1:,kl)
               dd2Wall => d2Wall(:,:,kl);   rrev => rev(1:,1:,kl)

           end select

           ! Loop over the owned faces of this subface. Therefore the
           ! nodal range of bcData must be used. The offset of +1 is
           ! present, because the starting index of the cell range is
           ! 1 larger than the starting index of the nodal range.

           do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
             do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

               ! Enforce k and epsilon in the 1st internal cell from
               ! the wall function table. There is an offset of -1 in
               ! the wall distance. Note that the offset compared to
               ! the current value must be stored. Also note that the
               ! curve fits contain the non-dimensional values.

               utau = viscSubface(nn)%utau(i,j)
               yp = ww(i,j,irho)*dd2Wall(i-1,j-1)*utau/rrlv(i,j)

             ! call curveTupYp(tup(itu5:itu5), yp, itu5, itu5)
             ! rrev(i,j) = tup(itu5)*rrlv(i,j)
             enddo
           enddo

         enddo bocos
       endif testWallFunctions

       end subroutine vfEddyViscosity

!      ==================================================================

       subroutine saEddyViscosity
!
!      ******************************************************************
!      *                                                                *
!      * saEddyViscosity computes the eddy-viscosity according to the   *
!      * Spalart-Allmaras model for the block given in blockPointers.   *
!      * This routine for both the original version as well as the      *
!      * modified version according to Edwards.                         *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use paramTurb
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
       real(kind=realType)   :: chi, chi3, fv1, rnuSA, cv13
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the cv1^3; cv1 is a constant of the Spalart-Allmaras model.

       cv13 = rsaCv1**3

       ! Loop over the cells of this block and compute the eddy viscosity.
       ! Do not include halo's.

       do k=2,kl
         do j=2,jl
           do i=2,il
             rnuSA      = w(i,j,k,itu1)*w(i,j,k,irho)
             chi        = rnuSA/rlv(i,j,k)
             chi3       = chi**3
             fv1        = chi3/(chi3+cv13)
             rev(i,j,k) = fv1*rnuSA
           enddo
         enddo
       enddo

       end subroutine saEddyViscosity

!      ==================================================================

       subroutine kwEddyViscosity
!
!      ******************************************************************
!      *                                                                *
!      * kwEddyViscosity computes the eddy viscosity according to the   *
!      * k-omega models (both the original Wilcox as well as the        *
!      * modified version) for the block given in blockPointers.        *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the cells of this block and compute the eddy viscosity.
       ! Do not include halo's.

       do k=2,kl
         do j=2,jl
           do i=2,il
             rev(i,j,k) = abs(w(i,j,k,irho)*w(i,j,k,itu1)/w(i,j,k,itu2))
           enddo
         enddo
       enddo

       end subroutine kwEddyViscosity

!      ==================================================================

       subroutine SSTEddyViscosity
!
!      ******************************************************************
!      *                                                                *
!      * SSTEddyViscosity computes the eddy viscosity according to      *
!      * menter's SST variant of the k-omega turbulence model for the   *
!      * block given in blockPointers.                                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use paramTurb
       use turbMod
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k

       real(kind=realType) :: t1, t2, arg2, f2, vortMag
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the vorticity squared in the cell centers. The reason
       ! for computing the vorticity squared is that a routine exists
       ! for it; for the actual eddy viscosity computation the vorticity
       ! itself is needed.

       prod  => dw(1:,1:,1:,iprod)
       vort  => prod
       call prodWmag2

       ! Loop over the cells of this block and compute the eddy viscosity.
       ! Do not include halo's.

       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the value of the function f2, which occurs in the
             ! eddy-viscosity computation.

             t1 = two*sqrt(w(i,j,k,itu1)) &
                / (0.09_realType*w(i,j,k,itu2)*d2Wall(i,j,k))
             t2 = 500.0_realType*rlv(i,j,k) &
                / (w(i,j,k,irho)*w(i,j,k,itu2)*d2Wall(i,j,k)**2)

             arg2 = max(t1,t2)
             f2   = tanh(arg2**2)

             ! And compute the eddy viscosity.

             vortMag    = sqrt(vort(i,j,k))
             rev(i,j,k) = w(i,j,k,irho)*rSSTA1*w(i,j,k,itu1) &
                        / max(rSSTA1*w(i,j,k,itu2), f2*vortMag)

           enddo
         enddo
       enddo

       end subroutine SSTEddyViscosity

!      ==================================================================

       subroutine ktEddyViscosity
!
!      ******************************************************************
!      *                                                                *
!      * ktEddyViscosity computes the eddy viscosity according to the   *
!      * k-tau turbulence model for the block given in blockPointers.   *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the cells of this block and compute the eddy viscosity.
       ! Do not include halo's.

       do k=2,kl
         do j=2,jl
           do i=2,il
             rev(i,j,k) = abs(w(i,j,k,irho)*w(i,j,k,itu1)*w(i,j,k,itu2))
           enddo
         enddo
       enddo

       end subroutine ktEddyViscosity
