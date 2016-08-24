!
!      ******************************************************************
!      *                                                                *
!      * File:          computeLamViscosity.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine computeLamViscosity
  !
  !      ******************************************************************
  !      *                                                                *
  !      * computeLamViscosity computes the laminar viscosity ratio in    *
  !      * the owned cell centers of the given block. Sutherland's law is *
  !      * used. It is assumed that the pointes already point to the      *
  !      * correct block before entering this subroutine.                 *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use constants
  use flowVarRefState
  use inputPhysics
  use iteration
  use utils, only : getCorrectForK
  implicit none
  !
  !      Local parameter.
  !
  real(kind=realType), parameter :: twoThird = two*third
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, ii
  real(kind=realType)   :: muSuth, TSuth, SSuth, T, pp
  logical               :: correctForK
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Return immediately if no laminar viscosity needs to be computed.

  if(.not. viscous ) return

  ! Determine whether or not the pressure must be corrected
  ! for the presence of the turbulent kinetic energy.

  correctForK = getCorrectForK()

  ! Compute the nonDimensional constants in sutherland's law.

  muSuth = muSuthDim/muRef
  TSuth  = TSuthDim/Tref
  SSuth  = SSuthDim/Tref

  ! Substract 2/3 rho k, which is a part of the normal turbulent
  ! stresses, in case the pressure must be corrected.

  if( correctForK ) then
#ifdef TAPENADE_FAST
     !$AD II-LOOP
     do ii=0,ie*je*ke - 1
        i = mod(ii, ie) + 1
        j = mod(ii/ie, je) + 1
        k = ii/(ie*je) + 1
#else
        do k=1,ke
           do j=1,je
              do i=1,ie
#endif             
                 pp = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)

                 T = pp/(RGas*w(i,j,k,irho))
                 rlv(i,j,k) = muSuth*((TSuth + SSuth)/(T + SSuth)) &
                      * ((T/TSuth)**1.5_realType)
#ifdef TAPENADE_FAST
              end do
#else
           enddo
        enddo
     enddo
#endif   
  else
     ! Loop over the owned cells *AND* first level halos of this
     ! block and compute the laminar viscosity ratio.
#ifdef TAPENADE_FAST
     !$AD II-LOOP
     do ii=0,ie*je*ke - 1
        i = mod(ii, ie) + 1
        j = mod(ii/ie, je) + 1
        k = ii/(ie*je) + 1
#else
        do k=1,ke
           do j=1,je
              do i=1,ie
#endif             
                 ! Compute the nonDimensional temperature and the
                 ! nonDimensional laminar viscosity.

                 T = p(i,j,k)/(RGas*w(i,j,k,irho))
                 rlv(i,j,k) = muSuth*((TSuth + SSuth)/(T + SSuth)) &
                      * ((T/TSuth)**1.5_realType)
#ifdef TAPENADE_FAST
              end do
#else
           enddo
        enddo
     enddo
#endif   
  end if
end subroutine computeLamViscosity

