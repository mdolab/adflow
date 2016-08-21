
!
!      ******************************************************************
!      *                                                                *
!      * File:          computeSpeedOfSoundSquared.F90                  *
!      * Author:        Gaetan K.W. Kenway                              *
!      * Starting date: 01-20-2014                                      *
!      * Last modified: 01-20-2014                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine computeSpeedOfSoundSquared

  !
  !      ******************************************************************
  !      *                                                                *
  !      * computeSpeedOfSoundSquared does what it says.                  *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants
  use blockPointers, only : ie, je, ke, w, p, aa, gamma
  implicit none
  !
  !      Local variables.
  !
  real(kind=realType), PARAMETER :: twothird=two*third
  integer(kind=intType) :: i, j, k, ii
  real(kind=realType) :: pp
  logical :: correctForK, getCorrectForK

  ! Determine if we need to correct for K
  correctForK = getCorrectForK()

  if (correctForK) then 
#ifdef TAPENADE_REVERSE
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
                 aa(i,j,k) = gamma(i,j,k)*pp/w(i,j,k,irho)
#ifdef TAPENADE_REVERSE
              end do
#else
           enddo
        enddo
     enddo
#endif   
  else
#ifdef TAPENADE_REVERSE
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
                 aa(i,j,k) = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
#ifdef TAPENADE_REVERSE
              end do
#else
           enddo
        enddo
     enddo
#endif   
  end if
end subroutine computeSpeedOfSoundSquared
