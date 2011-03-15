!
!      ******************************************************************
!      *                                                                *
!      * File:          computePressureAdj.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-19-2006                                      *
!      * Last modified: 03-20-2006                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine computePressureNKPC(wAdj, pAdj,nn,level,sps,sps2)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Simple routine to compute the pressure from the variables w.   *
  !      * A calorically perfect gas, i.e. constant gamma, is assumed.    *
  !      *                                                                *
  !      ******************************************************************
  !
  use flowVarRefState
  use inputPhysics
  use inputTimeSpectral !nIntervalTimespectral
  implicit none
  !
  !      Subroutine arguments
  !

  integer(kind=intType)::nn,level,sps,sps2
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) :: pAdj
  !
  !      Local variables
  !
  integer(kind=intType) :: i, j, k

  real(kind=realType) :: gm1, factK, v2

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  gm1 = gammaConstant - one

  ! Check the situation.

  if( kPresent ) then

     ! A separate equation for the turbulent kinetic energy is
     ! present. This variable must be taken into account.

     factK = five*third - gammaConstant

     do k=-2,2
        do j=-2,2
           do i=-2,2
              v2 = wAdj(i,j,k,ivx,sps2)**2 + wAdj(i,j,k,ivy,sps2)**2 &
                   + wAdj(i,j,k,ivz,sps2)**2

              pAdj(i,j,k,sps2) = gm1*(wAdj(i,j,k,irhoE,sps2)    &
                   - half*wAdj(i,j,k,irho,sps2)*v2) &
                   + factK*wAdj(i,j,k,irho,sps2)*wAdj(i,j,k,itu1,sps2)
           enddo
        enddo
     enddo

  else

     ! No separate equation for the turbulent kinetic enery.
     ! Use the standard formula.

     do k=-2,2
        do j=-2,2
           do i=-2,2
              v2 = wAdj(i,j,k,ivx,sps2)**2 + wAdj(i,j,k,ivy,sps2)**2 &
                   + wAdj(i,j,k,ivz,sps2)**2

              pAdj(i,j,k,sps2) = gm1*(wAdj(i,j,k,irhoE,sps2) &
                   - half*wAdj(i,j,k,irho,sps2)*v2)
           enddo
        enddo
     enddo

  endif

end subroutine computePressureNKPC
