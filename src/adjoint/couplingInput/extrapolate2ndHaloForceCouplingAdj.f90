!
!      ******************************************************************
!      *                                                                *
!      * File:          extrapolate2ndHaloForcesAdj.f90                 *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi,C.A.(Sandy) Mader                  *
!      * Starting date: 03-21-2006                                      *
!      * Last modified: 06-09-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine extrapolate2ndHaloForceCouplingAdj(nn,icBeg, icEnd, jcBeg, jcEnd, &
                                        wAdj0, wAdj1, wAdj2,           &
                                        pAdj0, pAdj1, pAdj2)
!
!      ******************************************************************
!      *                                                                *
!      * extrapolate2ndHaloAdj determines the states of the second      *
!      * layer halo cells of subface nn of the block to which the       *
!      * pointers in blockPointers currently point.                     *
!      *                                                                *
!      ******************************************************************
!
       use constants       !irhoE
       use flowVarRefState !gammaInf, kPresent
       use iteration       !nt1MG, nt2MG
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       integer(kind=intType), intent(in) :: icBeg, icEnd, jcBeg, jcEnd
       
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend,nw):: wAdj0, wAdj1
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend,nw):: wAdj2
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend) :: pAdj0, pAdj1
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend) :: pAdj2

!
!      Local parameter.
!
       real(kind=realType), parameter :: factor = 0.5_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l, ii, jj

       real(kind=realType) :: ovgm1, gm53, factK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Easier storage of variables involving gamma.

       ovgm1 = one/(gammaInf - one)
       gm53  =  gammaInf - five*third
       factK = -ovgm1*gm53

       ! Loop over the generic subface to set the state in the
       ! halo cells.

       do j=jcBeg, jcEnd
         do i=icBeg, icEnd

           ii = i 
           jj = j 

           ! Extrapolate the density, velocities and pressure.
           ! Make sure that a certain threshold is kept for the
           ! density and pressure.

           wAdj0(ii,jj,irho) = two*wAdj1(ii,jj,irho) - wAdj2(ii,jj,irho)
           wAdj0(ii,jj,irho) = max(factor*wAdj1(ii,jj,irho), &
                                          wAdj0(ii,jj,irho))

           wAdj0(ii,jj,ivx) = two*wAdj1(ii,jj,ivx) - wAdj2(ii,jj,ivx)
           wAdj0(ii,jj,ivy) = two*wAdj1(ii,jj,ivy) - wAdj2(ii,jj,ivy)
           wAdj0(ii,jj,ivz) = two*wAdj1(ii,jj,ivz) - wAdj2(ii,jj,ivz)

           !print *, 'HALO: ii, jj =', ii, jj
           pAdj0(ii,jj) = two*pAdj1(ii,jj) - pAdj2(ii,jj)
           pAdj0(ii,jj) = max(factor*pAdj1(ii,jj), pAdj0(ii,jj))

           ! Extrapolate the turbulent variables. Use constant
           ! extrapolation.

           do l=nt1MG,nt2MG
             wAdj0(i,j,l) = wAdj1(i,j,l)
           enddo

           ! Compute the total energy.

           wAdj0(ii,jj,irhoE) = ovgm1*pAdj0(ii,jj)       &
                              + half*wAdj0(ii,jj,irho)   &
                              *     (wAdj0(ii,jj,ivx)**2 &
                              +      wAdj0(ii,jj,ivy)**2 &
                              +      wAdj0(ii,jj,ivz)**2)

           if( kPresent )                            &
             wAdj0(ii,jj,irhoE) = wAdj0(ii,jj,irhoE) &
                                - factK*wAdj0(ii,jj,irho) &
                                *       wAdj0(ii,jj,itu1)
         enddo
       enddo

     end subroutine extrapolate2ndHaloForceCouplingAdj
