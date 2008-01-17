!
!      ******************************************************************
!      *                                                                *
!      * File:          unsteadyTurbSpectral.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-22-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine unsteadyTurbSpectral(ntu1, ntu2)
!
!      ******************************************************************
!      *                                                                *
!      * unsteadyTurbSpectral determines the spectral time derivative   *
!      * for all owned cells. This routine is called before the actual  *
!      * solve routines, such that the treatment is identical for all   *
!      * spectral solutions. The results is stored in the corresponding *
!      * entry in dw.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputPhysics
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: ntu1, ntu2
!
!      Local variables.
!
       integer(kind=intType) :: ii, mm, nn, sps, i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if not the time spectral equations are to
       ! be solved.

       if(equationMode /= timeSpectral) return

       ! Loop over the number of spectral modes and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, currentLevel, sps)

           ! Loop over the number of turbulent transport equations.

           nAdvLoop: do ii=ntu1, ntu2

             ! Initialize the time derivative to zero for the owned
             ! cell centers.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   dw(i,j,k,ii) = zero
                 enddo
               enddo
             enddo

             ! Loop over the number of terms which contribute to the
             ! time derivative.

             do mm=1,nTimeIntervalsSpectral

               ! Add the contribution to the time derivative for
               ! all owned cells.

               do k=2,kl
                 do j=2,jl
                   do i=2,il
                     dw(i,j,k,ii) = dw(i,j,k,ii)              &
                                  + dscalar(sectionID,sps,mm) &
                                  * flowDoms(nn,currentLevel,mm)%w(i,j,k,ii)
                   enddo
                 enddo
               enddo

             enddo

           enddo nAdvLoop
         enddo domains
       enddo spectralLoop

      end subroutine unsteadyTurbSpectral
