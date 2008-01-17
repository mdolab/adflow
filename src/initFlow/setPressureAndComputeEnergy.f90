!
!      ******************************************************************
!      *                                                                *
!      * File:          setPressureAndComputeEnergy.f90                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-07-2005                                      *
!      * Last modified: 09-13-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setPressureAndComputeEnergy(halosRead)
!
!      ******************************************************************
!      *                                                                *
!      * Due to the usage of the variable IOVar, which generalizes the  *
!      * IO and leads to reuse of code, currently the pressure is       *
!      * stored at the position of rhoE. In this routine that data is   *
!      * copied to the pressure array and the total energy is computed. *
!      * Note that this routine is only called when a restart is done.  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: halosRead
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, nHalo
       integer(kind=intType) :: i, j, k
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the value of nHalo, depending whether or not the halo cells
       ! have been read from the restart file.

       nHalo = 0
       if( halosRead ) nHalo = 1

       ! Loop over the number of time instances and the local blocks.
       ! As this routine is only called when a restart is performed,
       ! the MG start level is the finest level.

       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom

           ! Set the pointers to the correct block. As this routine is
           ! only called when a restart is performed, the MG start level
           ! is the finest level.

           call setPointers(nn,1_intType,sps)

           ! Determine the range for which the pressure must be computed.

           iBeg = 2 -nHalo; jBeg =  2-nHalo; kBeg =  2-nHalo
           iEnd = il+nHalo; jEnd = jl+nHalo; kEnd = kl+nHalo

           ! Copy the pressure for the required cells.

           do k=kBeg,kEnd
             do j=jBeg,jEnd
               do i=iBeg,iEnd
                 p(i,j,k) = w(i,j,k,irhoE)
               enddo
             enddo
           enddo

           ! Compute the total energy as well.

           call computeEtot(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd, kPresent)
         enddo
       enddo

       end subroutine setPressureAndComputeEnergy
