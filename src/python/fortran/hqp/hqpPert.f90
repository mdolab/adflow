!
!      ******************************************************************
!      *                                                                *
!      * File:          hqpPert.f90                                     *
!      * Author:        Jaina                                           *
!      * Starting date: 11-01-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine hqpPert(x,vel,nNodes,t)
!
!      ******************************************************************
!      *                                                                *
!      * Add grid velocities to test gust convection.                   *
!      *                                                                *
!      * nNodes:        Number of nodes at which data is requested.     *
!      * x(3,nNodes):   locations at which grid velocity is requested.  *
!      * vel(3,nNodes): velocities at those locations.                  *
!      * t:             time.                                           *
!      *                                                                *
!      ******************************************************************
!
       use precision
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nNodes
       real(kind=realType), dimension(3,nNodes),intent(in) :: x
       real(kind=realType), dimension(3,nNodes),intent(out) :: vel
       real(kind=realType), intent(in) :: t

!      Local variables

       integer(kind=intType) :: i
       real(kind=realType) :: angmax,fac,x0,vinf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       x0=-5.                   ! initialisation position of gust in chords
       angmax=0.0174            ! of magnitude 1 deg
       vinf=winf(ivx)

       do i=1,nNodes
 
          fac=(x(1,i)-x0-vinf*t/timeRef)

          if (fac > 0.) then
 
             vel(1,i)=0.
             vel(2,i)=0.
             vel(3,i)=0.
          else
 
             vel(1,i)=0.
             vel(2,i)=0.
             vel(3,i)=-angmax*vinf
          endif

       enddo

       return
       end
