!
!      ******************************************************************
!      *                                                                *
!      * File:          setUniformFlow.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-13-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setUniformFlow
!
!      ******************************************************************
!      *                                                                *
!      * setUniformFlow set the flow variables of all local blocks on   *
!      * the start level to the uniform flow field.                     *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use communication
       use flowVarRefState
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       use utils, only : setPointers
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, i, j, k, l

       real(kind=realType) :: tmp

       real(kind=realType), dimension(3) :: dirLoc, dirGlob
!
!      Interfaces
!
       interface
         subroutine velMagnAndDirectionSubface(vmag, dir, &
                                               BCData, mm)
           use block
           implicit none

           integer(kind=intType), intent(in) :: mm
           real(kind=realType), intent(out) :: vmag
           real(kind=realType), dimension(3), intent(inout) :: dir
           type(BCDataType), dimension(:), pointer :: BCData
         end subroutine velMagnAndDirectionSubface
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions and blocks.
       spectralLoop: do mm=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn,mgStartlevel,mm)

           ! Set the w-variables to the ones of the uniform flow field.
           do l=1,nw
             do k=0,kb
               do j=0,jb
                 do i=0,ib
                   w(i,j,k,l) = wInf(l)
                   dw(i,j,k,l) = zero
                 enddo
               enddo
             enddo
           enddo
           !set this here for a reinitialize flow to eliminate possible NAN's
           do l=1,nwf
             do k=0,kb
               do j=0,jb
                 do i=0,ib
                   fw(i,j,k,l) = zero
                 enddo
               enddo
             enddo
           enddo

           ! Set the pressure.

           p = pInfCorr

           ! Initialize the laminar and eddy viscosity, if appropriate,
           ! such that no uninitialized memory is present.

           if( viscous )   rlv = muInf
           if( eddyModel ) rev = eddyVisInfRatio*muInf

         enddo domains
       enddo spectralLoop

       ! Correct for the time spectral method in combination with an
       ! internal flow computation the velocity direction.
       ! It is possible that the prescribed direction is different
       ! for every time instance.

       testCorrection: if(equationMode == timeSpectral .and. &
                          flowType     == internalFlow) then

         ! Loop over the number of spectral solutions.

         spectralLoopCorr: do mm=1,nTimeIntervalsSpectral

           ! Initialize the local direction to zero. In the loop
           ! below this direction will be accumulated.

           dirLoc = zero

           ! Loop over the number of blocks and determine the average
           ! velocity direction prescribed at the inlets.

           domainLoop1: do nn=1,nDom

             ! Set the pointer for the BCData to make the code more
             ! readable.

             BCData => flowDoms(nn,mgStartlevel,mm)%BCData

             ! Loop over the number of boundary subfaces and update
             ! the local prescribed velocity direction.

             do l=1,flowDoms(nn,mgStartlevel,mm)%nBocos
               call velMagnAndDirectionSubface(tmp, dirLoc, BCData, l)
             enddo
           enddo domainLoop1

           ! Determine the sum of dirLoc and create a unit vector
           ! for the global direction.

           call mpi_allreduce(dirLoc, dirGlob, 3, sumb_real, mpi_sum, &
                              SUmb_comm_world, ierr)

           tmp         = one/max(eps,sqrt(dirGlob(1)**2 &
                       +                  dirGlob(2)**2 &
                       +                  dirGlob(3)**2))
           dirGlob(1) = tmp*dirGlob(1)
           dirGlob(2) = tmp*dirGlob(2)
           dirGlob(3) = tmp*dirGlob(3)

           ! Loop again over the local domains and correct the
           ! velocity direction.

           domainsLoop2: do nn=1,nDom

             ! Set the pointers for this block.

             call setPointers(nn,mgStartlevel,mm)

             do k=0,kb
               do j=0,jb
                 do i=0,ib
                   tmp          = sqrt(w(i,j,k,ivx)**2 &
                                +      w(i,j,k,ivy)**2 &
                                +      w(i,j,k,ivz)**2)
                   w(i,j,k,ivx) = tmp*dirGlob(1)
                   w(i,j,k,ivy) = tmp*dirGlob(2)
                   w(i,j,k,ivz) = tmp*dirGlob(3)
                 enddo
               enddo
             enddo

           enddo domainsLoop2
         enddo spectralLoopCorr

       endif testCorrection

       end subroutine setUniformFlow

!=================================================================

subroutine velMagnAndDirectionSubface(vmag, dir, BCData, mm)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * VelMagnAndDirectionSubface determines the maximum value    *
  !      * of the magnitude of the velocity as well as the sum of the     *
  !      * flow directions for the currently active subface.              *
  !      *                                                                *
  !      ******************************************************************
  !
  use block
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: mm

  real(kind=realType), intent(out) :: vmag
  real(kind=realType), dimension(3), intent(inout) :: dir

  type(BCDataType), dimension(:), pointer :: BCData
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j
  real(kind=realType)   :: vel
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Initialize vmag to -1.0.

  vmag = -one

  ! Check if the velocity is prescribed.

  if( associated(BCData(mm)%velx) .and. &
       associated(BCData(mm)%vely) .and. &
       associated(BCData(mm)%velz) ) then

     ! Loop over the owned faces of the subface. As the cell range
     ! may contain halo values, the nodal range is used.

     do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
        do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

           ! Compute the magnitude of the velocity and compare it
           ! with the current maximum. Store the maximum of the two.

           vel  = sqrt(BCData(mm)%velx(i,j)**2 &
                +      BCData(mm)%vely(i,j)**2 &
                +      BCData(mm)%velz(i,j)**2)
           vmag = max(vmag, vel)

           ! Compute the unit vector of the velocity and add it to dir.

           vel    = one/max(eps,vel)
           dir(1) = dir(1) + vel*BCData(mm)%velx(i,j)
           dir(2) = dir(2) + vel*BCData(mm)%vely(i,j)
           dir(3) = dir(3) + vel*BCData(mm)%velz(i,j)

        enddo
     enddo
  endif

  ! Check if the velocity direction is prescribed.

  if( associated(BCData(mm)%flowXdirInlet) .and. &
       associated(BCData(mm)%flowYdirInlet) .and. &
       associated(BCData(mm)%flowZdirInlet) ) then

     ! Add the unit vectors to dir by looping over the owned
     ! faces of the subfaces. Again the nodal range must be
     ! used for this.

     do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
        do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

           dir(1) = dir(1) + BCData(mm)%flowXdirInlet(i,j)
           dir(2) = dir(2) + BCData(mm)%flowYdirInlet(i,j)
           dir(3) = dir(3) + BCData(mm)%flowZdirInlet(i,j)

        enddo
     enddo

  endif

end subroutine velMagnAndDirectionSubface
