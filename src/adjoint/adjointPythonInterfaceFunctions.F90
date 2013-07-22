! This files contains several routines used for interfacing with the
! python level. The routines are as follows:

! 1. setAdjoint: Set the petsc adjoint vector from variables stored in Python
! 2. getAdjoint: Returns the variables in petsc adjoint vector to Python
! 3. getdrdwTVec: Multiply vec_in by dRdw^T to produce vec_out
! 4. getdRdaPsi: Multiply dRda^T*adjoint where adjoint is supplied, and output returned
! 5. getdRdxVPsi: Compute product dRdXv^T*psi and return result in dXv
! 6. getdFdxVec: Multiply vec_in by dFdx to produce vec_out
! 7. getdFdxTVec: Multiple vec_in by dFdx^T to produce vec_out
! 8. agumentRHS: Agument RHS of adjoint by dRdw^T*phi, where phi is supplied
! 9. getdFdwTVec: Return out_vec = dFdw^T*in_vec. 

subroutine setADjoint(nstate, adjoint)

#ifndef USE_NO_PETSC	
#define PETSC_AVOID_MPIF_
#include "finclude/petscdef.h"

  use ADjointPETSc, only : psi
  use petscvec
  use constants

  implicit none

  !  Input Variables
  integer(kind=intType),intent(in):: nstate
  real(kind=realType),dimension(nstate), intent(in) :: adjoint

  ! Local Variables
  integer(kind=intType) :: i, ierr
  real(kind=realType),pointer :: psi_pointer(:)

  ! Copy out adjoint vector:
  call VecGetArrayF90(psi, psi_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do a straight copy:
  do i=1,nstate
     psi_pointer(i) = adjoint(i)
  end do

  call VecRestoreArrayF90(psi, psi_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

end subroutine setADjoint

subroutine getADjoint(nstate, adjoint)

#ifndef USE_NO_PETSC	
#define PETSC_AVOID_MPIF_
#include "finclude/petscdef.h"

  use ADjointPETSc, only : psi
  use petscvec
  use constants

  implicit none

  ! Local Variables
  integer(kind=intType), intent(in):: nstate
  real(kind=realType), dimension(nstate), intent(out) :: adjoint

  ! Local Variables
  integer(kind=intType) :: i, ierr
  real(kind=realType), pointer :: psi_pointer(:)

  ! Copy out adjoint vector:
  call VecGetArrayF90(psi, psi_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do a straight copy:
  do i=1, nstate
     adjoint(i) = psi_pointer(i)
  end do

  call VecRestoreArrayF90(psi, psi_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

end subroutine getADjoint

subroutine getdRdwTVec(in_vec, out_vec, ndof)
#ifndef USE_NO_PETSC

  use ADjointPETSc
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in) :: in_vec(ndof)
  real(kind=realType), intent(inout) :: out_vec

  ! Working Variables
  integer(kind=intType) :: ierr

  ! We will use empty generic vectors w_like1, w_like2

  call VecPlaceArray(w_like1, in_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(w_like2, out_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call MatMult(dRdwT, w_like1, w_like2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(w_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(w_like2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

#endif
end subroutine getdRdwTVec

subroutine getdRdaPsi(output, ndv, adjoint, nstate)

#ifndef USE_NO_PETSC

  use communication
  use ADjointPETSc
  use ADjointVars
  use blockPointers
  use inputADjoint
  use section
  use inputTimeSpectral 
  use monitor 

  implicit none

  ! Input/Output Variables
  integer(kind=intType), intent(in) :: ndv, nstate
  real(kind=realType), intent(in) :: adjoint(nstate)
  real(kind=realType), intent(out) :: output(ndv)

  ! Local Variables
  integer(kind=intType) :: ierr, i

  ! put adjoint arry in w_like1
  call VecPlaceArray(w_like1, adjoint, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the result vector for dRda^T * psi
  call VecCreateMPI(SUMB_COMM_WORLD, PETSC_DECIDE, ndv, dRdaTPsi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do the Multiplication
  call MatMultTranspose(dRda, w_like1, dRdaTPsi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! This is a little wonkly, since dRdaTPsi only contains a handful of
  ! variables.
  ! Use a proper vec scatter here:
  call VecScatterCreateToAll(dRdaTPsi, dRdaTpsi_scatter, dRdaTPsi_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterBegin(dRdaTpsi_scatter, dRdaTPsi, dRdaTpsi_local, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecScatterEnd (dRdaTpsi_scatter, dRdaTPsi, dRdaTpsi_local, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now just pluck off the local values
  do i=1, nDv
     call VecGetValues(dRdaTPsi_local, 1, i-1, output(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do

  ! No longer need vectors/scatter
  call VecScatterDestroy(dRdaTpsi_scatter, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(dRdaTPsi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(dRdaTpsi_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(w_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine getdRdaPsi
!
subroutine getdRdXvPsi(dXv, ndof, adjoint, nstate)
#ifndef USE_NO_PETSC
 
#define PETSC_AVOID_MPIF_H
  use petscvec
  use ADjointPETSc, only: dRdx, xVec, w_like1
  use blockPointers
  use inputTimeSpectral 
  implicit none

  ! Input/Output Variables
  integer(kind=intType), intent(in) :: ndof, nstate
  real(kind=realType), intent(out)  :: dXv(ndof)
  real(kind=realType), intent(in)   :: adjoint(nstate)

  ! Local Variables
  integer(kind=intType) :: ierr, sps, i
   real(kind=realType), pointer :: xvec_pointer(:)

  ! Place adjoint in Vector
  call VecPlaceArray(w_like1, adjoint, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do the matMult with dRdx and put result into xVec. NOTE dRdx is
  ! already transposed and thus we just do a matMult NOT
  ! a matMultTranspose

  call MatMult(dRdx, w_like1, xVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Extract pointer for xVec
  call VecGetArrayF90(xVec, xvec_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Copy out the values to return
  do i=1, ndof
     dXv(i) = xvec_pointer(i)
  end do

  ! Reset the arrays
  call VecResetArray(w_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(xVec, xvec_pointer, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine getdRdXvPsi

subroutine spectralPrecscribedMotion(input, nin, dXv, nout)

  use blockPointers
  use section
  use inputTimeSpectral 
  use monitor 
  implicit none
  ! Input/Output Variables
  integer(kind=intType), intent(in) :: nin, nout
  real(kind=realType), intent(out)  :: dXv(nout)
  real(kind=realType), intent(in)   :: input(nin)

  ! Local Variables
  integer(kind=intType) :: ierr, sps, i, nn, mm, counter0, counter1
  integer(kind=intType) :: nodes_on_block, cum_nodes_on_block
  real(kind=realType), dimension(3)   :: rotationPoint, r
  real(kind=realType), dimension(3, 3) :: rotationMatrix  
  real(kind=realType) :: t(nSections), dt(nSections)
  real(kind=realType) :: tOld, tNew, pt(3)
  real(kind=realType), pointer :: xvec_pointer(:)
  real(kind=realType) :: time(3)
 
  !       For the TimeSpectral case, we need to include    *
  !     * the operation that rotates the base grid to each time instance *
  !     * This is basically the reverse of the operation that is done in *
  !     * setGrid.f90                                                    *
  !     *                                                                *
  !     * The operation in setGrid.f90 is the following                  *
  !     *                                                                *
  !     * X_sps = M(X - rotPoint) + rotPoint                             *
  !     *                                                                *
  !     * where                                                          *
  !     * X_sps is the set of coordinates at each time instance          *
  !     * M is the rotation matrix calculated by rotMatrixRigidBody      *
  !     * rotPoint is the point about which the motion takes place       *
  !     *                                                                *
  !     * It is easy to see dX_sps/dX = M                                *
  !     *                                                                *
  !     * What we are actually computing is the following:               *
  !     *                                                                *
  !     *            T          T                                        *
  !     *   /dX_sps \ /   dR   \                                         *
  !     *   |-------| |------- |  psi                                    *
  !     *   \  dX   / \ dX_sps /                                         *
  !     *                                                                *
  !     ******************************************************************
  
  ! Zero dXv for time spectral case since we add to array.
  dXv = zero
  
  ! Now we loop over the number of timeInstances and reduce xVec
  ! into GridVec
  
  do nn=1, nSections
     dt(nn) = sections(nn)%timePeriod &
          / real(nTimeIntervalsSpectral, realType)
  enddo
  
  timeUnsteady = zero
  counter0 = 0
  cum_nodes_on_block = 0
  ! The nDom loop followed by the sps loop is required to follow
  ! the globalNode ordering such that we can use the pointer from
  ! vecGetArrayF90

  do nn=1, nDom
     do sps = 1, nTimeIntervalsSpectral

        call setPointers(nn, 1, sps)
        nodes_on_block = il*jl*kl
        
        do mm=1, nSections
           t(mm) = (sps-1)*dt(mm)
        enddo
        
        ! Compute the displacements due to the rigid motion of the mesh.
        
        tNew = timeUnsteady + timeUnsteadyRestart
        tOld = tNew - t(1)

        call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)
        
        ! Take rotation Matrix Transpose
        rotationMatrix = transpose(rotationMatrix)
        
        counter1 = cum_nodes_on_block        
        
        ! Loop over the localally owned nodes:
        do i=1, nodes_on_block
           pt = (/input(3*counter0+1), &
                input(3*counter0+2), &
                input(3*counter0+3)/)
           
           dXv(3*counter1+1:3*counter1+3) = &
                dXv(3*counter1+1:3*counter1+3) + &
                matmul(rotationMatrix, pt)
           
           counter0 = counter0 + 1
           counter1 = counter1 + 1
        end do

     end do
     ! Increment the cumulative number of nodes by the nodes on the
     ! block we just did
     cum_nodes_on_block = cum_nodes_on_block + nodes_on_block
  end do

end subroutine spectralPrecscribedMotion

subroutine getdFdxVec(ndof, vec_in, vec_out)
#ifndef USE_NO_PETSC

  use communication
  use ADjointPETSc, only: dFdx, fVec1, fVec2, PETSC_VIEWER_STDOUT_SELF
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in)  :: vec_in(ndof)
  real(kind=realType), intent(out)  :: vec_out(ndof)
  integer(kind=intType) :: ierr
  vec_out(:) = 0.0

  call VecPlaceArray(fVec1, vec_in, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(fVec2, vec_out, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatMult(dFdx, fVec1, fVec2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(fVec1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(fVec2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

#endif
end subroutine getdFdxVec

subroutine getdFdxTVec(ndof, vec_in, vec_out)
#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply vec_in by dFdx to produce vec_out                     *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc, only: dFdx, fVec1, fVec2
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in)  :: vec_in(ndof)
  real(kind=realType), intent(out)  :: vec_out(ndof)
  integer(kind=intType) :: ierr
  vec_out(:) = 0.0
  
  call VecPlaceArray(fVec1, vec_in, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(fVec2, vec_out, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatMultTranspose(dFdx, fVec1, fVec2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(fVec1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(fVec2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

#endif
end subroutine getdFdxTVec

subroutine agumentRHS(ndof, phi)
#ifndef USE_NO_PETSC 

  use ADjointPETSc 
  use ADjointVars 
  use communication
  use inputADjoint

  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in) :: phi(ndof)
  integer(kind=intType) :: ierr

  call VecPlaceArray(fVec1, phi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! ------------- OLD Code using explit dFdw ---------
  ! ! Dump the result into adjointRHS
  ! call MatMultTranspose(dFdw, fVec1, adjointRHS, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)
  ! -------------------------------------------------

  ! New code using dFcdw computation from forward mode Assembly. This
  ! function requires the use of forward mode AD

  !w = x * y : VecPointwiseMult(Vec w, Vec x,Vec y)
  call VecPointwiseMult(fNode, fVec1, overArea, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatMultTranspose(dFndFc, fNode, fCell, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatMultTranspose(dFcdw, fCell, adjointRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call vecResetArray(fVec1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

#endif

end subroutine agumentRHS

subroutine getdFdwTVec(in_vec, in_dof, out_vec, out_dof)
#ifndef USE_NO_PETSC 

  use ADjointPETSc
  use communication

  implicit none

  ! Input/Ouput
  integer(kind=intType), intent(in) :: in_dof, out_dof
  real(kind=realType), intent(in) :: in_vec(in_dof)
  real(kind=realType), intent(inout) :: out_vec(out_dof)

  ! Working
  integer(kind=intType) :: ierr

  ! Put petsc wrapper around arrays
  call VecPlaceArray(fVec1, in_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecPlaceArray(w_like1, out_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Dump the result into adjointRHS since the we want to ADD this
  ! result to w_like1 below

  ! ------------ OldMethod
  ! call MatMultTranspose(dFdw, fVec1, adjointRHS, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  ! ------------ New Method
  call VecPointwiseMult(fNode, fVec1, overArea, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatMultTranspose(dFndFc, fNode, fCell, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatMultTranspose(dFcdw, fCell, adjointRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)
 
  ! do: w_like1 = w_like1 + adjointRHS
  call VecAxpy(w_like1, one, adjointRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call vecResetArray(fVec1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(w_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

end subroutine getdFdwTVec
