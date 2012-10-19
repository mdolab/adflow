!
!     ******************************************************************
!     *                                                                *
!     * File:          getdRdXvPsi.F90                                 *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 10-08-2010                                      *
!     * Last modified: 10-08-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getdRdXvPsi(dXv,ndof,adjoint,nstate)
#ifndef USE_NO_PETSC
  
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the current adjoint vector by dRdXv to get a vector   *
  !     * of length Xv. For the TimeSpectral case, we need to include    *
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
  !

#define PETSC_AVOID_MPIF_H
  use petscvec
  use ADjointPETSc, only: dRdx,xVec, w_like1
  use blockPointers
  use section
  use inputTimeSpectral 
  use monitor 
  implicit none

  integer(kind=intType), intent(in) :: ndof,nstate
  real(kind=realType), intent(out)  :: dXv(ndof)
  real(kind=realType), intent(in)   :: adjoint(nstate)

  integer(kind=intType) :: ierr,sps,i,nn,mm,counter0,counter1
  integer(kind=intType) :: nodes_on_block, cum_nodes_on_block
  real(kind=realType), dimension(3)   :: rotationPoint,r
  real(kind=realType), dimension(3,3) :: rotationMatrix  
  real(kind=realType) :: t(nSections),dt(nSections)
  real(kind=realType) :: tOld,tNew,pt(3)
  real(kind=realType),pointer :: xvec_pointer(:)
  real(kind=realType) :: time(3)

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

  ! If we only have 1 time instance (NOT TimeSpectral analysis, then
  ! xVec and gridVec would be the same, and dX_sps/dX would be the
  ! identity matrix and the calculation useless.

  if (nTimeIntervalsSpectral == 1) then 
     do i=1,ndof
        dXv(i) = xvec_pointer(i)
     end do
  else

     ! Zero dXv for time spectral case since we add to array.
     dXv = 0.0

     ! Now we loop over the number of timeInstances and reduce xVec
     ! into GridVec
  
     do nn=1,nSections
        dt(nn) = sections(nn)%timePeriod &
             / real(nTimeIntervalsSpectral,realType)
     enddo
  
     timeUnsteady = zero
     counter0 = 0
     cum_nodes_on_block = 0
     ! The nDom loop followed by the sps loop is required to follow
     ! the globalNode ordering such that we can use the pointer from
     ! vecGetArrayF90

     do nn=1,nDom
        do sps = 1,nTimeIntervalsSpectral

           call setPointers(nn,1,sps)
           nodes_on_block = il*jl*kl

           do mm=1,nSections
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
         

           do i=1,nodes_on_block

              pt = (/xvec_pointer(3*counter0+1),&
                     xvec_pointer(3*counter0+2),&
                     xvec_pointer(3*counter0+3)/)

              dXv(3*counter1+1:3*counter1+3) = &
                   dXv(3*counter1+1:3*counter1+3) + &
                   matmul(rotationMatrix,pt)
              
              counter0 = counter0 + 1
              counter1 = counter1 + 1
           end do

        end do
        ! Increment the cumulative number of nodes by the nodes on the
        ! block we just did
        cum_nodes_on_block = cum_nodes_on_block + nodes_on_block
     end do
  end if

  ! No longer need Vec and reset pointer
  call VecResetArray(w_like1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(xVec,xvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine getdRdXvPsi
