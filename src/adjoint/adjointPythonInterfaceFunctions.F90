! This files contains several routines used for interfacing with the
! python level. The routines are as follows:

! 3. getdrdwTVec: Multiply vec_in by dRdw^T to produce vec_out
! 6. getdFdxVec: Multiply vec_in by dFdx to produce vec_out

subroutine getdRdwTVec(in_vec, out_vec, ndof)

  use ADjointPETSc
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in) :: in_vec(ndof)
  real(kind=realType), intent(inout) :: out_vec(ndof)

  ! Working Variables
  integer(kind=intType) :: ierr

  ! We will use empty generic vectors psi_like1, psi_like2

  call VecPlaceArray(psi_like1, in_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(psi_like2, out_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call MatMult(dRdwT, psi_like1, psi_like2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(psi_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(psi_like2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine getdRdwTVec

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

