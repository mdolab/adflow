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
subroutine getdRdXvPsi(ndof,dXv)
  !
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
  use communication
  use ADjointPETSc, only: dRdx,psi,xVec,wVec
  use ADjointVars
  use blockPointers
  use warpingPETSC 
  use inputADjoint
  use section
  use inputTimeSpectral 
  use monitor 
 
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(out)  :: dXv(ndof)

  integer(kind=intType) :: ierr,sps,i,j,k,ind(3),nn,counter
  real(kind=realType), dimension(3)   :: rotationPoint,r
  real(kind=realType), dimension(3,3) :: rotationMatrix  
  real(kind=realType) :: t(nSections),dt(nSections)
  real(kind=realType) :: tOld,tNew,pt(3)

  ! Create a temporary vector which is the size of the grid*nTimeInstances
  call MatGetVecs(dRdx,xVec,wVec,ierr)
  call EChk(ierr,__file__,__line__)
  
  ! Do the matMultTranspose and put result into xVec
  call MatMultTranspose(dRdX,psi,xVec,ierr)
  call EChk(ierr,__file__,__line__)

  dXv = 0.0
  ! If we only have 1 time instance (NOT TimeSpectral analysis, then
  ! xVec and gridVec would be the same, and dX_sps/dX would be the
  ! identity matrix and the calculation useless.

  if (nTimeIntervalsSpectral == 1) then 
     counter = 0
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,1_intType)
        do k=1,kl
           do j=1,jl
              do i=1,il
                 ind = (/globalNode(i,j,k)*3  ,&
                        globalNode(i,j,k)*3+1,&
                        globalNode(i,j,k)*3+2/)
                 
                 call VecGetValues(xVec,3,ind,pt,ierr)
                 call EChk(ierr,__file__,__line__)
                 
                 dXv(counter*3+1:counter*3+3) = pt
                 
                 counter = counter + 1
                    
              end do
           end do
        end do
     end do ! Domain Loop
  else
  
     ! Now we loop over the number of timeInstances and reduce xVec
     ! into GridVec
  
     do nn=1,nSections
        dt(nn) = sections(nn)%timePeriod &
             / real(nTimeIntervalsSpectral,realType)
     enddo
  
     timeUnsteady = zero
    
     do sps = 1,nTimeIntervalsSpectral
        do nn=1,nSections
           t(nn) = (sps-1)*dt(nn)
        enddo
     
        ! Compute the displacements due to the rigid motion of the mesh.
     
        tNew = timeUnsteady + timeUnsteadyRestart
        tOld = tNew - t(1)

        call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

        ! Take rotation Matrix Transpose
        rotationMatrix = transpose(rotationMatrix)

        counter = 0
        do nn=1,nDom
           call setPointersAdj(nn,1_intType,sps)
           do k=1,kl
              do j=1,jl
                 do i=1,il
                    ind = (/globalNode(i,j,k)*3  ,&
                            globalNode(i,j,k)*3+1,&
                            globalNode(i,j,k)*3+2/)

                    call VecGetValues(xVec,3,ind,pt,ierr)
                    call EChk(ierr,__file__,__line__)

                    dXv(counter*3+1:counter*3+3) = dXv(counter*3+1:3*counter+3) + &
                         matmul(rotationMatrix,pt)

                    counter = counter + 1
                    
                 end do
              end do
           end do
        end do ! Domain Loop
     end do ! Sps loop
  end if

  ! No longer need xVec
  call VecDestroy(xVec,ierr)
  call EChk(ierr,__file__,__line__)

  call VecDestroy(wVec,ierr)
  call EChk(ierr,__file__,__line__)
  
end subroutine getdRdXvPsi
