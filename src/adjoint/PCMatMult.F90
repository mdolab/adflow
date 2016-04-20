subroutine testpc()

  use communication
  use adjointPETSc
  use adjointVars
  use flowvarrefstate
  use inputtimespectral
  use blockPointers
  implicit none
#ifndef USE_COMPLEX
  real(kind=realType) , dimension(:), allocatable :: v1, v2
  real(kind=realType) :: timeA, timeB, timeC, val
  integer(kind=intType) :: nDimw, ierr, i, nn, sps, j, k, ii, jj, kk, info
  logical :: useAD
  interface
     subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
          useObjective, frozenTurb, level, matrixTurb)
       use precision
       implicit none
#define PETSC_AVOID_MPIF_H
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

       Mat :: matrix
       Mat, optional :: matrixTurb
       ! Input Variables
       logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
       integer(kind=intType), intent(in) :: level
     end subroutine setupStateResidualMatrix
  end interface

  useAD = .True.
  call setupPCMatrix(useAD, .False., .False., 1)

  ! Now create all the factors
    call createPETscVars
  call setupStateResidualMatrix(drdwpret, usead, .True. , .False., &
       .False., .False., 1)

  nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  allocate(v1(nDimw), v2(nDimW))

  ! Set random into V1
  call random_number(v1)

  ! Dump psi into psi_like1 and RHS into psi_like2
  call VecPlaceArray(psi_like1, v1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(psi_like2, v2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call mpi_barrier(sumb_comm_world, ierr)
  timeA = mpi_wtime()
  call matMult(drdwpret, psi_like1, psi_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  timeB =mpi_wtime() 

  print *,'Petsc Time:', myid, timeB-timeA
  call vecNorm(psi_like2, NORM_2, val, ierr)
  print *,'PETsc Norm:', val


  call mpi_barrier(sumb_comm_world, ierr)
  timeA = mpi_wtime()
  call PCMatMult(drdwpret, psi_like1, psi_like2, ierr)
  timeB = mpi_wtime()

  print *,'My Time:', myid, timeB-timeA
  call vecNorm(psi_like2, NORM_2, val, ierr)
   print *,'My Norm:', val
  call VecResetArray(psi_like1, ierr)
  call VecResetArray(psi_like2, ierr)
#endif
end subroutine testpc

subroutine factorPCMatrix()

  use communication
  use adjointPETSc
  use adjointVars
  use flowvarrefstate
  use inputtimespectral
  use blockPointers
  implicit none

  integer(kind=intType) :: ierr, i, nn, sps, j, k, ii, jj, kk, info, timeA, timeB

  timeA = mpi_wtime()
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, 1, sps)
        
        ! ========= I-Lines ============
        do k=2, kl
           do j=2, jl
              
              ! Copy data from PCMat
              ii = 0
              jj = 0
              kk = 0
              do i=2, il
                 i_D_fact(:, nw*ii+1:nw*(ii+1), j, k) = PCMat(:,      1:1*nw, i, j, k)
                 ii = ii + 1

                 if (i > 2) then 
                    i_L_fact(:, nw*jj+1:nw*(jj+1), j, k) = PCMat(:, 1*nw+1:2*nw, i, j, k) 
                    jj = jj + 1
                 end if

                 if (i < il) then 
                    i_U_fact(:, nw*kk+1:nw*(kk+1), j, k) = PCMat(:, 2*nw+1:3*nw, i, j, k) 
                    kk = kk + 1
                 end if

              end do

              ! ! Perform factorization
              ! call dgeblttrf(nx, nw, i_D_fact(:, :, j, k), i_L_fact(:, :, j, k), &
              !      i_U_fact(:, :, j, k), i_U2_fact(:, :, j, k), i_ipiv(:, :, j, k), info)

           end do
        end do

        ! ========= J-Lines ============
        do k=2, kl
           do i=2, il
              
              ! Copy data from PCMat
              ii = 0
              jj = 0
              kk = 0
              do j=2, jl
                 j_D_fact(:, nw*ii+1:nw*(ii+1), i, k) = PCMat(:,      1:1*nw, i, j, k)
                 ii = ii + 1

                 if (j > 2) then 
                    j_L_fact(:, nw*jj+1:nw*(jj+1), i, k) = PCMat(:, 3*nw+1:4*nw, i, j, k) 
                    jj = jj + 1
                 end if

                 if (j < jl) then 
                    j_U_fact(:, nw*kk+1:nw*(kk+1), i, k) = PCMat(:, 4*nw+1:5*nw, i, j, k) 
                    kk =kk + 1
                 end if
              end do

              ! ! Perform factorization
              ! call dgeblttrf(ny, nw, j_D_fact(:, :, i, k), j_L_fact(:, :, i, k), &
              !      j_U_fact(:, :, i, k), j_U2_fact(:, :, i, k), j_ipiv(:, :, i, k), info)

           end do
        end do

        ! ========= k-Lines ============
        do j=2, jl
           do i=2, il
              
              ! Copy data from PCMat
              ii = 0 
              jj = 0
              kk = 0
              do k=2, kl
                 k_D_fact(:, nw*ii+1:nw*(ii+1), i, j) = PCMat(:,      1:1*nw, i, j, k)
                 ii = ii + 1

                 if (k > 2) then 
                    k_L_fact(:, nw*jj+1:nw*(jj+1), i, j) = PCMat(:, 5*nw+1:6*nw, i, j, k) 
                    jj = jj + 1
                 end if

                 if (k < kl) then 
                    k_U_fact(:, nw*kk+1:nw*(kk+1), i, j) = PCMat(:, 6*nw+1:7*nw, i, j, k) 
                    kk =kk + 1
                 end if
              end do

              ! ! Perform factorization
              ! call dgeblttrf(nz, nw, k_D_fact(:, :, i, j), k_L_fact(:, :, i,j), &
              !      k_U_fact(:, :, i, j), k_U2_fact(:, :, i, j), k_ipiv(:, :, i, j), info)

           end do
        end do
     end do
  end do

end subroutine factorPCMatrix

subroutine PCMatMult(A, vecX,  vecY, ierr)

  ! PETSc user-defied call back function for computing the product of
  ! PCMat with a vector.

  use constants
  use communication
  use blockPointers
  use iteration         
  use flowVarRefState     
  use inputAdjoint       
  use ADjointVars
  use inputTimeSpectral  

  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! PETSc Arguments
  Mat   A
  Vec   vecX, vecY 
  integer(kind=intType) ::ierr, i, j, k, l, nn, sps, ii
  real(kind=realType) :: sum1, sum2, sum3, sum4, sum5
  real(kind=realType) :: x1, x2, x3, x4, x5
  real(kind=realType), pointer :: yPtr(:)
  real(kind=realType) ::  RHS(nw*7)
  ! We first have to distribute xPtr into the halo cells for the local
  ! products.
  print *,'calling pcmatmult'
  call setPCVec(vecX)

  call VecGetArrayF90(vecY, yPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now we can compute the acutal matrix vector product.
  ii = 1
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, 1, sps)

        do k=2, kl
           do j=2, jl
              do i=2, il
    
                 ! Fill up the RHS

                 rhs(0*nw+1:1*nw) = PCVec1(:, i  , j, k)
                 rhs(1*nw+1:2*nw) = PCVec1(:, i-1, j, k)
                 rhs(2*nw+1:3*nw) = PCVec1(:, i+1, j, k)
                 rhs(3*nw+1:4*nw) = PCVec1(:, i, j-1, k)
                 rhs(4*nw+1:5*nw) = PCVec1(:, i, j+1, k)             
                 rhs(5*nw+1:6*nw) = PCVec1(:, i, j, k-1)
                 rhs(6*nw+1:7*nw) = PCVec1(:, i, j, k+1)

                 ! Call blass mat-vec. We can dump the result directly
                 ! into yPtr. There doesn't appear to be much
                 ! difference between the fortran matmul and blas for
                 ! these sized operations. 

                 !yPtr(ii:ii+nw) = matmul(PCMat(:, :, i, j, k), rhs)
                 call DGEMV('n', nw, nw*7, one, PCMat(:, :, i, j, k), nw, rhs, 1, zero, yPtr(ii), 1)
                 
                 ii = ii + nw

              end do
           end do
        end do
     end do
  end do

  call VecRestoreArrayF90(vecY, yPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ierr = 0

end subroutine PCMatMult

subroutine myShellPCApply(pc, vecX, vecY, ierr)

  use precision
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use communication
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! PETSc Arguments
  PC   pc
  Vec   vecX, vecY 
  integer(kind=intType) :: ierr, info, i, j, k, ii, nn, sps, ipiv(nw)
  real(kind=realType) :: blk(nw, nw), rhs(nw)
  real(kind=realType), pointer :: yPtr(:)


  ! First copy X to Y. This way we will continually transform vecY
  ! into the preconditioned vector.  
  call vecCopy(vecX, vecY, ierr)

  call VecGetArrayF90(vecY, yPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Do something useful here....

  call VecRestoreArrayF90(vecY, yPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine myShellPCApply

subroutine setPCVec(vecX)

  use constants
  use communication
  use blockPointers
  use inputTimeSpectral  
  use flowVarRefState
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif
  
  ! PETSc Arguments
  Vec   vecX
  integer(kind=intType) ::ierr, i, j, k, l, nn, mm, sps, nVar, size, index
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2, ii, jj, procID
  real(kind=realType), pointer :: xPtr(:)
  type(commType) :: commPattern
  type(internalCommType) :: internal
  integer, dimension(mpi_status_size) :: status
       
  call VecGetArrayReadF90(vecX, xPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! First set all the owned cells...this is basically just a straight
  ! copy in order
  ii = 0
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, 1, sps)
        do k=2, kl
           do j=2, jl
              do i=2, il
                 do l=1, nw
                    ii = ii + 1
                    flowDoms(nn, 1, sps)%PCVec1(l, i, j, k) = xPtr(ii)
                 end do
              end do
           end do
        end do
     end do
  end do

  ! Done with the vecX.
  call VecRestorearrayReadF90(vecX, xPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now we do a custom halo exchange. We can use the same pattern as
  ! whlo, but the ordering of the unknowns is different so we do our
  ! own.

  nVar = nw
  internal = internalCell_1st(1)
  commPattern = commPatternCell_1st(1)

  spectralModes: do mm=1,nTimeIntervalsSpectral

     ! Send the variables. The data is first copied into
     ! the send buffer after which the buffer is sent asap.

     ii = 1
     sends: do i=1,commPattern%nProcSend

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPattern%sendProc(i)
        size    = nVar*commPattern%nsend(i)
        
        ! Copy the data in the correct part of the send buffer.
        
        jj = ii
        do j=1,commPattern%nsend(i)

           ! Store the block id and the indices of the donor
           ! a bit easier.
           
           d1 = commPattern%sendList(i)%block(j)
           i1 = commPattern%sendList(i)%indices(j,1)
           j1 = commPattern%sendList(i)%indices(j,2)
           k1 = commPattern%sendList(i)%indices(j,3)
           
           ! Copy the given range of the working variables for
           ! this cell in the buffer. Update the counter jj.
           
           do k=1, nw
              sendBuffer(jj) = flowDoms(d1,1,mm)%PCVec1(k, i1, j1, k1)
              jj = jj + 1
           enddo

        enddo

        ! Send the data.
        
        call mpi_isend(sendBuffer(ii), size, sumb_real, procID,  &
             procID, SUmb_comm_world, sendRequests(i), &
             ierr)

        ! Set ii to jj for the next processor.

        ii = jj
        
     enddo sends

     ! Post the nonblocking receives.

     ii = 1
     receives: do i=1,commPattern%nProcRecv
        
        ! Store the processor id and the size of the message
        ! a bit easier.
        
        procID = commPattern%recvProc(i)
        size    = nVar*commPattern%nrecv(i)
        
        ! Post the receive.
        
        call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
             myID, SUmb_comm_world, recvRequests(i), ierr)
        
        ! And update ii.
        
        ii = ii + size
        
     enddo receives

     ! Copy the local data.
     
     localCopy: do i=1,internal%ncopy

        ! Store the block and the indices of the donor a bit easier.

        d1 = internal%donorBlock(i)
        i1 = internal%donorIndices(i,1)
        j1 = internal%donorIndices(i,2)
        k1 = internal%donorIndices(i,3)

        ! Idem for the halo's.

        d2 = internal%haloBlock(i)
        i2 = internal%haloIndices(i,1)
        j2 = internal%haloIndices(i,2)
        k2 = internal%haloIndices(i,3)

        ! Copy the given range of working variables.
        
        do k=1, nw
           flowDoms(d2,1,mm)%PCVec1(k, i2, j2, k2) = &
                flowDoms(d1,1,mm)%PCVec1(k, i1, j1, k1)
        enddo
        
     enddo localCopy

     ! Complete the nonblocking receives in an arbitrary sequence and
     ! copy the variables from the buffer into the halo's.
     
     size = commPattern%nProcRecv
     completeRecvs: do i=1,commPattern%nProcRecv
        
        ! Complete any of the requests.

        call mpi_waitany(size, recvRequests, index, status, ierr)
        
        ! Copy the data just arrived in the halo's.
        
        ii = index
        jj = nVar*commPattern%nrecvCum(ii-1)
        do j=1,commPattern%nrecv(ii)

           ! Store the block and the indices of the halo a bit easier.
           
           d2 = commPattern%recvList(ii)%block(j)
           i2 = commPattern%recvList(ii)%indices(j,1)
           j2 = commPattern%recvList(ii)%indices(j,2)
           k2 = commPattern%recvList(ii)%indices(j,3)

           do k=1, nw
              jj = jj + 1
              flowDoms(d2,1,mm)%PCVec1(k, i2, j2, k2) = recvBuffer(jj)
           enddo
           
        enddo

     enddo completeRecvs

     ! Complete the nonblocking sends.
     
     size = commPattern%nProcSend
     do i=1,commPattern%nProcSend
        call mpi_waitany(size, sendRequests, index, status, ierr)
     enddo
     
  enddo spectralModes

end subroutine setPCVec
