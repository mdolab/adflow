subroutine getCurrentResidual(rhoRes,totalRRes)
  use communication
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use iteration
  use inputIteration
  implicit none
  ! Compute the current resdiual of w
  real(kind=realType), intent(out) :: rhoRes,totalRRes
  real(kind=realType) :: ovv,r_sum,rho_sum
  integer(kind=intType) :: sps,nn,i,j,k,l,ierr
  currentLevel = 1
  groundLevel = 1
  rkStage = 0

  ! We have to explictly halo exchange the states since
  ! computeResidualNK no longer does this

  call whalo2(1_intType, 1_intType, nw, .False., &
       .False.,.False.)

  call computeResidualNK()

  r_sum = 0.0
  rho_sum = 0.0
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 ovv = 1/vol(i,j,k)
                 do l=1,nw
                    r_sum = r_sum + (dw(i,j,k,l)*ovv)**2
                 end do
                 rho_sum = rho_sum + (dw(i,j,k,irho)*ovv)**2
              end do
           end do
        end do
     end do
  end do

  totalRRes = zero ! Just to make valgrind happy...

  call mpi_allreduce(r_sum,totalRRes,1,sumb_real,mpi_sum,&
       SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call mpi_allreduce(rho_sum,rhoRes,1,sumb_real,mpi_sum,&
       SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ! curRes now has the inverse-volume weighted sum squared of the
  ! residuals, finally take the squareRoot

 totalRRes = sqrt(totalRRes)
 rhoRes = sqrt(rhoRes/nCellGlobal(currentLevel))

end subroutine getCurrentResidual
