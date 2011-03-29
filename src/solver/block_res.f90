! This is a super-combined function that combines the original
! functionality of: 
! Pressure Computation
! timeStep
! applyAllBCs
! initRes
! residual 

! The real difference between this and the original modules is that it
! it only operates on a single block at a time and as such the
! block/sps loop is outside the calculation. This routine is suitable
! for forward mode AD with Tapenade

subroutine block_res(nn,sps)
  
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use flowVarRefState     ! nw
  use inputPhysics 
  use iteration
  implicit none

  real(kind=realType) :: gm1,v2
  integer(kind=intType) :: nn,sps,i,j,k
  ! Compute the pressures
  gm1 = gammaConstant - one
  
  ! Compute P 
  do k=0,kb
     do j=0,jb
        do i=0,ib
           v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                + w(i,j,k,ivz)**2
           p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                - half*w(i,j,k,irho)*v2)
           p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
        enddo
     enddo
  enddo

  ! Apply all BC's
  call applyAllBC_block(.true.)
  
  ! Compute skin_friction Velocity
  call computeUtau_block

  ! Compute time step and spectral radius
  call timeStep_block(.false.)

  if( equations == RANSEquations ) then
     call initres_block(nt1MG, nMGVar,nn,sps) ! Initialize only the Turblent Variables
     call turbResidual_block
  endif

  ! Initialize the remaining flow variables
  call initRes_block(1,nwf,nn,sps)

  ! Actual residual calc
  call residual_block

end subroutine block_res
