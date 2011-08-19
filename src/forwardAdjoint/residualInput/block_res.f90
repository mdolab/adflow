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
  use inputTimeSpectral

  implicit none

  real(kind=realType) :: gm1,v2
  integer(kind=intType) :: nn,sps,i,j,k,sps2,mm,l
  logical :: correctForK

  ! Compute the pressures
  call setPointersOffTSInstance(nn,sps,sps)

  gm1 = gammaConstant - one
  correctForK = .False.
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

  !call computeEtot(0,ib,0,jb,0,kb,correctForK)
  
  call computeLamViscosity
  call computeEddyViscosity

  !  Apply all BC's
  call applyAllBC_block(.True.)
  
  ! Compute skin_friction Velocity
  call computeUtau_block
  
  ! Compute time step and spectral radius
  call timeStep_block(.false.)
  
  if( equations == RANSEquations ) then
     ! Initialize only the Turblent Variables
     call initres_block(nt1MG, nMGVar,nn,sps) 
     call turbResidual_block
  endif
  
  select case (equationMode)
  case (steady)
     ! Zero out just the flow variables
     dw(:,:,:,1:nwf) = 0.0
  case(timeSpectral)
     do sps2=1,nTimeIntervalsSpectral
        call setPointersOffTSInstance(nn,sps2,sps2)
        dw = 0.0
        do mm=1,nTimeIntervalsSpectral
           call setPointersoffTSInstance(nn,sps2,mm)
           call initRes_block_TS(1,nwf,nn,sps2,mm)
        end do
     end do
  end select
  
  ! Rest the pointers the the "on time instance"
  call setPointersOffTSInstance(nn,sps,sps)  
  
  ! Actual residual calc
  call residual_block

  ! Divide through by the volume
  do sps2 = 1,nTimeIntervalsSpectral
     call setPointersOffTSInstance(nn,sps2,sps2)
     do l=1,nwf
        do k=2,kl
           do j=2,jl
              do i=2,il
                 dw(i,j,k,l) = dw(i,j,k,l) / vol(i,j,k)
              end do
           end do
        end do
     end do
   
     ! Possibly scale turblent variables by 1e-3
     do l=nt1,nt2
        do k=2,kl
           do j=2,jl
              do i=2,il
                 dw(i,j,k,l) = dw(i,j,k,l) / vol(i,j,k)! * 1e-3
              end do
           end do
        end do
     end do
  end do

  call setPointersOffTSInstance(nn,sps,sps)
end subroutine block_res
