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

subroutine block_res_spatial(nn,sps)!,x_peturb)

  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use flowVarRefState     ! nw
  use inputPhysics 
  use iteration
  use inputTimeSpectral
  use section             !nsections
  use monitor             !timeunsteadyrestart

  implicit none


  real(kind=realType) :: gm1,v2
  integer(kind=intType) :: nn,sps,i,j,k,sps2,mm,l
  logical :: correctForK,useOldCoor=.false.
  !for grid velocities computation
  real(kind=realType), dimension(nSections) :: t

!  logical :: x_peturb(0:ie,0:je,0:ke,3)
  call setPointersOffTSInstance(nn,sps,sps)

  ! Do the spatial things first:
  call xhalo_block(1)!,x_peturb)
  call metric_block(nn,1,sps)

  ! Compute the time, which corresponds to this spectral solution.
  ! For steady and unsteady mode this is simply the restart time;
  ! for the spectral mode the periodic time must be taken into
  ! account, which can be different for every section.
  
  t = timeUnsteadyRestart
  
  if(equationMode == timeSpectral) then
     do nn=1,nSections
        t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
             /         real(nTimeIntervalsSpectral,realType)
     enddo
  endif
  call gridVelocitiesFineLevel_block(useOldCoor, t, sps) ! Required for TS
  call normalVelocities_block(sps) ! Required for TS
  !call slipVelocitiesFineLevel(.false., t, mm) !required for viscous

  ! Compute the pressures

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
  
  !  Apply all BC's
  call applyAllBC_block(.True.)
  
  ! Compute skin_friction Velocity
  call computeUtau_block
  
  ! Compute time step and spectral radius
  call timeStep_block(.false.)
  
  !   if( equations == RANSEquations ) then
  !      call initres_block(nt1MG, nMGVar,nn,sps) ! Initialize only the Turblent Variables
  !      call turbResidual_block
  !   endif
  
  select case (equationMode)
  case (steady)
     dw = 0.0
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
     do l=1,nw
        do k=2,kl
           do j=2,jl
              do i=2,il
                 dw(i,j,k,l) = dw(i,j,k,l) / vol(i,j,k)
              end do
           end do
        end do
     end do
  end do


end subroutine block_res_spatial
