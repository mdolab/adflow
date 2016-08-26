subroutine computeResidualANK()

  ! This is the residual evaluation driver for the ANK solver. It
  ! computes the residual for the mean flow but does not compute the
  ! turbulent residuals. 
  use constants
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use iteration
  use inputPhysics 
  use utils, only : setPointers
  use haloExchange, only : whalo2
  implicit none

  ! Local Variables
  integer(kind=intType) :: i, j, k, sps,nn
  logical secondHalo, correctForK
  real(kind=realType) :: gm1, factK, v2

  gm1 = gammaConstant - one
  rkStage = 0

  secondHalo = .false.
  correctForK = .false.
  if(currentLevel <= groundLevel) then
     secondHalo = .true.
     if (kPresent) then 
        correctForK = .True.
     end if
  end if

  ! Recompute pressure on ALL cells 
  spectralLoop: do sps=1, nTimeIntervalsSpectral
     domainsState: do nn=1, nDom
        ! Set the pointers to this block.
        call setPointers(nn, currentLevel, sps)
        factK = zero
        do k=0, kb
           do j=0, jb
              do i=0, ib
                 
                 gm1  = gamma(i, j, k) - one
                 factK = five*third - gamma(i, j ,k)
                 v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                      + w(i,j,k,ivz)**2

                 p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                      - half*w(i,j,k,irho)*v2) 

                 if( correctForK ) then 
                    p(i, j ,K) = p(i,j, k) + factK*w(i, j, k, irho) &
                         * w(i, j, k, itu1)
                 end if

                 ! Clip to make sure it is positive.
                 p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
              end do
           end do
        end do

        ! Compute Viscosities
        call computeLamViscosity  
        call computeEddyViscosity 
     end do domainsState
  end do spectralLoop

  ! Apply BCs
  call applyAllBC(secondHalo)

  ! Exchange halos
  call whalo2(currentLevel, 1_intType, nwf, .true., &
       .true., .true.)

  ! Compute time step (spectral radius is actually what we need)
  call timestep(.false.)

  ! Initialize Flow residuals
  call initres(1_intType, nwf)

  ! Actual Residual Calc
  call residual 

end subroutine computeResidualANK


! subroutine computeResidualANKTurb()

!   ! This is the residual evaluation driver for the ANK solver. It
!   ! computes the residual for the mean flow but does not compute the
!   ! turbulent residuals. 

!   use blockPointers
!   use inputTimeSpectral
!   use flowvarrefstate
!   use iteration
!   use inputPhysics 
!   implicit none

!   ! Local Variables
!   integer(kind=intType) :: nn, sps


!   call whalo2(currentLevel, nt1, nt2, .False., .False., .False.)
!   spectralLoop: do sps=1, nTimeIntervalsSpectral
!      domainsState: do nn=1, nDom
!         ! Set the pointers to this block.
!         call setPointers(nn, currentLevel, sps)
!         call computeEddyViscosity 
!      end do domainsState
!   end do spectralLoop
!   if (equations == RANSEquations) then 

!      call computeUTau
!      call initRes(nt1, nt2)
!      call turbResidual
!   end if

! end subroutine computeResidualANKTurb


