subroutine computeResidualNK()

  ! This is the residual evaluation driver for the NK solver. The
  ! actual function that is used for the matrix free matrix-vector
  ! products is FormFunction_mf (see formFunction.F90). This the
  ! routine that actually computes the residual. This works with
  ! Euler, Laminar and RANS equation modes. 

  ! This function uses the w that is currently stored in the flowDoms
  ! datastructure and leaves the resulting residual dw, in the same
  ! structure. setW() and setRVec() is used in formFunction to
  ! set/extract these values for communication with PETSc. 

  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use iteration
  use inputPhysics 
  use saModule
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

  if (equations == RANSequations) then 
     do nn=1,nDom
        do sps=1,nTimeIntervalsSpectral
           call setPointers(nn, currentLevel, sps)
           call bcTurbTreatment
           call applyAllTurbBCThisBLock(.True.)
        end do
     end do
  end if

  ! Exchange halos
  call whalo2(currentLevel, 1_intType, nw, .true., &
       .true., .true.)

  ! Compute time step (spectral radius is actually what we need)
  call timestep(.false.)

  ! Possible Turblent Equations
  if( equations == RANSEquations ) then
     ! Compute the skin-friction velocity (wall functions only)
     call computeUtau
     call initres(nt1, nt2) ! Initialize only the Turblent Variables
     call turbResidual
  endif

  ! Initialize Flow residuals
  call initres(1_intType, nwf)

  ! Actual Residual Calc
  call residual 

end subroutine computeResidualNK




