subroutine computeResidualNK()

  ! This is the residual evaluation driver for the NK solver. The
  ! actual function that is passed to petsc is FormFunction (see
  ! formFunction.F90). This the routine that actually computes the
  ! residual. This works with Euler,Laminar and NS equation
  ! modes. This function ONLY OPERATES ON THE FINEST GRID LEVEL. It
  ! does not coarser grid levels.

  ! This function uses the w that is currently stored in the flowDoms
  ! datastructure and leaves the resulting residual dw, in the same
  ! structure. setW() and setRVec() is used in formFunction to
  ! set/extract these values for communication with PETSc. 

  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use iteration
  use inputPhysics 
  implicit none

  ! Local Variables
  integer(kind=intType) :: ierr,i,j,k,l,sps,nn
  logical secondHalo ,correctForK
  real(kind=realType) :: gm1,v2,val

  secondHalo = .True. 

  currentLevel = 1_intType
  groundLevel = 1_intTYpe
  ! Next we need to compute the pressures
  gm1 = gammaConstant - one
  correctForK = .False.
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domainsState: do nn=1,nDom
        ! Set the pointers to this block.
        call setPointers(nn, currentLevel, sps)

        do k=2,kl
           do j=2,jl
              do i=2,il

                 v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                      + w(i,j,k,ivz)**2

                 p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                      - half*w(i,j,k,irho)*v2)
                 p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
              enddo
           enddo
        enddo

        call computeEtot(2_intType,il, 2_intType,jl, &
             2_intType,kl, correctForK)
        
     end do domainsState
  end do spectralLoop

  
  call computeLamViscosity
  call computeEddyViscosity

  !   Apply BCs
  call applyAllBC(secondHalo)

  ! Exchange solution -- always the fine level
  call whalo2(1_intType, 1_intType, nMGVar, .true., &
       .true., .true.)
  if (equations == RANSEquations) then
     call whalo2(1_intType, nt1, nt2, .false., .false., .true.)  
  end if

  ! Why does this need to be set?
  rkStage = 0
  
  ! Compute the skin-friction velocity
  call computeUtau

  ! Compute time step
  call timestep(.false.)

  ! Possible Turblent Equations
  if( equations == RANSEquations ) then
     call initres(nt1MG, nMGVar) ! Initialize only the Turblent Variables
     call turbResidual
  endif
  
  ! Initialize Flow residuals
  call initres(1_intType, nwf)

  ! Actual Residual Calc
  call residual 

end subroutine computeResidualNK
