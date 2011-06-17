subroutine computeResidualNK2()

  ! Debug copy of computeResidualNK to figure out where the nans are
  ! comming from
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
  integer(kind=intType) :: w_nan,p_nan,dw_nan
  secondHalo = .True. 

  currentLevel = 1_intType
  groundLevel = 1_intTYpe
  ! Next we need to compute the pressures
  gm1 = gammaConstant - one
  correctForK = .False.
  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check1:',w_nan,p_nan,dw_nan

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

  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check2:',w_nan,p_nan,dw_nan

  
  call computeLamViscosity
  call computeEddyViscosity

  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check3:',w_nan,p_nan,dw_nan
  

  !   Apply BCs
  call applyAllBC(secondHalo)

  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check4:',w_nan,p_nan,dw_nan

  ! Exchange solution -- always the fine level
  call whalo2(1_intType, 1_intType, nMGVar, .true., &
       .true., .true.)
  if (equations == RANSEquations) then
     call whalo2(1_intType, nt1, nt2, .false., .false., .true.)  
  end if

  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check5:',w_nan,p_nan,dw_nan


  ! Why does this need to be set?
  rkStage = 0
  
  ! Compute the skin-friction velocity
  call computeUtau

  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check6:',w_nan,p_nan,dw_nan


  ! Compute time step
  call timestep(.false.)

  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check7:',w_nan,p_nan,dw_nan

  ! Possible Turblent Equations
  if( equations == RANSEquations ) then
     call initres(nt1MG, nMGVar) ! Initialize only the Turblent Variables
     call turbResidual
  endif
  
  ! Initialize Flow residuals
  call initres(1_intType, nwf)
  
  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check8:',w_nan,p_nan,dw_nan

  ! Actual Residual Calc
  call residual2

  !call checkwforNan(w_nan)
  !call checkpforNan(p_nan)
  !call checkdwforNan(dw_nan)
  !print *, 'Check9:',w_nan,p_nan,dw_nan
  


end subroutine computeResidualNK2
