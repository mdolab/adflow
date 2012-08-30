subroutine computeResidualNK()

  ! This is the residual evaluation driver for the NK solver. The
  ! actual function that is used for the matrix free matrix-vector
  ! products is FormFunction_mf (see formFunction.F90). This the
  ! routine that actually computes the residual. This works with
  ! Euler,Laminar and NS equation modes. This function ONLY OPERATES
  ! ON THE FINEST GRID LEVEL. It does not coarser grid levels.

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
  real(kind=realType) :: gm1,v2

  secondHalo = .True. 
  currentLevel = 1_intType
  groundLevel = 1_intTYpe
  
  ! Next we need to compute the pressures
  gm1 = gammaConstant - one
  correctForK = .False.

  ! Why does this need to be set?
  rkStage = 0

  ! Recompute pressure on ALL cells 
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domainsState: do nn=1,nDom
        ! Set the pointers to this block.
        call setPointers(nn, currentLevel, sps)

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

        ! Compute Viscosity
        call computeLamViscosity  
        call computeEddyViscosity 

     end do domainsState
  end do spectralLoop

  !   Apply BCs
  call applyAllBC(secondHalo)

  call whalo2(1_intType, 1_intType, nw, .true., &
       .true., .true.)
  
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
subroutine computeResidualNK2()

  ! This is the residual evaluation driver for the NK solver. The
  ! actual function that is used for the matrix free matrix-vector
  ! products is FormFunction_mf (see formFunction.F90). This the
  ! routine that actually computes the residual. This works with
  ! Euler,Laminar and NS equation modes. This function ONLY OPERATES
  ! ON THE FINEST GRID LEVEL. It does not coarser grid levels.

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
  real(kind=realType) :: gm1,v2

  secondHalo = .True. 
  currentLevel = 1_intType
  groundLevel = 1_intTYpe
  
  ! Next we need to compute the pressures
  gm1 = gammaConstant - one
  correctForK = .False.

  ! Why does this need to be set?
  rkStage = 0
 

  ! Recompute pressure on ALL cells 
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domainsState: do nn=1,nDom
        ! Set the pointers to this block.
        call setPointers(nn, currentLevel, sps)

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

        ! Compute Viscosity
        call computeLamViscosity  
        call computeEddyViscosity 

     end do domainsState
  end do spectralLoop

  !Apply BCs
  call applyAllBC(secondHalo)

  call whalo2(1_intType, 1_intType, nMGVar, .true., &
       .true., .true.)
  

 
!   spectralLoop2: do sps=1,nTimeIntervalsSpectral
!      domains2: do nn=1,nDom
!         call setPointers(nn, 1, sps)
!         call bcSymm(secondHalo)
!         call bcNSWallAdiabatic2( secondHalo, correctForK)
!         call bcFarfield(secondHalo, correctForK)
!         call bcEulerWall(secondHalo, correctForK)
!      enddo domains2
!   enddo spectralLoop2
  
  ! Compute time step
  call timestep(.false.)
  
  !Compute the skin-friction velocity
  call computeUtau

  ! Possible Turblent Equations
  if( equations == RANSEquations ) then
     call initres(nt1MG, nMGVar) ! Initialize only the Turblent Variables
     call turbResidual
  endif
  
  ! Initialize Flow residuals
  call initres(1_intType, nwf)

  ! Actual Residual Calc
  call residual 


end subroutine computeResidualNK2
