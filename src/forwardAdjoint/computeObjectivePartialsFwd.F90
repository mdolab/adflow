subroutine getdIdw(costFunction, dIdw, nState)
#ifndef USE_COMPLEX
  ! Return the derivative: 'dIdw' of size 'nState' for the objective
  ! number 'costFunction'.

    use ADjointVars
  use ADjointPETSc
  use blockPointers   
  use communication  
  use costFunctions
  use inputTimeSpectral
  use inputPhysics
  use flowvarrefstate
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: costFunction, nState
  reaL(kind=realType), intent(out) :: dIdw(nState)
  
  ! Working Variables
  integer(kind=intType) :: i, ierr, sps
  real(kind=realtype) :: val
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, forceb
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: moment, momentb
  real(kind=realType), dimension(nTimeIntervalsSpectral) :: sepSensor, sepSensorb
  real(kind=realType), dimension(nTimeIntervalsSpectral) :: Cavitation, Cavitationb
  real(kind=realType) :: alpha, alphab, beta, betab
  real(kind=realType) :: objValue, objValueb
  integer(kind=intType) :: liftIndex, idim

  ! Compute the requiquired sensitivity of the objective with respect
  ! to the forces, moments and extra variables.

  do sps=1, nTimeIntervalsSpectral
     call getSolution(sps)
     force(1, sps) = functionValue(costFuncForceX)
     force(2, sps) = functionValue(costFuncForceY)
     force(3, sps) = functionValue(costFuncForceZ)
     moment(1, sps) = functionValue(costFuncMomX)
     moment(2, sps) = functionValue(costFuncMomY)
     moment(3, sps) = functionValue(costFuncMomZ)
     sepSensor(sps) = functionValue(costFuncSepSensor)
     Cavitation(sps) = functionValue(costFuncCavitation)
  end do

  objValueb = one
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  call getCostFunction_b(costFunction, force, forceb, moment, momentb, sepSensor,&
       sepSensorb, Cavitation, Cavitationb, alpha, alphab, beta, betab, liftIndex, objValue, objValueb)

  ! Set the supplied vector into a psi_like array
  call VecPlaceArray(psi_like1, dIdw, ierr)
  call EChk(ierr, __FILE__, __LINE__)  

  ! Zero Entries and multiply through by reverse-mode derivatives
  call VecZeroEntries(psi_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  do sps=1, nTimeIntervalsSpectral
     do i=1, 3
        call VecAXPY(psi_like1, forceb(i, sps), FMw(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAXPY(psi_like1, momentb(i, sps), FMw(i+3, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  if (costFunction == costFuncSepSensor) then 
     do sps=1,nTimeIntervalsSpectral
        call VecAXPY(psi_like1, sepSensorb(sps), FMw(iSepSensor, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end if
 
  if (costFunction == costFuncCavitation) then 
     do sps=1,nTimeIntervalsSpectral
        call VecAXPY(psi_like1, Cavitationb(sps), FMw(iCavitation, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end if

  ! Assemble dIdw (vector is called psi_like1)
  call VecAssemblyBegin(psi_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(psi_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Reset array
  call VecResetArray(psi_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine getdIdw

subroutine getdIdx(costFunction, dIdx, nSpatial)
#ifndef USE_COMPLEX
  ! Return the derivative: 'dIdx' of size 'nSpatial' for the objective
  ! number 'costFunction'.

  use ADjointVars
  use ADjointPETSc
  use blockPointers   
  use communication  
  use costFunctions
  use inputTimeSpectral
  use inputPhysics
  use flowvarrefstate
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: costFunction, nSpatial
  reaL(kind=realType), intent(out) :: dIdx(nSpatial)
  
  ! Working Variables
  integer(kind=intType) :: i, ierr, sps
  real(kind=realtype) :: val
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, forceb
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: moment, momentb
  real(kind=realType), dimension(nTimeIntervalsSpectral) :: sepSensor, sepSensorb
  real(kind=realType), dimension(nTimeIntervalsSpectral) :: Cavitation, Cavitationb
  real(kind=realType) :: alpha, alphab, beta, betab
  real(kind=realType) :: objValue, objValueb
  integer(kind=intType) :: liftIndex, idim

  ! Compute the requiquired sensitivity of the objective with respect
  ! to the forces, moments and extra variables.

  do sps=1, nTimeIntervalsSpectral
     call getSolution(sps)
     force(1, sps) = functionValue(costFuncForceX)
     force(2, sps) = functionValue(costFuncForceY)
     force(3, sps) = functionValue(costFuncForceZ)
     moment(1, sps) = functionValue(costFuncMomX)
     moment(2, sps) = functionValue(costFuncMomY)
     moment(3, sps) = functionValue(costFuncMomZ)
     sepSensor(sps) = functionValue(costFuncSepSensor)
     Cavitation(sps) = functionValue(costFuncCavitation)
  end do

  objValueb = one
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  call getCostFunction_b(costFunction, force, forceb, moment, momentb, sepSensor,&
       sepSensorb, Cavitation, Cavitationb, alpha, alphab, beta, betab, liftIndex, objValue, objValueb)

  call VecPlaceArray(x_like, dIdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Zero Entries and multiply by reverse-mode derivatives
  call VecZeroEntries(x_like, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do sps=1, nTimeIntervalsSpectral
     do i = 1, 3
        call VecAXPY(x_like, forceb(i, sps), FMx(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call VecAXPY(x_like, momentb(i, sps), FMx(i+3, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  ! Assemble x_like
  call VecAssemblyBegin(x_like, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(x_like, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (costFunction == costFuncSepSensor) then 
     do sps=1,nTimeIntervalsSpectral
        call VecAXPY(x_like, sepSensorb(sps), FMx(iSepSensor, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end if

  if (costFunction == costFuncCavitation) then 
     do sps=1,nTimeIntervalsSpectral
        call VecAXPY(x_like, Cavitationb(sps), FMx(iCavitation, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end if

  call VecResetArray(x_like, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine getdIdx

subroutine getdIda(costFunction)
#ifndef USE_COMPLEX
  ! Return the derivative: 'dIda' for the objective number 'costFunction'
  ! Result is stored in the dIda variable in the ADjointVars Module.

  use ADjointVars
  use ADjointPETSc
  use blockPointers   
  use communication  
  use costFunctions
  use inputTimeSpectral
  use inputPhysics
  use flowvarrefstate
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: costFunction
  
  ! Working Variables
  integer(kind=intType) :: i, ierr, sps
  real(kind=realtype) :: val
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, forceb
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: moment, momentb
  real(kind=realType), dimension(nTimeIntervalsSpectral) :: sepSensor, sepSensorb
  real(kind=realType), dimension(nTimeIntervalsSpectral) :: Cavitation, Cavitationb
  real(kind=realType) :: alpha, alphab, beta, betab
  real(kind=realType) :: objValue, objValueb
  integer(kind=intType) :: liftIndex, idim

  ! Compute the requiquired sensitivity of the objective with respect
  ! to the forces, moments and extra variables.
  
  do sps=1, nTimeIntervalsSpectral
     call getSolution(sps)
     force(1, sps) = functionValue(costFuncForceX)
     force(2, sps) = functionValue(costFuncForceY)
     force(3, sps) = functionValue(costFuncForceZ)
     moment(1, sps) = functionValue(costFuncMomX)
     moment(2, sps) = functionValue(costFuncMomY)
     moment(3, sps) = functionValue(costFuncMomZ)
     sepSensor(sps) = functionValue(costFuncSepSensor)
     Cavitation(sps) = functionValue(costFuncCavitation)
  end do

  objValueb = one
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  call getCostFunction_b(costFunction, force, forceb, moment, momentb, sepSensor,&
       sepSensorb, Cavitation, Cavitationb, alpha, alphab, beta, betab, liftIndex, objValue, objValueb)

  if (nDesignExtra > 0) then
     dIda = zero

     if (myid == 0) then
        if (nDesignAoA >= 0) then 
           dIda(nDesignAoA + 1) = dIda(nDesignAoA+1) + alphab
        end if
        if (nDesignSSA >= 0) then 
           dIda(nDesignSSA + 1) = dIda(nDesignSSA+1) + betab
        end if
        if (nDesignMach >= 0) then 
           dIda(nDesignMach + 1) = dIda(nDesignMach+1) + machcoefb + machb 
        end if
        
        if (nDesignMachGrid >= 0) then 
           dIda(nDesignMachGrid + 1) = dIda(nDesignMachGrid+1) + machgridb + machcoefb
        end if

        if (nDesignPressure >= 0) then 
           dIda(nDesignPressure + 1) = dIda(nDesignPressure+1) + Prefb
        end if

        if (nDesignTemperature >= 0) then 
           dIda(nDesignTemperature + 1) = dIda(nDesignTemperature+1) + tempfreestreamb
        end if

        if (nDesignReynolds >= 0) then 
           dIda(nDesignReynolds + 1) = dIda(nDesignTemperature+1) + reynoldsb
        end if

        if (nDesignLengthRef >= 0) then
           dIda(nDesignLengthRef + 1) = dIda(nDesignLengthRef+1) + lengthrefb
        end if

        ! Explict dependence on pointRef....Only on one proc since
        ! they are the same across all procs, and the sum would nProc
        ! too large if they were on all procs. 
        if (nDesignPointRefX >= 0) then
           dIda(nDesignPointRefX + 1) = dIda(nDesignPointRefX + 1) + pointrefb(1)
        end if

        if (nDesignPointRefY >= 0) then
           dIda(nDesignPointRefY + 1) = dIda(nDesignPointRefY + 1) + pointrefb(2)
        end if

        if (nDesignPointRefZ >= 0) then
           dIda(nDesignPointRefZ + 1) = dIda(nDesignPointRefZ + 1) + pointrefb(3)
        end if

     end if

     ! These three are a little different; The derivative wrt to the
     ! forces and moments for each spectral instance are computed when
     ! the extra residual matrix is computed. This is the only
     ! sensivitiy for pointRef. We also know the objective derivative
     ! wrt forces and moment so we chain-rule them together. Note that
     ! there is no dependence of 'force' on pointRef so it is not
     ! included here. Also these derivatives DO need to be summed over
     ! all procs.
     
     ! add missing dependence of mach
     if (nDesignMach > 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignMach+1) = dIda(nDesignMach+1) + &
                   dFMdExtra(idim, nDesignmach+1, sps)*forceb(idim, sps) 
              dIda(nDesignMach+1) = dIda(nDesignMach+1) + &
                   dFMdExtra(idim+3, nDesignmach+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignPointRefX >= 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignPointRefX+1) = dIda(nDesignPointRefX+1) + &
                   dFMdExtra(3+idim, nDesignPointRefX+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignPointRefY >= 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignPointRefY+1) = dIda(nDesignPointRefY+1) + &
                   dFMdExtra(3+idim, nDesignPointRefY+1, sps)*momentb(idim, sps)
           end do
        end do
     end if
     if (nDesignPointRefZ >= 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignPointRefZ+1) = dIda(nDesignPointRefZ+1) + &
                   dFMdExtra(3+idim, nDesignPointRefZ+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignLengthRef >= 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignLengthRef+1) = dIda(nDesignLengthRef+1) + &
                   dFMdExtra(idim, nDesignLengthRef+1, sps)*forceb(idim, sps)
              dIda(nDesignLengthRef+1) = dIda(nDesignLengthRef+1) + &
                   dFMdExtra(idim+3, nDesignLengthRef+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignPressure >= 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignPressure+1) = dIda(nDesignPressure+1) + &
                   dFMdExtra(idim, nDesignPressure+1, sps)*forceb(idim, sps)
              dIda(nDesignPressure+1) = dIda(nDesignPressure+1) + &
                   dFMdExtra(idim+3, nDesignPressure+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignTemperature >= 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignTemperature+1) = dIda(nDesignTemperature+1) + &
                   dFMdExtra(idim, nDesignTemperature+1, sps)*forceb(idim, sps)
              dIda(nDesignTemperature+1) = dIda(nDesignTemperature+1) + &
                   dFMdExtra(idim+3, nDesignTemperature+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

     if (nDesignReynolds >= 0) then
        do sps=1, nTimeIntervalsSpectral
           do idim=1,3
              dIda(nDesignReynolds+1) = dIda(nDesignReynolds+1) + &
                   dFMdExtra(idim, nDesignReynolds+1, sps)*forceb(idim, sps)
              dIda(nDesignReynolds+1) = dIda(nDesignReynolds+1) + &
                   dFMdExtra(idim+3, nDesignReynolds+1, sps)*momentb(idim, sps)
           end do
        end do
     end if

  end if
#endif
end subroutine getdIda


