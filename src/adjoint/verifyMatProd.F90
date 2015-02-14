subroutine verifyMatProd

  !     ******************************************************************
  !     *                                                                *
  !     * This is a verify routine that test the reverse mat prod        *
  !     * with the forward mode.                                         *
  !     *                                                                *
  !     ******************************************************************
  !
  use BCTypes
  use blockPointers_b
  use inputDiscretization 
  use inputTimeSpectral 
  use inputPhysics
  use iteration         
  use flowVarRefState     
  use inputAdjoint       
  use stencils
  use diffSizes
  use ADjointVars
  use inputDiscretization
  use cgnsGrid
  use inputMotion   
  implicit none

  ! Input Variables
  !logical, intent(in) :: firstRun, verifyState, verifySpatial, verifyExtra
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
  !Vec,  allocatable, dimension(:) :: vec1, vec2
  real(kind=realType),dimension(:),allocatable :: vec1, vec2 

  ! Local variables.
  integer(kind=intType) :: i, j, k, l, nn, ii, ierr, jj
  integer(kind=intType) :: nState, level, idxblk
   
  real(kind=realType) :: alpha, beta, force(3), moment(3), sepSensor, Cavitation
  real(kind=realType) :: alphab, betab, forceb(3), momentb(3), sepSensorb, Cavitationb
  real(kind=realType) :: fwdValue, revValue, ran
  real(kind=realType) :: time1, timeb, timed, time, xvbarsum1, xvbarsum2

  integer(kind=intType) :: liftIndex
  logical :: resetToRANS

#ifndef USE_COMPLEX

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif

  ! This routine will not use the extra variables to block_res or the
  ! extra outputs, so we must zero them here
   call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  
  ! Need to trick the residual evalution to use coupled (mean flow and
  ! turbulent) together.

  ! If we want to do the matrix on a coarser level, we must first
  ! restrict the fine grid solutions, since it is possible the
  ! NKsolver was used an the coarse grid solutions are (very!) out of
  ! date. 
  
  ! Assembling matrix on coarser levels is not entirely implemented yet. 
  level = 1
  currentLevel = level
  groundLevel = level

  ! Exchange data and call the residual to make sure its up to date
  ! withe current w
  call whalo2(1_intType, 1_intType, nw, .True., .True., .True.)
  call computeResidualNK ! This is the easiest way to do this

  ! If we are computing the jacobian for the RANS equations, we need
  ! to make block_res think that we are evauluating the residual in a
  ! fully coupled sense.  This is reset after this routine is
  ! finished.
  if (equations == RANSEquations) then
     nMGVar = nw
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if


  !Check state
  logicCheck1: if ( verifyState ) then
     nn = 1
     ! Set pointers to the first timeInstance...just to getSizes
     call setPointers(nn, level, 1)
     call setDiffSizes
     
     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     !call alloc_derivative_values( level)
     call alloc_derivative_values_bwd(level)
        
     ! Set pointers and derivative pointers
     call setPointers_b(nn, level, 1)
     !call setPointers_d(nn, level, 1)
     ! Reset All States and possibe AD seeds
     flowdomsb(1,1,1)%dw = zero 
     print *, ncellslocal(1)
     allocate(vec1(ncellslocal(1)*nState),vec2(ncellslocal(1)*nState))
     
     flowdomsb(1,1,1)%w = zero
     ii = 0
     do k=2, kl
        do j=2,jl
           do i=2,il
              do l = 1,nstate
                 call random_seed
                 call random_number(ran)
                 ii = ii + 1
                 vec1(ii) = ran
                 vec2(ii) = ran
                 flowdomsb(1,1,1)%dw(i, j, k, l) = ran
              end do
           end do
        end do
     end do
     call getdRdwTVec(vec1, vec2, ncellslocal(1)*nState)
     
     call BLOCK_RES_B(nn, 1, .False., alpha, alphab, beta, betab, &
          & liftindex, force, forceb, moment, momentb, sepsensor, sepsensorb, &
          & cavitation, cavitationb)
     

     ii = 0
     do k=2, kl
        do j=2,jl
           do i=2,il
              do l = 1,nstate
                 ii = ii + 1
                 if (abs(flowdomsb(1,1,1)%w(i,j,k,l) - vec2(ii))/abs(vec2(ii)) > 1e-3) then
                    print *,i,j,k,l,flowdomsb(1,1,1)%w(i, j, k, l)-vec2(ii)
                    print *, flowdomsb(1,1,1)%w(i, j, k, l), vec2(ii)
                 end if
              end do
           end do
        end do
     end do
     print *, 'state done'
     deallocate(vec1, vec2)
     !call dealloc_derivative_values(nn, level)
     call dealloc_derivative_values_bwd(level)

  end if logicCheck1

  !Check spatial
  logicCheck2: if ( verifySpatial ) then
     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     !call alloc_derivative_values(nn, level)
     call alloc_derivative_values_bwd(level)


     allocate(vec1(ncellslocal(1)*nState),vec2(nNodesLocal(1)*3))
     call random_seed

     xvbarsum1 = zero
     xvbarsum2 = zero
     ii = 0
     do nn=1,nDom
        call setPointers(nn, 1, 1)
        do k=2, kl
           do j=2,jl
              do i=2,il
                 do l = 1,5                
                    call random_number(ran)
                    ran = one
                    ii = ii + 1
                    vec1(ii) = ran
                 end do
              end do
           end do
        end do
     end do
     
     deallocate(vec1, vec2)
     !call dealloc_derivative_values(nn, level)
     call dealloc_derivative_values_bwd(level)
     print *,' xvbarsums:', xvbarsum1, xvbarsum2
  end if logicCheck2

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  ! Reset the paraters to use segrated turbulence solve. 
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nwf + 1
     nt2MG = nwf

     turbSegregated = .True.
     turbCoupled = .False.
     restrictEddyVis = .false.
     if( eddyModel ) restrictEddyVis = .true.
  end if

#endif
end subroutine verifyMatProd
