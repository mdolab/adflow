! ############### THIS ROUTINE IS OUT OF DATE! ###############

subroutine verifyAD
  !     ******************************************************************
  !     *                                                                *
  !     * This is a verify routine that runs compiled TGT debug mode.    *
  !     * It loop through all the inputs and check the output of AD      *
  !     * against finite difference:                                     *
  !     * firstRun: True, this is the first run that get piped into the  *
  !     *        second run.                                             *
  !     * usage: python runscript1.py | python runscript2.py             *
 
  !     ******************************************************************
  !
  use BCTypes
  use blockPointers_d      
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

  ! ! Input Variables
  ! !logical, intent(in) :: firstRun, verifyState, verifySpatial, verifyExtra

  ! ! Local variables.
  ! integer(kind=intType) :: i, j, k, l, nn
  ! integer(kind=intType) :: nState, level, idxblk
   
  ! real(kind=realType) :: alpha, beta, force(3), moment(3), sepSensor
  ! real(kind=realType) :: alphad, betad, forced(3), momentd(3), sepSensord

  ! integer(kind=intType) :: liftIndex
  ! logical :: resetToRANS

  print *,'verifyAD: THIS ROUTINE IS OUT OF DATE AND MUST BE REWRITTEN '
  stop


! #ifndef USE_COMPLEX

!   ! Setup number of state variable based on turbulence assumption
!   if ( frozenTurbulence ) then
!      nState = nwf
!   else
!      nState = nw
!   endif

!   ! This routine will not use the extra variables to block_res or the
!   ! extra outputs, so we must zero them here
!   alphad = zero
!   betad  = zero
!   machd  = zero
!   machGridd = zero
!   lengthRefd = zero
!   pointRefd  = zero
!   surfaceRefd = zero
!   call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  
!   ! Need to trick the residual evalution to use coupled (mean flow and
!   ! turbulent) together.

!   ! If we want to do the matrix on a coarser level, we must first
!   ! restrict the fine grid solutions, since it is possible the
!   ! NKsolver was used an the coarse grid solutions are (very!) out of
!   ! date. 
  
!   ! Assembling matrix on coarser levels is not entirely implemented yet. 
!   level = 1
!   currentLevel = level
!   groundLevel = level

!   ! Exchange data and call the residual to make sure its up to date
!   ! withe current w
!   call whalo2(1_intType, 1_intType, nw, .True., .True., .True.)
!   call computeResidualNK ! This is the easiest way to do this

!   ! If we are computing the jacobian for the RANS equations, we need
!   ! to make block_res think that we are evauluating the residual in a
!   ! fully coupled sense.  This is reset after this routine is
!   ! finished.
!   if (equations == RANSEquations) then
!      nMGVar = nw
!      nt1MG = nt1
!      nt2MG = nt2

!      turbSegregated = .False.
!      turbCoupled = .True.
!   end if

!   ! Determine if we want to use frozenTurbulent Adjoint
!   resetToRANS = .False. 
!   if (frozenTurbulence .and. equations == RANSEquations) then
!      equations = NSEquations 
!      resetToRANS = .True.
!   end if

!   !Check state
!   logicCheck1: if ( verifyState ) then
!      print *, "Verifying state ..."
!      ! Master Domain Loop
!      domainLoopAD1: do nn=1, nDom
        
!         ! Set pointers to the first timeInstance...just to getSizes
!         call setPointers(nn, level, 1)
!         call setDiffSizes
     
!         ! Allocate the memory we need for this block to do the forward
!         ! mode derivatives and copy reference values
!         call alloc_derivative_values(nn, level)
        
!         ! Set pointers and derivative pointers
!         call setPointers_d(nn, level, 1)
                  
!         ! verify block_res
!         do k=0, kb
!            do j=0, jb
!               do i=0, ib
!                  ! Master State Loop            
!                  stateLoop: do l=1, nState
                    
!                     flowDomsd(nn, 1, 1)%dw_deriv(:, :, :, :, :) = zero
!                     ! Reset All States and possibe AD seeds
!                     flowDoms(nn, level, 1)%w(:, :, :, :) =  flowDomsd(nn, 1, 1)%wtmp
!                     flowdomsd(nn, 1, 1)%w = zero ! This is actually w seed
                    
!                     ! Set the tagent direction
!                     flowDomsd(nn, 1, 1)%w(i, j, k, l) = 1.d0
                    
!                     ! initialize TGT debug
!                     if ( firstRun ) then
!                        call DEBUG_TGT_INIT1(1.d-6,1.d-18,.1d0)
!                     else
!                        call DEBUG_TGT_INIT2(1.d-6,1.d-18,.1d0)
!                     end if
                    
!                     ! let debugger know which one are purturbed
!                     call DEBUG_TGT_INITREAL8(flowDoms(nn, 1, 1)%w(i, j, k, l), &
!                     flowDomsd(nn, 1, 1)%w(i, j, k, l))
                    
!                     ! call block_res_t
!                     call DEBUG_TGT_CALL('block_res',.true.,.false.)
!                     call block_res_t(nn, 1, .True., &
!                          alpha, alphad, beta, betad, liftIndex, force, forced, moment, momentd, &
!                          sepSensor, sepSensord)
                    
!                     ! conclude debugger
!                     call DEBUG_TGT_EXIT()
!                     call DEBUG_TGT_CONCLUDEREAL8ARRAY('force', force, forced)
!                     call DEBUG_TGT_CONCLUDEREAL8ARRAY('moment', moment, momentd)
                    
!                  end do stateLoop
!               end do
!            end do
!         end do

!         call dealloc_derivative_values(nn, level)
!      end do domainLoopAD1
!   end if logicCheck1
      



!   !Check spatial
!   logicCheck2: if ( verifySpatial ) then
!      print *, "Verifying spatial ..."
!      ! Master Domain Loop
!      domainLoopAD2: do nn=1, nDom
        
!         ! Set pointers to the first timeInstance...just to getSizes
!         call setPointers(nn, level, 1)
!         call setDiffSizes
     
!         ! Allocate the memory we need for this block to do the forward
!         ! mode derivatives and copy reference values
!         call alloc_derivative_values(nn, level)
        
!         ! Set pointers and derivative pointers
!         call setPointers_d(nn, level, 1)
                    
!         ! verify block_res
!         do k=0, ke
!            do j=0, je
!               do i=0, ie
                 
!                  ! Master State Loop
!                  dofLoop: do l=1, 3
                    
!                     flowDomsd(nn, 1, 1)%dw_deriv(:, :, :, :, :) = zero
                    
!                     ! Reset All States and possibe AD seeds
!                     flowDoms(nn, level, 1)%x(:, :, :, :) =  flowDomsd(nn, 1, 1)%xtmp
!                     flowdomsd(nn, 1, 1)%x = zero ! This is actually w seed
                    
!                     ! Set the tagent direction
!                     flowDomsd(nn, 1, 1)%x(i, j, k, l) = 1.d0
                    
!                     ! initialize TGT debug
!                     if ( firstRun ) then
!                        call DEBUG_TGT_INIT1(1.d-6,1.d-18,.1d0)
!                     else
!                        call DEBUG_TGT_INIT2(1.d-6,1.d-18,.1d0)
!                     end if
                    
!                     ! let debugger know which one are purturbed
!                     call DEBUG_TGT_INITREAL8(flowDoms(nn, 1, 1)%x(i, j, k, l), &
!                     flowDomsd(nn, 1, 1)%x(i, j, k, l))
                    
!                     ! call block_res_t
!                     call DEBUG_TGT_CALL('block_res',.true.,.false.)
!                     call block_res_t(nn, 1, .True., &
!                          alpha, alphad, beta, betad, liftIndex, force, forced, moment, momentd,&
!                          sepSensor, sepSensord)

!                     ! conclude debugger
!                     call DEBUG_TGT_EXIT()
!                     call DEBUG_TGT_CONCLUDEREAL8ARRAY('force', force, forced)
!                     call DEBUG_TGT_CONCLUDEREAL8ARRAY('moment', moment, momentd)
               
!                  end do dofLoop
!               end do
!            end do
!         end do
        
!         call dealloc_derivative_values(nn, level)
!      end do domainLoopAD2
!   end if logicCheck2
  
  





!   !Check extra
!   logicCheck3: if ( verifyExtra ) then
!      print *, "Verifying Extra ..."
!      ! Master Domain Loop
!      domainLoopAD3: do nn=1, nDom
        
!         ! Set pointers to the first timeInstance...just to getSizes
!         call setPointers(nn, level, 1)
!         call setDiffSizes
     
!         ! Allocate the memory we need for this block to do the forward
!         ! mode derivatives and copy reference values
!         call alloc_derivative_values(nn, level)
        
!         ! Set pointers and derivative pointers
!         call setPointers_d(nn, level, 1)
!         idxblk = nbkGlobal

!         ! Set the seed
!         alphad = 0.0
!         betad = 0.0
!         machd = 0.0
!         machGridd = 0.0
!         machCoefd = 0.0
!         cgnsDomsd(idxblk)%rotrate(:) = 0.0
!         cgnsDomsd(idxblk)%rotcenter(:) = 0.0
!         rotpointd(:) = 0.0
!         pointrefd(:) = 0.0

!         ! zero out
!         flowDomsd(nn, 1, 1)%dw_deriv(:, :, :, :, :) = zero
        
!         ! Reset All x and possibe AD seeds
!         alphad = 1.d0

!         ! verify block_res
!         ! initialize TGT debug
!         if ( firstRun ) then
!            call DEBUG_TGT_INIT1(1.d-6,1.d-18,.1d0)
!         else
!            call DEBUG_TGT_INIT2(1.d-6,1.d-18,.1d0)
!         end if
                    
!         ! let debugger know which one are purturbed
!         call DEBUG_TGT_INITREAL8(alpha, alphad)
                    
!         ! call block_res_t
!         call DEBUG_TGT_CALL('block_res',.true.,.false.)
!         call block_res_t(nn, 1, .True., &
!              alpha, alphad, beta, betad, liftIndex, &
!              force, forced, moment, momentd, sepSensor, sepSensord)
                    
!         ! conclude debugger
!         call DEBUG_TGT_EXIT()
!         call DEBUG_TGT_CONCLUDEREAL8ARRAY('moment', moment, momentd)
        
!      end do domainLoopAD3
     
!      call dealloc_derivative_values(nn, level)
!   end if logicCheck3



!   ! Reset the correct equation parameters if we were useing the frozen
!   ! Turbulent 
!   if (resetToRANS) then
!      equations = RANSEquations
!   end if

!   ! Reset the paraters to use segrated turbulence solve. 
!   if (equations == RANSEquations) then
!      nMGVar = nwf
!      nt1MG = nwf + 1
!      nt2MG = nwf

!      turbSegregated = .True.
!      turbCoupled = .False.
!      restrictEddyVis = .false.
!      if( eddyModel ) restrictEddyVis = .true.
!   end if

! #endif
end subroutine verifyAD  
