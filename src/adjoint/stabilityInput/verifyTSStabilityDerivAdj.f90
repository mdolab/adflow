!
!     ******************************************************************
!     *                                                                *
!     * File:          verifyTSStabilityDerivAdj.f90                   *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 11-26-2009                                      *
!     * Last modified: 11-26-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifyTSStabilityDerivAdj(level)
!
!     ******************************************************************
!     *                                                                *
!     * Computes the Time spectral stability derivatives for the       *
!     * current configuration for the finest grid level over all time  *
!     * instances using the                                            *
!     * auxiliary routines modified for tapenade and compares them to  *
!     * the stability derivatives computed with the original code.     *
!     *                                                                *
!     * This is only executed in debug mode.                           *
!     *                                                                *
!     ******************************************************************
!
      use blockPointers
      use communication       ! myID
      use inputPhysics        !
      use flowVarRefState     !
      use inputDiscretization ! spaceDiscr
      use iteration           ! currentLevel
      use monitor             ! monLoc, MonGlob, nMonSum
      use inputTimeSpectral   ! nTimeInstancesMax
      use section             !nsections, sections%
      use bcTypes             ! EulerWall, NSWallAdiabatic, NSWallIsothermal
      use monitor             ! timeunsteady...
      use costFunctions
      
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level
!
!     Local variables.
!
      integer(kind=intType) :: discr, nHalo
      integer(kind=intType) :: mm, nn, sps,nnn

      logical :: fineGrid, correctForK, exchangeTurb

      integer(kind=intType)::liftIndex

      real(kind=realType), dimension(:,:,:,:), allocatable :: xAdj,wAdj
      real(kind=realType), dimension(:,:,:), allocatable :: pAdj
!      real(kind=realType), dimension(:,:,:,:), allocatable :: siAdj, sjAdj, skAdj

      real(kind=realType), dimension(3) ::rotRateAdj,rotCenterAdj
      
      integer(kind=intType) :: i2Beg, i2End, j2Beg, j2End
      integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd

      real(kind=realType), dimension(3) :: cFpAdj, cFvAdj
      real(kind=realType), dimension(3) :: cMpAdj, cMvAdj

      real(kind=realType), dimension(3) :: cFp, cFv
      real(kind=realType), dimension(3) :: cMp, cMv

      real(kind=realType), dimension(3) :: cFpFD, cFvFD
      real(kind=realType), dimension(3) :: cMpFD, cMvFD

      real(kind=realType), dimension(nTimeIntervalsSpectral) :: Cl,Cd,Cfx,Cfy,Cfz,Cmx,Cmy,Cmz
      real(kind=realType), dimension(nTimeIntervalsSpectral) :: ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj, &
                                                           CmxAdj,CmyAdj,CmzAdj
!      real(kind=realType), dimension(:,:,:),allocatable:: normAdj
      real(kind=realType), dimension(3) :: refPoint
      real(kind=realType) :: yplusMax

      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,machCoefAdj,machGridAdj
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj

      integer(kind=intType) :: nmonsum2 ,i
      real(kind=realType),  dimension(:), allocatable :: monLoc1, monGlob1
      real(kind=realType),  dimension(:), allocatable :: monLoc2, monGlob2

      real(kind=realType) :: fact!temporary

      logical :: contributeToForce, viscousSubface,secondhalo

      integer :: ierr
      
      real(kind=realType),dimension(nTimeIntervalsSpectral,8)::BaseCoef
      real(kind=realType),dimension(8)::dcdp,dcdpdot,dcdq,dcdqdot,dcdr,dcdrdot
      real(kind=realType),dimension(8)::dcdalpha,dcdalphadot,dcdbeta,dcdbetadot,dcdMach,dcdMachdot
      real(kind=realType),dimension(8)::Coef0,Coef0dot
      real(kind=realType), dimension(nCostFunction)::globalCFVals

!!$      real(kind=realType)::dcldp,dcldpdot,dcmzdp,dcmzdpdot         
!!$      real(kind=realType)::dcldq,dcldqdot,dcmzdq,dcmzdqdot
!!$      real(kind=realType)::dcldr,dcldrdot,dcmzdr,dcmzdrdot
!!$      real(kind=realType)::dcldalpha,dcldalphadot,dcmzdalpha,dcmzdalphadot
!!$      real(kind=realType)::dcldMach,dcldMachdot,dcmzdMach,dcmzdMachdot
!!$      real(kind=realType)::cl0,cl0dot,cmz0,cmz0dot

!!$      real(kind=realType)::cl0Adj,cmz0Adj,dcldalphaAdj,dcmzdalphaAdj

      real(kind=realType), dimension(nSections) :: t
!for debug
real(kind=realType), dimension(3) :: cfpadjout, cmpadjout
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      if(myID == 0) then
        write(*,*) "Running verifyTSStabilityDerivAdj..."
        !write(*,10) "CL","CD","Cfx","Cmx"
      endif

      ! Set the grid level of the current MG cycle, the value of the
      ! discretization and the logical fineGrid.

      currentLevel = level
      discr        = spaceDiscr
      fineGrid     = .true.

      ! Determine whether or not the total energy must be corrected
      ! for the presence of the turbulent kinetic energy and whether
      ! or not turbulence variables should be exchanged.

      correctForK  = .false.
      exchangeTurb = .false.
      secondhalo = .true.
!
!     ******************************************************************
!     *                                                                *
!     * Exchange halo data to make sure it is up-to-date.              *
!     * (originally called inside "rungeKuttaSmoother" subroutine).    *
!     *                                                                *
!     ******************************************************************

      ! allocate monLoc2, monGlob2

      nmonsum2 = 8
 !     print *,'allocating monsum'
      allocate(monLoc1(nmonsum2), monGlob1(nmonsum2))
      allocate(monLoc2(nmonsum2), monGlob2(nmonsum2))

!      print *,'exchanging halo data'
      ! Exchange the pressure if the pressure must be exchanged early.
      ! Only the first halo's are needed, thus whalo1 is called.
      ! Only on the fine grid.
      
      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
           .false., .false.)
      
      ! Apply all boundary conditions to all blocks on this level.
      
      call applyAllBC(secondHalo)
      
      ! Exchange the solution. Either whalo1 or whalo2
      ! must be called.
      
      if( secondHalo ) then
         call whalo2(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      else
         call whalo1(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      endif

      call mpi_barrier(SUmb_comm_world, ierr)      
!
!     ******************************************************************
!     *                                                                *
!     * Compute the stability derivatives using the original routine.  *
!     *                                                                *
!     ******************************************************************
!
      
      call computeTSDerivatives(coef0,dcdalpha,dcdalphadot,dcdq,dcdqdot)
      
      call mpi_barrier(SUmb_comm_world, ierr)

      do sps =1,nTimeIntervalsSpectral
         
         
         call computeAeroCoef(globalCFVals,sps)
         
         BaseCoef(sps,1) = globalCFVals(costFuncLiftCoef)
         BaseCoef(sps,2) = globalCFVals(costFuncDragCoef)
         BaseCoef(sps,3) = globalCFVals(costFuncForceXCoef)
         BaseCoef(sps,4) = globalCFVals(costFuncForceYCoef)
         BaseCoef(sps,5) = globalCFVals(costFuncForceZCoef)
         BaseCoef(sps,6) = globalCFVals(costFuncMomXCoef)
         BaseCoef(sps,7) = globalCFVals(costFuncMomYCoef)
         BaseCoef(sps,8) = globalCFVals(costFuncMomZCoef)
         
         
      end do

        
      !Now compute the stability derivatives

      call computeTSStabilityDerivAdj(BaseCoef,coef0,dcdalpha,&
           dcdalphadot,dcdq,dcdqdot)


      ! Root processor outputs results.
 !     print *,'printing results'
!!$      if(myID == 0) then
!!$         write(*,*) 'Stability Derivative verification results'
!!$         write(*,*) 'Cl0   cmz0   dcldalpha   dcmzdalpha'
!!$         write(*,20) "Original", cl0,cmz0,dcldalpha,dcmzdalpha
!!$         write(*,20) "Adjoint ", cl0Adj,cmz0Adj,dcldalphaAdj,dcmzdalphaAdj
!!$      endif
             
      !print *,'finished computing derivatives'
      ! Flush the output buffer and synchronize the processors.
       
      call f77flush()

       
      ! Output format.
       
10    format(1x,4a14)
20    format(1x,a,8(1x,e13.6))

    end subroutine verifyTSStabilityDerivAdj
     
