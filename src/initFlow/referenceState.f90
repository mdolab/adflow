!
!      ******************************************************************
!      *                                                                *
!      * File:          referenceState.f90                              *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 05-29-2003                                      *
!      * Last modified: 04-22-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine referenceState
!
!      ******************************************************************
!      *                                                                *
!      * referenceState computes the reference state values in case     *
!      * these have not been specified. A distinction is made between   *
!      * internal and external flows. In case nothing has been          *
!      * specified for the former a dimensional computation will be     *
!      * made. For the latter the reference state is set to an          *
!      * arbitrary state for an inviscid computation and computed for a *
!      * viscous computation. Furthermore for internal flows an average *
!      * velocity direction is computed from the boundary conditions,   *
!      * which is used for initialization.                              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use communication
       use constants
       use couplerParam
       use flowVarRefState
       use inputMotion
       use inputPhysics
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: sps, nn, mm

       real(kind=realType) :: gm1, ratio, tmp
       real(kind=realType) :: mx, my, mz, Re, v, TInfDim

       real(kind=realType), dimension(3) :: dirLoc, dirGlob
       real(kind=realType), dimension(5) :: valLoc, valGlob

       type(BCDataType), dimension(:), pointer :: BCData
!
!      Interfaces
!
       interface
         subroutine velMagnAndDirectionSubface(vmag, dir, &
                                               BCData, mm)
           use block
           implicit none

           integer(kind=intType), intent(in) :: mm
           real(kind=realType), intent(out) :: vmag
           real(kind=realType), dimension(3), intent(inout) :: dir
           type(BCDataType), dimension(:), pointer :: BCData
         end subroutine velMagnAndDirectionSubface
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the dimensional free stream temperature and pressure.
       ! From these values the density and viscosity is computed. For
       ! external viscous and internal computation this is corrected
       ! later on.

       pInfDim = pRef
       if(pRef <= zero) pInfDim = 101325.0_realType
       TInfDim = tempFreestream

       rhoInfDim = pInfDim/(RGasDim*TInfDim)
       muDim     = muSuthDim                                    &
                 * ((TSuthDim + SSuthDim)/(TInfDim + SSuthDim)) &
                 * ((TInfDim/TSuthDim)**1.5_realType)

       ! Check the flow type we are having here.

       flowTypeTest: if(flowType == internalFlow) then

         ! Internal flow computation. Initialize the array to store
         ! the local total temperature and pressure and the local static
         ! density pressure and velocity magnitude to -1. Also initialize
         ! the sum of the local flow direction to zero.

         valLoc = -one
         dirLoc = zero

         ! Loop over the number od spectral modes and local blocks
         ! to  determine the variables described above.

         do sps=1,nTimeIntervalsSpectral
           do nn=1,nDom

             ! Set the pointer for the boundary conditions to make the
             ! code more readable.

             BCData => flowDoms(nn,1,sps)%BCData

             ! Loop over the number of boundary faces of the
             ! computational block.

             do mm=1,flowDoms(nn,1,sps)%nBocos

               ! Determine the maximum value of the scalar quantities
               ! for this subface and store them if these are larger
               ! than the currently stored values.

               tmp = maxValueSubface(BCData(mm)%ptInlet)
               valLoc(1) = max(valLoc(1),tmp)

               tmp = maxValueSubface(BCData(mm)%ttInlet)
               valLoc(2) = max(valLoc(2),tmp)

               tmp = maxValueSubface(BCData(mm)%rho)
               valLoc(3) = max(valLoc(3),tmp)

               tmp = maxValueSubface(BCData(mm)%ps)
               valLoc(4) = max(valLoc(4),tmp)

               ! Determine the velocity magnitude and sum up the
               ! direction.

               call velMagnAndDirectionSubface(tmp, dirLoc, &
                                               BCData, mm)
               valLoc(5) = max(valLoc(5),tmp)

             enddo
           enddo
         enddo

         ! Determine the global maxima of valLoc and the sum
         ! of dirLoc.

         call mpi_allreduce(valLoc, valGlob, 5, sumb_real, mpi_max, &
                            SUmb_comm_world, ierr)

         call mpi_allreduce(dirLoc, dirGlob, 3, sumb_real, mpi_sum, &
                            SUmb_comm_world, ierr)

         ! Create a unit vector for the global direction.

         tmp         = one/max(eps,sqrt(dirGlob(1)**2 + dirGlob(2)**2 &
                     +                  dirGlob(3)**2))
         dirGlob(1) = tmp*dirGlob(1)
         dirGlob(2) = tmp*dirGlob(2)
         dirGlob(3) = tmp*dirGlob(3)

         ! Store this direction for the free stream; this will only be
         ! used for initialization.

         velDirFreestream = dirGlob

         ! Determine the situation we are having here.

         if(valGlob(1) > zero .and. valGlob(2) > zero) then

           ! Total conditions are present.
           ! Compute the value of gamma, which is gammaInf.
           ! This is not entirely correct, because there may be
           ! a difference between the static and total temperature,
           ! but it is only used for an initialization.

           call computeGamma(valGlob(2), gammaInf, 1_intType)
           gm1 = gammaInf - one

           ! Check if a static pressure is present. If so, estimate
           ! a much number. If not present, set Mach to 0.5.
           ! Limit the estimated Mach number to 0.5 for stability
           ! reasons.

           if(valGlob(4) > zero) then
             ratio = max((valGlob(1)/valGlob(4)), 1.007017518_realType)
       !     ratio = max((valGlob(1)/valGlob(4)), 1.0000007_realType)
             tmp   = ratio**(gm1/gammaInf)

             Mach = sqrt(two*(tmp-one)/gm1)
             Mach = min(Mach,half)
           else
             Mach = half
           endif

           ! Set a value of pInfDim and TInfDim. This is just for
           ! initialization. The final solution is independent of it.

           tmp = one + half*gm1*Mach*Mach
           pInfDim = valGlob(1)/(tmp**(gammaInf/gm1))
           TInfDim = valGlob(2)/tmp

           ! Compute the density.

           rhoInfDim = pInfDim/(RGasDim*TInfDim)

         else if(valGlob(3) > zero .and. valGlob(4) > zero .and. &
                 valGlob(5) > zero) then

           ! Density, pressure and velocity magnitude are present.
           ! Compute the dimensional temperature and the corresponding
           ! value of gamma.

           rhoInfDim = valGlob(3)
           pInfDim   = valGlob(4)
           TInfDim   = pInfDim/(RGasDim*rhoInfDim)

           call computeGamma(TInfDim, gammaInf, 1_intType)

           ! Compute the Mach number.

           Mach = valGlob(5)/sqrt(gammaInf*pInfDim/rhoInfDim)

         else

           ! Not enough boundary data is present for initialization.
           ! This typically occurs when running the code in coupled
           ! mode with another CFD code from which it gets the data.
           ! If the code is run in stand alone mode, terminate.

           if( standAloneMode ) then
             if(myID == 0)                      &
               call terminate("referenceState", &
                              "Not enough boundary data is present to &
                              &define a well posed problem for an &
                              &internal flow computation")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

           ! Multi-disciplinary mode.
           ! Use rhoIni, pIni, MachIni and velDirIni for initialization.
           ! Processor 0 prints a warning message.

           rhoInfDim        = rhoIni
           pInfDim          = pIni
           TInfDim          = pInfDim/(RGasDim*rhoInfDim)
           Mach             = MachIni
           velDirFreestream = velDirIni

           ! Compute the corresponding value of gamma.

           call computeGamma(TInfDim, gammaInf, 1_intType)

           if(myID == 0) then

             print "(a)", "#"
             print "(a)", "#*==================== !!! Warning !!! &
                          &======================"
             print "(a)", "# Not enough boundary data is present to &
                          &define a well posed problem for"
             print "(a)", "# an internal flow computation"
             print "(a)", "# It is assumed that the data is supplied &
                          &from a different code"
             print "(a)", "#*=====================================&
                          &======================"
             print "(a)", "#"

           endif

         endif

         ! Set MachCoef to Mach. Seen the previous lines this is quite
         ! arbitrary, but for an internal flow the coefficients are not
         ! so important anyway.

         MachCoef = Mach

         ! Compute the value of the molecular viscosity corresponding
         ! to the computed TInfDim.

         muDim = muSuthDim                                    &
               * ((TSuthDim + SSuthDim)/(TInfDim + SSuthDim)) &
               * ((TInfDim/TSuthDim)**1.5_realType)

         ! In case the reference pressure, density and temperature were
         ! not specified, set them to the 1.0, i.e. a dimensional
         ! computation is performed.

         if(pRef   <= zero) pRef   = one
         if(rhoRef <= zero) rhoRef = one
         if(TRef   <= zero) TRef   = one

       else flowTypeTest

         ! External flow. Compute the value of gammaInf.

         call computeGamma(tempFreestream, gammaInf, 1_intType)

         ! In case of a viscous problem, compute the
         ! dimensional free stream density and pressure.

         if(equations == NSEquations .or. &
            equations == RANSEquations) then

           ! Compute the x, y, and z-components of the Mach number
           ! relative to the body; i.e. the mesh velocity must be
           ! taken into account here.

           mx = MachCoef*velDirFreestream(1)
           my = MachCoef*velDirFreestream(2)
           mz = MachCoef*velDirFreestream(3)

           ! Reynolds number per meter, the viscosity using sutherland's
           ! law and the free stream velocity relative to the body.

           Re = Reynolds/ReynoldsLength
           muDim = muSuthDim                    &
                 * ((TSuthDim + SSuthDim)       &
                 /  (tempFreestream + SSuthDim)) &
                 * ((tempFreestream/TSuthDim)**1.5)
           v  = sqrt((mx*mx + my*my + mz*mz) &
              *      gammaInf*RGasDim*tempFreestream)

           ! Compute the free stream density and pressure.
           ! Set TInfDim to tempFreestream.

           rhoInfDim = Re*muDim/v
           pInfDim   = rhoInfDim*RGasDim*tempFreestream
           TInfDim   = tempFreestream

         endif

         ! In case the reference pressure, density and temperature were
         ! not specified, set them to the infinity values.

         if(pRef   <= zero) pRef   = pInfDim
         if(rhoRef <= zero) rhoRef = rhoInfDim
         if(TRef   <= zero) TRef   = TInfDim

       endif flowTypeTest

       ! Compute the value of muRef, such that the nonDimensional
       ! equations are identical to the dimensional ones.
       ! Note that in the non-dimensionalization of muRef there is
       ! a reference length. However this reference length is 1.0
       ! in this code, because the coordinates are converted to
       ! meters.

       muRef = sqrt(pRef*rhoRef)

       ! Compute timeRef for a correct nonDimensionalization of the
       ! unsteady equations. Some story as for the reference viscosity
       ! concerning the reference length.

       timeRef = sqrt(rhoRef/pRef)

       ! Compute the nonDimensional pressure, density, velocity,
       ! viscosity and gas constant.

       pInf   = pInfDim/pRef
       rhoInf = rhoInfDim/rhoRef
       uInf   = Mach*sqrt(gammaInf*pInf/rhoInf)
       RGas   = RGasDim*rhoRef*TRef/pRef
       muInf  = muDim/muRef

       !=================================================================

       contains

         !===============================================================

         function maxValueSubface(var)
!
!        ****************************************************************
!        *                                                              *
!        * maxValueSubface determines the maximum value of var for      *
!        * currently active subface. If var is not associated this      *
!        * function returs -1.0.                                        *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function type
!
         real(kind=realType) :: maxValueSubface
!
!        Function argument.
!
         real(kind=realType), dimension(:,:), pointer :: var
!
!        Local variables.
!
         integer(kind=intType) :: i, j
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Initialize the function to -1 and return immediately if
         ! var is not associated with data.

         maxValueSubface = -one
         if(.not. associated(var)) return

         ! Loop over the owned faces of the subface. As the cell range
         ! may contain halo values, the nodal range is used.

         do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
           do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd
             maxValueSubface = max(maxValueSubface,var(i,j))
           enddo
         enddo

         end function maxValueSubface

       end subroutine referenceState

       !=================================================================

       subroutine velMagnAndDirectionSubface(vmag, dir, BCData, mm)
!
!      ******************************************************************
!      *                                                                *
!      * VelMagnAndDirectionSubface determines the maximum value    *
!      * of the magnitude of the velocity as well as the sum of the     *
!      * flow directions for the currently active subface.              *
!      *                                                                *
!      ******************************************************************
!
       use block
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: mm

       real(kind=realType), intent(out) :: vmag
       real(kind=realType), dimension(3), intent(inout) :: dir

       type(BCDataType), dimension(:), pointer :: BCData
!
!      Local variables.
!
       integer(kind=intType) :: i, j
       real(kind=realType)   :: vel
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize vmag to -1.0.

       vmag = -one

       ! Check if the velocity is prescribed.

       if( associated(BCData(mm)%velx) .and. &
           associated(BCData(mm)%vely) .and. &
           associated(BCData(mm)%velz) ) then

         ! Loop over the owned faces of the subface. As the cell range
         ! may contain halo values, the nodal range is used.

         do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
           do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

             ! Compute the magnitude of the velocity and compare it
             ! with the current maximum. Store the maximum of the two.

             vel  = sqrt(BCData(mm)%velx(i,j)**2 &
                  +      BCData(mm)%vely(i,j)**2 &
                  +      BCData(mm)%velz(i,j)**2)
             vmag = max(vmag, vel)

             ! Compute the unit vector of the velocity and add it to dir.

             vel    = one/max(eps,vel)
             dir(1) = dir(1) + vel*BCData(mm)%velx(i,j)
             dir(2) = dir(2) + vel*BCData(mm)%vely(i,j)
             dir(3) = dir(3) + vel*BCData(mm)%velz(i,j)

           enddo
         enddo
       endif

       ! Check if the velocity direction is prescribed.

       if( associated(BCData(mm)%flowXdirInlet) .and. &
           associated(BCData(mm)%flowYdirInlet) .and. &
           associated(BCData(mm)%flowZdirInlet) ) then

         ! Add the unit vectors to dir by looping over the owned
         ! faces of the subfaces. Again the nodal range must be
         ! used for this.

         do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
           do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

             dir(1) = dir(1) + BCData(mm)%flowXdirInlet(i,j)
             dir(2) = dir(2) + BCData(mm)%flowYdirInlet(i,j)
             dir(3) = dir(3) + BCData(mm)%flowZdirInlet(i,j)

           enddo
         enddo

       endif

       end subroutine velMagnAndDirectionSubface
