!
!      ******************************************************************
!      *                                                                *
!      * File:          allocMemFlovar.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-06-2003                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine allocMemFlovarPart1(sps,level)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * allocMemFlovarPart1 allocates the memory for the flow          *
  !      * variables w and p for all the blocks on the given multigrid    *
  !      * level and spectral solution sps.                               *
  !      *                                                                *
  !      ******************************************************************
  !
  use block
  use constants
  use flowVarRefState
  use inputIteration
  use inputPhysics
  use inputUnsteady
  use inputphysics
  use iteration
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: sps, level
  !
  !      Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: nn
  integer(kind=intType) :: il, jl, kl, ie, je, ke, ib, jb, kb
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Loop over the domains.

  domains: do nn=1,nDom

     ! Store some dimensions a bit easier.

     il = flowDoms(nn,level,sps)%il
     jl = flowDoms(nn,level,sps)%jl
     kl = flowDoms(nn,level,sps)%kl

     ie = flowDoms(nn,level,sps)%ie
     je = flowDoms(nn,level,sps)%je
     ke = flowDoms(nn,level,sps)%ke

     ib = flowDoms(nn,level,sps)%ib
     jb = flowDoms(nn,level,sps)%jb
     kb = flowDoms(nn,level,sps)%kb

     ! Allocate the memory for the independent variables.
     ! Memory is allocated for the turbulent variables (if any) if
     ! the current level is smaller or equal to the multigrid start
     ! level or if the turbulent transport equations are solved in
     ! a coupled manner.

     if(level <= mgStartlevel .or. turbTreatment == coupled) then
        allocate(flowDoms(nn,level,sps)%w(0:ib,0:jb,0:kb,1:nw), &
             stat=ierr)
     else
        allocate(flowDoms(nn,level,sps)%w(0:ib,0:jb,0:kb,1:nMGVar), &
             stat=ierr)
     endif
     if(ierr /= 0)                           &
          call returnFail("allocMemFlovarPart1", &
          "Memory allocation failure for w")

     ! Alloc mem for dadi
     allocate(flowDoms(nn,level,sps)%dadidata(ie,je,ke,10), stat=ierr)

     ! Alloc mem for nodal gradients 
     allocate(flowDoms(nn,level,sps)%ux(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%uy(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%uz(il,jl,kl), stat=ierr)

     allocate(flowDoms(nn,level,sps)%vx(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%vy(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%vz(il,jl,kl), stat=ierr)

     allocate(flowDoms(nn,level,sps)%wx(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%wy(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%wz(il,jl,kl), stat=ierr)

     allocate(flowDoms(nn,level,sps)%qx(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%qy(il,jl,kl), stat=ierr)
     allocate(flowDoms(nn,level,sps)%qz(il,jl,kl), stat=ierr)

     ! Allocate memory for the pressure.
     allocate(flowDoms(nn,level,sps)%p(0:ib,0:jb,0:kb), stat=ierr)
     if(ierr /= 0)                           &
          call returnFail("allocMemFlovarPart1", &
          "Memory allocation failure for p")

     ! Allocate memory for the speed of sound squared
     allocate(flowDoms(nn,level,sps)%aa(0:ib,0:jb,0:kb), stat=ierr)
     if(ierr /= 0)                           &
          call returnFail("allocMemFlovarPart1", &
          "Memory allocation failure for p")

     ! The eddy viscosity for eddy viscosity models.
     ! Although a dependent variable, it is allocated on all grid
     ! levels, because the eddy viscosity might be frozen in the
     ! multigrid.

     ! Always allocate rev due to reverse mode - Peter Lyu 
     !if( eddyModel ) then
     allocate(flowDoms(nn,level,sps)%rev(0:ib,0:jb,0:kb), &
          stat=ierr)
     if(ierr /= 0)                           &
          call returnFail("allocMemFlovarPart1", &
          "Memory allocation failure for rev")
     !endif

     ! If this is the finest grid some more memory must be allocated.

     fineLevelTest: if(level == 1) then

        ! Allocate the memory for gamma and initialize it to
        ! the constant gamma value.

        allocate(flowDoms(nn,level,sps)%gamma(0:ib,0:jb,0:kb),    &
             stat=ierr)
        if(ierr /= 0)                           &
             call returnFail("allocMemFlovarPart1", &
             "Memory allocation failure for gamma.")

        flowDoms(nn,level,sps)%gamma = gammaConstant

        ! The laminar viscosity for viscous computations.
        ! Always allocate rlv due to reverse mode - Peter Lyu 
        !if( viscous ) then
        allocate(flowDoms(nn,level,sps)%rlv(0:ib,0:jb,0:kb), &
             stat=ierr)
        if(ierr /= 0)                           &
             call returnFail("allocMemFlovarPart1", &
             "Memory allocation failure for rlv")
        !endif

        ! The state vectors in the past for unsteady computations.

        if(equationMode          == unsteady .and. &
             timeIntegrationScheme == BDF) then
           allocate( &
                flowDoms(nn,level,sps)%wOld(nOldLevels,2:il,2:jl,2:kl,nw), &
                stat=ierr)
           if(ierr /= 0)                           &
                call returnFail("allocMemFlovarPart1", &
                "Memory allocation failure for wOld")

           ! Initialize wOld to zero, such that it is initialized.
           ! The actual values do not matter.

           flowDoms(nn,level,sps)%wOld = zero

        ! *******************************
        ! Added by HDN
        ! *******************************
        else if ( equationMode == unsteady .and. &
             timeIntegrationScheme == MD) then
           allocate( &
                flowDoms(nn,level,sps)%wOld(nOldLevels,2:il,2:jl,2:kl,nw), &
                stat=ierr)
           if(ierr /= 0)                           &
                call returnFail("allocMemFlovarPart1", &
                "Memory allocation failure for wOld")

           flowDoms(nn,level,sps)%wOld = zero
        endif

        ! If this is the 1st spectral solution (note that we are
        ! already on the finest grid) and the rans equations are
        ! solved, allocate the memory for the arrays used for the
        ! implicit boundary condition treatment. Normally this should
        ! only be allocated for RANS but the derivative calcs require
        ! these be allocated.

        sps1RansTest: if(sps == 1) then 
           allocate(flowDoms(nn,level,sps)%bmti1(je,ke,nt1:nt2,nt1:nt2), &
                flowDoms(nn,level,sps)%bmti2(je,ke,nt1:nt2,nt1:nt2), &
                flowDoms(nn,level,sps)%bmtj1(ie,ke,nt1:nt2,nt1:nt2), &
                flowDoms(nn,level,sps)%bmtj2(ie,ke,nt1:nt2,nt1:nt2), &
                flowDoms(nn,level,sps)%bmtk1(ie,je,nt1:nt2,nt1:nt2), &
                flowDoms(nn,level,sps)%bmtk2(ie,je,nt1:nt2,nt1:nt2), &
                flowDoms(nn,level,sps)%bvti1(je,ke,nt1:nt2), &
                flowDoms(nn,level,sps)%bvti2(je,ke,nt1:nt2), &
                flowDoms(nn,level,sps)%bvtj1(ie,ke,nt1:nt2), &
                flowDoms(nn,level,sps)%bvtj2(ie,ke,nt1:nt2), &
                flowDoms(nn,level,sps)%bvtk1(ie,je,nt1:nt2), &
                flowDoms(nn,level,sps)%bvtk2(ie,je,nt1:nt2), &
                stat=ierr)
           if(ierr /= 0)                           &
                call returnFail("allocMemFlovarPart1", &
                "Memory allocation failure for bmti1, etc")

        endif sps1RansTest

     endif fineLevelTest

  enddo domains

end subroutine allocMemFlovarPart1

!      ==================================================================

subroutine allocMemFlovarPart2(sps, level)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * AllocMemFlovarPart2 allocates the memory for the dependent  *
  !      * flow variables and iteration variables for all the blocks on   *
  !      * the given multigrid level and spectral solution sps. Some      *
  !      * variables are only allocated on the coarser grids, e.g. the    *
  !      * multigrid forcing terms and the state vector upon entrance on  *
  !      * the mg level. Other variables are only allocated on the finest *
  !      * mesh. These are typically dependent variables like laminar     *
  !      * viscosity, or residuals, time step, etc. Exceptions are        *
  !      * pressure and eddy viscosity. Although these are dependent      *
  !      * variables, they are allocated on all grid levels.              *
  !      *                                                                *
  !      ******************************************************************
  !
  use block
  use constants
  use flowVarRefState
  use inputPhysics
  use inputDiscretization
  use inputIteration
  use inputUnsteady
  use iteration
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: sps, level
  !
  !      Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: nn, mm
  integer(kind=intType) :: il, jl, kl, ie, je, ke, ib, jb, kb
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Loop over the domains.

  domains: do nn=1,nDom

     ! Store some dimensions a bit easier.

     il = flowDoms(nn,level,sps)%il
     jl = flowDoms(nn,level,sps)%jl
     kl = flowDoms(nn,level,sps)%kl

     ie = flowDoms(nn,level,sps)%ie
     je = flowDoms(nn,level,sps)%je
     ke = flowDoms(nn,level,sps)%ke

     ib = flowDoms(nn,level,sps)%ib
     jb = flowDoms(nn,level,sps)%jb
     kb = flowDoms(nn,level,sps)%kb

     ! Block is moving. Allocate the memory for s, sFaceI,
     ! sFaceJ and sFaceK.

     allocate(flowDoms(nn,level,sps)%s(ie,je,ke,3),  &
          flowDoms(nn,level,sps)%sFaceI(0:ie,je,ke), &
          flowDoms(nn,level,sps)%sFaceJ(ie,0:je,ke), &
          flowDoms(nn,level,sps)%sFaceK(ie,je,0:ke), stat=ierr)
     if(ierr /= 0)                              &
          call returnFail("allocMemFlovarPart2", &
          "Memory allocation failure for s, &
          &sFaceI, sFaceJ and sFaceK.")

     ! Extra face velocities for ALE
     if (equationMode == unSteady .and. useALE) then 
        allocate( &
             flowDoms(nn,level,sps)%sVeloIALE(0:ie,je,ke,3), &
             flowDoms(nn,level,sps)%sVeloJALE(ie,0:je,ke,3), &
             flowDoms(nn,level,sps)%sVeloKALE(ie,je,0:ke,3), &
             flowDoms(nn,level,sps)%sFaceIALE(0:nALEsteps,0:ie,je,ke), &
             flowDoms(nn,level,sps)%sFaceJALE(0:nALEsteps,ie,0:je,ke), &
             flowDoms(nn,level,sps)%sFaceKALE(0:nALEsteps,ie,je,0:ke), stat=ierr)
        if(ierr /= 0)                              &
             call returnFail("allocMemFlovarPart2", &
             "Memory allocation failure for &
             sVeloIALE, sVeloJALE and sVeloKALE; &
             sFaceIALE, sFaceJALE and sFaceKALE.")
     end if


     ! Test if we are on the finest mesh.

     fineLevelTest: if(level == 1) then

        ! Allocate the memory that must always be allocated.

        allocate( &
             flowDoms(nn,level,sps)%dw(0:ib,0:jb,0:kb,1:nw),  &
             flowDoms(nn,level,sps)%fw(0:ib,0:jb,0:kb,1:nwf), &
             flowDoms(nn,level,sps)%dtl(1:ie,1:je,1:ke),      &
             flowDoms(nn,level,sps)%radI(1:ie,1:je,1:ke),     &
             flowDoms(nn,level,sps)%radJ(1:ie,1:je,1:ke),     &
             flowDoms(nn,level,sps)%radK(1:ie,1:je,1:ke),     &
             flowDoms(nn,level,sps)%scratch(0:ib,0:jb,0:kb,5), &
             stat=ierr)
        if(ierr /= 0)                              &
             call returnFail("allocMemFlovarPart2", &
             "Memory allocation failure for dw, fw, dwOld, fwOld, &
             &gamma, dtl and the spectral radii.")

        ! Initialize dw and fw to zero.

        flowDoms(nn,level,sps)%dw = zero
        flowDoms(nn,level,sps)%fw = zero


        ! Extra variables for ALE
        if (equationMode == unSteady .and. useALE) then 
           allocate( &
                flowDoms(nn,level,sps)%dwALE(0:nALEsteps,0:ib,0:jb,0:kb,1:nw),  &
                flowDoms(nn,level,sps)%fwALE(0:nALEsteps,0:ib,0:jb,0:kb,1:nwf), &
                stat=ierr)
           if(ierr /= 0)                              &
                call returnFail("allocMemFlovarPart2", &
                "Memory allocation failure for dwALE, fwALE.")
           
           flowDoms(nn,level,sps)%dwALE = zero
           flowDoms(nn,level,sps)%fwALE = zero
        end if

        ! Allocate the memory for the zeroth runge kutta stage
        allocate(flowDoms(nn,level,sps)%wn(2:il,2:jl,2:kl,1:nMGVar), &
             flowDoms(nn,level,sps)%pn(2:il,2:jl,2:kl), stat=ierr)
        if(ierr /= 0)                              &
             call returnFail("allocMemFlovarPart2", &
             "Memory allocation failure for wn and pn")

        ! For unsteady mode using Runge-Kutta schemes allocate the
        ! memory for dwOldRK.

        if(equationMode          == unsteady .and. &
             timeIntegrationScheme == explicitRK) then

           mm = nRKStagesUnsteady - 1
           allocate(flowDoms(nn,level,sps)%dwOldRK(mm,il,jl,kl,nw), &
                stat=ierr)
           if(ierr /= 0) &
                call returnFail("allocMemFlovarPart2", &
                "Memory allocation failure for dwOldRK.")
        endif

     else fineLevelTest

        ! Coarser level. Allocate the memory for the multigrid
        ! forcing term and the state variables upon entry.

        allocate(flowDoms(nn,level,sps)%p1(1:ie,1:je,1:ke),          &
             flowDoms(nn,level,sps)%w1(1:ie,1:je,1:ke,1:nMGVar), &
             flowDoms(nn,level,sps)%wr(2:il,2:jl,2:kl,1:nMGVar), &
             stat=ierr)
        if(ierr /= 0)                              &
             call returnFail("allocMemFlovarPart2", &
             "Memory allocation failure for p1, w1 &
             &and wr")

        ! Initialize w1 and p1 to zero, just that they
        ! are initialized.

        flowDoms(nn,level,sps)%p1 = zero
        flowDoms(nn,level,sps)%w1 = zero

     endif fineLevelTest

  enddo domains

end subroutine allocMemFlovarPart2
