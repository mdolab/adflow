!
!     ******************************************************************
!     *                                                                *
!     * File:          assertionsADjoint.f90                           *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine assertionsADjoint(level)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Perform all the necessary assertions before running the        *
  !     * discrete adjoint solver. This takes into account the physical  *
  !     * models and boundary conditions that are currently supported.   *
  !     * Since the BCs are identical for all grid levels and time       *
  !     * instances, only the finest grid and 1st time instance are      *
  !     * tested.                                                        *
  !     *                                                                *
  !     ******************************************************************
  !
  use BCTypes
  use blockpointers
  use flowVarRefState     ! viscous
  use inputDiscretization ! orderFlow
  use inputPhysics        ! equationMode
  use cgnsgrid            ! oversetPresent
  use section             ! secID
  use inputtimespectral   ! nTimeSpectralInterval
  implicit none
  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level
  !
  !     Local variables.
  !
  integer(kind=intType) :: nn, mm, boundary
  character(len=2*maxStringLen) :: errorMessage

  integer(kind=intType) ::  sps, secID, nTime
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Physical model assertions.                                     *
  !     *                                                                *
  !     ******************************************************************
  !

  ! Assert that it is a steady problem.

  if( equationMode == unsteady ) &
       call terminate("assertionsADjoint", &
       "Cannot handle unsteady yet.")

  ! Assert that it is an inviscid problem.

  if( viscous ) &
       call terminate("assertionsADjoint", &
       "Cannot handle viscous terms yet.")

  ! Assert that there are no overset grids.

  if( oversetPresent ) &
       call terminate("assertionsADjoint", &
       "Cannot handle overset grids yet.")

  ! Loop over the local domains.

  domainLoop: do nn=1,nDom
     nTime = nTimeIntervalsSpectral

     spectralLoop: do sps = 1, nTime

        !set the correct block

        call setPointers(nn,level,sps)

        ! Loop over the subfaces.

        subfaceLoop: do mm=1,nSubface

           ! Verify ADjoint support for the boundary condition type
           ! (check adjoint/residualInput/ApplyBCsAdj.f90)

           boundary = BCType(mm)

           if( boundary == SymmPolar        .or. &
                boundary == NSWallAdiabatic  .or. &
                boundary == NSWallIsothermal .or. &
                boundary == SubsonicInflow   .or. &
                boundary == SubsonicOutflow  .or. &
                boundary == MassBleedInflow  .or. &
                boundary == MassBleedOutflow .or. &
                boundary == mDot             .or. &
                boundary == Thrust           .or. &
                boundary == Extrap           .or. &
                boundary == DomainInterfaceAll  .or. &
                boundary == SupersonicInflow ) then
              write(errorMessage,99) "Cannot handle specified BC type (",&
                   boundary, &
                   ") yet."
              call terminate("assertionsADjoint", errorMessage)
           endif

        enddo subfaceLoop

     enddo spectralLoop

  enddo domainLoop
         
99       format(a,i3,a)

end subroutine assertionsADjoint

