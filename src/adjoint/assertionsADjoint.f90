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
      use block
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
      print *,'checking problem'
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

      print *,'problem ok'
!
!     ******************************************************************
!     *                                                                *
!     * Boundary condition assertions.                                 *
!     *                                                                *
!     ******************************************************************
!
      ! Relevant variables in module block
      !
      ! nDom:                  number of local computational blocks.
      ! flowDoms(:,:,:):       array of blocks. Dimensions are
      !                        (nDom,nLevels,nTimeIntervalsSpectral).
      !                        Type blockType.
      !
      ! type blockType
       ! nSubface:             Number of subfaces on this block.
       ! nViscBocos:           Number of viscous wall boundary subfaces.
       !                       These are numbered first.
       ! nInvBocos:            Number of inviscid wall boundary
       !                       subfaces. These follow the viscous ones.
       ! nBocos:               Total number of physical boundary
       !                       subfaces.
       ! n1to1:                Number of 1 to 1 block boundaries.
       ! nNonMatch:            Number of non-matching block boundaries.
       ! BCType(:):            Boundary condition type for each
       !                       subface. See the module BCTypes for
       !                       the possibilities.
       ! BCFaceID(:):          Block face location of each subface.
       !                       possible values are: iMin, iMax, jMin,
       !                       jMax, kMin, kMax. see also module
       !                       BCTypes.
      !
     
     
      print *,'checking boundary conditions'
      ! Loop over the local domains.

      domainLoop: do nn=1,nDom
         print *,'indomain',nn,level
         secID = flowDoms(nn,level,1)%sectionID
         nTime = nTimeIntervalsSpectral!sections(secID)%nTimeInstances

         spectralLoop: do sps = 1, nTime
            print *,'timeinstance',sps,ntime            
!!$            ! Loop over the subfaces.
!!$            
!!$            subfaceLoop: do mm=1,flowDoms(nn,level,sps)%nSubface
!!$               
!!$               ! Verify ADjoint support for the boundary condition type
!!$               ! (check adjoint/residualInput/penaltyStateBCsAdj.f90)
!!$               print *,'subface',mm,flowDoms(nn,level,sps)%BCType!(mm)
!!$!               boundary = flowDoms(nn,level,sps)%BCType(mm)
!!$!               print *,'boundary',boundary
!!$               if( boundary == SymmPolar        .or. &
!!$                    boundary == NSWallAdiabatic  .or. &
!!$                    boundary == NSWallIsothermal .or. &
!!$                    boundary == SubsonicInflow   .or. &
!!$                    boundary == SubsonicOutflow  .or. &
!!$                    boundary == MassBleedInflow  .or. &
!!$                    boundary == MassBleedOutflow .or. &
!!$                    boundary == mDot             .or. &
!!$                    boundary == Thrust           .or. &
!!$                    boundary == Extrap           .or. &
!!$                    boundary == DomainInterfaceAll  .or. &
!!$                    boundary == SupersonicInflow ) then
!!$                  write(errorMessage,99) "Cannot handle specified BC type (",&
!!$                       boundary, &
!!$                       ") yet."
!!$                  call terminate("assertionsADjoint", errorMessage)
!!$               endif
!!$               
!!$               
!!$            enddo subfaceLoop
            
         enddo spectralLoop
       
      enddo domainLoop

         ! Output format.
         
99       format(a,i3,a)
         
       end subroutine assertionsADjoint
       
