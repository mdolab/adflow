!
!      ******************************************************************
!      *                                                                *
!      * File:          applyAllBCAdj.f90                               *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      *                C.A.(Sandy) Mader                               *
!      * Starting date: 04-16-2008                                      *
!      * Last modified: 04-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine applyAllBCAdj(wInfAdj,pInfCorrAdj,wAdj, pAdj, &
                              siAdj, sjAdj, skAdj, volAdj, normAdj, &
                              iCell, jCell, kCell,secondHalo)
!
!      ******************************************************************
!      *                                                                *
!      * applyAllBCAdj applies the possible boundary conditions for the *
!      * halo cells adjacent to the cell for which the residual needs   *
!      * to be computed.                                                *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers!, only : ie, ib, je, jb, ke, kb, nBocos, &
                         !         BCFaceID, BCType, BCData,p,w
       use flowVarRefState
       use inputDiscretization !precond,choimerkle, etc...
       implicit none

!!$       use blockPointers
!!$       use flowVarRefState
!!$       use inputDiscretization
!!$       use inputTimeSpectral
!!$       use iteration
!!$       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo

       integer(kind=intType) :: iCell, jCell, kCell

       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
                                                   intent(inout) :: wAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2),    &
                                                   intent(inout) :: pAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,3), &
                                                   intent(in) :: siAdj, sjAdj, skAdj
       !real(kind=realType), dimension(0:0,0:0,0:0), intent(in) :: volAdj
       real(kind=realType), intent(in) :: volAdj,pInfCorrAdj
       real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(in) :: normAdj
       real(kind=realType), dimension(nw),intent(in)::wInfAdj

!
!      Local variables.
!
       integer(kind=intType) :: nn, sps
       integer(kind=intType)::i,j,k,ii,jj,kk,l
       integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

       logical :: correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       
!moved outside
!!$       ! Determine whether or not the total energy must be corrected
!!$       ! for the presence of the turbulent kinetic energy.
!!$
!!$       if( kPresent ) then
!!$         if((currentLevel <= groundLevel) .or. turbCoupled) then
!!$           correctForK = .true.
!!$         else
!!$           correctForK = .false.
!!$         endif
!!$       else
!!$         correctForK = .false.
!!$       endif

!!$       ! Loop over the number of spectral solutions.
!!$
!!$       spectralLoop: do sps=1,nTimeIntervalsSpectral
!!$
!!$         ! Loop over the number of blocks.
!!$
!!$         domains: do nn=1,nDom
!!$
!!$           ! Set the pointers for this block.
!!$
!!$           call setPointers(nn, currentLevel, sps)
!!$
!!$           ! Apply all the boundary conditions. The order is important.

           ! The symmetry boundary conditions.
       
!*************************
!       call bcSymmAdj(wAdj,pAdj,normAdj,iCell,jCell,kCell,secondHalo)
!**************************
!###       call bcSymmPolar(secondHalo)

!!$       ! call bcEulerWall(secondHalo, correctForK)
!!$
!!$       ! The viscous wall boundary conditions.
!!$
!!$       call bcNSWallAdiabatic( secondHalo, correctForK)
!!$       call bcNSWallIsothermal(secondHalo, correctForK)
!!$
!!$       ! The farfield is a special case, because the treatment
!!$       ! differs when preconditioning is used. Make that distinction
!!$       ! and call the appropriate routine.
!!$       

!!$!*******************************
       select case (precond)
          
       case (noPrecond)

          call bcFarfieldAdj(secondHalo,wInfAdj,pInfCorrAdj, wAdj,pAdj,      &
               siAdj, sjAdj, skAdj, normAdj,iCell,jCell,kCell)

       case (Turkel)
          call terminate("applyAllBC", "Farfield boundary conditions for Turkel preconditioner not implemented")
          
       case (ChoiMerkle)
          call terminate("applyAllBC", "Farfield boundary conditions for Choi and Merkle preconditioner not implemented")

       end select
!!$!******************************8

       
!!$
!!$       ! Subsonic outflow and bleed outflow boundaries.
!!$       
!!$       call bcSubsonicOutflow(secondHalo, correctForK)
!!$       
!!$       ! Subsonic inflow boundary.
!!$       
!!$       call bcSubsonicInflow(secondHalo, correctForK)
!!$       
!!$       ! Bleed inflow regions.
!!$       
!!$       call bcBleedInflow( secondHalo, correctForK)
!!$       
!!$       ! Engine boundary conditions. Not implemented yet.
!!$
!!$       call bcMdot(secondHalo, correctForK)
!!$       call bcThrust(secondHalo, correctForK)
!!$       
!!$       ! Extrapolation boundary conditions; this also includes
!!$       ! the supersonic outflow boundary conditions. The difference
!!$       ! between the two is that the extrap boundary conditions
!!$       ! correspond to singular lines and supersonic outflow
!!$       ! boundaries to physical boundaries. The treatment however
!!$       ! is identical.
!!$       
!!$       call bcExtrap(secondHalo, correctForK)
!!$       
       ! Inviscid wall boundary conditions.
 
       call bcEulerWallAdj(secondHalo, wAdj,pAdj,      &
            siAdj, sjAdj, skAdj, normAdj,iCell,jCell,kCell)
       
!!$       ! Domain-interface boundary conditions,
!!$       ! when coupled with other solvers.
!!$       
!!$       call bcDomainInterface(secondHalo, correctForK)
!!$       
!!$       ! Supersonic inflow boundary conditions.
!!$       
!!$       call bcSupersonicInflow(secondHalo, correctForK)

!!$         enddo domains
!!$       enddo spectralLoop

     end subroutine applyAllBCAdj
