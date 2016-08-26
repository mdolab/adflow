!
!       File:          numberOfIsoVariables.f90                        
!       Author:        Gaetan K. W. Kenway                             
!       Starting date: 07-21-2013                                      
!       Last modified: 07-21-2013                                      
!
subroutine numberOfIsoSurfVariables(nIsoSolVar)
  !
  !       numberOfVolSolVariables determines the number of variables     
  !       to be written on the isosurface. These are similar to the      
  !       volume variables.                                              
  !
  use flowVarRefState
  use inputPhysics
  use extraOutput
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(out) :: nIsoSolVar
  !
  !       Begin execution                                                
  !
  nIsoSolvar   = 0

  ! Check whether or not some additional solution variables must
  ! be written.
  if (isoWriteRho)           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteVx )           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteVy )           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteVz )           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteP )            nIsoSolvar = nIsoSolvar + 1
  if( isoWriteTurb )         nIsoSolvar  = nIsoSolvar + (nw - nwf)
  if( isoWriteMx )           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteMy )           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteMz )           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteRVx )          nIsoSolvar = nIsoSolvar + 1
  if( isoWriteRVy )          nIsoSolvar = nIsoSolvar + 1
  if( isoWriteRVz )          nIsoSolvar = nIsoSolvar + 1
  if( isoWriteRhoe )         nIsoSolvar = nIsoSolvar + 1
  if( isoWriteTemp )         nIsoSolvar = nIsoSolvar + 1
  if( isoWriteCp )           nIsoSolvar = nIsoSolvar + 1
  if( isoWriteMach )         nIsoSolvar = nIsoSolvar + 1
  if( isoWriteRMach )        nIsoSolvar = nIsoSolvar + 1
  if( isoWriteMachTurb )     nIsoSolvar = nIsoSolvar + 1
  if( isoWriteEddyVis )      nIsoSolvar = nIsoSolvar + 1
  if( isoWriteRatioEddyVis ) nIsoSolvar = nIsoSolvar + 1
  if( isoWriteDist )         nIsoSolvar = nIsoSolvar + 1
  if( isoWriteVort )         nIsoSolvar = nIsoSolvar + 1
  if( isoWriteVortx )        nIsoSolvar = nIsoSolvar + 1
  if( isoWriteVorty )        nIsoSolvar = nIsoSolvar + 1
  if( isoWriteVortz )        nIsoSolvar = nIsoSolvar + 1
  if( isoWritePtotloss )     nIsoSolvar = nIsoSolvar + 1
  if( isoWriteShock )        nIsoSolVar = nIsoSolVar + 1
  if (isoWriteFilteredShock) nIsoSolVar = nIsoSolVar + 1

  ! Check the discrete variables.

  if( isoWriteResRho )  nIsoSolvar  = nIsoSolvar + 1
  if( isoWriteResMom )  nIsoSolvar  = nIsoSolvar + 3
  if( isoWriteResRhoe ) nIsoSolvar  = nIsoSolvar + 1
  if( isoWriteResTurb ) nIsoSolvar  = nIsoSolvar   &
       + (nw - nwf)

  if( isoWriteBlank ) nIsoSolvar  = nIsoSolvar + 1

end subroutine numberOfIsoSurfVariables
