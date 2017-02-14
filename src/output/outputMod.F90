module outputMod
  !
  !       This local module contains variables used when writing the     
  !       grid and solution files.                                       
  !
  use constants, only : intType, maxStringLen
  implicit none

  ! nblocksCGNSblock(0:cgnsNDom): The number of local blocks per
  !                               cgns block in cumulative storage
  !                               format.
  ! blocksCGNSblock(nDom):        The corresponding local block ID's.

  integer(kind=intType), dimension(:), allocatable :: nblocksCGNSblock
  integer(kind=intType), dimension(:), allocatable :: blocksCGNSblock

  ! nDomPerProc(0:nProc):  The number of flow domains for each 
  !                        processor in cumulative storage format.
  ! IDsBegOrAllDoms(4,..): The CGNS block numbers and the beginning
  !                        indices for all of the flow domains on
  !                        every processor.

  integer(kind=intType), dimension(:),   allocatable :: nDomPerProc
  integer(kind=intType), dimension(:,:), allocatable :: IDsBegOrAllDoms

  ! nGridsToWrite:   Number of grid files to write.
  ! nVolSolToWrite:  Number of volume solution files to write.
  !                  For CGNS nVolSolToWrite == nGridsToWrite.
  ! nSurfSolToWrite: Number of surface solution files to write.

  integer(kind=intType) :: nGridsToWrite
  integer(kind=intType) :: nVolSolToWrite
  integer(kind=intType) :: nSurfSolToWrite

  ! gridFileNames(nGridsToWrite):      Names of the grid files to
  !                                    write.
  ! volSolFileNames(nVolSolToWrite):   Names of the volume solution
  !                                    files to write.
  ! surfSolFileNames(nSurfSolToWrite): Names of the surface solution
  !                                    files to write.
  ! fileIDs(nFilesToWrite):            Identifiers for the files to
  !                                    write. As the grids, volume
  !                                    solution and surface solution
  !                                    files are written one after
  !                                    the other, only one set is
  !                                    needed.
  ! cgnsBases(nFilesToWrite):          The CGNS base IDs of the
  !                                    files.

  character(len=maxStringLen), dimension(:), allocatable :: &
       gridFileNames, volSolFileNames, surfSolFileNames

  integer, dimension(:), allocatable :: fileIDs
  integer, dimension(:), allocatable :: cgnsBases
  integer, dimension(:), allocatable :: cgnsIsoSurfBases
  integer, dimension(:), allocatable :: cgnsLiftDistBases

  ! useLinksInCGNS:     Whether or not to use links in CGNS between
  !                     the grid and volume solution files. If not,
  !                     the grid and solution are written in the
  !                     same file.

  logical :: useLinksInCGNS

contains 

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

  subroutine numberOfSurfSolVariables(nSolVar)
    !
    !       numberOfSurfSolVariables determines the number of surface      
    !       variables to be written to the surface solution file.          
    !
    use precision
    use extraOutput
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(out) :: nSolVar

    ! Initialize the number of solution variables to zero.

    nSolVar = 0

    ! Determine the number of surface variables to be written.

    if( surfWriteRho )      nSolVar = nSolVar +1
    if( surfWriteP   )      nSolVar = nSolVar +1
    if( surfWriteTemp )     nSolVar = nSolVar +1
    if( surfWriteVx )       nSolVar = nSolVar +1
    if( surfWriteVy )       nSolVar = nSolVar +1
    if( surfWriteVz )       nSolVar = nSolVar +1
    if( surfWriteRVx )      nSolVar = nSolVar +1
    if( surfWriteRVy )      nSolVar = nSolVar +1
    if( surfWriteRVz )      nSolVar = nSolVar +1
    if( surfWriteCp )       nSolVar = nSolVar +1
    if( surfWritePtotloss ) nSolVar = nSolVar +1
    if( surfWriteMach )     nSolVar = nSolVar +1
    if( surfWriteRMach )    nSolVar = nSolVar +1
    if( surfWriteCf )       nSolVar = nSolVar +1
    if( surfWriteCh )       nSolVar = nSolVar +1
    if( surfWriteYplus )    nSolVar = nSolVar +1
    if( surfWriteCfx )      nSolVar = nSolVar +1
    if( surfWriteCfy )      nSolVar = nSolVar +1
    if( surfWriteCfz )      nSolVar = nSolVar +1
    if( surfWriteBlank )    nSolVar = nSolVar +1
    if( surfWriteSepSensor )      nSolVar = nSolVar +1
    if( surfWriteCavitation )     nsolVar = nsolVar +1
    if( surfWriteGC )             nsolVar = nsolVar +1

  end subroutine numberOfSurfSolVariables

  subroutine numberOfVolSolVariables(nVolSolvar, nVolDiscrVar)
    !
    !       numberOfVolSolVariables determines the number of volume        
    !       variables to be written to the solution file. A distinction is 
    !       made between solution variables and discrete variables. The    
    !       former discribes the actual solution, the latter is additional 
    !       info such as equation residuals.                               
    !
    use flowVarRefState
    use inputPhysics
    use extraOutput
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(out) :: nVolSolvar, nVolDiscrVar

    ! Initialize the number of solution variables to the number of
    ! independent variables and the number of discrete variables to 0.

    nVolSolvar   = nw
    nVolDiscrVar = 0

    ! Check whether or not some additional solution variables must
    ! be written.

    if( volWriteMx )           nVolSolvar = nVolSolvar + 1
    if( volWriteMy )           nVolSolvar = nVolSolvar + 1
    if( volWriteMz )           nVolSolvar = nVolSolvar + 1
    if( volWriteRVx )          nVolSolvar = nVolSolvar + 1
    if( volWriteRVy )          nVolSolvar = nVolSolvar + 1
    if( volWriteRVz )          nVolSolvar = nVolSolvar + 1
    if( volWriteRhoe )         nVolSolvar = nVolSolvar + 1
    if( volWriteTemp )         nVolSolvar = nVolSolvar + 1
    if( volWriteCp )           nVolSolvar = nVolSolvar + 1
    if( volWriteMach )         nVolSolvar = nVolSolvar + 1
    if( volWriteRMach )        nVolSolvar = nVolSolvar + 1
    if( volWriteMachTurb )     nVolSolvar = nVolSolvar + 1
    if( volWriteEddyVis )      nVolSolvar = nVolSolvar + 1
    if( volWriteRatioEddyVis ) nVolSolvar = nVolSolvar + 1
    if( volWriteDist )         nVolSolvar = nVolSolvar + 1
    if( volWriteVort )         nVolSolvar = nVolSolvar + 1
    if( volWriteVortx )        nVolSolvar = nVolSolvar + 1
    if( volWriteVorty )        nVolSolvar = nVolSolvar + 1
    if( volWriteVortz )        nVolSolvar = nVolSolvar + 1
    if( volWritePtotloss )     nVolSolvar = nVolSolvar + 1
    if( volWriteShock )        nVolSolvar = nVolSolvar + 1
    if( volWriteFilteredShock) nVolSolvar = nVolSolvar + 1
    if( volWriteGC  )          nVolSolvar = nVolSolvar + 1
    if( volWriteStatus )       nVolSolvar = nVolSolvar + 1
    ! Check the discrete variables.

    if( volWriteResRho )  nVolDiscrVar  = nVolDiscrVar + 1
    if( volWriteResMom )  nVolDiscrVar  = nVolDiscrVar + 3
    if( volWriteResRhoe ) nVolDiscrVar  = nVolDiscrVar + 1
    if( volWriteResTurb ) nVolDiscrVar  = nVolDiscrVar   &
         + (nw - nwf)

    if( volWriteBlank ) nVolDiscrVar  = nVolDiscrVar + 1

  end subroutine numberOfVolSolVariables

  subroutine copyDataBufSinglePrecision(val, buffer,      &
       iBeg, jBeg, kBeg, &
       iEnd, jEnd, kEnd, subRange)
    !
    !       copyDataBufSinglePrecision stores the given 1D buffer into the 
    !       subrange of the 3D single precision val array.                 
    !
    use precision
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: iBeg, jBeg, kBeg
    integer(kind=intType), intent(in) :: iEnd, jEnd, kEnd
    integer(kind=intType), dimension(3,2), intent(in) :: subRange

    real(kind=realType), dimension(*), intent(in) :: buffer
    real(kind=4), dimension(iBeg:iEnd,jBeg:jEnd,kBeg:kEnd), &
         intent(inout) :: val
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ll

    ! Copy the subrange into val.

    ll = 0
    do k=subRange(3,1), subRange(3,2)
       do j=subRange(2,1), subRange(2,2)
          do i=subRange(1,1), subRange(1,2)
             ll = ll + 1
             val(i,j,k) = real(buffer(ll), singleType)
          enddo
       enddo
    enddo

  end subroutine copyDataBufSinglePrecision

  !      ==================================================================

  subroutine copyDataBufDoublePrecision(val, buffer,      &
       iBeg, jBeg, kBeg, &
       iEnd, jEnd, kEnd, subRange)
    !
    !       copyDataBufDoublePrecision stores the given 1D buffer into the 
    !       subrange of the 3D double precision val array.                 
    !
    use precision
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: iBeg, jBeg, kBeg
    integer(kind=intType), intent(in) :: iEnd, jEnd, kEnd
    integer(kind=intType), dimension(3,2), intent(in) :: subRange

    real(kind=realType), dimension(*), intent(in) :: buffer
    real(kind=8), dimension(iBeg:iEnd,jBeg:jEnd,kBeg:kEnd), &
         intent(inout) :: val
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ll

    ! Copy the subrange into val.

    ll = 0
    do k=subRange(3,1), subRange(3,2)
       do j=subRange(2,1), subRange(2,2)
          do i=subRange(1,1), subRange(1,2)
             ll = ll + 1
             val(i,j,k) = real(buffer(ll), doubleType)
          enddo
       enddo
    enddo

  end subroutine copyDataBufDoublePrecision

  subroutine volSolNames(solNames)
    !
    !       volSolNames sets the names for the volume variables to be      
    !       written to the volume solution file. Sids convention names are 
    !       used as much as possible.                                      
    !
    use constants
    use cgnsNames
    use inputPhysics
    use flowVarRefState
    use extraOutput
    implicit none
    !
    !      Subroutine argument.
    !
    character(len=*), dimension(*), intent(out) :: solNames
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn
    !
    ! First store the names of the independent flow variables.

    solNames(1) = cgnsDensity
    solNames(2) = cgnsVelx
    solNames(3) = cgnsVely
    solNames(4) = cgnsVelz
    solNames(5) = cgnsPressure

    ! The turbulent variables if the RANS equations are solved.
    ! Note that these are the primitive variables and not the
    ! conservative ones. The reason is that the sids conventions only
    ! defines these names and not the conservative ones.

    if(equations == RANSEquations) then

       select case(turbModel)

       case(spalartAllmaras, spalartAllmarasEdwards)
          solNames(itu1) = cgnsTurbSaNu

       case(komegaWilcox, komegaModified, menterSST)
          solNames(itu1) = cgnsTurbK
          solNames(itu2) = cgnsTurbOmega

       case(ktau)
          solNames(itu1) = cgnsTurbK
          solNames(itu2) = cgnsTurbTau

       case(v2f)
          solNames(itu1) = cgnsTurbK
          solNames(itu2) = cgnsTurbEpsilon
          solNames(itu3) = cgnsTurbV2
          solNames(itu4) = cgnsTurbF

       end select

    endif

    ! Initialize nn to the number of independent variables.

    nn = nw

    ! Check the additional variables to be written.

    if( volWriteMx )  then
       nn = nn + 1
       solNames(nn) = cgnsMomx
    endif

    if( volWriteMy )  then
       nn = nn + 1
       solNames(nn) = cgnsMomy
    endif

    if( volWriteMz )  then
       nn = nn + 1
       solNames(nn) = cgnsMomz
    endif

    if( volWriteRVx )  then
       nn = nn + 1
       solNames(nn) = cgnsRelVelx
    endif

    if( volWriteRVy )  then
       nn = nn + 1
       solNames(nn) = cgnsRelVely
    endif

    if( volWriteRVz )  then
       nn = nn + 1
       solNames(nn) = cgnsRelVelz
    endif

    if( volWriteRhoe ) then
       nn = nn + 1
       solNames(nn) = cgnsEnergy
    endif

    if( volWriteTemp ) then
       nn = nn + 1
       solNames(nn) = cgnsTemp
    endif

    if( volWriteCp ) then
       nn = nn + 1
       solNames(nn) = cgnsCp
    endif

    if( volWriteMach ) then
       nn = nn + 1
       solNames(nn) = cgnsMach
    endif

    if( volWriteRMach ) then
       nn = nn + 1
       solNames(nn) = cgnsRelMach
    endif

    if( volWriteMachTurb ) then
       nn = nn + 1
       solNames(nn) = cgnsMachTurb
    endif

    if( volWriteEddyVis ) then
       nn = nn + 1
       solNames(nn) = cgnsEddy
    endif

    if( volWriteRatioEddyVis ) then
       nn = nn + 1
       solNames(nn) = cgnsEddyRatio
    endif

    if( volWriteDist ) then
       nn = nn + 1
       solNames(nn) = cgNSWallDist
    endif

    if( volWriteVort ) then
       nn = nn + 1
       solNames(nn) = cgnsVortMagn
    endif

    if( volWriteVortx ) then
       nn = nn + 1
       solNames(nn) = cgnsVortx
    endif

    if( volWriteVorty ) then
       nn = nn + 1
       solNames(nn) = cgnsVorty
    endif

    if( volWriteVortz ) then
       nn = nn + 1
       solNames(nn) = cgnsVortz
    endif

    if( volWritePtotloss ) then
       nn = nn + 1
       solNames(nn) = cgnsPtotloss
    endif

    if( volWriteResRho ) then
       nn = nn + 1
       solNames(nn) = cgnsResRho
    endif

    if( volWriteResMom ) then
       nn = nn + 1
       solNames(nn) = cgnsResMomx

       nn = nn + 1
       solNames(nn) = cgnsResMomy

       nn = nn + 1
       solNames(nn) = cgnsResMomz
    endif

    if( volWriteResRhoe) then
       nn = nn + 1
       solNames(nn) = cgnsResRhoe
    endif

    if( volWriteResTurb ) then

       select case(turbModel)

       case(spalartAllmaras, spalartAllmarasEdwards)
          nn = nn + 1
          solNames(nn) = cgnsResNu

       case(komegaWilcox, komegaModified, menterSST)
          nn = nn + 1
          solNames(nn) = cgnsResK

          nn = nn + 1
          solNames(nn) = cgnsResOmega

       case(ktau)
          nn = nn + 1
          solNames(nn) = cgnsResK

          nn = nn + 1
          solNames(nn) = cgnsResTau

       case(v2f)
          nn = nn + 1
          solNames(nn) = cgnsResK

          nn = nn + 1
          solNames(nn) = cgnsResEpsilon

          nn = nn + 1
          solNames(nn) = cgnsResV2

          nn = nn + 1
          solNames(nn) = cgnsResF

       end select

    endif

    if ( volWriteShock) then
       nn = nn + 1
       solNames(nn) = cgnsShock
    end if

    if ( volWriteFilteredShock) then
       nn = nn + 1
       solNames(nn) = cgnsFilteredShock
    end if

    if( volWriteBlank) then
       nn = nn + 1
       solNames(nn) = cgnsBlank
    endif

    if( volWriteGC) then
       nn = nn + 1
       solNames(nn) = cgnsGC
    endif

    if( volWriteStatus) then
       nn = nn + 1
       solNames(nn) = cgnsStatus
    endif

  end subroutine volSolNames

  subroutine surfSolNames(solNames)
    !
    !       surfSolNames sets the names for the surface variables to be    
    !       written to the surface solution file. Sids convention names    
    !       are used as much as possible.                                  
    !
    use cgnsNames
    use extraOutput
    implicit none
    !
    !      Subroutine argument.
    !
    character(len=*), dimension(*), intent(out) :: solNames
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    ! Initialize nn to 0.

    nn = 0

    ! Check which surfaces variables must be written and set
    ! solNames accordingly.

    if( surfWriteRho ) then
       nn = nn + 1
       solNames(nn) = cgnsDensity
    endif

    if( surfWriteP   ) then
       nn = nn + 1
       solNames(nn) = cgnsPressure
    endif

    if( surfWriteTemp ) then
       nn = nn + 1
       solNames(nn) = cgnsTemp
    endif

    if( surfWriteVx ) then
       nn = nn + 1
       solNames(nn) = cgnsVelx
    endif

    if( surfWriteVy ) then
       nn = nn + 1
       solNames(nn) = cgnsVely
    endif

    if( surfWriteVz ) then
       nn = nn + 1
       solNames(nn) = cgnsVelz
    endif

    if( surfWriteRVx ) then
       nn = nn + 1
       solNames(nn) = cgnsRelVelx
    endif

    if( surfWriteRVy ) then
       nn = nn + 1
       solNames(nn) = cgnsRelVely
    endif

    if( surfWriteRVz ) then
       nn = nn + 1
       solNames(nn) = cgnsRelVelz
    endif

    if( surfWriteCp ) then
       nn = nn + 1
       solNames(nn) = cgnsCp
    endif

    if( surfWritePtotloss ) then
       nn = nn + 1
       solNames(nn) = cgnsPtotloss
    endif

    if( surfWriteMach ) then
       nn = nn + 1
       solNames(nn) = cgnsMach
    endif

    if( surfWriteRMach ) then
       nn = nn + 1
       solNames(nn) = cgnsRelMach
    endif

    if( surfWriteCf ) then
       nn = nn + 1
       solNames(nn) = cgnsSkinFmag
    endif

    if( surfWriteCh ) then
       nn = nn + 1
       solNames(nn) = cgnsStanton
    endif

    if( surfWriteYplus ) then
       nn = nn + 1
       solNames(nn) = cgnsYplus
    endif

    if( surfWriteCfx )  then
       nn = nn + 1
       solNames(nn) = cgnsSkinFx
    endif

    if( surfWriteCfy ) then
       nn = nn + 1
       solNames(nn) = cgnsSkinFy
    endif

    if( surfWriteCfz ) then
       nn = nn + 1
       solNames(nn) = cgnsSkinFz
    endif

    if( surfWriteBlank ) then
       nn = nn + 1
       solNames(nn) = cgnsBlank
    endif

    if (surfWriteSepSensor) then
       nn = nn + 1
       solNames(nn) = cgnsSepSensor
    end if

    if (surfWriteCavitation) then
       nn = nn + 1
       solNames(nn) = cgnsCavitation
    end if

    if (surfWriteGC) then
       nn = nn + 1
       solNames(nn) = cgnsGC
    end if

  end subroutine surfSolNames

  subroutine storeSolInBuffer(buffer, copyInBuffer, solName, &
       iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    !
    !       StoreSolInBuffer stores the given range of the variable        
    !       indicated by solName in IOVar and copies it into buffer if     
    !       desired. It is assumed that the variables in blockPointers     
    !       already point to the correct block.                            
    !
    use constants
    use blockPointers
    use cgnsGrid
    use cgnsNames
    use flowVarRefState
    use inputPhysics
    use IOModule
    use flowUtils, only : computePTot
    use utils, only : terminate
    use overset, only : oversetPresent
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
    integer(kind=intType), intent(in) :: kBeg, kEnd

    real(kind=realType), dimension(*), intent(out) :: buffer
    character(len=*), intent(in)                   :: solName

    logical, intent(in) :: copyInBuffer
    !
    !      Local parameters
    !
    real(kind=realType), parameter :: plim   = 0.001_realType
    real(kind=realType), parameter :: rholim = 0.001_realType
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, jj, kk, nn

    real(kind=realType) :: uuy, uuz, vvx, vvz, wwx, wwy, tmp
    real(kind=realType) :: vortx, vorty, vortz, a2, ptotInf, ptot
    real(kind=realType) :: a, UovA(3), gradP(3)

    real(kind=realType), dimension(:,:,:,:), pointer :: wIO

    ! Set the pointer to the correct entry of IOVar. I'm cheating a
    ! bit here, because I know that only memory has been allocated
    ! for the first solution ID of IOVar.

    wIO => IOVar(nbkLocal,1)%w

    ! Determine the variable to be stored, compute it and store
    ! it in the 1D array buffer.

    select case(solName)

    case (cgnsDensity)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,irho)
             enddo
          enddo
       enddo

    case (cgnsMomx)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,irho)*w(i,j,k,ivx)
             enddo
          enddo
       enddo

    case (cgnsMomy)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,irho)*w(i,j,k,ivy)
             enddo
          enddo
       enddo

    case (cgnsMomz)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,irho)*w(i,j,k,ivz)
             enddo
          enddo
       enddo

    case (cgnsEnergy)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,irhoE)
             enddo
          enddo
       enddo

    case (cgnsTurbSaNu,cgnsTurbK)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,itu1)
             enddo
          enddo
       enddo

    case (cgnsTurbOmega,cgnsTurbTau,cgnsTurbEpsilon)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,itu2)
             enddo
          enddo
       enddo

    case (cgnsTurbV2)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,itu3)
             enddo
          enddo
       enddo

    case (cgnsTurbF)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,itu4)
             enddo
          enddo
       enddo

    case (cgnsVelx)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,ivx)
             enddo
          enddo
       enddo

    case (cgnsVely)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,ivy)
             enddo
          enddo
       enddo

    case (cgnsVelz)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,ivz)
             enddo
          enddo
       enddo

    case (cgnsRelVelx)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,ivx)-s(i,j,k,1)
             enddo
          enddo
       enddo

    case (cgnsRelVely)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,ivy)-s(i,j,k,2)
             enddo
          enddo
       enddo

    case (cgnsRelVelz)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = w(i,j,k,ivz)-s(i,j,k,3)
             enddo
          enddo
       enddo

    case (cgnsPressure)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = p(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsTemp)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = p(i,j,k)/(RGas*w(i,j,k,irho))
             enddo
          enddo
       enddo

    case (cgnsCp)
       tmp = two/(gammaInf*pInf*MachCoef*MachCoef)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = tmp*(p(i,j,k) - pInf)
             enddo
          enddo
       enddo

    case (cgnsMach)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                     / max(w(i,j,k,irho),rholim)
                tmp = (w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                     +  w(i,j,k,ivz)**2)/a2
                wIO(i,j,k,1) = sqrt(max(zero,tmp))
             enddo
          enddo
       enddo

    case (cgnsRelMach)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                     / max(w(i,j,k,irho),rholim)
                tmp = ((w(i,j,k,ivx)-s(i,j,k,1))**2 +&
                     (w(i,j,k,ivy)-s(i,j,k,2))**2 &
                     +(w(i,j,k,ivz)-s(i,j,k,3))**2)/a2
                wIO(i,j,k,1) = sqrt(max(zero,tmp))
             enddo
          enddo
       enddo


    case (cgnsMachTurb)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                tmp = w(i,j,k,irho)*w(i,j,k,itu1) &
                     / (gamma(i,j,k)*max(p(i,j,k),plim))
                wIO(i,j,k,1) = sqrt(max(zero,tmp))
             enddo
          enddo
       enddo

    case (cgnsEddy)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = rev(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsEddyRatio)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = rev(i,j,k)/rlv(i,j,k)
             enddo
          enddo
       enddo

    case (cgNSWallDist)
       do k=kBeg,kEnd
          kk = max(2_intType,k); kk = min(kl,kk)
          do j=jBeg,jEnd
             jj = max(2_intType,j); jj = min(jl,jj)
             do i=iBeg,iEnd
                ii = max(2_intType,i); ii = min(il,ii)
                wIO(i,j,k,1) = d2Wall(ii,jj,kk)
             enddo
          enddo
       enddo

    case (cgnsVortMagn)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                tmp = half/vol(i,j,k)
                uuy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                     - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                     + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                     - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                     + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                     - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                uuz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                     - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                     + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                     - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                     + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                     - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                vvx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                     - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                     + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                     - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                     + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                     - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                vvz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                     - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                     + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                     - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                     + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                     - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                wwx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                     - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                     + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                     - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                     + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                     - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                wwy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                     - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                     + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                     - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                     + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                     - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                vortx = wwy - vvz; vorty = uuz - wwx; vortz = vvx - uuy

                wIO(i,j,k,1) = tmp*sqrt(vortx**2 + vorty**2 + vortz**2)
             enddo
          enddo
       enddo

    case (cgnsVortx)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                tmp = half/vol(i,j,k)
                vvz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                     - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                     + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                     - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                     + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                     - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                wwy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                     - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                     + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                     - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                     + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                     - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                wIO(i,j,k,1) = tmp*(wwy - vvz)
             enddo
          enddo
       enddo

    case (cgnsVorty)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                tmp = half/vol(i,j,k)
                uuz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                     - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                     + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                     - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                     + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                     - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                wwx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                     - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                     + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                     - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                     + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                     - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                wIO(i,j,k,1) = tmp*(uuz - wwx)
             enddo
          enddo
       enddo

    case (cgnsVortz)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                tmp = half/vol(i,j,k)
                uuy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                     - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                     + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                     - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                     + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                     - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                vvx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                     - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                     + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                     - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                     + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                     - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                wIO(i,j,k,1) = tmp*(vvx - uuy)
             enddo
          enddo
       enddo

    case (cgnsPtotloss)

       ! Compute the free stream total pressure.

       call computePtot(rhoInf, uInf, zero, zero, &
            pInf, ptotInf)
       ptotInf = one/ptotInf

       ! Loop over the cell centers and compute the
       ! total pressure loss.

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                call computePtot(w(i,j,k,irho), w(i,j,k,ivx), &
                     w(i,j,k,ivy),  w(i,j,k,ivz), &
                     p(i,j,k),      ptot)

                wIO(i,j,k,1) = one - ptot*ptotInf
             enddo
          enddo
       enddo

    case (cgnsResRho)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,irho)
                wIO(i,j,k,1) = dw(i,j,k,irho)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResMomx)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,imx)
                wIO(i,j,k,1) = dw(i,j,k,imx)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResMomy)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,imy)
                wIO(i,j,k,1) = dw(i,j,k,imy)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResMomz)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,imz)
                wIO(i,j,k,1) = dw(i,j,k,imz)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResRhoE)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,irhoE)
                wIO(i,j,k,1) = dw(i,j,k,irhoE)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResNu,cgnsResK)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,itu1)
                wIO(i,j,k,1) = dw(i,j,k,itu1)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResOmega,cgnsResTau,cgnsResEpsilon)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,itu2)
                wIO(i,j,k,1) = dw(i,j,k,itu2)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResV2)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,itu3)
                wIO(i,j,k,1) = dw(i,j,k,itu3)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsResF)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                !   wIO(i,j,k,1) = dw(i,j,k,itu4)
                wIO(i,j,k,1) = dw(i,j,k,itu4)/vol(i,j,k)
             enddo
          enddo
       enddo

    case (cgnsBlank)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = real(iblank(i,j,k),realType)
             enddo
          enddo
       enddo

    case (cgnsGC)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                wIO(i,j,k,1) = real(globalcell(i,j,k),realType)
             enddo
          enddo
       enddo

    case (cgnsstatus)
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                if (oversetPresent) then 
                   wIO(i,j,k,1) = real(status(i, j, k))
                else
                   wIO(i,j,k,1) = 0
                end if
             enddo
          enddo
       enddo
       
    case (cgnsShock)

       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd

                ! Here we compute U/a <dot> grad P / ||grad P||
                ! Whre U is the velocity vector, a is the speed of
                ! sound and P is the pressure. 

                ! U / a
                a  = sqrt(gamma(i,j,k)*max(p(i,j,k),plim) &
                     / max(w(i,j,k,irho),rholim))
                UovA = (/w(i,j,k,ivx), w(i,j,k,ivy), w(i,j,k,ivz)/)/a

                ! grad P / ||grad P||

                gradP(1) = si(i,  j,k,1)*P(i+1,j,k) &
                     - si(i-1,j,k,1)*P(i-1,j,k) &
                     + sj(i,j,  k,1)*P(i,j+1,k) &
                     - sj(i,j-1,k,1)*P(i,j-1,k) &
                     + sk(i,j,k,  1)*P(i,j,k+1) &
                     - sk(i,j,k-1,1)*P(i,j,k-1)

                gradP(2) = si(i,  j,k,2)*P(i+1,j,k) &
                     - si(i-1,j,k,2)*P(i-1,j,k) &
                     + sj(i,j,  k,2)*P(i,j+1,k) &
                     - sj(i,j-1,k,2)*P(i,j-1,k) &
                     + sk(i,j,k,  2)*P(i,j,k+1) &
                     - sk(i,j,k-1,2)*P(i,j,k-1)

                gradP(3) = si(i,  j,k,3)*P(i+1,j,k) &
                     - si(i-1,j,k,3)*P(i-1,j,k) &
                     + sj(i,j,  k,3)*P(i,j+1,k) &
                     - sj(i,j-1,k,3)*P(i,j-1,k) &
                     + sk(i,j,k,  3)*P(i,j,k+1) &
                     - sk(i,j,k-1,3)*P(i,j,k-1)

                gradP = gradP / sqrt(gradP(1)**2 + gradP(2)**2 + gradP(3)**2 )
                ! Dot product
                wIO(i,j,k,1) = UovA(1)*gradP(1) + UovA(2)*gradP(2) + UovA(3)*gradP(3)
             end do
          end do
       end do


    case default
       call terminate("storeSolInBuffer", &
            "This should not happen")

    end select

    ! Copy the data in the 1D buffer, if desired.

    if( copyInBuffer ) then
       nn = 0
       do k=kBeg,kEnd
          do j=jBeg,jEnd
             do i=iBeg,iEnd
                nn = nn + 1
                buffer(nn) = wIO(i,j,k,1)
             enddo
          enddo
       enddo
    endif

  end subroutine storeSolInBuffer

  subroutine storeSurfsolInBuffer(sps, buffer, nn, blockID,   &
       faceID, cellRange, solName, &
       viscousSubface, useRindLayer)
    !
    !       storeSurfsolInBuffer stores the variable indicated by          
    !       solName of the given block ID in the buffer. As the solution   
    !       must be stored in the center of the boundary face the average  
    !       value of the first internal cell and its corresponding halo is 
    !       computed. The counter nn is updated in this routine. However   
    !       it is not initialized, because multiple contributions may be   
    !       stored in buffer.                                              
    !
    use blockPointers
    use cgnsNames
    use constants
    use flowVarRefState
    use inputPhysics
    use inputIO
    use communication 
    use utils, only : setPointers
    use flowUtils, only : computePtot
    use inputCostFunctions
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)    :: sps, blockID, faceID
    integer(kind=intType), intent(inout) :: nn
    integer(kind=intType), dimension(3,2), intent(in) :: cellRange
    real(kind=realType), dimension(*), intent(out) :: buffer
    character(len=*), intent(in) :: solName
    logical, intent(in) :: viscousSubface, useRindLayer
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, sigma
    integer(kind=intType) :: ii, jj, mm, iiMax, jjMax, offVis

    integer(kind=intType), dimension(2,2) :: rangeFace
    integer(kind=intType), dimension(3,2) :: rangeCell

    integer(kind=intType), dimension(:,:), pointer :: viscPointer
    integer(kind=intType), dimension(:,:), pointer :: iblank2

    real(kind=realType) :: fact, gm1, ptotInf, ptot, psurf, rsurf
    real(kind=realType) :: usurf, vsurf, wsurf, m2surf, musurf
    real(kind=realType) :: fx, fy, fz, fn, a2Tot, a2, qw
    real(kind=realType) :: tauxx, tauyy, tauzz
    real(kind=realType) :: tauxy, tauxz, tauyz
    real(kind=realType) :: pm1, a, sensor, plocal, sensor1
    real(kind=realType), dimension(3) :: norm, V

    real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
    real(kind=realType), dimension(:,:,:), pointer :: ss1, ss2, ss
    real(kind=realType), dimension(:,:),   pointer :: pp1, pp2

    real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
    real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
    real(kind=realType), dimension(:,:),   pointer :: dd2Wall

    ! Set the pointers to this block.

    call setPointers(blockID, 1_intType, sps)

    ! Set the offset for the viscous data, such that the range is
    ! limited to the actual physical face. Viscous data, like skin
    ! friction, need gradIent information, which is not available
    ! in the halo's.

    offVis = 0
    if(  useRindLayer ) offVis = 1

    ! CellRange contains the range of the current block in the
    ! original cgns block. Substract the offset and store the local
    ! range in rangeCell.

    rangeCell(1,1) = cellRange(1,1) - iBegor + 1
    rangeCell(1,2) = cellRange(1,2) - iBegor + 1

    rangeCell(2,1) = cellRange(2,1) - jBegor + 1
    rangeCell(2,2) = cellRange(2,2) - jBegor + 1

    rangeCell(3,1) = cellRange(3,1) - kBegor + 1
    rangeCell(3,2) = cellRange(3,2) - kBegor + 1
    !
    !                Viscous variables for a non-viscous wall.             
    !                Simply set the variables to zero and return.          
    !
    if(.not. viscousSubface) then

       select case (solName)

       case (cgnsSkinFmag, cgnsStanton, cgnsYplus, &
            cgnsSkinFx, cgnsSkinFy, cgnsSkinFz)

          ! Update the counter and set this entry of buffer to 0.

          do k=rangeCell(3,1),rangeCell(3,2)
             do j=rangeCell(2,1),rangeCell(2,2)
                do i=rangeCell(1,1),rangeCell(1,2)
                   nn = nn + 1
                   buffer(nn) = zero
                enddo
             enddo
          enddo

          ! Work has been done for this variable. So return.

          return

       end select
    endif
    !
    !       Determine the face on which the subface is located and set     
    !       a couple of variables accordingly. In this way a generic       
    !       treatment is possible and there is no need to repeat the code  
    !       for each of the six block faces.                               
    !       Note that for dd2Wall a slightly different notation must be    
    !       used. Reason is that d2Wall starts at index 2, rather than 0.  
    !

    select case (faceID)

    case (iMin)
       rangeFace(1,1:2) = rangeCell(2,1:2)
       rangeFace(2,1:2) = rangeCell(3,1:2)
       iiMax = jl; jjMax = kl

       ww1    => w(1,1:,1:,:);   ww2    => w(2,1:,1:,:)
       pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)
       ss => si(1,:,:,:) ; fact = -one

       pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)           
       gamma1 => gamma(1,1:,1:); gamma2 => gamma(2,1:,1:)

       if( blockIsMoving)then
          ss1    => s(1,1:,1:,:);   ss2    => s(2,1:,1:,:)
       end if

       iblank2 => iblank(2,1:,1:)
       viscPointer => viscIminPointer

       if( viscous ) then
          rlv1 => rlv(1,1:,1:); rlv2 => rlv(2,1:,1:)
       endif

       if(equations == RANSEquations) dd2Wall => d2Wall(2,:,:)

       !===============================================================

    case (iMax)
       rangeFace(1,1:2) = rangeCell(2,1:2)
       rangeFace(2,1:2) = rangeCell(3,1:2)
       iiMax = jl; jjMax = kl

       ww1    => w(ie,1:,1:,:);   ww2    => w(il,1:,1:,:)
       ss => si(il,:,:,:) ; fact = one

       pp1    => p(ie,1:,1:);     pp2    => p(il,1:,1:)
       gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)

       if( blockIsMoving)then
          ss1    => s(ie-1,1:,1:,:);   ss2    => s(ie,1:,1:,:)
       end if

       iblank2 => iblank(il,1:,1:)
       viscPointer => viscImaxPointer

       if( viscous ) then
          rlv1 => rlv(ie,1:,1:); rlv2 => rlv(il,1:,1:)
       endif

       if(equations == RANSEquations) dd2Wall => d2Wall(il,:,:)

       !===============================================================

    case (jMin)
       rangeFace(1,1:2) = rangeCell(1,1:2)
       rangeFace(2,1:2) = rangeCell(3,1:2)
       iiMax = il; jjMax = kl

       ww1    => w(1:,1,1:,:);   ww2    => w(1:,2,1:,:)
       ss => sj(:,1,:,:) ; fact = -one

       pp1    => p(1:,1,1:);     pp2    => p(1:,2,1:)
       gamma1 => gamma(1:,1,1:); gamma2 => gamma(1:,2,1:)

       if( blockIsMoving)then
          ss1    => s(1:,1,1:,:);   ss2    => s(1:,2,1:,:)
       end if

       iblank2 => iblank(1:,2,1:)
       viscPointer => viscJminPointer

       if( viscous ) then
          rlv1 => rlv(1:,1,1:); rlv2 => rlv(1:,2,1:)
       endif

       if(equations == RANSEquations) dd2Wall => d2Wall(:,2,:)

       !===============================================================

    case (jMax)
       rangeFace(1,1:2) = rangeCell(1,1:2)
       rangeFace(2,1:2) = rangeCell(3,1:2)
       iiMax = il; jjMax = kl

       ww1    => w(1:,je,1:,:);   ww2    => w(1:,jl,1:,:)
       ss => sj(:,jl,:,:); fact = one

       pp1    => p(1:,je,1:);     pp2    => p(1:,jl,1:)
       gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)

       if( blockIsMoving)then
          ss1    => s(1:,je-1,1:,:);   ss2    => s(1:,je,1:,:)
       end if

       iblank2 => iblank(1:,jl,1:)
       viscPointer => viscJmaxPointer

       if( viscous ) then
          rlv1 => rlv(1:,je,1:); rlv2 => rlv(1:,jl,1:)
       endif

       if(equations == RANSEquations) dd2Wall => d2Wall(:,jl,:)

       !===============================================================

    case (kMin)
       rangeFace(1,1:2) = rangeCell(1,1:2)
       rangeFace(2,1:2) = rangeCell(2,1:2)
       iiMax = il; jjMax = jl

       ww1    => w(1:,1:,1,:);   ww2    => w(1:,1:,2,:)
       ss => sk(:,:,1,:);  fact = -one

       pp1    => p(1:,1:,1);     pp2    => p(1:,1:,2)
       gamma1 => gamma(1:,1:,1); gamma2 => gamma(1:,1:,2)

       if( blockIsMoving)then
          ss1    => s(1:,1:,1,:);   ss2    => s(1:,1:,2,:)
       end if

       iblank2 => iblank(1:,1:,2)
       viscPointer => viscKminPointer

       if( viscous ) then
          rlv1 => rlv(1:,1:,1); rlv2 => rlv(1:,1:,2)
       endif

       if(equations == RANSEquations) dd2Wall => d2Wall(:,:,2)

       !===============================================================

    case (kMax)
       rangeFace(1,1:2) = rangeCell(1,1:2)
       rangeFace(2,1:2) = rangeCell(2,1:2)
       iiMax = il; jjMax = jl

       ww1    => w(1:,1:,ke,:);   ww2    => w(1:,1:,kl,:)
       ss => sk(:,:,kl,:);  fact = one

       pp1    => p(1:,1:,ke);     pp2    => p(1:,1:,kl)
       gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)

       if( blockIsMoving)then
          ss1    => s(1:,1:,ke-1,:);   ss2    => s(1:,1:,ke,:)
       end if

       iblank2 => iblank(1:,1:,kl)
       viscPointer => viscKmaxPointer

       if( viscous ) then
          rlv1 => rlv(1:,1:,ke); rlv2 => rlv(1:,1:,kl)
       endif

       if(equations == RANSEquations) dd2Wall => d2Wall(:,:,kl)

    end select
    !
    !       The actual part for storing the data. Determine the variable   
    !       to be written and loop over the boundary faces of the subface. 
    !
    ! Determine the variable to be written.

    varName: select case (solName)

    case (cgnsDensity)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             buffer(nn) = half*(ww1(i,j,irho) + ww2(i,j,irho))
          enddo
       enddo

       !===============================================================

    case (cgnsPressure)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             buffer(nn) = half*(pp1(i,j) + pp2(i,j))
          enddo
       enddo

       !===============================================================

    case (cgnsTemp)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             buffer(nn) = (pp1(i,j) + pp2(i,j)) &
                  / (RGas*(ww1(i,j,irho) + ww2(i,j,irho)))
          enddo
       enddo

       !===============================================================

    case (cgnsVelx)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             if (viscousSurfaceVelocities .and. viscous) then
                buffer(nn) = ww2(i,j,ivx)
             else
                buffer(nn) = half*(ww1(i,j,ivx) + ww2(i,j,ivx))
             end if
          enddo
       enddo

       !===============================================================

    case (cgnsVely)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             if (viscousSurfaceVelocities .and. viscous) then
                buffer(nn) = ww2(i,j,ivy)
             else
                buffer(nn) = half*(ww1(i,j,ivy) + ww2(i,j,ivy))
             end if
          enddo
       enddo

       !===============================================================

    case (cgnsVelz)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             if (viscousSurfaceVelocities .and. viscous) then
                buffer(nn) = ww2(i,j,ivz)
             else
                buffer(nn) = half*(ww1(i,j,ivz) + ww2(i,j,ivz))
             end if

          enddo
       enddo

       !===============================================================

    case (cgnsRelVelx)
       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             if (viscousSurfaceVelocities .and. viscous) then
                buffer(nn) = ww2(i,j,ivx) - ss2(i,j,1)
             else
                buffer(nn) = half*(ww1(i,j,ivx) + ww2(i,j,ivx))-half*(ss1(i,j,1) + ss2(i,j,1))
             end if
          enddo
       enddo

       !===============================================================

    case (cgnsRelVely)
       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             if (viscousSurfaceVelocities .and. viscous) then
                buffer(nn) = ww2(i,j,ivy) - ss2(i,j,2)
             else
                buffer(nn) = half*(ww1(i,j,ivy) + ww2(i,j,ivy))-half*(ss1(i,j,2) + ss2(i,j,2))
             end if
          enddo
       enddo

       !===============================================================

    case (cgnsRelVelz)
       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             if (viscousSurfaceVelocities .and. viscous) then
                buffer(nn) = ww2(i,j,ivz) - ss2(i,j,3)
             else
                buffer(nn) = half*(ww1(i,j,ivz) + ww2(i,j,ivz))-half*(ss1(i,j,3) + ss2(i,j,3))
             end if
          enddo
       enddo


       !================================================================

    case (cgnsCp)

       ! Factor multiplying p-pInf

       fact = two/(gammaInf*pInf*MachCoef*MachCoef)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             buffer(nn) = fact*(half*(pp1(i,j) + pp2(i,j)) - pInf)
          enddo
       enddo

       !===============================================================

    case (cgnsPtotloss)

       ! First compute the total pressure of the free stream.

       call computePtot(rhoInf, uInf, zero, zero, &
            pInf, ptotInf)
       ptotInf = one/ptotInf

       ! Loop over the faces and compute the total pressure loss.

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)

             psurf = half*(pp1(i,j) + pp2(i,j))
             rsurf = half*(ww1(i,j,irho) + ww2(i,j,irho))
             usurf = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
             vsurf = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
             wsurf = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))

             call computePtot(rsurf, usurf, vsurf, wsurf, &
                  psurf, ptot)

             nn = nn + 1
             buffer(nn) = one - ptot*ptotInf
          enddo
       enddo

       !===============================================================

    case (cgnsMach)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)

             psurf  = half*(pp1(i,j) + pp2(i,j))
             rsurf  = half*(ww1(i,j,irho) + ww2(i,j,irho))
             usurf  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
             vsurf  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
             wsurf  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))
             m2surf = rsurf*(usurf**2 + vsurf**2 + wsurf**2) &
                  / (half*(gamma1(i,j) + gamma2(i,j))*psurf)

             nn = nn + 1
             buffer(nn) = sqrt(m2surf)
          enddo
       enddo


       !===============================================================

    case (cgnsRelMach)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)

             psurf  = half*(pp1(i,j) + pp2(i,j))
             rsurf  = half*(ww1(i,j,irho) + ww2(i,j,irho))
             usurf  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))-half*(ss1(i,j,1) + ss2(i,j,1))
             vsurf  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))-half*(ss1(i,j,2) + ss2(i,j,2))
             wsurf  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))-half*(ss1(i,j,3) + ss2(i,j,3))
             m2surf = rsurf*(usurf**2 + vsurf**2 + wsurf**2) &
                  / (half*(gamma1(i,j) + gamma2(i,j))*psurf)

             nn = nn + 1
             buffer(nn) = sqrt(m2surf)
          enddo
       enddo

       !        ================================================================

    case (cgnsSkinFmag, cgnsYplus, &
         cgnsSkinFx, cgnsSkinFy, cgnsSkinFz)

       ! To avoid a lot of code duplication these 5 variables are
       ! treated together.

       ! Multiplication factor to obtain the skin friction from
       ! the wall shear stress.

       fact = two/(gammaInf*pInf*MachCoef*MachCoef)

       ! Loop over the given range of faces. As the viscous data is
       ! only present in the owned faces, the values of the halo's
       ! are set equal to the nearest physical face. Therefore the
       ! working indices are ii and jj.

       do j=rangeFace(2,1), rangeFace(2,2)
          if(j == rangeFace(2,1)) then
             jj = min(j + offVis, rangeFace(2,2))
          else if(j == rangeFace(2,2)) then
             jj = max(j - offVis, rangeFace(2, 1))
          else
             jj = j
          endif

          do i=rangeFace(1,1), rangeFace(1,2)
             if(i == rangeFace(1,1)) then
                ii = min(i + offVis, rangeFace(1,2))
             else if(i == rangeFace(1,2)) then
                ii = max(i - offVis, rangeFace(1,1))
             else
                ii = i
             endif

             ! Determine the viscous subface on which this
             ! face is located.

             mm = viscPointer(ii,jj)

             ! Store the 6 components of the viscous stress tensor
             ! a bit easier.

             tauxx = viscSubface(mm)%tau(ii,jj,1)
             tauyy = viscSubface(mm)%tau(ii,jj,2)
             tauzz = viscSubface(mm)%tau(ii,jj,3)
             tauxy = viscSubface(mm)%tau(ii,jj,4)
             tauxz = viscSubface(mm)%tau(ii,jj,5)
             tauyz = viscSubface(mm)%tau(ii,jj,6)

             ! Compute the "unit" force on this face. The unit normal
             ! is outward pointing per definition. A minus sign is
             ! present, because of the definition of the viscous
             ! stress tensor. Note that in the normal the indices i
             ! and j could be used. However this is not done.

             norm(1) = BCData(mm)%norm(ii,jj,1)
             norm(2) = BCData(mm)%norm(ii,jj,2)
             norm(3) = BCData(mm)%norm(ii,jj,3)

             fx = -(tauxx*norm(1) + tauxy*norm(2) + tauxz*norm(3))
             fy = -(tauxy*norm(1) + tauyy*norm(2) + tauyz*norm(3))
             fz = -(tauxz*norm(1) + tauyz*norm(2) + tauzz*norm(3))

             fn = fx*norm(1) + fy*norm(2) + fz*norm(3)

             fx = fx - fn*norm(1)
             fy = fy - fn*norm(2)
             fz = fz - fn*norm(3)

             ! Determine the variable to be stored and compute it.
             ! Note that an offset of -1 must be used in dd2Wall,
             ! because the original array, d2Wall, starts at 2.
             ! First update the counter nn.

             nn = nn + 1

             select case (solName)
             case (cgnsSkinFmag)
                buffer(nn) = fact*sqrt(fx*fx + fy*fy + fz*fz)

             case (cgnsSkinFx)
                buffer(nn) = fact*fx

             case (cgnsSkinFy)
                buffer(nn) = fact*fy

             case (cgnsSkinFz)
                buffer(nn) = fact*fz

             case (cgnsYplus)
                rsurf      = half*(ww1(ii,jj,irho) + ww2(ii,jj,irho))
                musurf     = half*(rlv1(ii,jj)     + rlv2(ii,jj))
                buffer(nn) = sqrt(rsurf*sqrt(fx*fx + fy*fy + fz*fz)) &
                     * dd2Wall(ii-1,jj-1)/musurf
             end select

          enddo
       enddo

       !        ================================================================

    case (cgnsStanton)

       ! Some constants needed to compute the stanton number.

       gm1   = gammaInf - one
       a2Tot = gammaInf*pInf*(one + half*gm1*MachCoef*MachCoef) &
            / rhoInf
       fact   = MachCoef*sqrt(gammaInf*pInf*rhoInf)/gm1

       ! Loop over the given range of faces. As the viscous data is
       ! only present in the owned faces, the values of the halo's
       ! are set equal to the nearest physical face. Therefore the
       ! working indices are ii and jj.
       do j=rangeFace(2,1), rangeFace(2,2)
          if(j == rangeFace(2,1)) then
             jj = min(j + offVis, rangeFace(2,2))
          else if(j == rangeFace(2,2)) then
             jj = max(j - offVis, rangeFace(2, 1))
          else
             jj = j
          endif

          do i=rangeFace(1,1), rangeFace(1,2)
             if(i == rangeFace(1,1)) then
                ii = min(i + offVis, rangeFace(1,2))
             else if(i == rangeFace(1,2)) then
                ii = max(i - offVis, rangeFace(1,1))
             else
                ii = i
             endif

             ! Determine the viscous subface on which this
             ! face is located.

             mm = viscPointer(ii,jj)

             ! Compute the heat flux. Multipy with the sign of the
             ! normal to obtain the correct value.

             qw = viscSubface(mm)%q(ii,jj,1)*BCData(mm)%norm(ii,jj,1) &
                  + viscSubface(mm)%q(ii,jj,2)*BCData(mm)%norm(ii,jj,2) &
                  + viscSubface(mm)%q(ii,jj,3)*BCData(mm)%norm(ii,jj,3)

             ! Compute the speed of sound squared at the wall and
             ! the stanton number, which is stored in buffer.

             a2 = half*(gamma1(ii,jj)   + gamma2(ii,jj)) &
                  *      (pp1(ii,jj)      + pp2(ii,jj))    &
                  /      (ww1(ii,jj,irho) + ww2(ii,jj,irho))

             nn = nn + 1
             buffer(nn) = qw/(fact*(a2Tot-a2))

          enddo
       enddo

       !        ================================================================

    case (cgnsBlank)

       ! Loop over the given range of faces. Since iblanks are set
       ! to 2 for boundary conditions and >= 10 for the boundary,
       ! take the minimum of the value and 1, so that cells with
       ! valid data always have an iblank of 1.

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1
             buffer(nn) = real(min(iblank2(i,j), 1_intType), realType)
          enddo
       enddo

    case (cgnsSepSensor)

       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)
             nn = nn + 1

             ! Get normalized surface velocity:
             v(1) = ww2(i, j, ivx)
             v(2) = ww2(i, j, ivy)
             v(3) = ww2(i, j, ivz)

             ! Normalize
             v = v / (sqrt(v(1)**2 + v(2)**2 + v(3)**2) + 1e-16)

             ! Dot product with free stream
             sensor = -dot_product(v, velDirFreeStream)

             !Now run through a smooth heaviside function:
             sensor = one/(one + exp(-2*sepSensorSharpness*(sensor - sepSensorOffset)))
             buffer(nn) = sensor
          enddo
       enddo

    case (cgnsCavitation)
       fact = two/(gammaInf*pInf*MachCoef*MachCoef)
       do j=rangeFace(2,1), rangeFace(2,2)
          do i=rangeFace(1,1), rangeFace(1,2)

             nn = nn + 1
             ! Get local pressure
             plocal = half*(pp1(i,j) + pp2(i,j))

             sigma = 1.4 
             sensor1 = (-(fact)*(plocal-pInf))- sigma
             sensor1 = one/(one + exp(-2*10*sensor1))
             buffer(nn) = sensor1
             !print*, sensor
          enddo
       enddo
    end select varName

  end subroutine storeSurfsolInBuffer

  subroutine storeOldSolInBuffer(buffer, ind, wID, &
       iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    !
    !       storeOldSolInBuffer stores the given range of the wID'th       
    !       conservative variable of an old solution in buffer. Needed for 
    !       a time accurate restart. It is assumed that the variables in   
    !       blockPointers already point to the correct block.              
    !
    use blockPointers
    use constants
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: ind, wID
    integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
    integer(kind=intType), intent(in) :: kBeg, kEnd

    real(kind=realType), dimension(*), intent(out) :: buffer
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, nOld, nn

    ! Store the index in wOld a bit easier.

    nOld = ind - 1

    ! Loop over the cell range of the block and copy the wID'th
    ! variable in buffer.

    nn = 0

    do k=kBeg,kEnd
       do j=jBeg,jEnd
          do i=iBeg,iEnd
             nn = nn + 1; buffer(nn) = wOld(nOld,i,j,k,wID)
          enddo
       enddo
    enddo

  end subroutine storeOldSolInBuffer

  subroutine describeScheme(string)
    !
    !       describeScheme gives a short description about the scheme      
    !       used to obtain the solution. The description is stored in the  
    !       character array string.                                        
    !
    use constants
    use inputDiscretization
    use inputPhysics
    use flowVarRefState
    implicit none
    !
    !      Subroutine arguments.
    !
    character(len=*), intent(out) :: string

    ! Write the basic scheme info.

    select case(spaceDiscr)
    case (dissScalar)
       write(string,100) "Scalar dissipation scheme", vis2, vis4
    case (dissMatrix)
       write(string,100) "Matrix dissipation scheme", vis2, vis4
    case (dissCusp)
       write(string,100) "CUSP dissipation scheme", vis2, vis4
    case (upwind)
       select case (limiter)
       case (firstOrder)
          write(string,110) "First order upwind scheme."
       case (noLimiter)
          write(string,111) kappaCoef
       case (vanAlbeda)
          write(string,112) kappaCoef
       case (minmod)
          write(string,113) kappaCoef
       end select

       select case (riemann)
       case (Roe)
          write(string,130) trim(string), &
               "Roe's approximate Riemann Solver."
       case (vanLeer)
          write(string,130) trim(string), &
               "Van Leer flux vector splitting."
       case (ausmdv)
          write(string,130) trim(string), &
               "ausmdv flux vector splitting."
       end select

    end select

100 format(1X, A,", k2 = ", e12.5, ", k4 = ", e12.5,".")
110 format(1X, A)
111 format(1X, "Second order upwind scheme using linear reconstruction, &
         &i.e. no limiter, kappa =", 1X, f7.3,".")
112 format(1X, "Second order upwind scheme with Van Albeda limiter, &
         &kappa =", 1X, f7.3,".")
113 format(1X, "Second order upwind scheme with Minmod limiter, &
         &kappa =", 1X, f7.3,".")
130 format(1X, A, 1X, A)

    ! In case of the scalar dissipation scheme, write whether or not
    ! directional scaling has been applied.

    if(spaceDiscr == dissScalar) then
       if( dirScaling ) then
          write(string,200) trim(string), adis
       else
          write(string,210) trim(string)
       endif
    endif

200 format(1X, A, 1X, "Directional scaling of dissipation with exponent", &
         1X,e12.5, ".")
210 format(1X, A, 1X, "No directional scaling of dissipation.")

    ! For the Euler equations, write the inviscid wall boundary
    ! condition treatment.

    if(equations == EulerEquations) then
       select case (eulerWallBcTreatment)
       case (constantPressure)
          write(string,300) trim(string), &
               "Zero normal pressure gradIent"
       case (linExtrapolPressure)
          write(string,300) trim(string), &
               "Linear extrapolation of normal &
               &pressure gradIent"
       case (quadExtrapolPressure)
          write(string,300) trim(string), &
               "Quadratic extrapolation of normal &
               &pressure gradIent"
       case (normalMomentum)
          write(string,300) trim(string), &
               "Normal momentum equation used to &
               &determine pressure gradient"
       end select
    endif
300 format(1X, A, 1X, A, 1X,"for inviscid wall boundary conditions.")

    ! If preconditioning is used, write the preconditioner.

    select case(precond)
    case (Turkel)
       write(string,400) trim(string), &
            "Turkel preconditioner for inviscid fluxes."
    case (ChoiMerkle)
       write(string,400) trim(string), &
            "Choi Merkle preconditioner for inviscid fluxes."
    end select
400 format(1X, A, 1X, A)

    ! For a viscous computation write that a central discretization
    ! is used for the viscous fluxes.

    if( viscous ) then
       write(string,400) trim(string), &
            "Central discretization for viscous fluxes."
    endif

  end subroutine describeScheme

  subroutine isoSurfNames(solNames)
    !
    !       isoNames sets the names for the volume variables to be         
    !       written to the isosurfaces. Sids convention names are          
    !       used as much as possible.                                      
    !
    use constants
    use cgnsNames
    use inputPhysics
    use flowVarRefState
    use extraOutput
    implicit none
    !
    !      Subroutine argument.
    !
    character(len=*), dimension(*), intent(out) :: solNames
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn


    ! Check the additional variables to be written -- there are no
    ! default variables already written
    nn = 0
    if (isoWriteRho) then
       nn = nn + 1
       solNames(nn) = cgnsDensity
    end if

    if (isoWriteVx) then
       nn = nn + 1
       solNames(nn) = cgnsVelx
    end if

    if (isoWriteVy) then
       nn = nn + 1
       solNames(nn) = cgnsVely
    end if

    if (isoWriteVz) then
       nn = nn + 1
       solNames(nn) = cgnsVelz
    end if

    if (isoWriteP) then
       nn = nn + 1
       solNames(nn) = cgnsPressure
    end if

    if( isoWriteTurb ) then

       select case(turbModel)

       case(spalartAllmaras, spalartAllmarasEdwards)
          nn = nn + 1
          solNames(nn) = cgnsTurbSaNu

       case(komegaWilcox, komegaModified, menterSST)
          nn = nn + 1
          solNames(nn) = cgnsTurbK
          nn = nn + 1
          solNames(nn) = cgnsTurbOmega

       case(ktau)
          nn = nn + 1
          solNames(nn) = cgnsTurbK
          nn = nn + 1
          solNames(nn) = cgnsTurbTau

       case(v2f)
          nn = nn + 1
          solNames(nn) = cgnsTurbK
          nn = nn + 1
          solNames(nn) = cgnsTurbEpsilon
          nn = nn + 1
          solNames(nn) = cgnsTurbV2
          nn = nn + 1
          solNames(nn) = cgnsTurbF

       end select

    endif

    if( isoWriteMx )  then
       nn = nn + 1
       solNames(nn) = cgnsMomx
    endif

    if( isoWriteMy )  then
       nn = nn + 1
       solNames(nn) = cgnsMomy
    endif

    if( isoWriteMz )  then
       nn = nn + 1
       solNames(nn) = cgnsMomz
    endif

    if( isoWriteRVx )  then
       nn = nn + 1
       solNames(nn) = cgnsRelVelx
    endif

    if( isoWriteRVy )  then
       nn = nn + 1
       solNames(nn) = cgnsRelVely
    endif

    if( isoWriteRVz )  then
       nn = nn + 1
       solNames(nn) = cgnsRelVelz
    endif

    if( isoWriteRhoe ) then
       nn = nn + 1
       solNames(nn) = cgnsEnergy
    endif

    if( isoWriteTemp ) then
       nn = nn + 1
       solNames(nn) = cgnsTemp
    endif

    if( isoWriteCp ) then
       nn = nn + 1
       solNames(nn) = cgnsCp
    endif

    if( isoWriteMach ) then
       nn = nn + 1
       solNames(nn) = cgnsMach
    endif

    if( isoWriteRMach ) then
       nn = nn + 1
       solNames(nn) = cgnsRelMach
    endif

    if( isoWriteMachTurb ) then
       nn = nn + 1
       solNames(nn) = cgnsMachTurb
    endif

    if( isoWriteEddyVis ) then
       nn = nn + 1
       solNames(nn) = cgnsEddy
    endif

    if( isoWriteRatioEddyVis ) then
       nn = nn + 1
       solNames(nn) = cgnsEddyRatio
    endif

    if( isoWriteDist ) then
       nn = nn + 1
       solNames(nn) = cgNSWallDist
    endif

    if( isoWriteVort ) then
       nn = nn + 1
       solNames(nn) = cgnsVortMagn
    endif

    if( isoWriteVortx ) then
       nn = nn + 1
       solNames(nn) = cgnsVortx
    endif

    if( isoWriteVorty ) then
       nn = nn + 1
       solNames(nn) = cgnsVorty
    endif

    if( isoWriteVortz ) then
       nn = nn + 1
       solNames(nn) = cgnsVortz
    endif

    if( isoWritePtotloss ) then
       nn = nn + 1
       solNames(nn) = cgnsPtotloss
    endif

    if( isoWriteResRho ) then
       nn = nn + 1
       solNames(nn) = cgnsResRho
    endif

    if( isoWriteResMom ) then
       nn = nn + 1
       solNames(nn) = cgnsResMomx

       nn = nn + 1
       solNames(nn) = cgnsResMomy

       nn = nn + 1
       solNames(nn) = cgnsResMomz
    endif

    if( isoWriteResRhoe) then
       nn = nn + 1
       solNames(nn) = cgnsResRhoe
    endif

    if( isoWriteResTurb ) then

       select case(turbModel)

       case(spalartAllmaras, spalartAllmarasEdwards)
          nn = nn + 1
          solNames(nn) = cgnsResNu

       case(komegaWilcox, komegaModified, menterSST)
          nn = nn + 1
          solNames(nn) = cgnsResK

          nn = nn + 1
          solNames(nn) = cgnsResOmega

       case(ktau)
          nn = nn + 1
          solNames(nn) = cgnsResK

          nn = nn + 1
          solNames(nn) = cgnsResTau

       case(v2f)
          nn = nn + 1
          solNames(nn) = cgnsResK

          nn = nn + 1
          solNames(nn) = cgnsResEpsilon

          nn = nn + 1
          solNames(nn) = cgnsResV2

          nn = nn + 1
          solNames(nn) = cgnsResF

       end select

    endif

    if (isoWriteShock) then
       nn = nn + 1
       solNames(nn) = cgnsShock
    end if

    if (isoWriteFilteredShock) then
       nn = nn + 1
       solNames(nn) = cgnsFilteredShock
    end if

    if( isoWriteBlank) then
       nn = nn + 1
       solNames(nn) = cgnsBlank
    endif

  end subroutine isoSurfNames

  subroutine setHelpVariablesWriting
    !
    !       setHelpVariablesWriting determines the variables, which are    
    !       needed to write the CGNS files.                                
    !
    use block
    use cgnsGrid
    use communication
    use monitor
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr, nSend
    integer, dimension(nProc) :: recvCounts, displs

    integer(kind=intType) :: i, nn

    integer(kind=intType), dimension(cgnsNDom) :: tmp
    integer(kind=intType), dimension(4,nDom)   :: buffer

    ! Determine for each CGNS block how many (sub) blocks are stored
    ! on this processor. Note that this info is the same for all
    ! spectral solutions, so the 1st is fine.

    allocate(nBlocksCGNSblock(0:cgnsNDom), blocksCGNSblock(nDom), &
         stat=ierr)
    if(ierr /= 0)                               &
         call terminate("setHelpVariablesWriting", &
         "Memory allocation failure for &
         &nBlocksCGNSblock and blocksCGNSblock.")

    nBlocksCGNSblock = 0
    do nn=1,nDom
       i = flowDoms(nn,1,1)%cgnsBlockID
       nBlocksCGNSblock(i) = nBlocksCGNSblock(i) + 1
    enddo

    ! Put nBlocksCGNSblock in cumulative storage format.
    ! Store this accumulated value in tmp, which serves as
    ! a counter later on.

    do i=1,cgnsNDom
       tmp(i)              = nBlocksCGNSblock(i-1)
       nBlocksCGNSblock(i) = nBlocksCGNSblock(i) + tmp(i)
    enddo

    ! Determine the values for blocksCGNSblock.

    do nn=1,nDom
       i = flowDoms(nn,1,1)%cgnsBlockID
       tmp(i) = tmp(i) + 1
       blocksCGNSblock(tmp(i)) = nn
    enddo

  end subroutine setHelpVariablesWriting

  subroutine releaseHelpVariablesWriting
    !
    !       releaseHelpVariablesWriting releases the memory of the         
    !       variables, which were needed to write the CGNS files.          
    !
    use cgnsGrid
    use monitor
    use utils, only : terminate
    implicit none
    !
    !      Local variables
    !
    integer :: ierr

    ! Release the memory of the allocatable arrays in outputMod.

    deallocate(nBlocksCGNSblock, blocksCGNSblock, stat=ierr)
    if(ierr /= 0)                                   &
         call terminate("releaseHelpVariablesWriting", &
         "Deallocation failure for nBlocksCGNSblock, &
         &etc.")

  end subroutine releaseHelpVariablesWriting
  subroutine writeCGNSHeader(cgnsInd, base)
    !
    !       writeCGNSHeader writes a descriptive header to the given base  
    !       of the given CGNS file. Only processor 0 performs this task.   
    !
    use constants
    use cgnsGrid
    use cgnsNames
    use flowVarRefState
    use su_cgns
    use inputPhysics
    use inputTimeSpectral
    use monitor
    use utils, only : terminate, setCGNSRealType

    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in) :: cgnsInd, base
    !
    !      Local variables.
    !
    integer :: ierr, realTypeCGNS

    real(kind=cgnsRealType) :: val

    character(len=2048) :: message
    character(len=7)    :: integerString
    character(len=12)   :: realString

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Go to the correct position in the CGNS file.

    call cg_goto_f(cgnsInd, base, ierr, "end")
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling cg_goto_f")

    ! Create a data class type node to indicate that nonDimensional
    ! solution data is written for which the reference state
    ! is known.

    call cg_dataclass_write_f(NormalizedByDimensional,ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling &
         &cg_dataclass_write_f")

    ! Write the info about the solver used.

    call cg_descriptor_write_f("SolverInfo", &
         "ADflow multiblock code", ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling &
         &cg_descriptor_write_f")

    ! Write the info about the scheme used; message is used as
    ! storage for the string containing the scheme description.

    call describeScheme(message)
    call cg_descriptor_write_f("DiscretizationScheme", message, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling &
         &cg_descriptor_write_f")

    ! Write the similation type to the CGNS file.

    select case (equationMode)

    case (steady)

       ! Steady mode. Just write this info.

       call cg_simulation_type_write_f(cgnsInd, base, &
            nonTimeaccurate, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_simulation_type_write_f")

       !===============================================================

    case (unsteady)

       ! Unsteady mode. First write the simulation type.

       call cg_simulation_type_write_f(cgnsInd, base, &
            timeaccurate, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_simulation_type_write_f")

       ! Write some additional stuff, like time step and
       ! physical time. First store it in the big string message.

       write(integerString,"(i7)") timeStepUnsteady + &
            nTimeStepsRestart
       write(realString,"(e12.5)") timeUnsteady + &
            timeUnsteadyRestart

       integerString = adjustl(integerString)
       realString    = adjustl(realString)

       write(message,101) trim(integerString), trim(realString)
101    format("Unsteady time step ",a,", physical time ",a, &
            " seconds")

       ! And write the info.

       call cg_descriptor_write_f("UnsteadyInfo", message, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_descriptor_write_f")

       !===============================================================

    case (timeSpectral)

       ! Time spectral mode. This is not a predefined mode in CGNS
       ! and therefore use userDefined.

       call cg_simulation_type_write_f(cgnsInd, base, &
            UserDefined, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_simulation_type_write_f")

       ! Write some info to the string message.

       write(integerString,"(i7)") nTimeIntervalsSpectral
       integerString = adjustl(integerString)

       write(message,102) trim(integerString)
102    format("Time spectral mode for periodic problems; ",a,  &
            " spectral solutions have been used to model the &
            &problem.")

       ! And write the info.

       call cg_descriptor_write_f("PeriodicInfo", message, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_descriptor_write_f")
    end select

    ! Go back to the given base in the cgns file.

    call cg_goto_f(cgnsInd, base, ierr, "end")
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling cg_goto_f")

    ! Create a flow equation set.

    call cg_equationset_write_f(cgnsPhysDim, ierr)
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling &
         &cg_equationset_write_f")

    ! Write the rest of the physical model under the flow
    ! equation set just created.

    call cg_goto_f(cgnsInd, base, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling cg_goto_f")

    ! Write the governing equations solved.

    select case (equations)
    case (EulerEquations)
       call cg_governing_write_f(Euler, ierr)

    case (NSEquations)
       call cg_governing_write_f(nsLaminar, ierr)

    case (RANSEquations)
       call cg_governing_write_f(nsTurbulent, ierr)
    end select

    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSHeader", &
         "Something wrong when calling &
         &cg_governing_write_f")

    ! Write the information about the gas model used.
    ! Determine the cp model used in the computation.

    select case (cpModel)

    case (cpConstant)

       ! Constant cp and thus constant gamma.

       call cg_model_write_f("GasModel_t", Ideal, ierr)
       if(ierr /= CG_OK)                &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_model_write_f")

       ! Write the actual value of gamma; this must be done under
       ! gas model type, which explains the goto statement.

       call cg_goto_f(cgnsInd, base, ierr, "FlowEquationSet_t", &
            1, "GasModel_t", 1, "end")
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling cg_goto_f")

       val = gammaConstant
       call cg_array_write_f(cgnsHeatRatio, realTypeCGNS, &
            1, 1, val, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_array_write_f")

       ! And create a data class under SpecificHeatRatio to tell that
       ! this is a nonDimensional parameter.

       call cg_goto_f(cgnsInd, base, ierr,       &
            "FlowEquationSet_t", 1,    &
            "GasModel_t", 1, "DataArray_t", 1,"end")
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling cg_goto_f")

       call cg_dataclass_write_f(NonDimensionalParameter,ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_dataclass_write_f")

       !===============================================================

    case (cpTempCurveFits)

       ! Cp as function of the temperature is given via curve fits.

       call cg_model_write_f("GasModel_t", ThermallyPerfect, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_model_write_f")

    end select

    ! The rest of physical model description is only
    ! for viscous flows.

    viscousTest: if( viscous ) then

       ! Write the info of the viscosity model. Under the flow
       ! equation set.

       call cg_goto_f(cgnsInd, base, ierr, &
            "FlowEquationSet_t", 1, "end")
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling cg_goto_f")

       call cg_model_write_f("ViscosityModel_t", sutherlandlaw, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_model_write_f")

       ! Write the info about the thermal conductivity, i.e.
       ! Constant Prandtl number. Write the used value as well.

       call cg_model_write_f("ThermalConductivityModel_t", &
            constantPrandtl, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_model_write_f")

       call cg_goto_f(cgnsInd, base, ierr, "FlowEquationSet_t", 1,&
            "ThermalConductivityModel_t", 1, "end")
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling cg_goto_f")

       val = prandtl
       call cg_array_write_f(cgnsPrandtl, realTypeCGNS, 1, 1, &
            val, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_array_write_f")

       ! And create a data class under Prandtl number to tell that
       ! this is a nonDimensional parameter.

       call cg_goto_f(cgnsInd, base, ierr, "FlowEquationSet_t", &
            1, "ThermalConductivityModel_t", 1,       &
            "DataArray_t", 1,"end")
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling cg_goto_f")

       call cg_dataclass_write_f(NonDimensionalParameter,ierr)
       if(ierr /= CG_OK)                  &
            call terminate("writeCGNSHeader", &
            "Something wrong when calling &
            &cg_dataclass_write_f")

       ! The rest of the physical model description is only for the
       ! RANS equations.

       turbulentTest: if(equations == RANSEquations) then

          select case(turbModel)

          case (spalartAllmaras)
             call writeCGNSSaInfo(cgnsInd, base)

          case (spalartAllmarasEdwards)
             call writeCGNSSaeInfo(cgnsInd, base)

          case (komegaWilcox)
             call writeCGNSKomegaWilcoxInfo(cgnsInd, base)

          case (komegaModified)
             call writeCGNSKomegaModifiedInfo(cgnsInd, base)

          case (ktau)
             call writeCGNSKtauInfo(cgnsInd, base)

          case (menterSST)
             call writeCGNSMenterSSTInfo(cgnsInd, base)

          case (v2f)
             call writeCGNSV2fInfo(cgnsInd, base)
          end select

       endif turbulentTest

    endif viscousTest

    ! Write the reference state.

    call writeCGNSReferenceState(cgnsInd, base)
  end subroutine writeCGNSHeader

  subroutine writeCGNSKomegaModifiedInfo(cgnsInd, cgnsBase)
    !
    !       writeCGNSKomegaModifiedInfo writes information about the       
    !       modified k-omega turbulence model to the cgns file.            
    !
    use inputPhysics
    use cgnsNames
    use su_cgns
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, ierr

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Write the info of the turbulence model under the flow equation
    ! set. So move to this location first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                              &
         call terminate("writeCGNSKomegaModifiedInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write that the k-omega model is used.

    call cg_model_write_f("TurbulenceModel_t", &
         TwoEquation_Wilcox, ierr)
    if(ierr /= CG_OK)                              &
         call terminate("writeCGNSKomegaModifiedInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the turbulent closure type.

    call cg_model_write_f("TurbulenceClosure_t", EddyViscosity, ierr)
    if(ierr /= CG_OK)                              &
         call terminate("writeCGNSKomegaModifiedInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the details of the turbulence model under the turbulent
    ! closure type.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1, &
         "TurbulenceClosure_t", 1, "end")
    if(ierr /= CG_OK)                              &
         call terminate("writeCGNSKomegaModifiedInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write the value of the turbulent prandtl number.

    val = prandtlTurb
    call cg_array_write_f(cgnsPrandtlTurb, realTypeCGNS, &
         1, 1, val, ierr)
    if(ierr /= CG_OK)                              &
         call terminate("writeCGNSKomegaModifiedInfo", &
         "Something wrong when calling cg_array_write_f")

    ! Indicate that this is a nonDimensional parameter.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "DataArray_t", 1,"end")
    if(ierr /= CG_OK)                                  &
         call terminate("writeCGNSKomegaModifiedInfo", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter,ierr)
    if(ierr /= CG_OK)                              &
         call terminate("writeCGNSKomegaModifiedInfo", &
         "Something wrong when calling &
         &cg_dataclass_write_f")
  end subroutine writeCGNSKomegaModifiedInfo

  subroutine writeCGNSKomegaWilcoxInfo(cgnsInd, cgnsBase)
    !
    !       writeCGNSKomegaWilcoxInfo writes information about the         
    !       standard Wilcox k-omega turbulence model to the cgns file.     
    !
    use inputPhysics
    use cgnsNames
    use su_cgns
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, ierr

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Write the info of the turbulence model under the flow equation
    ! set. So move to this location first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                            &
         call terminate("writeCGNSKomegaWilcoxInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write that the k-omega model is used.

    call cg_model_write_f("TurbulenceModel_t", &
         TwoEquation_Wilcox, ierr)
    if(ierr /= CG_OK)                            &
         call terminate("writeCGNSKomegaWilcoxInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the turbulent closure type.

    call cg_model_write_f("TurbulenceClosure_t", EddyViscosity, ierr)
    if(ierr /= CG_OK)                            &
         call terminate("writeCGNSKomegaWilcoxInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the details of the turbulence model under the turbulent
    ! closure type.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "end")
    if(ierr /= CG_OK)                            &
         call terminate("writeCGNSKomegaWilcoxInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write the value of the turbulent prandtl number.

    val = prandtlTurb
    call cg_array_write_f(cgnsPrandtlTurb, realTypeCGNS, &
         1, 1, val, ierr)
    if(ierr /= CG_OK)                            &
         call terminate("writeCGNSKomegaWilcoxInfo", &
         "Something wrong when calling cg_array_write_f")

    ! Indicate that this is a nonDimensional parameter.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "DataArray_t", 1,"end")
    if(ierr /= CG_OK)                            &
         call terminate("writeCGNSKomegaWilcoxInfo", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter,IERR)
    if(ierr /= CG_OK)                            &
         call terminate("writeCGNSKomegaWilcoxInfo", &
         "Something wrong when calling &
         &cg_dataclass_write_f")
  end subroutine writeCGNSKomegaWilcoxInfo

  subroutine writeCGNSKtauInfo(cgnsInd, cgnsBase)
    !
    !       WriteCGNSKtauInfo writes information about the k-tau           
    !       turbulence model to the cgns file.                             
    !
    use inputPhysics
    use cgnsNames
    use su_cgns
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, ierr

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Write the info of the turbulence model under the flow equation
    ! set. So move to this location first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSKtauInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write that user defined model is used; k-tau is not
    ! supported by cgns.

    call cg_model_write_f("TurbulenceModel_t", &
         UserDefined, ierr)
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSKtauInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the turbulent closure type.

    call cg_model_write_f("TurbulenceClosure_t", EddyViscosity, ierr)
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSKtauInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the details of the turbulence model under the turbulent
    ! closure type.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "end")
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSKtauInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write the value of the turbulent prandtl number.

    val = prandtlTurb
    call cg_array_write_f(cgnsPrandtlTurb, realTypeCGNS, &
         1, 1, val, ierr)
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSKtauInfo", &
         "Something wrong when calling cg_array_write_f")

    ! Indicate that this is a nonDimensional parameter.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "DataArray_t", 1,"end")
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSKtauInfo", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter,ierr)
    if(ierr /= CG_OK)                    &
         call terminate("writeCGNSKtauInfo", &
         "Something wrong when calling &
         &cg_dataclass_write_f")
  end subroutine writeCGNSKtauInfo

  subroutine writeCGNSMenterSSTInfo(cgnsInd, cgnsBase)
    !
    !       WriteCGNSMenterSSTInfo writes information about menter's       
    !       SST turbulence model to the cgns file.                         
    !
    use inputPhysics
    use cgnsNames
    use su_cgns
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, ierr

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.
    ! Note that this info is only written to the 1st base.

    realTypeCGNS = setCGNSRealType()

    ! Write the info of the turbulence model under the flow equation
    ! set. So move to this location first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                         &
         call terminate("writeCGNSMenterSSTInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write that the SST variant of the kOmega model is used.

    call cg_model_write_f("TurbulenceModel_t", &
         TwoEquation_MenterSST, ierr)
    if(ierr /= CG_OK)                         &
         call terminate("writeCGNSMenterSSTInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the turbulent closure type.

    call cg_model_write_f("TurbulenceClosure_t", EddyViscosity, ierr)
    if(ierr /= CG_OK)                         &
         call terminate("writeCGNSMenterSSTInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the details of the turbulence model under the turbulent
    ! closure type.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "end")
    if(ierr /= CG_OK)                         &
         call terminate("writeCGNSMenterSSTInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write the value of the turbulent prandtl number.

    val = prandtlTurb
    call cg_array_write_f(cgnsPrandtlTurb, realTypeCGNS, &
         1, 1, val, ierr)
    if(ierr /= CG_OK)                         &
         call terminate("writeCGNSMenterSSTInfo", &
         "Something wrong when calling cg_array_write_f")

    ! Indicate that this is a nonDimensional parameter.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "DataArray_t", 1,"end")
    if(ierr /= CG_OK)                         &
         call terminate("writeCGNSMenterSSTInfo", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter,ierr)
    if(ierr /= CG_OK)                         &
         call terminate("writeCGNSMenterSSTInfo", &
         "Something wrong when calling &
         &cg_dataclass_write_f")
  end subroutine writeCGNSMenterSSTInfo

  subroutine writeCGNSReferenceState(cgnsInd, cgnsBase)
    !
    !       writeCGNSReferenceState writes the reference state to the      
    !       cgns file. Enough info is specified such that a restart can be 
    !       performed by a different solver, which uses a different        
    !       nonDimensionalization.                                         
    !
    use constants
    use cgnsNames
    use su_cgns
    use inputPhysics
    use flowVarRefState
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: ierr, realTypeCGNS, ii

    integer(kind=intType) :: i

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Go to the base.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "end")
    if(ierr /= CG_OK)                      &
         call terminate("writeReferenceState", &
         "Something wrong when calling cg_goto_f")

    ! Create the reference state node with a nice description.

    call cg_state_write_f("Reference state variables for &
         &nonDimensional data. Variables are &
         &nonDimensionalized using the reference &
         &density, pressure and temperature.", ierr)
    if(ierr /= CG_OK)                      &
         call terminate("writeReferenceState", &
         "Something wrong when calling cg_state_write_f")

    ! The actual data should be written below the reference state
    ! node. So go there first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "ReferenceState_t", 1, "end")
    if(ierr /= CG_OK)                      &
         call terminate("writeReferenceState", &
         "Something wrong when calling cg_goto_f")

    ! Write the Mach number and indicate that it is a nonDimensional
    ! parameter

    val = Mach
    call cg_array_write_f(cgnsMach, realTypeCGNS, 1, 1, val, ierr)
    if(ierr /= CG_OK)                      &
         call terminate("writeReferenceState", &
         "Something wrong when calling cg_array_write_f")

    ii = 1
    call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
         "DataArray_t", ii, "end")
    if(ierr /= CG_OK)                      &
         call terminate("writeReferenceState", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter, ierr)
    if(ierr /= CG_OK)                      &
         call terminate("writeReferenceState", &
         "Something wrong when calling &
         &cg_dataclass_write_f")

    ! Write the 3 flow angles. The units are degrees.

    velocityDir: do i=1,3

       ! Go to the reference state node.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, &
            "ReferenceState_t", 1, "end")
       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling cg_goto_f")

       ! Store component i of the direction in val.

       val = velDirFreestream(i)

       select case(i)
       case (1_intType)

          call cg_array_write_f(cgnsVelVecX, realTypeCGNS, &
               1, 1, val, ierr)

       case (2_intType)
          call cg_array_write_f(cgnsVelVecY, realTypeCGNS, &
               1, 1, val, ierr)

       case (3_intType)
          call cg_array_write_f(cgnsVelVecZ, realTypeCGNS, &
               1, 1, val, ierr)
       end select

       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling &
            &cg_array_write_f")

       ! Write the info that the unit vector is nondimensional.

       ii = ii + 1
       call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
            "DataArray_t", ii, "end")
       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling cg_goto_f")

       call cg_dataclass_write_f(NonDimensionalParameter, ierr)
       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling &
            &cg_dataclass_write_f")

    enddo velocityDir


    ! Write some reference values of the density, pressure, temperature,
    ! velocity and length.

    refLoop: do i=1,5

       ! Go to the reference state node.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, &
            "ReferenceState_t", 1, "end")
       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling cg_goto_f")

       ! Write a value, depending on i.

       select case(i)
       case (1_intType)
          val = rhoref
          call cg_array_write_f(cgnsDensity, realTypeCGNS, &
               1, 1, val, ierr)
       case (2_intType)
          val = pref
          call cg_array_write_f(cgnsPressure, realTypeCGNS, &
               1, 1, val, ierr)

       case (3_intType)
          val = Tref
          call cg_array_write_f(cgnsTemp, realTypeCGNS, &
               1, 1, val, ierr)

       case (4_intType)
          val = sqrt(pref/rhoref)
          call cg_array_write_f(cgnsVelocity, realTypeCGNS, &
               1, 1, val, ierr)

       case (5_intType)
          val = one
          call cg_array_write_f(cgnsLength, realTypeCGNS, &
               1, 1, val, ierr)
       end select

       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling &
            &cg_array_write_f")

       ! Write the info that the this reference value is dimensional
       ! and based on si units.

       ii = ii + 1
       call cg_goto_f(cgnsInd, cgnsBase, ierr, "ReferenceState_t", 1, &
            "DataArray_t", ii, "end")
       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling cg_goto_f")

       call cg_dataclass_write_f(Dimensional, ierr)
       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling &
            &cg_dataclass_write_f")

       call cg_units_write_f(Kilogram, Meter, Second, Kelvin, &
            Null, ierr)
       if(ierr /= CG_OK)                      &
            call terminate("writeReferenceState", &
            "Something wrong when calling &
            &cg_units_write_f")

    enddo refLoop
  end subroutine writeCGNSReferenceState

  subroutine writeCGNSSaInfo(cgnsInd, cgnsBase)
    !
    !       WriteCGNSSaInfo writes information about the Spalart           
    !       Allmaras turbulence model to the cgns file.                    
    !
    use inputPhysics
    use cgnsNames
    use su_cgns
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, ierr

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Write the info of the turbulence model under the flow equation
    ! set. So move to this location first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write that the spalart-allmaras model is used.

    call cg_model_write_f("TurbulenceModel_t", &
         OneEquation_SpalartAllmaras, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the turbulent closure type.

    call cg_model_write_f("TurbulenceClosure_t", EddyViscosity, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the details of the turbulence model under the turbulent
    ! closure type.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "end")
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write the value of the turbulent prandtl number.

    val = prandtlTurb
    call cg_array_write_f(cgnsPrandtlTurb, realTypeCGNS, &
         1, 1, val, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_array_write_f")

    ! Indicate that this is a nonDimensional parameter.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "DataArray_t", 1,"end")
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter,ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling &
         &cg_dataclass_write_f")
  end subroutine writeCGNSSaInfo

  subroutine writeCGNSSaeInfo(cgnsInd, cgnsBase)
    !
    !       WriteCGNSSaeInfo writes information about the Spalart          
    !       Allmaras turbulence model using the Edwards modification to    
    !       the cgns file.                                                 
    !
    use inputPhysics
    use cgnsNames
    use su_cgns
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, ierr

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Write the info of the turbulence model under the flow equation
    ! set. So move to this location first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write that the spalart-allmaras model is used.

    call cg_model_write_f("TurbulenceModel_t", &
         OneEquation_SpalartAllmaras, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the turbulent closure type.

    call cg_model_write_f("TurbulenceClosure_t", EddyViscosity, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the details of the turbulence model under the turbulent
    ! closure type.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "end")
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write the value of the turbulent prandtl number.

    val = prandtlTurb
    call cg_array_write_f(cgnsPrandtlTurb, realTypeCGNS, &
         1, 1, val, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_array_write_f")

    ! Indicate that this is a nonDimensional parameter.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "DataArray_t", 1,"end")
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter,ierr)
    if(ierr /= CG_OK)                  &
         call terminate("writeCGNSSaInfo", &
         "Something wrong when calling &
         &cg_dataclass_write_f")
  end subroutine writeCGNSSaeInfo

  subroutine writeCGNSV2fInfo(cgnsInd, cgnsBase)
    !
    !       WriteCGNSV2fInfo writes information about Durbin's v2f         
    !       turbulence model to the cgns file.                             
    !
    use inputPhysics
    use cgnsNames
    use su_cgns
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: realTypeCGNS, ierr

    real(kind=cgnsRealType) :: val

    ! Set the cgns real type.

    realTypeCGNS = setCGNSRealType()

    ! Write the info of the turbulence model under the flow equation
    ! set. So move to this location first.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, &
         "FlowEquationSet_t", 1, "end")
    if(ierr /= CG_OK)                   &
         call terminate("writeCGNSV2fInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write that user defined model is used; v2-f is not
    ! supported by cgns.

    call cg_model_write_f("TurbulenceModel_t", UserDefined, ierr)

    if(ierr /= CG_OK)                   &
         call terminate("writeCGNSV2fInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the turbulent closure type.

    call cg_model_write_f("TurbulenceClosure_t", EddyViscosity, ierr)
    if(ierr /= CG_OK)                   &
         call terminate("writeCGNSV2fInfo", &
         "Something wrong when calling cg_model_write_f")

    ! Write the details of the turbulence model under the turbulent
    ! closure type.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "end")
    if(ierr /= CG_OK)                   &
         call terminate("writeCGNSV2fInfo", &
         "Something wrong when calling cg_goto_f")

    ! Write the value of the turbulent prandtl number.

    val = prandtlTurb
    call cg_array_write_f(cgnsPrandtlTurb, realTypeCGNS, &
         1, 1, val, ierr)
    if(ierr /= CG_OK)                   &
         call terminate("writeCGNSV2fInfo", &
         "Something wrong when calling cg_array_write_f")

    ! Indicate that this is a nonDimensional parameter.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "FlowEquationSet_t", 1,&
         "TurbulenceClosure_t", 1, "DataArray_t", 1,"end")
    if(ierr /= CG_OK)                   &
         call terminate("writeCGNSV2fInfo", &
         "Something wrong when calling cg_goto_f")

    call cg_dataclass_write_f(NonDimensionalParameter,ierr)
    if(ierr /= CG_OK)                   &
         call terminate("writeCGNSV2fInfo", &
         "Something wrong when calling &
         &cg_dataclass_write_f")
  end subroutine writeCGNSV2fInfo

end module outputMod
