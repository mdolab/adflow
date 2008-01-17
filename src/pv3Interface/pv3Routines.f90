!
!      ******************************************************************
!      *                                                                *
!      * File:          pv3Routines.f90                                 *
!      * Author:        Juan J. Alonso, Edwin van der Weide,            *
!      *                Steve Repsher, Seonghyeon Hahn                  *
!      * Starting date: 04-19-2003                                      *
!      * Last modified: 11-04-2005                                      *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      *   SUmb PV3 multiblock visualization interface                  *
!      *                                                                *
!      *   All provided subroutines for PV3 are based on pV3 Rev 2.00   *
!      *   hence a pV3Server with version > 2.00 should be used.        *
!      *                                                                *
!      ******************************************************************
!
       subroutine initializePV3
!
!      ******************************************************************
!      *                                                                *
!      * InitializePV3 is called at the beginning of program execution  *
!      * to setup all the necessary data structures for real-time       *
!      * visualization using PV3 (R.Haimes-MIT, http://raphael.mit.edu) *
!      *                                                                *
!      * Note that this routine can only be called once the basic data  *
!      * structures in SUmb, including domain decomposition, have       *
!      * been allocated and setup.                                      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use constants
       use flowVarRefState
       use inputPhysics
       use inputVisualization
       implicit none
!
!      Local parameters.
!
       integer(kind=intPV3Type), parameter :: nKeys = 26, nCuts = 5
!
!      Local variables.
!
       integer(kind=intPV3Type) :: mKeys
       integer(kind=intPV3Type) :: clientID, mirror, maxblock
       integer(kind=intPV3Type) :: iopt, PV3Status
       integer(kind=intPV3Type), dimension(nKeys) :: keyType, keyBind

       real(kind=realPV3Type), dimension(2,nKeys) :: fLims
       real(kind=realPV3Type), dimension(4,4)     :: repMat

       character (len=10) :: string
       character (len=20) :: clientName, discipline
       character (len=32) :: PV3Title, cutplaneName(nCuts)
       character (len=32) :: keyTitle(nKeys)
!
!      ******************************************************************
!      *                                                                *
!      * Initialize variables                                           *
!      *                                                                *
!      ******************************************************************
!
       !
       ! Set title to be printed out to the screen
       !

       PV3Title = "# SUmb (pV3 Rev. 2.0) client"
       if (myID == 0) print "(a)", PV3Title

       !
       ! Set names of programmed cutting planes
       !

       cutplaneName(1) = 'x coordinate'
       cutplaneName(2) = 'y coordinate'
       cutplaneName(3) = 'z coordinate'
       cutplaneName(4) = 'r coordinate'
       cutplaneName(5) = 'theta coordinate'

       !
       ! Key titles
       !

       ! Scalar fields

       keyTitle(:)  = '-- not defined --               '
       keyTitle(1)  = 'density                         '
       keyTitle(2)  = 'x-momentum                      '
       keyTitle(3)  = 'y-momentum                      '
       keyTitle(4)  = 'z-momentum                      '
       keyTitle(5)  = 'total energy                    '
       keyTitle(6)  = 'Mach number                     '
       keyTitle(7)  = 'velocity magnitude              '
       keyTitle(8)  = 'pressure                        '
       keyTitle(9)  = 'total pressure                  '
       keyTitle(10) = 'temperature                     '
       keyTitle(11) = 'total temperature               '
       keyTitle(12) = 'enthalpy                        '
       keyTitle(13) = 'total enthalpy                  '
       keyTitle(14) = 'entropy                         '
       keyTitle(15) = 'Cp                              '
       keyTitle(16) = 'processor number                '
       keyTitle(17) = 'block number                    '
       keyTitle(18) = 'overset regions                 '

       ! Vector fields

       keyTitle(19) = 'velocity vector                 '
       keyTitle(20) = 'momentum vector                 '

       !
       ! Key types, 1 - scalar, 2 - vector
       !

       keyType(:)      = 1
       keyType(1:18)   = 1
       keyType(19:20)  = 2

       !
       ! Key bindings (in ascii code)
       !

       keyBind(:)  = 0
       keyBind(1)  = 100      ! -- 'd'
       keyBind(2)  = 120      ! -- 'x'
       keyBind(3)  = 121      ! -- 'y'
       keyBind(4)  = 122      ! -- 'z'
       keyBind(5)  = 69       ! -- 'E'
       keyBind(6)  = 109      ! -- 'm'
       keyBind(7)  = 113      ! -- 'q'
       keyBind(8)  = 112      ! -- 'p'
       keyBind(9)  = 80       ! -- 'P'
       keyBind(10) = 116      ! -- 't'
       keyBind(11) = 84       ! -- 'T'
       keyBind(12) = 104      ! -- 'h'
       keyBind(13) = 72       ! -- 'H'
       keyBind(14) = 115      ! -- 's'
       keyBind(15) = 67       ! -- 'C'
       keyBind(16) = 35       ! -- '#'
       keyBind(17) = 82       ! -- 'R'
       keyBind(18) = 105      ! -- 'i'
       keyBind(19) = 86       ! -- 'V'
       keyBind(20) = 77       ! -- 'M'

       !
       ! Key information for turbulence related variables.
       ! Furthermore set the value of mKeys, the number of pre-defined
       ! keys. This number depends on the equations solved and the
       ! turbulence model (for RANS).
       !

       equationModel: select case (equations)
          case (EulerEquations)
            mKeys = 20

          !==============================================================

          case (NSEquations)

             keyTitle(21) = 'laminar viscosity               '
             keyBind (21) = 108      ! ("l")
             keyType (21) = 1

             mKeys = 21

          case (RANSEquations)

             ! Determine the turbulence model and set the parameters
             ! accordingly.

             select case (turbModel)

               case (spalartAllmaras, spalartAllmarasEdwards)

                 keyTitle(21) = 'eddy viscosity                  '
                 keyTitle(22) = 'eddy/laminar viscosity ratio    '
                 keyTitle(23) = 'nu tilde (working variable)     '
                 keyTitle(24) = 'chi                             '

                 keyBind(21)  = 101      ! ("e")
                 keyBind(22)  = 114      ! ("r")
                 keyBind(23)  = 110      ! ("n")
                 keyBind(24)  = 99       ! ("c")

                 keyType(21)  = 1
                 keyType(22)  = 1
                 keyType(23)  = 1
                 keyType(24)  = 1

                 mKeys = 24

               !=========================================================

               case (komegaWilcox, komegaModified, menterSST)

                 keyTitle(21) = 'eddy viscosity                  '
                 keyTitle(22) = 'eddy/laminar viscosity ratio    '
                 keyTitle(23) = 'turbulent kinetic energy(k)     '
                 keyTitle(24) = 'turbulent dissip. rate(omega)   '

                 keyBind(21)  = 101      ! ("e")
                 keyBind(22)  = 114      ! ("r")
                 keyBind(23)  = 107      ! ("k")
                 keyBind(24)  = 111      ! ("o")

                 keyType(21)  = 1
                 keyType(22)  = 1
                 keyType(23)  = 1
                 keyType(24)  = 1

                 mKeys = 24

               !=========================================================

               case (ktau)

                 keyTitle(21) = 'eddy viscosity                  '
                 keyTitle(22) = 'eddy/laminar viscosity ratio    '
                 keyTitle(23) = 'turbulent kinetic energy(k)     '
                 keyTitle(24) = 'inverse turb. dissip. rate(tau) '

                 keyBind(21)  = 101      ! ("e")
                 keyBind(22)  = 114      ! ("r")
                 keyBind(23)  = 107      ! ("k")
                 keyBind(24)  = 111      ! ("o"), t and T are
                                          !        Already taken

                 keyType(21)  = 1
                 keyType(22)  = 1
                 keyType(23)  = 1
                 keyType(24)  = 1

                 mKeys = 24

               !=========================================================

               case (v2f)

                 keyTitle(21) = 'eddy viscosity                  '
                 keyTitle(22) = 'eddy/laminar viscosity ratio    '
                 keyTitle(23) = 'turbulent kinetic energy(k)     '
                 keyTitle(24) = 'turbulent dissipation (epsilon) '
                 keyTitle(25) = 'turbulent scalar v2             '
                 keyTitle(26) = 'turbulent scalar f              '

                 keyBind(21)  = 101      ! ("e")
                 keyBind(22)  = 114      ! ("r")
                 keyBind(23)  = 107      ! ("k")
                 keyBind(24)  = 111      ! ("o"), e and E are
                                         !        already taken
                 keyBind(25)  = 118      ! ("v")
                 keyBind(26)  = 102      ! ("f")

                 keyType(21)  = 1
                 keyType(22)  = 1
                 keyType(23)  = 1
                 keyType(24)  = 1
                 keyType(25)  = 1
                 keyType(26)  = 1

                 mKeys = 26

             end select

       end select equationModel

       !
       ! Scalar field limits (set equal to each other for automatic
       !  limits being set by PV3 to min and max of the whole field)
       !

       fLims(:,:) = 1.0_realPV3Type
!
!      ******************************************************************
!      *                                                                *
!      * Setup PV3 constants                                            *
!      *                                                                *
!      * PV3Title      Title (up to 80 characters)                      *
!      * cid           Unique integer client id                         *
!      * clientname    Unique character client name                     *
!      * discipline    Discipline name to which this proc belongs       *
!      *                                                                *
!      * Iopt   = -3   Structure unsteady with connectivity supplied    *
!      * iopt   = -2   Unsteady grid/data with connectivity supplied    *
!      * iopt   = -1   Steady grid and unsteady data with connectivity  *
!      *               supplied                                         *
!      * iopt   =  0   Steady grid and data                             *
!      * iopt   =  1   Steady grid and unsteady data                    *
!      * iopt   =  2   Unsteady grid and data                           *
!      *                                                                *
!      * nCuts         Number of programmed cuts                        *
!      * cutplaneName  Names of programmed cutting planes               *
!      * mKeys         Number of pre-defined keys for fields            *
!      * keyBind       Ascii key bindings for fields                    *
!      * keyTitle      Key titles                                       *
!      * keyType       Key types (scalar - 1, vector - 2)               *
!      * fLims(2,*)    Upper and lower limits for fields                *
!      *                                                                *
!      * mirr   < 0    nRep replicate the data -nrep times using repMat *
!      * mirr   = 0    No mirroring                                     *
!      * mirr   = 1    Mirroring about x=0.0 plane                      *
!      * mirr   = 2    Mirroring about y=0.0 plane                      *
!      * mirr   = 3    Mirroring about z=0.0 plane                      *
!      *                                                                *
!      * repMat        Replication matrix                               *
!      *                                                                *
!      * maxblock      Max number of structured blocks in this client   *
!      *                                                                *
!      * PV3Status     PV3 initialization input and return status       *
!      *                                                                *
!      ******************************************************************
!
!      Client ID [1,nProc], client name, and discipline.
!
       clientID = myID +1

       write(string,'(i7)') clientID
       clientName = 'SUmb proc'//trim(adjustl(string))

       discipline = "SUmb"
!
!      Visualization option (*** check with Edwin to generalize ***)
!
       iopt  = -3
       if( PV3VisOnly ) iopt = 0
!
!      Mirror option - default to no mirroring
!
       mirror = 0
!
!      Repetition matrix (only used if mirror < 0, for replication
!                         with rotational symmetry)
!      Default to no replication
!
       repMat = 0.0_realPV3Type
!
!      Maximum number of blocks in this client
!
       maxblock = nDom
!
!      Initialize PV3Status according to PV3 manual.
!
!      Input value of PV3Status is addition of:
!
!      0 - do not wait for server to be started
!          do not terminate when server is stopped
!      1 - wait for the server to startup to proceed with calculation
!      2 - terminte computation when server terminates
!      4 - non time-accurate mode
!
       PV3Status = 0
       if( PV3VisOnly ) PV3Status = 3
!
!      Initialize PV3
!
       call pv_init(PV3Title, clientID, clientName, discipline, iopt,   &
                    nCuts, cutplaneName, mKeys, keyBind, keyTitle,      &
                    keyType, fLims, mirror, repMat, maxblock, PV3Status)
!
!      Report the status of PV3 initialization
!
       if (myID == 0) then
          print "(a,1x,i5)", "# PV3 initialization status =", &
                PV3Status
          print "(a)", "#"
       end if
!
       return
       end subroutine initializePV3

!      ==================================================================

       subroutine pvstruc(knode,kequiv,kcel1,kcel2,kcel3,kcel4,knptet,  &
                          kptet,knblock,blocks,kphedra,ksurf,knsurf,hint)
!
!      ******************************************************************
!      *                                                                *
!      * pvstruc provides size information so that PV3 can allocate     *
!      *  its own visualization data structures.  It is called only     *
!      *  once, except for structure unsteady cases (where the number   *
!      *  of blocks, their sizes, connectivity, etc. Change during the  *
!      *  calculation, as is possible during grid sequencing)           *
!      *                                                                *
!      * knode         Number of non-structured block nodes             *
!      *                                                                *
!      * kequiv        Number of equivalence pairs                      *
!      *                                                                *
!      * kcel1         Number of tetrahedra                             *
!      * kcel2         Number of pyramids                               *
!      * kcel3         Number of prism                                  *
!      * kcel4         Number of hexahedra                              *
!      * knptet        Number of poly-tetrahedral strips                *
!      * kptet         Number of tetrahedral cells in all strips        *
!      *                                                                *
!      * knblock       Number of structured blocks                      *
!      * iblocks(3,*)  Block dimensions                                 *
!      *                                                                *
!      * ksurf         Number of surface faces                          *
!      * knsurf        Number of surface groups                         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type) :: knode, kequiv, kcel1, kcel2, kcel3
       integer(kind=intPV3Type) :: kcel4, knptet, kptet, knblock, kphedra
       integer(kind=intPV3Type) :: ksurf, knsurf, hint
       integer(kind=intPV3Type), dimension (3,*) :: blocks
       integer(kind=intPV3Type), dimension (nBCs+1) :: bcPresent
!
!      Local variables.
!
       integer(kind=intType) :: kk, nn, ibc, i, idim, jdim, kdim
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      Number of non-block nodes
!
       knode = 0
!
!      Number of node equivalency pairs
!
       kequiv = 0
!
!      Number of tets, pyramids, prisms, hexs, polyhedra
!
       kcel1   = 0
       kcel2   = 0
       kcel3   = 0
       kcel4   = 0
       kphedra = 0
!
!      Number of poly-tetrahedral strips and poly-tetrahedral cells
!
       knptet = 0
       kptet  = 0
!
!      Number of structured blocks and their sizes. Note knblock
!      is negative to indicate that blanking is supplied, only if
!      overset grids are present.
!
       knblock = nDom
       if (oversetPresent) knblock = -nDom

       domains: do nn=1,nDom

         ! The size for all spectral modes is the same.

         blocks(1,nn) = flowDoms(nn,groundLevel,1)%il
         blocks(2,nn) = flowDoms(nn,groundLevel,1)%jl
         blocks(3,nn) = flowDoms(nn,groundLevel,1)%kl

       end do domains
!
!      Number of surface cells and number of surface cell group.
!       *** modify later when actual group definitions are in place ***
!
!      For the time being, lump all subfaces with the same kind of boundary
!       condition onto a group named after that boundary condition
!
!      Count existing types of boundary conditions by looping over the
!       subfaces of all blocks
!
       bcPresent = 0
       ksurf     = 0
       knsurf    = 0
 
       blockLoop: do nn=1,nDom

         ! Set the pointers for this block. The information for all
         ! spectral modes is identical, so it is okay to set the
         ! pointers to the 1st spectral mode.

         call setPointers(nn, groundLevel, 1_intType)

         ! Loop over all boundary condition subfaces

         subfaces: do i=1,nSubface

            ibc = 0

            select case (BCType(i))
               case (BCNull)
                  ibc = 1
               case(Symm)
                  ibc = 2
               case(SymmPolar)
                  ibc = 3
               case(NSWallAdiabatic)
                  ibc = 4
               case(NSWallIsothermal)
                  ibc = 5
               case(EulerWall)
                  ibc = 6
               case(FarField)
                  ibc = 7
               case(SupersonicInflow)
                  ibc = 8
               case(SubsonicInflow)
                  ibc = 9
               case(SupersonicOutflow)
                  ibc = 10
               case(SubsonicOutflow)
                  ibc = 11
               case(MassBleedInflow)
                  ibc = 12
               case(MassBleedOutflow)
                  ibc = 13
               case(mDot)
                  ibc = 14
               case(Thrust)
                  ibc = 15
               case(Extrap)
                  ibc = 16
               case(SlidingInterface)
                  ibc = 17
               case(OversetOuterBound)
                  ibc = 18
               case(DomainInterfaceAll)
                  ibc = 19
               case(DomainInterfaceRhoUVW)
                  ibc = 20
               case(DomainInterfaceP)
                  ibc = 21
               case(DomainInterfaceRho)
                  ibc = 22
               case(DomainInterfaceTotal)
                  ibc = 23
               case(B2BMatch, B2BMismatch)

                  ! Check for periodicity.

                  kk = cgnsSubface(i)
                  if(kk > 0) then
                    if( cgnsDoms(nbkGlobal)%conn1to1(kk)%periodic ) &
                      ibc = 24
                  endif
            end select

            if (ibc /= 0) then
            if (bcPresent(ibc) == 0) then
               knsurf         = knsurf +1
               bcPresent(ibc) = 1
            end if
            end if

            if (ibc > 0) then
               idim  = max(1_intType,(inEnd(i)-inBeg(i)))
               jdim  = max(1_intType,(jnEnd(i)-jnBeg(i)))
               kdim  = max(1_intType,(knEnd(i)-knBeg(i)))
               ksurf = ksurf + idim*jdim*kdim
            end if
         end do subfaces

       end do blockLoop
!
!      Hint for PV3 particle locations.  Since we will not be providing
!      block connectivity information (pvlocate and pvconnect) the
!      streamlines will not go through block boundaries.  This can be
!      fixed if pvlocate and pvconnect are written.
!
       hint = 1

       return
       end subroutine pvstruc

!      ==================================================================

       subroutine pvgrid(xyz)
!
!      ******************************************************************
!      *                                                                *
!      * pvgrid provides the node locations in 3D space                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realPV3Type), dimension (3,*) :: xyz
!
!      Local variables.
!
       integer(kind=intPV3Type) :: nn, i, j, k, nNode
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nNode = 0

       domains: do nn=1,nDom

         ! Set the pointers for this block. Consider only the 1st
         ! spectral solution for now.

         call setPointers(nn, groundLevel, 1_intType)
 
         do k=1,kl
           do j=1,jl
              do i=1,il
                nNode = nNode +1
                xyz(1,nNode) = x(i,j,k,1)
                xyz(2,nNode) = x(i,j,k,2)
                xyz(3,nNode) = x(i,j,k,3)
              end do
            end do
         end do

       end do domains

       return
       end subroutine pvgrid

!      ==================================================================

       subroutine pvblank(pblnk, tbcon)
!
!      ******************************************************************
!      *                                                                *
!      * Pvblank provides PV3 with the iblanking data for the mesh. The *
!      * arguments are:                                                 *
!      *                                                                *
!      *   Pblnk =  0 - hole                                            *
!      *            1 - field                                           *
!      *           -# - continue streamline on local block #            *
!      *                                                                *
!      *   Tbcon = clientID of continuing block for pblnk < 0           *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use communication
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type), dimension (*) :: pblnk, tbcon
!
!      Local variables.
!
       integer(kind=intPV3Type) :: nn, i, j, k, nNode, ibmax, ss
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nNode = 0

       domains: do nn=1,nDom

         ! Set the pointers for this block. Consider only the 1st
         ! spectral solution for now. Then change the iblanks on the
         ! fringe to make sure they aren't blanked.

         call setPointers(nn, groundLevel, 1_intType)
         call changeIblanks(.false., 1_intType)
 
         do k=1,kl
           do j=1,jl
             do i=1,il
 
               nNode = nNode + 1
 
               if (any(iblank(i:i+1,j:j+1,k:k+1) <= 0)) then
                 pblnk(nNode) = 0
                 tbcon(nNode) = 0
               else
                 ibmax = maxval(iblank(i:i+1,j:j+1,k:k+1))
                 if (ibmax >= 10) then
                   ss = ibmax/10
                   pblnk(nNode) = -neighBlockOver(ss)
                   if (neighProcOver(ss) == myId) then
                     tbcon(nNode) = 0
                   else
                     tbcon(nNode) = neighProcOver(ss) + 1
                   end if
                 else
                   pblnk(nNode) = 1
                   tbcon(nNode) = 0
                 end if
               end if

             end do
           end do
         end do

         ! Re-zero the iblanks on the fringe to blank them in the solver.

         call changeIblanks(.false., 0_intType)

       end do domains

       return
       end subroutine pvblank

!      ==================================================================

       subroutine pvsurface(nsurf,scon,scel,tsurf)
!
!      ******************************************************************
!      *                                                                *
!      * pvsurface supplies PV3 with the surface data structures        *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type) :: nSurf(3,*), scon(*), scel(4,*)
       character(len=20)        :: tsurf(*)
!
!      Local variables.
!
       integer(kind=intType) :: knSurf, ksurf, ibc, i, j, k, l
       integer(kind=intType) :: ibase, nn, kk, kn, PV3BC
       logical               :: bcFound, storeSubface
       character(len=20)     :: bcName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      Fill in domain surface data
!
       ksurf  = 0
       knSurf = 0
!
!      Display all surfaces with boundary conditions different from
!       flow through. Note that if more boundary conditions are supported
!       additional cases will need to be added here
!      *** This will change when groups are properly supported ***
!
!      For the time being, lump all subfaces with the same kind of boundary
!       condition onto a group named after that boundary condition.
!
       bcs: do ibc=1,nBCs+1
         ibase   = 0
         bcFound = .false.
!
!        Process one boundary condition at a time
!
         PV3BC = bcNotValid

         select case (ibc)
           case(1_intType)
             PV3BC = BCNull
           case(2_intType)
             PV3BC = Symm
           case(3_intType)
             PV3BC = SymmPolar
           case(4_intType)
             PV3BC = NSWallAdiabatic
           case(5_intType)
             PV3BC = NSWallIsothermal
           case(6_intType)
             PV3BC = EulerWall
           case(7_intType)
             PV3BC = FarField
           case(8_intType)
             PV3BC = SupersonicInflow
           case(9_intType)
             PV3BC = SubsonicInflow
           case(10_intType)
             PV3BC = SupersonicOutflow
           case(11_intType)
             PV3BC = SubsonicOutflow
           case(12_intType)
             PV3BC = MassBleedInflow
           case(13_intType)
             PV3BC = MassBleedOutflow
           case(14_intType)
             PV3BC = mDot
           case(15_intType)
             PV3BC = Thrust
           case(16_intType)
             PV3BC = Extrap
           case(17_intType)
             PV3BC = SlidingInterface
           case(18_intType)
             PV3BC = OversetOuterBound
           case(19_intType)
             PV3BC = DomainInterfaceAll
           case(20_intType)
             PV3BC = DomainInterfaceRhoUVW
           case(21_intType)
             PV3BC = DomainInterfaceP
           case(22_intType)
             PV3BC = DomainInterfaceRho
           case(23_intType)
             PV3BC = DomainInterfaceTotal
           case(24_intType)
             PV3BC = B2BMatch   ! Just to indicate an internal block
                                ! boundary for periodic.
         end select
!
!        Loop over blocks and subfaces
!
         if (PV3BC /= BCNotValid) then

         blocks: do nn=1,nDom

           ! Set the pointers for this block. The information for all
           ! spectral modes is identical, so it is okay to set the
           ! pointers to the 1st spectral mode.

           call setPointers(nn, groundLevel, 1_intType)

           ! Loop over all boundary condition subfaces

           subfaces: do l=1,nSubface

              ! Distinguish between physical boundaries and internal
              ! boundaries; the latter may be periodic.

              storeSubface = .false.

              if(l <= nBocos) then
                if(BCType(l) == PV3BC) then
                  bcFound      = .true.
                  storeSubface = .true.
                endif
              else if(PV3BC == b2bMatch) then  ! Flag for periodic
                kk = cgnsSubface(l)
                if(kk > 0) then
                  if(cgnsDoms(nbkGlobal)%conn1to1(kk)%periodic) then
                    bcFound      = .true.
                    storeSubface = .true.
                  endif
                endif
              endif

              if( storeSubface ) then

                 select case (BCFaceID(l))
                    case (iMin)
                       do k=knBeg(l),knEnd(l)-1
                          do j=jnBeg(l),jnEnd(l)-1
                             ksurf = ksurf +1
                             kn    = ibase +1 +(j-1)*il +(k-1)*il*jl
                             scel(1,ksurf) = kn
                             scel(2,ksurf) = kn + il
                             scel(3,ksurf) = kn + il*jl +il
                             scel(4,ksurf) = kn + il*jl
                             scon(ksurf)   = 0
                          end do
                       end do

                    case (iMax)
                       do k=knBeg(l),knEnd(l)-1
                          do j=jnBeg(l),jnEnd(l)-1
                             ksurf = ksurf +1
                             kn    = ibase +il +(j-1)*il +(k-1)*il*jl
                             scel(1,ksurf) = kn
                             scel(2,ksurf) = kn + il
                             scel(3,ksurf) = kn + il*jl +il
                             scel(4,ksurf) = kn + il*jl
                             scon(ksurf)   = 0
                          end do
                       end do

                    case (jMin)
                       do k=knBeg(l),knEnd(l)-1
                          do i=inBeg(l),inEnd(l)-1
                             ksurf = ksurf +1
                             kn    = ibase +i +(k-1)*il*jl
                             scel(1,ksurf) = kn
                             scel(2,ksurf) = kn + 1
                             scel(3,ksurf) = kn + il*jl +1
                             scel(4,ksurf) = kn + il*jl
                             scon(ksurf)   = 0
                          end do
                       end do
 
                    case (jMax)
                       do k=knBeg(l),knEnd(l)-1
                          do i=inBeg(l),inEnd(l)-1
                             ksurf = ksurf +1
                             kn    = ibase +i +(jl-1)*il +(k-1)*il*jl
                             scel(1,ksurf) = kn
                             scel(2,ksurf) = kn + 1
                             scel(3,ksurf) = kn + il*jl +1
                             scel(4,ksurf) = kn + il*jl
                             scon(ksurf)   = 0
                          end do
                       end do

                    case (kMin)
                       do j=jnBeg(l),jnEnd(l)-1
                          do i=inBeg(l),inEnd(l)-1
                             ksurf = ksurf +1
                             kn    = ibase +i +(j-1)*il
                             scel(1,ksurf) = kn
                             scel(2,ksurf) = kn + 1
                             scel(3,ksurf) = kn + il +1
                             scel(4,ksurf) = kn + il
                             scon(ksurf)   = 0
                          end do
                       end do
 
                    case (kMax)
                       do j=jnBeg(l),jnEnd(l)-1
                          do i=inBeg(l),inEnd(l)-1
                             ksurf = ksurf +1
                             kn    = ibase +i +(j-1)*il +(kl-1)*il*jl
                             scel(1,ksurf) = kn
                             scel(2,ksurf) = kn + 1
                             scel(3,ksurf) = kn + il +1
                             scel(4,ksurf) = kn + il
                             scon(ksurf)   = 0
                          end do
                       end do

                 end select
              end if

           end do subfaces

           ibase = ibase +il*jl*kl

         end do blocks
         end if

         if ( bcFound ) then
            knSurf = knSurf +1
            nSurf(1,knSurf) = ksurf
            nSurf(2,knSurf) = 1
            nSurf(3,knSurf) = ibc

            ! Set human-readable bc name

            select case (ibc)
              case(1_intType)
                 bcName = "BCNull"
              case(2_intType)
                 bcName = "Symm"
              case(3_intType)
                 bcName = "SymmPolar"
              case(4_intType)
                 bcName = "NSWallAdiabatic"
              case(5_intType)
                 bcName = "NSWallIsothermal"
              case(6_intType)
                 bcName = "EulerWall"
              case(7_intType)
                 bcName = "FarField"
              case(8_intType)
                 bcName = "SupersonicInflow"
              case(9_intType)
                 bcName = "SubsonicInflow"
              case(10_intType)
                 bcName = "SupersonicOutflow"
              case(11_intType)
                 bcName = "SubsonicOutflow"
              case(12_intType)
                 bcName = "MassBleedInflow"
              case(13_intType)
                 bcName = "MassBleedOutflow"
              case(14_intType)
                 bcName = "mDot"
              case(15_intType)
                 bcName = "Thrust"
              case(16_intType)
                 bcName = "Extrap"
              case(17_intType)
                 bcName = "SlidingInterface"
              case(18_intType)
                 bcName = "OversetOuterBound"
              case(19_intType)
                 bcName = "DomainInterfaceAll"
              case(20_intType)
                 bcName = "DomainInterfaceRhoUVW"
              case(21_intType)
                 bcName = "DomainInterfaceP"
              case(22_intType)
                 bcName = "DomainInterfaceRho"
              case(23_intType)
                 bcName = "DomainInterfaceTotal"
              case(24_intType)
                 bcName = "Periodic"
            end select
            tsurf(knSurf) = bcName
         end if

       end do bcs

       return
       end subroutine pvsurface

!      ==================================================================

       subroutine pvscal(key,v)
!
!      ******************************************************************
!      *                                                                *
!      * pvscal provides the values of the components of scalar fields. *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use communication
       use flowVarRefState
       use inputPhysics
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type), intent(in)              :: key
       real(kind=realPV3Type), dimension(*), intent(out) :: v
!
!      Local parameter.
!
       real(kind=realType), parameter :: twoThird = two*third
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: knode, nn, i, j, k, l, i1(3), i2(3), bb

       integer(kind=intType), dimension(0:nDom) :: offset

       real(kind=realType) :: dummyK, PV3Gamma, v2, pl, tl, fact

       real(kind=realType), allocatable, dimension(:,:,:,:) :: wtemp
       real(kind=realType), allocatable, dimension(:,:,:)   :: lv, ev
       real(kind=realType), allocatable, dimension(:,:,:)   :: tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Displaying the overset regions is a special case so check if
       ! that is the key.

       displayOverset: if (key == 18) then

         ! Initialize the scalar array using the iblank variable. This
         ! will set the nodes of every boundary cell to 1.0 and zero
         ! everywhere else.

         knode = 0

         domains0: do nn=1,nDom

           ! Set the pointers for this block. Consider only the 1st
           ! spectral solution for now.

           call setPointers(nn, groundLevel, 1_intType)
           call changeIblanks(.false., 1_intType)

           do k=1,kl
             do j=1,jl
                do i=1,il

                  knode = knode +1
                  if (any(iblank(i:i+1,j:j+1,k:k+1) >= 9)) then
                    v(knode) = one
                  else
                    v(knode) = zero
                  end if

                end do
              end do
           end do

           call changeIblanks(.false., 0_intType)

         end do domains0

         ! Compute the displacements for each block into the scalar
         ! array being returned.

         offset(0) = 0
         do nn = 1,nDom
           offset(nn) = offset(nn-1) + flowDoms(nn,groundLevel,1)%il &
                                     * flowDoms(nn,groundLevel,1)%jl &
                                     * flowDoms(nn,groundLevel,1)%kl
         end do

         ! Loop over the communication sending list and set the value
         ! to 2.0 for the whole stencil (i.e. all the nodes of the
         ! cells making up the stencil.

         do nn = 1,commPatternOverset(groundLevel,1)%nProcSend
           do l = 1,commPatternOverset(groundLevel,1)%nsend(nn)

             bb = commPatternOverset(groundLevel,1)%sendList(nn)%block(l)
             i1 = commPatternOverset(groundLevel,1)%sendList(nn)%indices(l,:)

             ! Set the beginning and ending indices of the stencil.
             ! Be aware it may straddle a boundary so ignore that part.

             i2(1) = min(i1(1)+1, flowDoms(bb,groundLevel,1)%il)
             i2(2) = min(i1(2)+1, flowDoms(bb,groundLevel,1)%jl)
             i2(3) = min(i1(3)+1, flowDoms(bb,groundLevel,1)%kl)

             i1 = max(i1-1, 1_intType)

             ! Loop over the stencil and set the scalar to 2.

             do k = i1(3), i2(3)
               do j = i1(2), i2(2)
                 do i = i1(1), i2(1)

                   knode = offset(bb-1) + &
                           flowDoms(bb,groundLevel,1)%jl * &
                           flowDoms(bb,groundLevel,1)%il * (k - 1) + &
                           flowDoms(bb,groundLevel,1)%il * (j - 1) + i
                   v(knode) = two

                 end do
               end do
             end do

           end do
         end do

         ! Repeat the process for the internal list.

         do l = 1,internalOverset(groundLevel,1)%ncopy

           bb = internalOverset(groundLevel,1)%donorBlock(l)
           i1 = internalOverset(groundLevel,1)%donorIndices(l,:)

           ! Set the beginning and ending indices of the stencil.
           ! Be aware it may straddle a boundary so ignore that part.

           i2(1) = min(i1(1)+1, flowDoms(bb,groundLevel,1)%il)
           i2(2) = min(i1(2)+1, flowDoms(bb,groundLevel,1)%jl)
           i2(3) = min(i1(3)+1, flowDoms(bb,groundLevel,1)%kl)

           i1 = max(i1-1, 1_intType)

           ! Loop over the stencil and set the scalar to 2.

           do k = i1(3), i2(3)
             do j = i1(2), i2(2)
               do i = i1(1), i2(1)

                 knode = offset(bb-1) + &
                         flowDoms(bb,groundLevel,1)%jl * &
                         flowDoms(bb,groundLevel,1)%il * (k - 1) + &
                         flowDoms(bb,groundLevel,1)%il * (j - 1) + i
                 v(knode) = two

               end do
             end do
           end do

         end do

         ! Nothing further needs to be done.

         return

       end if displayOverset

       ! Carry on.

       knode = 0

       domains: do nn=1,nDom

         ! Set the pointers for this block. Consider only the 1st
         ! spectral solution for now.

         call setPointers(nn, groundLevel, 1_intType)

         ! Transfer flow field values from cell centers to nodes
         ! note that SUmb stores rho, u, v, w, and (rho E). However
         ! the entry for (rho E) will be used to store the pressure.
         ! This is easier for the non-constant gamma case.

         allocate(wtemp(il,jl,kl, nw), tmp(il,jl,kl), stat=ierr)
         if(ierr /= 0)              &
           call terminate("pvscal", &
                          "Memory allocation failure for wtemp &
                          &and tmp")

         do l=1,nw
           if(l == irhoE) then

             do k=1,kl
             do j=1,jl
             do i=1,il
               wtemp(i,j,k,l) = eighth*(p(i,j  ,k)   + p(i+1,j  ,k)   &
                              +         p(i,j+1,k)   + p(i+1,j+1,k)   &
                              +         p(i,j  ,k+1) + p(i+1,j  ,k+1) &
                              +         p(i,j+1,k+1) + p(i+1,j+1,k+1))
             enddo
             enddo
             enddo

           else

             do k=1,kl
             do j=1,jl
             do i=1,il
               wtemp(i,j,k,l) = eighth*(w(i  ,j  ,k  ,l) &
                              +         w(i+1,j  ,k  ,l) &
                              +         w(i  ,j+1,k  ,l) &
                              +         w(i+1,j+1,k  ,l) &
                              +         w(i  ,j  ,k+1,l) &
                              +         w(i+1,j  ,k+1,l) &
                              +         w(i  ,j+1,k+1,l) &
                              +         w(i+1,j+1,k+1,l))
             end do
             end do
             end do

           endif
         end do

         ! Allocate and compute the laminar viscosity for a viscous
         ! computation.

         if(equations ==   NSEquations .or. &
            equations == RANSEquations) then

           allocate(lv(il,jl,kl), stat=ierr)
           if(ierr /= 0)              &
             call terminate("pvscal", &
                            "Memory allocation failure for lv")

           do k=1,kl
           do j=1,jl
           do i=1,il
             lv(i,j,k) = eighth*(rlv(i,j  ,k)   + rlv(i+1,j  ,k)   &
                       +         rlv(i,j+1,k)   + rlv(i+1,j+1,k)   &
                       +         rlv(i,j  ,k+1) + rlv(i+1,j  ,k+1) &
                       +         rlv(i,j+1,k+1) + rlv(i+1,j+1,k+1))
           enddo
           enddo
           enddo

         endif

         ! Allocate and compute the eddy viscosity if an eddy model
         ! is used.

         if( eddyModel ) then

           allocate(ev(il,jl,kl), stat=ierr)
           if(ierr /= 0)              &
             call terminate("pvscal", &
                            "Memory allocation failure for ev")

           do k=1,kl
           do j=1,jl
           do i=1,il
             ev(i,j,k) = eighth*(rev(i,j  ,k)   + rev(i+1,j  ,k)   &
                       +         rev(i,j+1,k)   + rev(i+1,j+1,k)   &
                       +         rev(i,j  ,k+1) + rev(i+1,j  ,k+1) &
                       +         rev(i,j+1,k+1) + rev(i+1,j+1,k+1))
           enddo
           enddo
           enddo

         endif
!
!        ****************************************************************
!        *                                                              *
!        * Set return scalar fields                                     *
!        *                                                              *
!        ****************************************************************
!
         scalarField: select case(key)

           case (1_intPV3Type)        ! Density
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,irho)
             end do
             end do
             end do

           !=============================================================

           case (2_intPV3Type)        ! X-momentum
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,irho) * wtemp(i,j,k,ivx)
             end do
             end do
             end do

           !=============================================================

           case (3_intPV3Type)        ! Y-momentum
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,irho) * wtemp(i,j,k,ivy)
             end do
             end do
             end do

           !=============================================================

           case (4_intPV3Type)       ! Z-momentum
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,irho) * wtemp(i,j,k,ivz)
             end do
             end do
             end do

           !=============================================================

           case (5_intPV3Type)       ! Total energy

             ! Store the total number of nodes in j.

             j = il*jl*kl

             ! Make a distinction whether or not k is present.
             ! If not, just send a dummy in its place; it is not
             ! used in etotArray.

             if( kPresent ) then
               call etotArray(wtemp(1,1,1,irho),  wtemp(1,1,1,ivx),  &
                              wtemp(1,1,1,ivy),   wtemp(1,1,1,ivz),  &
                              wtemp(1,1,1,irhoE), wtemp(1,1,1,itu1), &
                              tmp, kPresent, j)
             else
               dummyK = zero
               call etotArray(wtemp(1,1,1,irho),  wtemp(1,1,1,ivx), &
                              wtemp(1,1,1,ivy),   wtemp(1,1,1,ivz), &
                              wtemp(1,1,1,irhoE), dummyK,           &
                              tmp, kPresent, j)
             endif

             ! Copy the values of tmp in v; v cannot be used directly
             ! in the call to etotArray, because of different accuracy.

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode = knode + 1
               v(knode) = tmp(i,j,k)
             enddo
             enddo
             enddo

           !=============================================================

           case(6_intPV3Type)        ! Mach number
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode = knode +1
               v2    = wtemp(i,j,k,ivx)**2 + wtemp(i,j,k,ivy)**2 &
                     + wtemp(i,j,k,ivz)**2
               pl    = wtemp(i,j,k,irhoE)
               if( kPresent ) pl = pl - twoThird &
                                 * wtemp(i,j,k,irho)*wtemp(i,j,k,itu1)
               Tl    = Tref*pl/(RGas*wtemp(i,j,k,irho))

               call computeGamma(Tl, PV3Gamma, 1_intType)

               v(knode) = sqrt(wtemp(i,j,k,irho)*v2 &
                        /      (PV3Gamma*wtemp(i,j,k,irhoE)))
             end do
             end do
             end do

           !=============================================================

           case(7_intPV3Type)        ! Velocity magnitude
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = sqrt(wtemp(i,j,k,ivx)**2 + wtemp(i,j,k,ivy)**2 &
                        +      wtemp(i,j,k,ivz)**2)
             end do
             end do
             end do

           !=============================================================

           case(8_intPV3Type)        ! Pressure
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,irhoE)
             end do
             end do
             end do

           !=============================================================

           case(9_intPV3Type)        ! Total pressure

             ! Store the total number of nodes in j.

             j = il*jl*kl

             call computePtot(wtemp(1,1,1,irho),  wtemp(1,1,1,ivx), &
                              wtemp(1,1,1,ivy),   wtemp(1,1,1,ivz), &
                              wtemp(1,1,1,irhoE), tmp, j)

             ! Copy the values of tmp in v; v cannot be used directly
             ! in the call to computePtot, because of different
             ! accuracy.

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode = knode + 1
               v(knode) = tmp(i,j,k)
             enddo
             enddo
             enddo

           !=============================================================

           case(10_intPV3Type)        ! Temperature
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode = knode +1
               pl    = wtemp(i,j,k,irhoE)
               if( kPresent ) pl = pl - twoThird &
                                 * wtemp(i,j,k,irho)*wtemp(i,j,k,itu1)

               v(knode) = pl/(RGas*wtemp(i,j,k,irho))
             end do
             end do
             end do

           !=============================================================

           case(11_intPV3Type)        ! Total temperature

             ! Store the total number of nodes in j.

             j = il*jl*kl

             call computeTtot(wtemp(1,1,1,irho),  wtemp(1,1,1,ivx), &
                              wtemp(1,1,1,ivy),   wtemp(1,1,1,ivz), &
                              wtemp(1,1,1,irhoE), tmp, j)

             ! Copy the values of tmp in v; v cannot be used directly
             ! in the call to computeTtot, because of different
             ! accuracy.

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode = knode + 1
               v(knode) = tmp(i,j,k)
             enddo
             enddo
             enddo

           !=============================================================

           case(12_intPV3Type)        ! Enthalpy

             ! Store the total number of nodes in j.

             j = il*jl*kl

             ! First compute the internal energy per unit mass.
             ! Make a distinction whether or not k is present.
             ! If not, just send a dummy in its place; it is not
             ! used in eintArray.

             if( kPresent ) then
               call eintArray(wtemp(1,1,1,irho), wtemp(1,1,1,irhoE), &
                              wtemp(1,1,1,itu1), tmp, kPresent, j)
             else
               dummyK = zero
               call eintArray(wtemp(1,1,1,irho), wtemp(1,1,1,irhoE), &
                              dummyK, tmp, kPresent, j)
             endif

             ! Add the p/rho term.

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = tmp(i,j,k) &
                        + wtemp(i,j,k,irhoE)/wtemp(i,j,k,irho)
             end do
             end do
             end do

           !=============================================================

           case(13_intPV3Type)        ! Total enthalpy

             ! Store the total number of nodes in j.

             j = il*jl*kl

             ! First compute the internal energy per unit mass.
             ! Make a distinction whether or not k is present.
             ! If not, just send a dummy in its place; it is not
             ! used in eintArray.

             if( kPresent ) then
               call eintArray(wtemp(1,1,1,irho), wtemp(1,1,1,irhoE), &
                              wtemp(1,1,1,itu1), tmp, kPresent, j)
             else
               dummyK = zero
               call eintArray(wtemp(1,1,1,irho), wtemp(1,1,1,irhoE), &
                              dummyK, tmp, kPresent, j)
             endif

             ! Add the p/rho term and the kinetic energy.

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v2       = wtemp(i,j,k,ivx)**2 + wtemp(i,j,k,ivy)**2 &
                        + wtemp(i,j,k,ivz)**2
               v(knode) = tmp(i,j,k) + half*v2                      &
                        + wtemp(i,j,k,irhoE)/wtemp(i,j,k,irho)
             end do
             end do
             end do

           !=============================================================

           case(14_intPV3Type)        ! Entropy. formula is probably
             do k=1,kl                ! not valid for variable gamma.
             do j=1,jl
             do i=1,il
               knode = knode +1
               pl    = wtemp(i,j,k,irhoE)
               if( kPresent ) pl = pl - twoThird &
                                 * wtemp(i,j,k,irho)*wtemp(i,j,k,itu1)
               Tl    = Tref*pl/(RGas*wtemp(i,j,k,irho))

               call computeGamma(Tl, PV3Gamma, 1_intType)

               v(knode) = log(wtemp(i,j,k,irhoE)) &
                        - PV3Gamma*log(wtemp(i,j,k,irho))
             end do
             end do
             end do

           !=============================================================

           case(15_intPV3Type)        ! Cp

             fact = two/(gammaInf*pInf*MachCoef*MachCoef)

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = fact*(wtemp(i,j,k,irhoE) -pInf)
             end do
             end do
             end do

           !=============================================================

           case(16_intPV3Type)        ! Proc#
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = myID +1
             end do
             end do
             end do

           !=============================================================

           case(17_intPV3Type)        ! Global block number
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
           !   V(knode) = real(nn,realPV3Type)
               v(knode) = real(nbkGlobal,realPV3Type)
             end do
             end do
             end do

           !=============================================================

           case(21_intPV3Type)        ! Laminar or eddy viscosity.

             if(equations == NSEquations) then

               do k=1,kl
               do j=1,jl
               do i=1,il
                 knode    = knode +1
                 v(knode) = lv(i,j,k)
               end do
               end do
               end do

             else if( eddyModel ) then

               do k=1,kl
               do j=1,jl
               do i=1,il
                 knode    = knode +1
                 v(knode) = ev(i,j,k)
               end do
               end do
               end do

             endif

           !=============================================================

           case(22_intPV3Type)   ! eddy viscosity ratio

             if( eddyModel ) then

               do k=1,kl
               do j=1,jl
               do i=1,il
                 knode    = knode +1
                 v(knode) = ev(i,j,k)/lv(i,j,k)
               end do
               end do
               end do

             endif

           !=============================================================

           case(23_intPV3Type)   ! 1st turbulent transport variable.
                                 ! Either nu tilde or k.

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,itu1)
             end do
             end do
             end do

           !=============================================================

           case(24_intPV3Type)   ! Either chi (spalart-allmaras type)
                                 ! or the 2nd transport variable.

             select case (turbModel)
               case (spalartAllmaras, spalartAllmarasEdwards)

                 do k=1,kl
                 do j=1,jl
                 do i=1,il
                   knode    = knode +1
                   v(knode) = wtemp(i,j,k,irho)*wtemp(i,j,k,itu1) &
                            / lv(i,j,k)
                 end do
                 end do
                 end do

               !=========================================================

               case default

                 do k=1,kl
                 do j=1,jl
                 do i=1,il
                   knode    = knode +1
                   v(knode) = wtemp(i,j,k,itu2)
                 end do
                 end do
                 end do

             end select

           !=============================================================

           case(25_intPV3Type)   ! The 3rd turbulent transport variable

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,itu3)
             end do
             end do
             end do

           !=============================================================

           case(26_intPV3Type)   ! The 4th turbulent transport variable

             do k=1,kl
             do j=1,jl
             do i=1,il
               knode    = knode +1
               v(knode) = wtemp(i,j,k,itu4)
             end do
             end do
             end do

         end select scalarField

         ! Deallocate the temp memory for this block.

         deallocate(wtemp, tmp, stat=ierr)
         if(ierr /= 0) call terminate("pvscal", &
                                      "Deallocation failure for wtemp &
                                      &and tmp")

         if(equations ==   NSEquations .or. &
            equations == RANSEquations) then

           deallocate(lv, stat=ierr)
           if(ierr /= 0)              &
             call terminate("pvscal", &
                            "Deallocation failure for lv")
         endif

         if( eddyModel ) then

           deallocate(ev, stat=ierr)
           if(ierr /= 0)              &
             call terminate("pvscal", &
                            "Deallocation failure for ev")
         endif

       end do domains

       return
       end subroutine pvscal

!      ==================================================================

       subroutine pvvect(key,v)
!
!      ******************************************************************
!      *                                                                *
!      * pvvect provides the values of the components of vector fields  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type)                :: key
       real(kind=realPV3Type), dimension (3,*) :: v
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intPV3Type) :: nn, i, j, k, knode

       real(kind=realType), allocatable, dimension(:,:,:) :: rho, vx
       real(kind=realType), allocatable, dimension(:,:,:) :: vy,  vz
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       knode = 0

       domains: do nn=1,nDom

         ! Set the pointers for this block. Consider only the 1st
         ! spectral solution for now.

         call setPointers(nn, groundLevel, 1_intType)

         ! Allocate the memory for the density and velocities in the
         ! nodes and compute them by averaging.

         allocate(vx(il,jl,kl),  vy(il,jl,kl), vz(il,jl,kl), &
                  rho(il,jl,kl), stat=ierr)
         if(ierr /= 0) &
           call terminate("pvvect", &
                          "Memory allocation failure for rho, vx, &
                          &vy and vz")

         do k=1,kl
         do j=1,jl
         do i=1,il
           rho(i,j,k) = eighth*(w(i,j  ,k  ,irho) + w(i+1,j  ,k  ,irho) &
                      +         w(i,j+1,k  ,irho) + w(i+1,j+1,k  ,irho) &
                      +         w(i,j  ,k+1,irho) + w(i+1,j  ,k+1,irho) &
                      +         w(i,j+1,k+1,irho) + w(i+1,j+1,k+1,irho))

           vx(i,j,k)  = eighth*(w(i,j  ,k  ,ivx)  + w(i+1,j  ,k  ,ivx) &
                      +         w(i,j+1,k  ,ivx)  + w(i+1,j+1,k  ,ivx) &
                      +         w(i,j  ,k+1,ivx)  + w(i+1,j  ,k+1,ivx) &
                      +         w(i,j+1,k+1,ivx)  + w(i+1,j+1,k+1,ivx))

           vy(i,j,k)  = eighth*(w(i,j  ,k  ,ivy)  + w(i+1,j  ,k  ,ivy) &
                      +         w(i,j+1,k  ,ivy)  + w(i+1,j+1,k  ,ivy) &
                      +         w(i,j  ,k+1,ivy)  + w(i+1,j  ,k+1,ivy) &
                      +         w(i,j+1,k+1,ivy)  + w(i+1,j+1,k+1,ivy))

           vz(i,j,k)  = eighth*(w(i,j  ,k  ,ivz)  + w(i+1,j  ,k  ,ivz) &
                      +         w(i,j+1,k  ,ivz)  + w(i+1,j+1,k  ,ivz) &
                      +         w(i,j  ,k+1,ivz)  + w(i+1,j  ,k+1,ivz) &
                      +         w(i,j+1,k+1,ivz)  + w(i+1,j+1,k+1,ivz))
         end do
         end do
         end do
!
!        ****************************************************************
!        *                                                              *
!        * Set return vector fields                                     *
!        *                                                              *
!        ****************************************************************
!
         vectorField: select case(key)

           case(19_intPV3Type)   ! Velocity vector.
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode = knode +1
               v(1,knode) = vx(i,j,k)
               v(2,knode) = vy(i,j,k)
               v(3,knode) = vz(i,j,k)
             end do
             end do
             end do

           !=============================================================

           case(20_intPV3Type)   ! Momentum vector.
             do k=1,kl
             do j=1,jl
             do i=1,il
               knode = knode +1
               v(1,knode) = rho(i,j,k)*vx(i,j,k)
               v(2,knode) = rho(i,j,k)*vy(i,j,k)
               v(3,knode) = rho(i,j,k)*vz(i,j,k)
             end do
             end do
             end do

         end select vectorField

         ! Deallocate the memory of the nodal variables again.

         deallocate(vx, vy, vz, rho, stat=ierr)
         if(ierr /= 0)              &
           call terminate("pvvect", &
                          "Deallocation failure for rho, vx, &
                          &vy and vz")

       end do domains

       return
       end subroutine pvvect

!      ==================================================================

       subroutine pvthres(key,xyz,t)
!
!      ******************************************************************
!      *                                                                *
!      * Return scalar thresholding values for the 3D field             *
!      *                                                                *
!      *   where:  key     - Index for the selected scalar function     *
!      *           xyz     - Coordinate triads for all the 3-d nodes.   *
!      *                     This is the same data as set in pvgrid.    *
!      *                                                                *
!      * Note: Key and xyz must not be modified!                        *
!      *                                                                *
!      * The following paramters must be returned by this routine:      *
!      *                                                                *
!      *       t       - the threshold value for all the 3D nodes       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type) :: key
       real(kind=realPV3Type), dimension (*)   :: t
       real(kind=realPV3Type), dimension (3,*) :: xyz
!
!      Local variables.
!
       integer(kind=intType) :: nn, i, j, k, nNode
       integer(kind=intType) :: il, jl, kl
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nNode = 0

       domains: do nn=1,nDom

         ! Abbreviate the upper boundaries a bit easier.

         il = flowDoms(nn,groundLevel,1)%il
         jl = flowDoms(nn,groundLevel,1)%jl
         kl = flowDoms(nn,groundLevel,1)%kl

         do k=1,kl
           do j=1,jl
              do i=1,il
                nNode = nNode +1
                t(nNode) = sqrt(xyz(3,nNode)*xyz(3,nNode) +  &
                                xyz(2,nNode)*xyz(2,nNode))
              end do
            end do
         end do

       end do domains

       return
       end subroutine pvthres

!      ==================================================================

       subroutine pvzprime(idCut,xyz,nNode,zp,zprime,xpc,ypc,halfw)
!
!      ******************************************************************
!      *                                                                *
!      * Return z-prime for programmer defined cutting surfaces         *
!      *                                                                *
!      *  This routine is called when the programmed cutting plane      *
!      *  function key is hit to set up the cutting surface data.       *
!      *                                                                *
!      *  xyz     - Coordinate triads for all the 3-D nodes.            *
!      *  nNode   - The number of 3-D nodes.                            *
!      *                                                                *
!      *  note: xyz and nNode must not be modified!                     *
!      *                                                                *
!      *  The following values are returned by this routine:            *
!      *                                                                *
!      *  zp      - The calculated z-prime values for all 3-D nodes.    *
!      *            This vector must be completely filled.              *
!      *  zprime  - Starting z-prime value                              *
!      *  xpc     - Starting x-prime center value for the 2-D window    *
!      *  ypc     - Starting y-prime center value for the 2-D window    *
!      *  halfw   - Starting 1/2 width for the 2-D window in            *
!      *            x-prime/y-prime space.                              *
!      *                                                                *
!      ******************************************************************
!
       use PV3state
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type) :: idCut, nNode
       real(kind=realPV3Type)                  :: zprime,xpc,ypc,halfw
       real(kind=realPV3Type), dimension (*)   :: zp
       real(kind=realPV3Type), dimension (3,*) :: xyz
!
!      Local variables.
!
       integer(kind=intPV3Type) :: i
       real(kind=realPV3Type)   :: xAve,yAve,zAve,rAve,tAve
       real(kind=realPV3Type)   :: xDev,yDev,zDev,rDev,tDev
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      * Note that this routine has been coded with turbo machinery     *
!      *  geometries in mind (for r and theta cutting planes) and       *
!      *  assumes that in those cases, the axis of the turbo machine is *
!      *  along the x-axis                                              *
!      *                                                                *
!      ******************************************************************
!
!      Store current value of requested programmed cut so that pvxyprime
!       may use it to decided on the proper mapping
!
       cutSurface = idCut

!      Only return zprime, xpc, ypc, halfw if idCut > 0.

       if (idCut > 0) then

          ! Calculate averages for x, y, z, r, and theta

          rAve = 0.0_realPV3Type
          tAve = 0.0_realPV3Type
          xAve = 0.0_realPV3Type
          yAve = 0.0_realPV3Type
          zAve = 0.0_realPV3Type

          do i=1,nNode
             rAve = rAve + sqrt(xyz(2,i)**2 + xyz(3,i)**2)
             tAve = tAve + atan2(xyz(2,i), xyz(3,i))
             xAve = xAve + xyz(1,i)
             yAve = yAve + xyz(2,i)
             zAve = zAve + xyz(3,i)
          end do

          rAve = rAve/real(nNode,realPV3Type)
          tAve = tAve/real(nNode,realPV3Type)
          xAve = xAve/real(nNode,realPV3Type)
          yAve = yAve/real(nNode,realPV3Type)
          zAve = zAve/real(nNode,realPV3Type)

          ! Calculate standard deviations for x, y, z, r, and theta

          rDev = 0.0_realPV3Type
          tDev = 0.0_realPV3Type
          xDev = 0.0_realPV3Type
          yDev = 0.0_realPV3Type
          zDev = 0.0_realPV3Type

          do i=1,nNode
             rDev = rDev + abs(sqrt(xyz(2,i)**2 +xyz(3,i)**2) -rAve)
             tDev = tDev + abs(atan2(xyz(2,i),xyz(3,i))-tAve)
             xDev = xDev + abs(xyz(1,i)-xAve)
             yDev = yDev + abs(xyz(2,i)-yAve)
             zDev = zDev + abs(xyz(3,i)-zAve)
          end do

          rDev = rDev/real(nNode,realPV3Type)
          tDev = tDev/real(nNode,realPV3Type)
          xDev = xDev/real(nNode,realPV3Type)
          yDev = yDev/real(nNode,realPV3Type)
          zDev = zDev/real(nNode,realPV3Type)
 
       endif
!
!      ******************************************************************
!      *                                                                *
!      *  Setup coordinates for the various cutting surfaces.  PV3 will *
!      *   create a cutting surface along isosurfaces of the scalar     *
!      *   fields provided here (x, y, z, r, theta)                     *
!      *                                                                *
!      ******************************************************************
!
       select case (idCut)

        case (1,-1)

           ! X - plane

           do i = 1, nNode
              zp(i)  = xyz(1,i)
           end do

           if (idCut > 0) then
              zprime = xAve
              xpc    = yAve
              ypc    = zAve
              halfw  = yDev +zDev
           end if

        case (2,-2)

           ! Y - plane

           do i = 1, nNode
              zp(i)  = xyz(2,i)
           end do

           if (idCut > 0) then
              zprime = yAve
              xpc    = xAve
              ypc    = zAve
              halfw  = xDev +zDev
           end if

        case (3,-3)

           ! Z - plane

           do i = 1, nNode
              zp(i)  = xyz(3,i)
           end do

           if (idCut > 0) then
              zprime = zAve
              xpc    = xAve
              ypc    = yAve
              halfw  = xDev +yDev
           end if

        case (4,-4)

           ! R - plane

           do i = 1, nNode
              zp(i)  = sqrt(xyz(2,i)**2 +xyz(3,i)**2)
           end do

           if (idCut > 0) then
              zprime = rAve
              xpc    = xAve
              ypc    = rAve*tAve
              halfw  = rAve*tDev +xDev
           end if

        case (5,-5)

           ! Theta - plane

           do i = 1, nNode
              zp(i) = atan2(xyz(2,i),xyz(3,i))
           end do

           if (idCut > 0) then
              zprime = tAve
              xpc    = xAve
              ypc    = rAve
              halfw  = xDev +rDev
           endif

       end select

       return
       end subroutine pvzprime

!      ==================================================================

       subroutine pvxyprime(zprime, kn, xyz, n, xyp)
!
!      ******************************************************************
!      *                                                                *
!      * Return x-prime and y-prime for programmer defined cutting      *
!      *  surfaces                                                      *
!      *                                                                *
!      * This routine is called during the data-collection phase        *
!      *  of the program defined surfaces.  This routine may be         *
!      *  called many times during this phase.                          *
!      *  There is no ordering of the nodes in kn, infact kn            *
!      *  may contain the same node more than once in a single call.    *
!      *                                                                *
!      * where:  zprime  - The current z-prime value                    *
!      *         kn      - A vector of indices for the nodes to be      *
!      *                   transformed                                  *
!      *         xyz     - Coordinate triads for 3-d nodes.             *
!      *         n       - The number of nodes to be transformed.       *
!      *                                                                *
!      *         Note: zprime, kn, xyz and n must not be modified!      *
!      *                                                                *
!      *  The following array is returned by this routine:              *
!      *                                                                *
!      *         xyp     - Calculated x-prime and y-prime values for    *
!      *                   the 3-D nodes                                *
!      *                                                                *
!      ******************************************************************
!
       use PV3state
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intPV3Type)                :: n
       integer(kind=intPV3Type), dimension(n)  :: kn
       real(kind=realPV3Type)                  :: zprime
       real(kind=realPV3Type), dimension (2,*) :: xyp
       real(kind=realPV3Type), dimension (3,*) :: xyz
!
!      Local variables.
!
       integer(kind=intPV3Type) :: i, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       select case (cutSurface)

       case (1,-1)

          ! x - plane

          do i = 1, n
             k = kn(i)
             xyp(1,i) = xyz(2,k)
             xyp(2,i) = xyz(3,k)
          end do

       case (2,-2)

          ! y - plane

          do i = 1, n
             k = kn(i)
             xyp(1,i) = xyz(1,k)
             xyp(2,i) = xyz(3,k)
          end do

       case (3,-3)

          ! z - plane

          do i = 1, n
             k = kn(i)
             xyp(1,i) = xyz(1,k)
             xyp(2,i) = xyz(2,k)
          end do

       case (4,-4)

          ! R - plane

          do i = 1, n
             k = kn(i)
             xyp(1,i) = xyz(1,k)
             xyp(2,i) = zprime*atan2(xyz(2,k),xyz(3,k))
          end do

       case (5,-5)

          ! Theta - plane

          do i = 1, n
             k = kn(i)
             xyp(1,i) = xyz(1,k)
             xyp(2,i) = sqrt(xyz(2,k)**2 +xyz(3,k)**2)
          end do

       end select

      return
      end subroutine pvxyprime
