! @File    :   cudaBlock.f90
! @Desc    :   This module contains the block data type and provides subroutines for copying data between the CPU and GPU
!              For variable names, follow the same convention as 'block.F90'

module cudaBlock
    use constants, only: realType, intType, porType, maxCGNSNameLen
    implicit none
    save
    ! =========================================================
    !           Types
    ! =========================================================
    type cudaBCDataType
        ! integer(kind=intType), pointer, device :: inBeg, inEnd, jnBeg, jnEnd
        ! integer(kind=intType), pointer, device :: icBeg, icEnd, jcBeg, jcEnd

        real(kind=realType), dimension(:, :, :), pointer, device :: norm
        real(kind=realType), dimension(:, :), pointer, device :: rface
        real(kind=realType), dimension(:, :, :), pointer, device :: F, Fv, Fp
        real(kind=realType), dimension(:, :, :), pointer, device :: T, Tv, Tp
        real(kind=realType), dimension(:, :), pointer, device :: area
        real(kind=realType), dimension(:, :), pointer, device :: CpTarget
        integer(kind=intType), dimension(:, :), pointer, device :: surfIndex

        real(kind=realType), dimension(:, :), pointer, device :: nodeVal
        real(kind=realType), dimension(:, :), pointer, device :: cellVal
        real(kind=realType), dimension(:), pointer, device :: symNorm

        logical :: symNormSet

        integer(kind=intType), pointer, device :: subsonicInletTreatment

        real(kind=realType), dimension(:, :, :), pointer, device :: uSlip
        real(kind=realType), dimension(:, :), pointer, device :: TNS_Wall

        character(maxCGNSNameLen), pointer, device :: family
        integer(kind=intType), pointer, device :: famID

        real(kind=realType), dimension(:, :, :, :), pointer, device :: normALE
        real(kind=realType), dimension(:, :, :), pointer, device :: rFaceALE
        real(kind=realType), dimension(:, :, :, :), pointer, device :: uSlipALE
        real(kind=realType), dimension(:, :), pointer, device :: nodeHeatFlux
        real(kind=realType), dimension(:, :), pointer, device :: cellHeatFlux

        real(kind=realType), dimension(:, :), pointer, device :: ptInlet, ttInlet, htInlet
        real(kind=realType), dimension(:, :), pointer, device :: flowXDirInlet, flowYDirInlet, flowZDirInlet

        real(kind=realType), dimension(:, :, :), pointer, device :: turbInlet

        real(kind=realType), dimension(:, :), pointer, device :: rho
        real(kind=realType), dimension(:, :), pointer, device :: velX, velY, velZ
        real(kind=realType), dimension(:, :), pointer, device :: ps

        integer(kind=intType), dimension(:, :), pointer, device :: iblank

        real(Kind=realType), dimension(:, :), pointer, device :: delta
        real(Kind=realType), dimension(:, :), pointer, device :: deltaNode
    end type cudaBCDataType

    type cudablockType
        integer(kind=intType) :: nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb

        logical :: rightHanded

        ! Current indices into the original block
        integer(kind=intType) :: ii, jj, kk

        ! Double halos
        real(kind=realType), dimension(:,:,:,:),pointer,device :: w
        real(kind=realType), dimension(:,:,:),pointer,device :: P,gamma,ss

        ! Single halos
        real(kind=realType), dimension(:,:,:,:), pointer,device :: x
        real(kind=realType), dimension(:,:,:),device,pointer :: rlv, rev, vol, aa
        real(kind=realType), dimension(:,:,:),device,pointer :: radI, radJ, radK, dtl
        real(kind=realType), dimension(:,:,:, :),device,pointer :: dss ! Shock sensor

        ! No halos
        real(kind=realType), dimension(:,:,:),pointer,device :: volRef, d2wall
        integer(kind=intType), dimension(:,:,:),pointer,device :: iblank

        ! Face Porosities
        integer(kind=porType), dimension(:,:,:),pointer,device :: porI
        integer(kind=porType), dimension(:,:,:),pointer,device :: porJ
        integer(kind=porType), dimension(:,:,:),pointer,device :: porK

        ! Single halos (only owned cells significant)
        real(kind=realType), dimension(:,:,:, :),pointer,device :: fw
        real(kind=realType), dimension(:,:,:, :),pointer,device :: dw

        ! Face projected areas
        real(kind=realType), dimension(:,:,:, :),pointer,device :: sI
        real(kind=realType), dimension(:,:,:, :),pointer,device :: sJ
        real(kind=realType), dimension(:,:,:, :),pointer,device :: sK

        ! Face velocities
        real(kind=realType), dimension(:,:,:), pointer,device :: sFaceI
        real(kind=realType), dimension(:,:,:), pointer,device :: sFaceJ
        real(kind=realType), dimension(:,:,:), pointer,device :: sFaceK

        ! Nodal gradients
        real(kind=realType), dimension(:,:,:),pointer,device :: ux, uy, uz
        real(kind=realType), dimension(:,:,:),pointer,device :: vx, vy, vz
        real(kind=realType), dimension(:,:,:),pointer,device :: wx, wy, wz
        real(kind=realType), dimension(:,:,:),pointer,device :: qx, qy, qz

        ! --- BCData ---
        integer(kind=intType), pointer, device :: nSubface, n1to1, nBocos, nViscBocos

        integer(kind=intType), dimension(:), pointer, device :: BCType
        integer(kind=intType), dimension(:), pointer, device :: BCFaceID

        integer(kind=intType), dimension(:), pointer, device :: cgnsSubface

        integer(kind=intType), dimension(:), pointer, device :: inBeg, inEnd
        integer(kind=intType), dimension(:), pointer, device :: jnBeg, jnEnd
        integer(kind=intType), dimension(:), pointer, device :: knBeg, knEnd

        integer(kind=intType), dimension(:), pointer, device :: dinBeg, dinEnd
        integer(kind=intType), dimension(:), pointer, device :: djnBeg, djnEnd
        integer(kind=intType), dimension(:), pointer, device :: dknBeg, dknEnd

        integer(kind=intType), dimension(:), pointer, device :: icBeg, icEnd
        integer(kind=intType), dimension(:), pointer, device :: jcBeg, jcEnd
        integer(kind=intType), dimension(:), pointer, device :: kcBeg, kcEnd

        integer(kind=intType), dimension(:), pointer, device :: neighBlock
        integer(kind=intType), dimension(:), pointer, device :: neighProc
        integer(kind=intType), dimension(:), pointer, device :: l1, l2, l3
        integer(kind=intType), dimension(:), pointer, device :: groupNum

    end type cudaBlockType

    ! Device (GPU) and host (CPU) domains
    type(cudaBlockType), device, allocatable,dimension(:,:) :: d_cudaDoms(:,:)
    type(cudaBlockType), allocatable,dimension(:,:) :: h_cudaDoms(:,:)
    type(cudaBCDataType), device, allocatable,dimension(:) :: d_BCData(:)
    type(cudaBCDataType), allocatable,dimension(:) :: h_BCData(:)

    contains 
    
    subroutine copyCudaBlock
        ! subroutine to allocate cuda domains and copy data from host to device
        ! probably use this stack exchange post to help allocate
        ! https://stackoverflow.com/questions/44680150/how-to-allocate-arrays-of-arrays-in-structure-with-cuda-fortran
        use constants, only: zero, RANSEquations, EulerEquations, Symm, SymmPolar, NSWallAdiabatic, FarField
        use blockPointers, only: nDom, &
            bnx=>nx, bny => ny, bnz => nz, &
            bil => il, bjl => jl, bkl => kl, &
            bie => ie, bje => je, bke => ke, &
            bib => ib, bjb => jb, bkb => kb, &
            bw => w, bp => p, bgamma => gamma, &
            bradi => radi, bradj => radj, bradk => radk, &
            bux => ux, buy => uy, buz => uz, &
            bvx => vx, bvy => vy, bvz => vz, &
            bwx => wx, bwy => wy, bwz => wz, &
            bqx => qx, bqy => qy, bqz => qz, &
            bx => x, brlv => rlv, brev => rev, bvol => vol, bVolRef => volRef, bd2wall => d2wall, &
            biblank => iblank, bPorI => porI, bPorJ => porJ, bPorK => porK, bdw => dw, bfw => fw, &
            bShockSensor => shockSensor, &
            bsi => si, bsj => sj, bsk => sk, &
            bsFaceI => sFaceI, bsFaceJ => sFaceJ, bsFaceK => sFaceK, &
            bdtl => dtl, baa => aa, &
            addGridVelocities, brightHanded => rightHanded, &
            ! --- BC ---
            bBCData => BCData, &
            ! Stuff that wee need but aren't renaming in this module
            nBocos, BCType
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowVarRefState, only: nw,nwf,nwt
        use inputPhysics, only: equations
        use utils, only: setPointers
        
        implicit none
        
        integer(kind=intType) :: nn, sps, mm    ! domain and time spectral interval indices
        integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iNodeBeg, iNodeEnd, jNodeBeg, jNodeEnd


        ! Allocate device 
        allocate(d_cudaDoms(nDom, nTimeIntervalsSpectral))
        ! Allocate host
        allocate(h_cudaDoms(nDom, nTimeIntervalsSpectral))

        
        domainsLoop: do nn = 1, nDom ! loop over domains
            spectralLoop:do sps = 1, nTimeIntervalsSpectral ! loop over time spectral intervals
                

                call setPointers(nn, 1, sps)
                allocate(d_BCData(nBocos))
                allocate(h_BCData(nBocos))

                ! Indices
                h_cudaDoms(nn, sps)%nx = bnx
                h_cudaDoms(nn, sps)%ny = bny
                h_cudaDoms(nn, sps)%nz = bnz
                h_cudaDoms(nn, sps)%il = bil
                h_cudaDoms(nn, sps)%jl = bjl
                h_cudaDoms(nn, sps)%kl = bkl
                h_cudaDoms(nn, sps)%ie = bie
                h_cudaDoms(nn, sps)%je = bje
                h_cudaDoms(nn, sps)%ke = bke
                h_cudaDoms(nn, sps)%ib = bib
                h_cudaDoms(nn, sps)%jb = bjb
                h_cudaDoms(nn, sps)%kb = bkb

                h_cudaDoms(nn, sps)%rightHanded = brightHanded

                ! allocate double halos
                allocate(h_cudaDoms(nn,sps)%w(0:bib, 0:bjb, 0:bkb, 1:nw)) ! flowfield states
                allocate(h_cudaDoms(nn,sps)%p(0:bib, 0:bjb, 0:bkb)) ! static pressures
                allocate(h_cudaDoms(nn,sps)%gamma(0:bib, 0:bjb, 0:bkb)) ! spec heat ratio
                allocate(h_cudaDoms(nn,sps)%ss(0:bib, 0:bjb, 0:bkb)) ! shock sensors
                ! and then set them
                h_cudaDoms(nn, sps)%w = bw
                h_cudaDoms(nn, sps)%p = bp
                h_cudaDoms(nn, sps)%gamma = bgamma
                h_cudaDoms(nn, sps)%ss = bShockSensor

                ! allocate single halos
                allocate(h_cudaDoms(nn, sps)%x(0:bie, 0:bje, 0:bke, 1:3))
                allocate(h_cudaDoms(nn, sps)%rlv(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%rev(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%vol(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%aa(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%radI(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%radJ(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%radK(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%dtl(1:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%dss(1:bie, 1:bje, 1:bke, 1:3))
                ! and set them
                ! NOTE: these need to be computed 
                h_cudaDoms(nn, sps)%x = bx
                h_cudaDoms(nn, sps)%rlv = brlv
                h_cudaDoms(nn, sps)%rev = brev
                h_cudaDoms(nn, sps)%vol = bvol(1:bie, 1:bje, 1:bke)
                ! initialize to zero
                h_cudaDoms(nn, sps)%aa = zero
                h_cudaDoms(nn, sps)%radI = zero
                h_cudaDoms(nn, sps)%radJ = zero
                h_cudaDoms(nn, sps)%radK = zero
                h_cudaDoms(nn, sps)%dtl = zero

                h_cudaDoms(nn, sps)%dss = zero

                ! no halos
                if (equations == RANSEquations) then 
                    allocate(h_cudaDoms(nn, sps)%d2wall(2:bil,2:bjl,2:bkl))
                    h_cudaDoms(nn, sps)%d2wall = bd2wall
                end if 
                allocate(h_cudaDoms(nn, sps)%volRef(2:bil,2:bjl,2:bkl))
                allocate(h_cudaDoms(nn, sps)%iblank(2:bil,2:bjl,2:bkl))
                h_cudaDoms(nn, sps)%volRef = bvolRef
                h_cudaDoms(nn, sps)%iblank  = biblank

                ! face porosities
                allocate(h_cudaDoms(nn, sps)%porI(1:bil,2:bjl,2:bkl))
                allocate(h_cudaDoms(nn, sps)%porJ(2:bil,1:bjl,2:bkl))
                allocate(h_cudaDoms(nn, sps)%porK(2:bil,2:bjl,1:bkl))
                h_cudaDoms(nn,sps)%porI = bPorI
                h_cudaDoms(nn,sps)%porJ = bPorJ
                h_cudaDoms(nn,sps)%porK = bPorK

                ! single halos (only owned cells significant)
                allocate(h_cudaDoms(nn, sps)%fw(1:bie, 1:bje, 1:bke, 1:nwf))
                allocate(h_cudaDoms(nn, sps)%dw(1:bie, 1:bje, 1:bke, 1:nw))
                h_cudaDoms(nn, sps)%fw = zero
                h_cudaDoms(nn, sps)%dw = zero

                ! face projected areas
                allocate(h_cudaDoms(nn,sps)%sI(0:bie, 1:bje, 1:bke,3))
                allocate(h_cudaDoms(nn,sps)%sJ(1:bie, 0:bje, 1:bke,3))
                allocate(h_cudaDoms(nn,sps)%sK(1:bie, 1:bje, 0:bke,3))
                h_cudaDoms(nn, sps)%sI = bsi
                h_cudaDoms(nn, sps)%sJ = bsj
                h_cudaDoms(nn, sps)%sK = bsk

                ! Face velocities
                allocate(h_cudaDoms(nn, sps)%sFaceI(0:bie, 1:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%sFaceJ(1:bie, 0:bje, 1:bke))
                allocate(h_cudaDoms(nn, sps)%sFaceK(1:bie, 1:bje, 0:bke))
                if (addGridVelocities) then
                    h_cudaDoms(nn, sps)%sFaceI = bsFaceI
                    h_cudaDoms(nn, sps)%sFaceJ = bsFaceJ
                    h_cudaDoms(nn, sps)%sFaceK = bsFaceK
                else
                    h_cudaDoms(nn, sps)%sFaceI = zero
                    h_cudaDoms(nn, sps)%sFaceJ = zero
                    h_cudaDoms(nn, sps)%sFaceK = zero
                end if

                ! Nodal gradients
                allocate(h_cudaDoms(nn, sps)%ux(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%uy(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%uz(1:bil, 1:bjl, 1:bkl))
                h_cudaDoms(nn, sps)%ux = zero; h_cudaDoms(nn, sps)%uy = zero; h_cudaDoms(nn, sps)%uz = zero

                allocate(h_cudaDoms(nn, sps)%vx(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%vy(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%vz(1:bil, 1:bjl, 1:bkl))
                h_cudaDoms(nn, sps)%vx = zero; h_cudaDoms(nn, sps)%vy = zero; h_cudaDoms(nn, sps)%vz = zero

                allocate(h_cudaDoms(nn, sps)%wx(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%wy(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%wz(1:bil, 1:bjl, 1:bkl))
                h_cudaDoms(nn, sps)%wx = zero; h_cudaDoms(nn, sps)%wy = zero; h_cudaDoms(nn, sps)%wz = zero

                allocate(h_cudaDoms(nn, sps)%qx(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%qy(1:bil, 1:bjl, 1:bkl))
                allocate(h_cudaDoms(nn, sps)%qz(1:bil, 1:bjl, 1:bkl))
                h_cudaDoms(nn, sps)%qx = zero; h_cudaDoms(nn, sps)%qy = zero; h_cudaDoms(nn, sps)%qz = zero

                ! ! --- BCData allocation within host ---
                ! ! cudaDoms is an array of blocks
                ! ! BC data is a pointer of type cudaBCDataType
                ! ! each block has this
                
                ! The following is similar to the allocMemBCData() subroutine in BCData.F90
                bocoLoop: do mm = 1, nBocos
                
                    iBeg = bBCData(mm)%icbeg; iEnd = bBCData(mm)%icend
                    jBeg = bBCData(mm)%jcbeg; jEnd = bBCData(mm)%jcend

                    inodeBeg = bBCData(mm)%inbeg; inodeEnd = bBCData(mm)%inend
                    jnodeBeg = bBCData(mm)%jnbeg; jnodeEnd = bBCData(mm)%jnend

                    select case (BCType(mm))
                        
                    ! case (NSWallAdiabatic)
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%uSlip(iBeg:iEnd, jBeg:jEnd, 3))
                    !     ! allocate(h_cudaDoms(nn, sps)%BCData(mm)%uSlipALE(0:nALEsteps, iBeg:iEnd, jBeg:jEnd, 3))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%F(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd, 3))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%T(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd, 3))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%Tp(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd, 3))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%Tv(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd, 3))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%Fp(iNodeBeg + 1:iNodeEnd, jNodeBeg + 1:jNodeEnd, 3))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%Fv(iNodeBeg + 1:iNodeEnd, jNodeBeg + 1:jNodeEnd, 3))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%area(iNodeBeg + 1:iNodeEnd, jNodeBeg + 1:jNodeEnd))
                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%CpTarget(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd))
                    
                    !     h_cudaDoms(nn, sps)%BCData(mm)%uslip = bBCData(mm)%uslip
                    !     ! h_cudaDoms(nn, sps)%BCData(mm)%uSlipALE = bBCData(mm)%uSlipALE
                    !     h_cudaDoms(nn, sps)%BCData(mm)%F = bBCData(mm)%F
                    !     h_cudaDoms(nn, sps)%BCData(mm)%Fv = bBCData(mm)%Fv
                    !     h_cudaDoms(nn, sps)%BCData(mm)%Fp = bBCData(mm)%Tp
                    !     h_cudaDoms(nn, sps)%BCData(mm)%T = bBCData(mm)%T
                    !     h_cudaDoms(nn, sps)%BCData(mm)%Tv = bBCData(mm)%Tv
                    !     h_cudaDoms(nn, sps)%BCData(mm)%Tp = bBCData(mm)%Tp
                    !     h_cudaDoms(nn, sps)%BCData(mm)%area = bBCData(mm)%area
                    !     h_cudaDoms(nn, sps)%BCData(mm)%CpTarget = zero
                        
                    ! case (farField)

                    !     allocate(h_cudaDoms(nn, sps)%BCData(mm)%rface(iBeg:iEnd, jBeg:jEnd))
                    !     ! allocate(h_cudaDoms(nn, sps)%BCData(mm)%rfaceALE(0:nALEsteps, iBeg:iEnd, jBeg:jEnd))
                    !     h_cudaDoms(nn, sps)%BCData(mm)%rface = bBCData(mm)%rface
                    !     ! h_cudaDoms(nn, sps)%BCData(mm)%rfaceALE = bBCData(mm)%rfaceALE
                    
                    case (symm, symmPolar)

                        allocate(h_BCData(mm)%rface(iBeg:iEnd, jBeg:jEnd))
                        ! allocate(h_BCData(mm)%rfaceALE(0:nALEsteps, iBeg:iEnd, jBeg:jEnd))
                        h_BCData(mm)%rface = bBCData(mm)%rface
                        ! h_BCData(mm)%rfaceALE = bBCData(mm)%rfaceALE

                    ! case (SupersonicInflow, DomainInterfaceAll) ! TODO:
                    ! case (SupersonicOutflow)
                    ! case (SupersonicInflow)
                    ! case (SupersonicOutflow, MassBleedOutflow, DomainInterfaceP)
                    ! case (DomainInterfaceRhoUVW)
                    ! case (DomainInterfaceTotal)
                    ! case (domainInterfaceRho)

                    end select

                end do bocoLoop
                ! TODO GN: pickup here with debugging and allocating and setting vars btwn cpu and gpu

            end do spectralLoop
        end do domainsLoop

        ! --- Copy data from host back to device---
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                d_cudaDoms(nn, sps) = h_cudaDoms(nn, sps)
            end do
        end do

    end subroutine copyCudaBlock
end module cudaBlock