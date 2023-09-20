! @File    :   cudaBlock.f90
! @Desc    :   This module contains the block data type and provides subroutines for copying data between the CPU and GPU
!              For variable names, follow the same convention as 'block.F90'

module cudaBlock
    use constants, only: realType, intType,porType
    implicit none
    save

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
    end type cudaBlockType

    ! Device (GPU) and host (CPU) domains
    type(cudaBlockType), device, allocatable,dimension(:,:) :: d_cudaDoms(:,:)
    type(cudaBlockType), allocatable,dimension(:,:) :: h_cudaDoms(:,:)

    contains 
    
    subroutine copyCudaBlock
        ! subroutine to allocate cuda domains and copy data from host to device
        ! probably use this stack exchange post to help allocate
        ! https://stackoverflow.com/questions/44680150/how-to-allocate-arrays-of-arrays-in-structure-with-cuda-fortran
        use constants, only: zero, RANSEquations, EulerEquations
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
            addGridVelocities, brightHanded => rightHanded
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowVarRefState, only: nw,nwf,nwt
        use inputPhysics, only: equations
        use utils, only: setPointers
        
        implicit none
        
        integer(kind=intType) :: nn, sps    ! domain and time spectral interval indices

        ! Allocate device 
        allocate(d_cudaDoms(nDom, nTimeIntervalsSpectral))
        ! Allocate host
        allocate(h_cudaDoms(nDom, nTimeIntervalsSpectral))
        
        do nn = 1, nDom ! loop over domains
            do sps = 1, nTimeIntervalsSpectral ! loop over time spectral intervals
                
                call setPointers(nn,1, sps)
                
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


            end do
        end do

        ! --- Copy data from host back to device---
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                d_cudaDoms(nn, sps) = h_cudaDoms(nn, sps)
            end do
        end do

    end subroutine copyCudaBlock
end module cudaBlock