! @File    :   cudaResidual.f90
! @Desc    :   This module contains the routines for residual evaluation using CUDA FORTRAN
!              They generally have the same subroutine names as the CPU version

module cudaResidual
    ! --- ADflow modules to use ---
    use precision, only: realType, intType
    use cudaBlock, only: cudaDoms=>d_cudaDoms, h_cudaDoms
    use constants
    use cudaBCPointers
    ! --- CUDA FORTRAN module ---
    use cudafor
    implicit none

    ! Actual dimensions to execute
    !  nx, ny, nz - Block integer dimensions for no halo cell based
    !               quantities.
    !  il, jl, kl - Block integer dimensions for no halo node based
    !               quantities.
    !  ie, je, ke - Block integer dimensions for single halo
    !               cell-centered quantities.
    !  ib, jb, kb - Block integer dimensions for double halo
    !               cell-centered quantities.
    integer(kind=intType) :: h_nx, h_ny, h_nz, h_il, h_jl, h_kl, h_ie, h_je, h_ke, h_ib, h_jb, h_kb
    integer(kind=intType),constant :: nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb

    logical, device :: rightHanded

    ! Current indices into the original block
    integer(kind=intType) :: ii, jj, kk

    ! Double halos
    real(kind=realType), dimension(:,:,:,:),allocatable,device :: w
    real(kind=realType), dimension(:,:,:),allocatable,device :: P,gamma,ss

    ! Single halos
    real(kind=realType), dimension(:,:,:,:), allocatable,device :: x
    real(kind=realType), dimension(:,:,:),device,allocatable :: rlv, rev, vol, aa
    real(kind=realType), dimension(:,:,:),device,allocatable :: radI, radJ, radK, dtl
    real(kind=realType), dimension(:,:,:, :),device,allocatable :: dss ! Shock sensor

    ! No halos
    real(kind=realType), dimension(:,:,:),allocatable,device :: volRef, d2wall
    integer(kind=intType), dimension(:,:,:),allocatable,device :: iblank

    ! Face Porosities
    integer(kind=porType), dimension(:,:,:),allocatable,device :: porI
    integer(kind=porType), dimension(:,:,:),allocatable,device :: porJ
    integer(kind=porType), dimension(:,:,:),allocatable,device :: porK

    ! Single halos (only owned cells significant)
    real(kind=realType), dimension(:,:,:, :),allocatable,device :: fw
    real(kind=realType), dimension(:,:,:, :),allocatable,device :: dw

    ! Face projected areas
    real(kind=realType), dimension(:,:,:, :),allocatable,device :: sI
    real(kind=realType), dimension(:,:,:, :),allocatable,device :: sJ
    real(kind=realType), dimension(:,:,:, :),allocatable,device :: sK

    ! Face velocities
    real(kind=realType), dimension(:,:,:), allocatable,device :: sFaceI
    real(kind=realType), dimension(:,:,:), allocatable,device :: sFaceJ
    real(kind=realType), dimension(:,:,:), allocatable,device :: sFaceK

    ! Nodal gradients
    real(kind=realType), dimension(:,:,:),allocatable,device :: ux, uy, uz
    real(kind=realType), dimension(:,:,:),allocatable,device :: vx, vy, vz
    real(kind=realType), dimension(:,:,:),allocatable,device :: wx, wy, wz
    real(kind=realType), dimension(:,:,:),allocatable,device :: qx, qy, qz

    contains

    subroutine copydata
        use constants, only: zero, RANSEquations, EulerEquations
        use blockPointers, only: &
            bnx => nx, bny => ny, bnz => nz, &
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
        use flowVarRefState, only: nw,nwf,nwt
        use inputPhysics, only: equations
        implicit none

            ! Copy block pointers dimensions
            h_nx = bnx; 
            h_ny = bny; 
            h_nz = bnz
            h_il = bil; 
            h_jl = bjl; 
            h_kl = bkl
            h_ie = bie; 
            h_je = bje; 
            h_ke = bke
            h_ib = bib; 
            h_jb = bjb; 
            h_kb = bkb

            nx = bnx; ny = bny; nz = bnz
            il = bil; jl = bjl; kl = bkl
            ie = bie; je = bje; ke = bke
            ib = bib; jb = bjb; kb = bkb
            rightHanded = brightHanded

            !allocate double halos
            allocate(w(0:h_ib,0:h_jb,0:h_kb,1:nw))
            allocate(p(0:h_ib,0:h_jb,0:h_kb))
            allocate(gamma(0:h_ib,0:h_jb,0:h_kb))
            allocate(ss(0:h_ib,0:h_jb,0:h_kb))

            w = bw; p = bp; gamma = bgamma; ss = bShockSensor

            !allocate single halos
            allocate(x(0:h_ie,0:h_je,0:h_ke,1:3))
            allocate(rlv(1:h_ie,1:h_je,1:h_ke))
            allocate(rev(1:h_ie,1:h_je,1:h_ke))
            allocate(vol(1:h_ie,1:h_je,1:h_ke))
            allocate(aa(1:h_ie,1:h_je,1:h_ke))
            allocate(radI(1:h_ie,1:h_je,1:h_ke))
            allocate(radJ(1:h_ie,1:h_je,1:h_ke))
            allocate(radK(1:h_ie,1:h_je,1:h_ke))
            allocate(dtl(1:h_ie,1:h_je,1:h_ke))
            allocate(dss(1:h_ie,1:h_je,1:h_ke,1:3))

            x = bx
            rlv = brlv
            rev = brev
            vol = bvol(1:bie, 1:bje, 1:bke)
            !these need to be comptued 

            aa = zero
            radI = zero
            radJ = zero
            radK = zero
            dtl = zero

            !set this to zero
            dss = zero

            !no halos
            allocate(volRef(2:h_il,2:h_jl,2:h_kl))
            if (equations == RANSEquations) then 
                allocate(d2wall(2:h_il,2:h_jl,2:h_kl))
                d2wall = bd2wall
            end if 
            allocate(iblank(2:h_il,2:h_jl,2:h_kl))
            volRef = bvolRef
            iblank  = biblank

            !face porosities
            allocate(porI(1:h_il,2:h_jl,2:h_kl))
            allocate(porJ(2:h_il,1:h_jl,2:h_kl))
            allocate(porK(2:h_il,2:h_jl,1:h_kl))
            porI = bPorI
            porJ = bPorJ
            porK = bPorK

            !single halos (only owned cells significant)
            allocate(fw(1:h_ie,1:h_je,1:h_ke,1:nwf))
            allocate(dw(1:h_ie,1:h_je,1:h_ke,1:nw))
            fw = zero
            dw = zero

            !face projected areas
            allocate(sI(0:h_ie, 1:h_je, 1:h_ke,3))
            allocate(sJ(1:h_ie, 0:h_je, 1:h_ke,3))
            allocate(sK(1:h_ie, 1:h_je, 0:h_ke,3))
            sI = bsi
            sJ = bsj
            sK = bsk
            ! Face velocities
            allocate(sFaceI(0:h_ie, 1:h_je, 1:h_ke))
            allocate(sFaceJ(1:h_ie, 0:h_je, 1:h_ke))
            allocate(sFaceK(1:h_ie, 1:h_je, 0:h_ke))
            if (addGridVelocities) then
                sFaceI = bsFaceI
                sFaceJ = bsFaceJ
                sFaceK = bsFaceK
            else
                sFaceI = zero
                sFaceJ = zero
                sFaceK = zero

            end if
            ! Nodal gradients
            allocate(ux(1:h_il, 1:h_jl, 1:h_kl))
            allocate(uy(1:h_il, 1:h_jl, 1:h_kl))
            allocate(uz(1:h_il, 1:h_jl, 1:h_kl))
            ux = zero; uy = zero; uz = zero

            allocate(vx(1:h_il, 1:h_jl, 1:h_kl))
            allocate(vy(1:h_il, 1:h_jl, 1:h_kl))
            allocate(vz(1:h_il, 1:h_jl, 1:h_kl))
            vx = zero; vy = zero; vz = zero

            allocate(wx(1:h_il, 1:h_jl, 1:h_kl))
            allocate(wy(1:h_il, 1:h_jl, 1:h_kl))
            allocate(wz(1:h_il, 1:h_jl, 1:h_kl))
            wx = zero; wy = zero; wz = zero

            allocate(qx(1:h_il, 1:h_jl, 1:h_kl))
            allocate(qy(1:h_il, 1:h_jl, 1:h_kl))
            allocate(qz(1:h_il, 1:h_jl, 1:h_kl))
            qx = zero; qy = zero; qz = zero
    end subroutine copydata




    ! ========================================================================================================
    ! OLD COMMIT SUBROUTINES THAT WORKED
    ! ========================================================================================================
    !all non flux subroutines 
    attributes(global) subroutine computeSpeedOfSoundSquared_v1
        use constants, only: two, third,irho,itu1
        implicit none
        integer :: i,j,k
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z
        
        if (i <= ie .AND. j <= je .AND. k <= ke) then
            !TODO this will not work if TKE is present check 
            !computeSpeedOfSoundSquared in blockette.F90 need
            !check to correct for K  
            aa(i,j,k) = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
        end if
    end subroutine computeSpeedOfSoundSquared_v1

    attributes(global) subroutine allNodalGradients_v1
        use constants, only: zero,one,fourth,ivx,ivy,ivz
        use precision, only: realType
        implicit none

        real(kind=realType) :: a2, oVol, uBar,vBar,wBar,sx,sy,sz
        real(kind=realType) :: tmp

        integer(kind=intType) :: i, j, k
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z

        ! First part. Contribution in the k-direction.
        ! The contribution is scattered to both the left and right node
        ! in k-direction.
        
        !loop k 1 to ke, j 1 to jl and i 1 to il
        if (i <= il .AND. j <= jl .AND. k <= ke) then
            ! Compute 8 times the average normal for this part of
            ! the control volume. The factor 8 is taken care of later
            ! on when the division by the volume takes place.
            sx = sk(i, j, k - 1, 1) + sk(i + 1, j, k - 1, 1) &
                + sk(i, j + 1, k - 1, 1) + sk(i + 1, j + 1, k - 1, 1) &
                + sk(i, j, k, 1) + sk(i + 1, j, k, 1) &
                + sk(i, j + 1, k, 1) + sk(i + 1, j + 1, k, 1)
            sy = sk(i, j, k - 1, 2) + sk(i + 1, j, k - 1, 2) &
                + sk(i, j + 1, k - 1, 2) + sk(i + 1, j + 1, k - 1, 2) &
                + sk(i, j, k, 2) + sk(i + 1, j, k, 2) &
                + sk(i, j + 1, k, 2) + sk(i + 1, j + 1, k, 2)
            sz = sk(i, j, k - 1, 3) + sk(i + 1, j, k - 1, 3) &
                + sk(i, j + 1, k - 1, 3) + sk(i + 1, j + 1, k - 1, 3) &
                + sk(i, j, k, 3) + sk(i + 1, j, k, 3) &
                + sk(i, j + 1, k, 3) + sk(i + 1, j + 1, k, 3)
            ! Compute the average velocities and speed of sound squared
            ! for this integration point. Node that these variables are
            ! stored in w(ivx), w(ivy), w(ivz) and p.

            ubar = fourth * (w(i, j, k, ivx) + w(i + 1, j, k, ivx) &
                            + w(i, j + 1, k, ivx) + w(i + 1, j + 1, k, ivx))
            vbar = fourth * (w(i, j, k, ivy) + w(i + 1, j, k, ivy) &
                            + w(i, j + 1, k, ivy) + w(i + 1, j + 1, k, ivy))
            wbar = fourth * (w(i, j, k, ivz) + w(i + 1, j, k, ivz) &
                            + w(i, j + 1, k, ivz) + w(i + 1, j + 1, k, ivz))

            a2 = fourth * (aa(i, j, k) + aa(i + 1, j, k) + aa(i, j + 1, k) + aa(i + 1, j + 1, k))
            
            ! Add the contributions to the surface integral to the node
            ! j-1 and substract it from the node j. For the heat flux it
            ! is reversed, because the negative of the gradient of the
            ! speed of sound must be computed.
            if (k > 1) then

                tmp = atomicadd(ux(i, j, k - 1), ubar * sx)
                tmp = atomicadd(uy(i, j, k - 1), ubar * sy)
                tmp = atomicadd(uz(i, j, k - 1), ubar * sz)
                tmp = atomicadd(vx(i, j, k - 1), vbar * sx)
                tmp = atomicadd(vy(i, j, k - 1), vbar * sy)
                tmp = atomicadd(vz(i, j, k - 1), vbar * sz)

                tmp = atomicadd(wx(i, j, k - 1), wbar * sx)
                tmp = atomicadd(wy(i, j, k - 1), wbar * sy)
                tmp = atomicadd(wz(i, j, k - 1), wbar * sz)
                
                tmp = atomicsub(qx(i, j, k - 1), a2 * sx)
                tmp = atomicsub(qy(i, j, k - 1), a2 * sy)
                tmp = atomicsub(qz(i, j, k - 1), a2 * sz)

                ! ux(i, j, k - 1) = ux(i, j, k - 1) + ubar * sx
                ! uy(i, j, k - 1) = uy(i, j, k - 1) + ubar * sy
                ! uz(i, j, k - 1) = uz(i, j, k - 1) + ubar * sz

                ! vx(i, j, k - 1) = vx(i, j, k - 1) + vbar * sx
                ! vy(i, j, k - 1) = vy(i, j, k - 1) + vbar * sy
                ! vz(i, j, k - 1) = vz(i, j, k - 1) + vbar * sz

                ! wx(i, j, k - 1) = wx(i, j, k - 1) + wbar * sx
                ! wy(i, j, k - 1) = wy(i, j, k - 1) + wbar * sy
                ! wz(i, j, k - 1) = wz(i, j, k - 1) + wbar * sz

                ! qx(i, j, k - 1) = qx(i, j, k - 1) - a2 * sx
                ! qy(i, j, k - 1) = qy(i, j, k - 1) - a2 * sy
                ! qz(i, j, k - 1) = qz(i, j, k - 1) - a2 * sz

            end if 

            if (k < ke) then

                tmp = atomicsub(ux(i, j, k), ubar * sx)
                tmp = atomicsub(uy(i, j, k), ubar * sy)
                tmp = atomicsub(uz(i, j, k), ubar * sz)

                tmp = atomicsub(vx(i, j, k), vbar * sx)
                tmp = atomicsub(vy(i, j, k), vbar * sy)
                tmp = atomicsub(vz(i, j, k), vbar * sz)

                tmp = atomicsub(wx(i, j, k), wbar * sx)
                tmp = atomicsub(wy(i, j, k), wbar * sy)
                tmp = atomicsub(wz(i, j, k), wbar * sz)
                
                tmp = atomicadd(qx(i, j, k), a2 * sx) 
                tmp = atomicadd(qy(i, j, k), a2 * sy)
                tmp = atomicadd(qz(i, j, k), a2 * sz)

                ! ux(i, j, k) = ux(i, j, k) - ubar * sx
                ! uy(i, j, k) = uy(i, j, k) - ubar * sy
                ! uz(i, j, k) = uz(i, j, k) - ubar * sz

                ! vx(i, j, k) = vx(i, j, k) - vbar * sx
                ! vy(i, j, k) = vy(i, j, k) - vbar * sy
                ! vz(i, j, k) = vz(i, j, k) - vbar * sz

                ! wx(i, j, k) = wx(i, j, k) - wbar * sx
                ! wy(i, j, k) = wy(i, j, k) - wbar * sy
                ! wz(i, j, k) = wz(i, j, k) - wbar * sz

                ! qx(i, j, k) = qx(i, j, k) + a2 * sx
                ! qy(i, j, k) = qy(i, j, k) + a2 * sy
                ! qz(i, j, k) = qz(i, j, k) + a2 * sz

            end if 
        end if
        ! Second part. Contribution in the j-direction.
        ! The contribution is scattered to both the left and right node
        ! in j-direction.
        !loop k 1 to kl j 1 to je and i 1 to il
        if (i <= il .AND. j <= je .AND. k <= kl) then
            ! Compute 8 times the average normal for this part of
            ! the control volume. The factor 8 is taken care of later
            ! on when the division by the volume takes place.
            sx = sj(i, j - 1, k, 1) + sj(i + 1, j - 1, k, 1) &
                    + sj(i, j - 1, k + 1, 1) + sj(i + 1, j - 1, k + 1, 1) &
                    + sj(i, j, k, 1) + sj(i + 1, j, k, 1) &
                    + sj(i, j, k + 1, 1) + sj(i + 1, j, k + 1, 1)
            sy = sj(i, j - 1, k, 2) + sj(i + 1, j - 1, k, 2) &
                    + sj(i, j - 1, k + 1, 2) + sj(i + 1, j - 1, k + 1, 2) &
                    + sj(i, j, k, 2) + sj(i + 1, j, k, 2) &
                    + sj(i, j, k + 1, 2) + sj(i + 1, j, k + 1, 2)
            sz = sj(i, j - 1, k, 3) + sj(i + 1, j - 1, k, 3) &
                    + sj(i, j - 1, k + 1, 3) + sj(i + 1, j - 1, k + 1, 3) &
                    + sj(i, j, k, 3) + sj(i + 1, j, k, 3) &
                    + sj(i, j, k + 1, 3) + sj(i + 1, j, k + 1, 3)
            ! Compute the average velocities and speed of sound squared
            ! for this integration point. Node that these variables are
            ! stored in w(ivx), w(ivy), w(ivz) and p.
            ubar = fourth * (w(i, j, k, ivx) + w(i + 1, j, k, ivx) &
                    + w(i, j, k + 1, ivx) + w(i + 1, j, k + 1, ivx))
            vbar = fourth * (w(i, j, k, ivy) + w(i + 1, j, k, ivy) &
                                + w(i, j, k + 1, ivy) + w(i + 1, j, k + 1, ivy))
            wbar = fourth * (w(i, j, k, ivz) + w(i + 1, j, k, ivz) &
                                + w(i, j, k + 1, ivz) + w(i + 1, j, k + 1, ivz))

            a2 = fourth * (aa(i, j, k) + aa(i + 1, j, k) + aa(i, j, k + 1) + aa(i + 1, j, k + 1))
            
            ! Add the contributions to the surface integral to the node
            ! j-1 and substract it from the node j. For the heat flux it
            ! is reversed, because the negative of the gradient of the
            ! speed of sound must be computed.
    
            if (j > 1) then

                tmp = atomicadd(ux(i, j - 1, k), ubar * sx)
                tmp = atomicadd(uy(i, j - 1, k), ubar * sy)
                tmp = atomicadd(uz(i, j - 1, k), ubar * sz)

                tmp = atomicadd(vx(i, j - 1, k), vbar * sx)
                tmp = atomicadd(vy(i, j - 1, k), vbar * sy)
                tmp = atomicadd(vz(i, j - 1, k), vbar * sz)

                tmp = atomicadd(wx(i, j - 1, k), wbar * sx)
                tmp = atomicadd(wy(i, j - 1, k), wbar * sy)
                tmp = atomicadd(wz(i, j - 1, k), wbar * sz)
                
                tmp = atomicsub(qx(i, j - 1, k), a2 * sx)
                tmp = atomicsub(qy(i, j - 1, k), a2 * sy)
                tmp = atomicsub(qz(i, j - 1, k), a2 * sz)

                ! ux(i, j - 1, k) = ux(i, j - 1, k) + ubar * sx
                ! uy(i, j - 1, k) = uy(i, j - 1, k) + ubar * sy
                ! uz(i, j - 1, k) = uz(i, j - 1, k) + ubar * sz

                ! vx(i, j - 1, k) = vx(i, j - 1, k) + vbar * sx
                ! vy(i, j - 1, k) = vy(i, j - 1, k) + vbar * sy
                ! vz(i, j - 1, k) = vz(i, j - 1, k) + vbar * sz

                ! wx(i, j - 1, k) = wx(i, j - 1, k) + wbar * sx
                ! wy(i, j - 1, k) = wy(i, j - 1, k) + wbar * sy
                ! wz(i, j - 1, k) = wz(i, j - 1, k) + wbar * sz

                ! qx(i, j - 1, k) = qx(i, j - 1, k) - a2 * sx
                ! qy(i, j - 1, k) = qy(i, j - 1, k) - a2 * sy
                ! qz(i, j - 1, k) = qz(i, j - 1, k) - a2 * sz

            end if

            if (j < je) then

                tmp = atomicsub(ux(i, j, k), ubar * sx) 
                tmp = atomicsub(uy(i, j, k), ubar * sy)
                tmp = atomicsub(uz(i, j, k), ubar * sz)

                tmp = atomicsub(vx(i, j, k), vbar * sx)
                tmp = atomicsub(vy(i, j, k), vbar * sy)
                tmp = atomicsub(vz(i, j, k), vbar * sz)

                tmp = atomicsub(wx(i, j, k), wbar * sx)
                tmp = atomicsub(wy(i, j, k), wbar * sy)
                tmp = atomicsub(wz(i, j, k), wbar * sz)
                
                tmp = atomicadd(qx(i, j, k), a2 * sx)
                tmp = atomicadd(qy(i, j, k), a2 * sy)
                tmp = atomicadd(qz(i, j, k), a2 * sz)

                ! ux(i, j, k) = ux(i, j, k) - ubar * sx
                ! uy(i, j, k) = uy(i, j, k) - ubar * sy
                ! uz(i, j, k) = uz(i, j, k) - ubar * sz

                ! vx(i, j, k) = vx(i, j, k) - vbar * sx
                ! vy(i, j, k) = vy(i, j, k) - vbar * sy
                ! vz(i, j, k) = vz(i, j, k) - vbar * sz

                ! wx(i, j, k) = wx(i, j, k) - wbar * sx
                ! wy(i, j, k) = wy(i, j, k) - wbar * sy
                ! wz(i, j, k) = wz(i, j, k) - wbar * sz

                ! qx(i, j, k) = qx(i, j, k) + a2 * sx
                ! qy(i, j, k) = qy(i, j, k) + a2 * sy
                ! qz(i, j, k) = qz(i, j, k) + a2 * sz

            end if 
        end if
        
        ! Third part. Contribution in the i-direction.
        ! The contribution is scattered to both the left and right node
        ! in i-direction.
        
        ! loop k 1 to kl j 1 to jl and i 1 to ie
        if (i <= ie .AND. j <= jl .AND. k <= kl) then
            ! Compute 8 times the average normal for this part of
            ! the control volume. The factor 8 is taken care of later
            ! on when the division by the volume takes place.
            sx = si(i - 1, j, k, 1) + si(i - 1, j + 1, k, 1) &
                    + si(i - 1, j, k + 1, 1) + si(i - 1, j + 1, k + 1, 1) &
                    + si(i, j, k, 1) + si(i, j + 1, k, 1) &
                    + si(i, j, k + 1, 1) + si(i, j + 1, k + 1, 1)
            sy = si(i - 1, j, k, 2) + si(i - 1, j + 1, k, 2) &
                    + si(i - 1, j, k + 1, 2) + si(i - 1, j + 1, k + 1, 2) &
                    + si(i, j, k, 2) + si(i, j + 1, k, 2) &
                    + si(i, j, k + 1, 2) + si(i, j + 1, k + 1, 2)
            sz = si(i - 1, j, k, 3) + si(i - 1, j + 1, k, 3) &
                    + si(i - 1, j, k + 1, 3) + si(i - 1, j + 1, k + 1, 3) &
                    + si(i, j, k, 3) + si(i, j + 1, k, 3) &
                    + si(i, j, k + 1, 3) + si(i, j + 1, k + 1, 3)
            ! Compute the average velocities and speed of sound squared
            ! for this integration point. Node that these variables are
            ! stored in w(ivx), w(ivy), w(ivz) and p.
            ubar = fourth * (w(i, j, k, ivx) + w(i, j + 1, k, ivx) &
                                + w(i, j, k + 1, ivx) + w(i, j + 1, k + 1, ivx))
            vbar = fourth * (w(i, j, k, ivy) + w(i, j + 1, k, ivy) &
                                + w(i, j, k + 1, ivy) + w(i, j + 1, k + 1, ivy))
            wbar = fourth * (w(i, j, k, ivz) + w(i, j + 1, k, ivz) &
                                + w(i, j, k + 1, ivz) + w(i, j + 1, k + 1, ivz))

            a2 = fourth * (aa(i, j, k) + aa(i, j + 1, k) + aa(i, j, k + 1) + aa(i, j + 1, k + 1))
            ! Add the contributions to the surface integral to the node
            ! j-1 and substract it from the node j. For the heat flux it
            ! is reversed, because the negative of the gradient of the
            ! speed of sound must be computed.
            if (i > 1) then

                tmp = atomicadd(ux(i - 1, j, k), ubar * sx) 
                tmp = atomicadd(uy(i - 1, j, k), ubar * sy)
                tmp = atomicadd(uz(i - 1, j, k), ubar * sz)

                tmp = atomicadd(vx(i - 1, j, k), vbar * sx)
                tmp = atomicadd(vy(i - 1, j, k), vbar * sy)
                tmp = atomicadd(vz(i - 1, j, k), vbar * sz)

                tmp = atomicadd(wx(i - 1, j, k), wbar * sx)
                tmp = atomicadd(wy(i - 1, j, k), wbar * sy)
                tmp = atomicadd(wz(i - 1, j, k), wbar * sz)
                
                tmp = atomicsub(qx(i - 1, j, k), a2 * sx)
                tmp = atomicsub(qy(i - 1, j, k), a2 * sy)
                tmp = atomicsub(qz(i - 1, j, k), a2 * sz)

                ! ux(i - 1, j, k) = ux(i - 1, j, k) + ubar * sx
                ! uy(i - 1, j, k) = uy(i - 1, j, k) + ubar * sy
                ! uz(i - 1, j, k) = uz(i - 1, j, k) + ubar * sz

                ! vx(i - 1, j, k) = vx(i - 1, j, k) + vbar * sx
                ! vy(i - 1, j, k) = vy(i - 1, j, k) + vbar * sy
                ! vz(i - 1, j, k) = vz(i - 1, j, k) + vbar * sz

                ! wx(i - 1, j, k) = wx(i - 1, j, k) + wbar * sx
                ! wy(i - 1, j, k) = wy(i - 1, j, k) + wbar * sy
                ! wz(i - 1, j, k) = wz(i - 1, j, k) + wbar * sz

                ! qx(i - 1, j, k) = qx(i - 1, j, k) - a2 * sx
                ! qy(i - 1, j, k) = qy(i - 1, j, k) - a2 * sy
                ! qz(i - 1, j, k) = qz(i - 1, j, k) - a2 * sz

            end if

            if (i < ie) then

                tmp = atomicsub(ux(i, j, k), ubar * sx)
                tmp = atomicsub(uy(i, j, k), ubar * sy)
                tmp = atomicsub(uz(i, j, k), ubar * sz)

                tmp = atomicsub(vx(i, j, k), vbar * sx)
                tmp = atomicsub(vy(i, j, k), vbar * sy)
                tmp = atomicsub(vz(i, j, k), vbar * sz)

                tmp = atomicsub(wx(i, j, k), wbar * sx)
                tmp = atomicsub(wy(i, j, k), wbar * sy)
                tmp = atomicsub(wz(i, j, k), wbar * sz)
                
                tmp = atomicadd(qx(i, j, k), a2 * sx)
                tmp = atomicadd(qy(i, j, k), a2 * sy)
                tmp = atomicadd(qz(i, j, k), a2 * sz)

                ! ux(i, j, k) = ux(i, j, k) - ubar * sx
                ! uy(i, j, k) = uy(i, j, k) - ubar * sy
                ! uz(i, j, k) = uz(i, j, k) - ubar * sz

                ! vx(i, j, k) = vx(i, j, k) - vbar * sx
                ! vy(i, j, k) = vy(i, j, k) - vbar * sy
                ! vz(i, j, k) = vz(i, j, k) - vbar * sz

                ! wx(i, j, k) = wx(i, j, k) - wbar * sx
                ! wy(i, j, k) = wy(i, j, k) - wbar * sy
                ! wz(i, j, k) = wz(i, j, k) - wbar * sz

                ! qx(i, j, k) = qx(i, j, k) + a2 * sx
                ! qy(i, j, k) = qy(i, j, k) + a2 * sy
                ! qz(i, j, k) = qz(i, j, k) + a2 * sz

            end if
        end if 

        !need to get rid of this to avoid race condition
        !each cell only scales one node 
        !however need to wait for each nodes accumulation to be done
        !if thread on edge of boundary need to wait for another block to be done
        !therefore split this part into a separate kernel

        ! ! Divide by 8 times the volume to obtain the correct gradients.
        ! if (i <= il .AND. j <= jl .AND. k <= kl) then
        !     ! Compute the inverse of 8 times the volume for this node.
        !     oVol = one / (vol(i, j, k) + vol(i, j, k + 1) &
        !                           + vol(i + 1, j, k) + vol(i + 1, j, k + 1) &
        !                           + vol(i, j + 1, k) + vol(i, j + 1, k + 1) &
        !                           + vol(i + 1, j + 1, k) + vol(i + 1, j + 1, k + 1))
        !     ! Compute the correct velocity gradients and "unit" heat
        !     ! fluxes. The velocity gradients are stored in ux, etc.
        !     ux(i, j, k) = ux(i, j, k) * oVol
        !     uy(i, j, k) = uy(i, j, k) * oVol
        !     uz(i, j, k) = uz(i, j, k) * oVol

        !     vx(i, j, k) = vx(i, j, k) * oVol
        !     vy(i, j, k) = vy(i, j, k) * oVol
        !     vz(i, j, k) = vz(i, j, k) * oVol

        !     wx(i, j, k) = wx(i, j, k) * oVol
        !     wy(i, j, k) = wy(i, j, k) * oVol
        !     wz(i, j, k) = wz(i, j, k) * oVol

        !     qx(i, j, k) = qx(i, j, k) * oVol
        !     qy(i, j, k) = qy(i, j, k) * oVol
        !     qz(i, j, k) = qz(i, j, k) * oVol
        ! end if 
    end subroutine allNodalGradients_v1

    attributes(global) subroutine scaleNodalGradients_v1
        use constants, only: zero,fourth,ivx,ivy,ivz,one
        use precision, only: realType
        implicit none

        real(kind=realType) :: a2, oVol, uBar,vBar,wBar,sx,sy,sz

        integer(kind=intType) :: i, j, k
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z
        if (i <= il .AND. j <= jl .AND. k <= kl) then
            ! Compute the inverse of 8 times the volume for this node.
            oVol = one / (vol(i, j, k) + vol(i, j, k + 1) &
                                  + vol(i + 1, j, k) + vol(i + 1, j, k + 1) &
                                  + vol(i, j + 1, k) + vol(i, j + 1, k + 1) &
                                  + vol(i + 1, j + 1, k) + vol(i + 1, j + 1, k + 1))
            ! Compute the correct velocity gradients and "unit" heat
            ! fluxes. The velocity gradients are stored in ux, etc.
            ux(i, j, k) = ux(i, j, k) * oVol
            uy(i, j, k) = uy(i, j, k) * oVol
            uz(i, j, k) = uz(i, j, k) * oVol

            vx(i, j, k) = vx(i, j, k) * oVol
            vy(i, j, k) = vy(i, j, k) * oVol
            vz(i, j, k) = vz(i, j, k) * oVol

            wx(i, j, k) = wx(i, j, k) * oVol
            wy(i, j, k) = wy(i, j, k) * oVol
            wz(i, j, k) = wz(i, j, k) * oVol

            qx(i, j, k) = qx(i, j, k) * oVol
            qy(i, j, k) = qy(i, j, k) * oVol
            qz(i, j, k) = qz(i, j, k) * oVol
        end if 

    end subroutine scaleNodalGradients_v1
    
    attributes(global) subroutine viscousFlux_v1
        use precision, only: realType, intType
        use constants, only: half, zero,one,two, third,fourth,eighth,ivx,ivy,ivz,irhoE,irho,itu1,imx,imy,imz,noFlux
        use cudaInputPhysics, only: useQCR, prandtl, prandtlturb
        use cudaFlowvarRefState, only: eddyModel
        use cudaIteration, only: rFil
        implicit none

        ! Variables for viscous flux
        real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
        real(kind=realType) :: gm1, factLamHeat, factTurbHeat
        real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
        real(kind=realType) :: q_x, q_y, q_z
        real(kind=realType) :: corr, ssx, ssy, ssz, fracDiv, snrm
        real(kind=realType) :: tauxx, tauyy, tauzz
        real(kind=realType) :: tauxy, tauxz, tauyz
        real(kind=realType) :: tauxxS, tauyyS, tauzzS
        real(kind=realType) :: tauxyS, tauxzS, tauyzS
        real(kind=realType) :: ubar, vbar, wbar
        real(kind=realType) :: exx, eyy, ezz
        real(kind=realType) :: exy, exz, eyz
        real(kind=realType) :: Wxx, Wyy, Wzz
        real(kind=realType) :: Wxy, Wxz, Wyz, Wyx, Wzx, Wzy
        real(kind=realType) :: den, Ccr1
        real(kind=realType) :: fmx, fmy, fmz, frhoE, fact
        integer(kind=intType) :: i, j, k, io, jo, ko
        real(kind=realType), parameter :: xminn = 1.e-10_realType
        real(kind=realType), parameter :: twoThird = two * third
        real(kind=realType) :: tmp
        !thread indices start at 1
        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x 
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y 
        k = (blockIdx%z - 1) * blockDim%z + threadIdx%z 

        ! Set QCR parameters
        Ccr1 = 0.3_realType
        rFilv = rFil

        ! The diagonals of the vorticity tensor components are always zero
        Wxx = zero
        Wyy = zero
        Wzz = zero

        !
        !         viscous fluxes in the k-direction.
        !
        mue = zero
        !loop k 1 to kl, j 2 to jl, i 2 to il
        if (k <= kl .and. j <= jl .and. i <= il .and. j>=2 .and. i>=2 .and. k>=1) then
            ! Set the value of the porosity. If not zero, it is set
            ! to average the eddy-viscosity and to take the factor
            ! rFilv into account.
            por = half * rFilv
            if (porK(i, j, k) == noFlux) por = zero

            ! Compute the laminar and (if present) the eddy viscosities
            ! multiplied by the porosity. Compute the factor in front of
            ! the gradients of the speed of sound squared for the heat
            ! flux.
            mul = por * (rlv(i, j, k) + rlv(i, j, k + 1))
            mue = por * (rev(i, j, k) + rev(i, j, k + 1))
            mut = mul + mue

            gm1 = half * (gamma(i, j, k) + gamma(i, j, k + 1)) - one
            factLamHeat = one / (prandtl * gm1)
            factTurbHeat = one / (prandtlTurb * gm1)

            heatCoef = mul * factLamHeat + mue * factTurbHeat
            ! Compute the gradients at the face by averaging the four
                    ! nodal values.

            u_x = fourth * (ux(i - 1, j - 1, k) + ux(i, j - 1, k) &
                        + ux(i - 1, j, k) + ux(i, j, k))
            u_y = fourth * (uy(i - 1, j - 1, k) + uy(i, j - 1, k) &
                        + uy(i - 1, j, k) + uy(i, j, k))
            u_z = fourth * (uz(i - 1, j - 1, k) + uz(i, j - 1, k) &
                        + uz(i - 1, j, k) + uz(i, j, k))

            v_x = fourth * (vx(i - 1, j - 1, k) + vx(i, j - 1, k) &
                        + vx(i - 1, j, k) + vx(i, j, k))
            v_y = fourth * (vy(i - 1, j - 1, k) + vy(i, j - 1, k) &
                        + vy(i - 1, j, k) + vy(i, j, k))
            v_z = fourth * (vz(i - 1, j - 1, k) + vz(i, j - 1, k) &
                        + vz(i - 1, j, k) + vz(i, j, k))

            w_x = fourth * (wx(i - 1, j - 1, k) + wx(i, j - 1, k) &
                        + wx(i - 1, j, k) + wx(i, j, k))
            w_y = fourth * (wy(i - 1, j - 1, k) + wy(i, j - 1, k) &
                        + wy(i - 1, j, k) + wy(i, j, k))
            w_z = fourth * (wz(i - 1, j - 1, k) + wz(i, j - 1, k) &
                        + wz(i - 1, j, k) + wz(i, j, k))

            q_x = fourth * (qx(i - 1, j - 1, k) + qx(i, j - 1, k) &
                        + qx(i - 1, j, k) + qx(i, j, k))
            q_y = fourth * (qy(i - 1, j - 1, k) + qy(i, j - 1, k) &
                        + qy(i - 1, j, k) + qy(i, j, k))
            q_z = fourth * (qz(i - 1, j - 1, k) + qz(i, j - 1, k) &
                        + qz(i - 1, j, k) + qz(i, j, k))
            ! The gradients in the normal direction are corrected, such
            ! that no averaging takes places here.
            ! First determine the vector in the direction from the
            ! cell center k to cell center k+1.

            ssx = eighth * (x(i - 1, j - 1, k + 1, 1) - x(i - 1, j - 1, k - 1, 1) &
                        + x(i - 1, j, k + 1, 1) - x(i - 1, j, k - 1, 1) &
                        + x(i, j - 1, k + 1, 1) - x(i, j - 1, k - 1, 1) &
                        + x(i, j, k + 1, 1) - x(i, j, k - 1, 1))
            ssy = eighth * (x(i - 1, j - 1, k + 1, 2) - x(i - 1, j - 1, k - 1, 2) &
                        + x(i - 1, j, k + 1, 2) - x(i - 1, j, k - 1, 2) &
                        + x(i, j - 1, k + 1, 2) - x(i, j - 1, k - 1, 2) &
                        + x(i, j, k + 1, 2) - x(i, j, k - 1, 2))
            ssz = eighth * (x(i - 1, j - 1, k + 1, 3) - x(i - 1, j - 1, k - 1, 3) &
                        + x(i - 1, j, k + 1, 3) - x(i - 1, j, k - 1, 3) &
                        + x(i, j - 1, k + 1, 3) - x(i, j - 1, k - 1, 3) &
                        + x(i, j, k + 1, 3) - x(i, j, k - 1, 3))

            ! Determine the length of this vector and create the
            ! unit normal.

            snrm = one / sqrt(ssx * ssx + ssy * ssy + ssz * ssz)
            ssx = snrm * ssx
            ssy = snrm * ssy
            ssz = snrm * ssz
    
            ! Correct the gradients.

            corr = u_x * ssx + u_y * ssy + u_z * ssz &
                    - (w(i, j, k + 1, ivx) - w(i, j, k, ivx)) * snrm
            u_x = u_x - corr * ssx
            u_y = u_y - corr * ssy
            u_z = u_z - corr * ssz

            corr = v_x * ssx + v_y * ssy + v_z * ssz &
                    - (w(i, j, k + 1, ivy) - w(i, j, k, ivy)) * snrm
            v_x = v_x - corr * ssx
            v_y = v_y - corr * ssy
            v_z = v_z - corr * ssz

            corr = w_x * ssx + w_y * ssy + w_z * ssz &
                    - (w(i, j, k + 1, ivz) - w(i, j, k, ivz)) * snrm
            w_x = w_x - corr * ssx
            w_y = w_y - corr * ssy
            w_z = w_z - corr * ssz

            corr = q_x * ssx + q_y * ssy + q_z * ssz &
                    + (aa(i, j, k + 1) - aa(i, j, k)) * snrm
            q_x = q_x - corr * ssx
            q_y = q_y - corr * ssy
            q_z = q_z - corr * ssz
            
            ! Compute the stress tensor and the heat flux vector.
            ! We remove the viscosity from the stress tensor (tau)
            ! to define tauS since we still need to separate between
            ! laminar and turbulent stress for QCR.
            ! Therefore, laminar tau = mue*tauS, turbulent
            ! tau = mue*tauS, and total tau = mut*tauS.

            fracDiv = twoThird * (u_x + v_y + w_z)
            tauxxS = two * u_x - fracDiv
            tauyyS = two * v_y - fracDiv
            tauzzS = two * w_z - fracDiv

            tauxyS = u_y + v_x
            tauxzS = u_z + w_x
            tauyzS = v_z + w_y

            q_x = heatCoef * q_x
            q_y = heatCoef * q_y
            q_z = heatCoef * q_z
            
            ! Add QCR corrections if necessary
            if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x * u_x + u_y * u_y + u_z * u_z + &
                           v_x * v_x + v_y * v_y + v_z * v_z + &
                           w_x * w_x + w_y * w_y + w_z * w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue * Ccr1 / den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact * (Wxy * tauxyS + Wxz * tauxzS) * two
                eyy = fact * (Wyx * tauxyS + Wyz * tauyzS) * two
                ezz = fact * (Wzx * tauxzS + Wzy * tauyzS) * two

                exy = fact * (Wxy * tauyyS + Wxz * tauyzS + &
                              Wyx * tauxxS + Wyz * tauxzS)
                exz = fact * (Wxy * tauyzS + Wxz * tauzzS + &
                              Wzx * tauxxS + Wzy * tauxyS)
                eyz = fact * (Wyx * tauxzS + Wyz * tauzzS + &
                              Wzx * tauxyS + Wzy * tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut * tauxxS - exx
                tauyy = mut * tauyyS - eyy
                tauzz = mut * tauzzS - ezz
                tauxy = mut * tauxyS - exy
                tauxz = mut * tauxzS - exz
                tauyz = mut * tauyzS - eyz

            else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut * tauxxS
                tauyy = mut * tauyyS
                tauzz = mut * tauzzS
                tauxy = mut * tauxyS
                tauxz = mut * tauxzS
                tauyz = mut * tauyzS

            end if

            ! Compute the average velocities for the face. Remember that
            ! the velocities are stored and not the momentum.
            ubar = half * (w(i, j, k, ivx) + w(i, j, k + 1, ivx))
            vbar = half * (w(i, j, k, ivy) + w(i, j, k + 1, ivy))
            wbar = half * (w(i, j, k, ivz) + w(i, j, k + 1, ivz))

            ! Compute the viscous fluxes for this k-face.
            fmx = tauxx * sk(i, j, k, 1) + tauxy * sk(i, j, k, 2) &
                    + tauxz * sk(i, j, k, 3)
            fmy = tauxy * sk(i, j, k, 1) + tauyy * sk(i, j, k, 2) &
                    + tauyz * sk(i, j, k, 3)
            fmz = tauxz * sk(i, j, k, 1) + tauyz * sk(i, j, k, 2) &
                    + tauzz * sk(i, j, k, 3)
            frhoE = (ubar * tauxx + vbar * tauxy + wbar * tauxz) * sk(i, j, k, 1)
            frhoE = frhoE + (ubar * tauxy + vbar * tauyy + wbar * tauyz) * sk(i, j, k, 2)
            frhoE = frhoE + (ubar * tauxz + vbar * tauyz + wbar * tauzz) * sk(i, j, k, 3)
            frhoE = frhoE - q_x * sk(i, j, k, 1) - q_y * sk(i, j, k, 2) - q_z * sk(i, j, k, 3)

            ! Update the residuals of cell k and k+1.
            
            tmp = atomicsub(fw(i, j, k, imx), fmx)
            tmp = atomicsub(fw(i, j, k, imy), fmy)
            tmp = atomicsub(fw(i, j, k, imz), fmz)
            tmp = atomicsub(fw(i, j, k, irhoE), frhoE)
            ! fw(i, j, k, imx) = fw(i, j, k, imx) - fmx
            ! fw(i, j, k, imy) = fw(i, j, k, imy) - fmy
            ! fw(i, j, k, imz) = fw(i, j, k, imz) - fmz
            ! fw(i, j, k, irhoE) = fw(i, j, k, irhoE) - frhoE

            tmp = atomicadd(fw(i, j, k + 1, imx), fmx)
            tmp = atomicadd(fw(i, j, k + 1, imy), fmy)
            tmp = atomicadd(fw(i, j, k + 1, imz), fmz)
            tmp = atomicadd(fw(i, j, k + 1, irhoE), frhoE)
            ! fw(i, j, k + 1, imx) = fw(i, j, k + 1, imx) + fmx
            ! fw(i, j, k + 1, imy) = fw(i, j, k + 1, imy) + fmy
            ! fw(i, j, k + 1, imz) = fw(i, j, k + 1, imz) + fmz
            ! fw(i, j, k + 1, irhoE) = fw(i, j, k + 1, irhoE) + frhoE

            ! Temporarily store the shear stress and heat flux, even
            ! if we won't need it. This can still vectorize
        end if

        !
        !         Viscous fluxes in the j-direction.
        !
        !loop k 2 to kl, j 1 to jl, i 2 to il
        if (k <= kl .and. j <= jl .and. i <= il .and. k>=2 .and. i>=2 .and. j>=1) then
            ! Set the value of the porosity. If not zero, it is set
                    ! to average the eddy-viscosity and to take the factor
                    ! rFilv into account.

            por = half * rFilv
            if (porJ(i, j, k) == noFlux) por = zero

            ! Compute the laminar and (if present) the eddy viscosities
            ! multiplied by the porosity. Compute the factor in front of
            ! the gradients of the speed of sound squared for the heat
            ! flux.

            mul = por * (rlv(i, j, k) + rlv(i, j + 1, k))
            mue = por * (rev(i, j, k) + rev(i, j + 1, k))
            mut = mul + mue

            gm1 = half * (gamma(i, j, k) + gamma(i, j + 1, k)) - one
            factLamHeat = one / (prandtl * gm1)
            factTurbHeat = one / (prandtlTurb * gm1)

            heatCoef = mul * factLamHeat + mue * factTurbHeat

            ! Compute the gradients at the face by averaging the four
            ! nodal values.

            u_x = fourth * (ux(i - 1, j, k - 1) + ux(i, j, k - 1) &
                            + ux(i - 1, j, k) + ux(i, j, k))
            u_y = fourth * (uy(i - 1, j, k - 1) + uy(i, j, k - 1) &
                            + uy(i - 1, j, k) + uy(i, j, k))
            u_z = fourth * (uz(i - 1, j, k - 1) + uz(i, j, k - 1) &
                            + uz(i - 1, j, k) + uz(i, j, k))

            v_x = fourth * (vx(i - 1, j, k - 1) + vx(i, j, k - 1) &
                            + vx(i - 1, j, k) + vx(i, j, k))
            v_y = fourth * (vy(i - 1, j, k - 1) + vy(i, j, k - 1) &
                            + vy(i - 1, j, k) + vy(i, j, k))
            v_z = fourth * (vz(i - 1, j, k - 1) + vz(i, j, k - 1) &
                            + vz(i - 1, j, k) + vz(i, j, k))

            w_x = fourth * (wx(i - 1, j, k - 1) + wx(i, j, k - 1) &
                            + wx(i - 1, j, k) + wx(i, j, k))
            w_y = fourth * (wy(i - 1, j, k - 1) + wy(i, j, k - 1) &
                            + wy(i - 1, j, k) + wy(i, j, k))
            w_z = fourth * (wz(i - 1, j, k - 1) + wz(i, j, k - 1) &
                            + wz(i - 1, j, k) + wz(i, j, k))

            q_x = fourth * (qx(i - 1, j, k - 1) + qx(i, j, k - 1) &
                            + qx(i - 1, j, k) + qx(i, j, k))
            q_y = fourth * (qy(i - 1, j, k - 1) + qy(i, j, k - 1) &
                            + qy(i - 1, j, k) + qy(i, j, k))
            q_z = fourth * (qz(i - 1, j, k - 1) + qz(i, j, k - 1) &
                            + qz(i - 1, j, k) + qz(i, j, k))

            ! The gradients in the normal direction are corrected, such
            ! that no averaging takes places here.
            ! First determine the vector in the direction from the
            ! cell center j to cell center j+1.

            ssx = eighth * (x(i - 1, j + 1, k - 1, 1) - x(i - 1, j - 1, k - 1, 1) &
                            + x(i - 1, j + 1, k, 1) - x(i - 1, j - 1, k, 1) &
                            + x(i, j + 1, k - 1, 1) - x(i, j - 1, k - 1, 1) &
                            + x(i, j + 1, k, 1) - x(i, j - 1, k, 1))
            ssy = eighth * (x(i - 1, j + 1, k - 1, 2) - x(i - 1, j - 1, k - 1, 2) &
                            + x(i - 1, j + 1, k, 2) - x(i - 1, j - 1, k, 2) &
                            + x(i, j + 1, k - 1, 2) - x(i, j - 1, k - 1, 2) &
                            + x(i, j + 1, k, 2) - x(i, j - 1, k, 2))
            ssz = eighth * (x(i - 1, j + 1, k - 1, 3) - x(i - 1, j - 1, k - 1, 3) &
                            + x(i - 1, j + 1, k, 3) - x(i - 1, j - 1, k, 3) &
                            + x(i, j + 1, k - 1, 3) - x(i, j - 1, k - 1, 3) &
                            + x(i, j + 1, k, 3) - x(i, j - 1, k, 3))

            ! Determine the length of this vector and create the
            ! unit normal.

            snrm = one / sqrt(ssx * ssx + ssy * ssy + ssz * ssz)
            ssx = snrm * ssx
            ssy = snrm * ssy
            ssz = snrm * ssz

            ! Correct the gradients.

            corr = u_x * ssx + u_y * ssy + u_z * ssz &
                    - (w(i, j + 1, k, ivx) - w(i, j, k, ivx)) * snrm
            u_x = u_x - corr * ssx
            u_y = u_y - corr * ssy
            u_z = u_z - corr * ssz

            corr = v_x * ssx + v_y * ssy + v_z * ssz &
                    - (w(i, j + 1, k, ivy) - w(i, j, k, ivy)) * snrm
            v_x = v_x - corr * ssx
            v_y = v_y - corr * ssy
            v_z = v_z - corr * ssz

            corr = w_x * ssx + w_y * ssy + w_z * ssz &
                    - (w(i, j + 1, k, ivz) - w(i, j, k, ivz)) * snrm
            w_x = w_x - corr * ssx
            w_y = w_y - corr * ssy
            w_z = w_z - corr * ssz

            corr = q_x * ssx + q_y * ssy + q_z * ssz &
                    + (aa(i, j + 1, k) - aa(i, j, k)) * snrm
            q_x = q_x - corr * ssx
            q_y = q_y - corr * ssy
            q_z = q_z - corr * ssz

            ! Compute the stress tensor and the heat flux vector.
            ! We remove the viscosity from the stress tensor (tau)
            ! to define tauS since we still need to separate between
            ! laminar and turbulent stress for QCR.
            ! Therefore, laminar tau = mue*tauS, turbulent
            ! tau = mue*tauS, and total tau = mut*tauS.

            fracDiv = twoThird * (u_x + v_y + w_z)

            tauxxS = two * u_x - fracDiv
            tauyyS = two * v_y - fracDiv
            tauzzS = two * w_z - fracDiv

            tauxyS = u_y + v_x
            tauxzS = u_z + w_x
            tauyzS = v_z + w_y

            q_x = heatCoef * q_x
            q_y = heatCoef * q_y
            q_z = heatCoef * q_z

            ! Add QCR corrections if necessary
            if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x * u_x + u_y * u_y + u_z * u_z + &
                            v_x * v_x + v_y * v_y + v_z * v_z + &
                            w_x * w_x + w_y * w_y + w_z * w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue * Ccr1 / den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact * (Wxy * tauxyS + Wxz * tauxzS) * two
                eyy = fact * (Wyx * tauxyS + Wyz * tauyzS) * two
                ezz = fact * (Wzx * tauxzS + Wzy * tauyzS) * two

                exy = fact * (Wxy * tauyyS + Wxz * tauyzS + &
                                Wyx * tauxxS + Wyz * tauxzS)
                exz = fact * (Wxy * tauyzS + Wxz * tauzzS + &
                                Wzx * tauxxS + Wzy * tauxyS)
                eyz = fact * (Wyx * tauxzS + Wyz * tauzzS + &
                                Wzx * tauxyS + Wzy * tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut * tauxxS - exx
                tauyy = mut * tauyyS - eyy
                tauzz = mut * tauzzS - ezz
                tauxy = mut * tauxyS - exy
                tauxz = mut * tauxzS - exz
                tauyz = mut * tauyzS - eyz

            else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut * tauxxS
                tauyy = mut * tauyyS
                tauzz = mut * tauzzS
                tauxy = mut * tauxyS
                tauxz = mut * tauxzS
                tauyz = mut * tauyzS

            end if

            ! Compute the average velocities for the face. Remember that
            ! the velocities are stored and not the momentum.

            ubar = half * (w(i, j, k, ivx) + w(i, j + 1, k, ivx))
            vbar = half * (w(i, j, k, ivy) + w(i, j + 1, k, ivy))
            wbar = half * (w(i, j, k, ivz) + w(i, j + 1, k, ivz))

            ! Compute the viscous fluxes for this j-face.

            fmx = tauxx * sj(i, j, k, 1) + tauxy * sj(i, j, k, 2) &
                    + tauxz * sj(i, j, k, 3)
            fmy = tauxy * sj(i, j, k, 1) + tauyy * sj(i, j, k, 2) &
                    + tauyz * sj(i, j, k, 3)
            fmz = tauxz * sj(i, j, k, 1) + tauyz * sj(i, j, k, 2) &
                    + tauzz * sj(i, j, k, 3)
            frhoE = (ubar * tauxx + vbar * tauxy + wbar * tauxz) * sj(i, j, k, 1) &
                    + (ubar * tauxy + vbar * tauyy + wbar * tauyz) * sj(i, j, k, 2) &
                    + (ubar * tauxz + vbar * tauyz + wbar * tauzz) * sj(i, j, k, 3) &
                    - q_x * sj(i, j, k, 1) - q_y * sj(i, j, k, 2) - q_z * sj(i, j, k, 3)

            ! Update the residuals of cell j and j+1.

            tmp = atomicsub(fw(i, j, k, imx), fmx)
            tmp = atomicsub(fw(i, j, k, imy), fmy)
            tmp = atomicsub(fw(i, j, k, imz), fmz)
            tmp = atomicsub(fw(i, j, k, irhoE), frhoE)
            ! fw(i, j, k, imx) = fw(i, j, k, imx) - fmx
            ! fw(i, j, k, imy) = fw(i, j, k, imy) - fmy
            ! fw(i, j, k, imz) = fw(i, j, k, imz) - fmz
            ! fw(i, j, k, irhoE) = fw(i, j, k, irhoE) - frhoE

            tmp = atomicadd(fw(i, j + 1, k, imx), fmx)
            tmp = atomicadd(fw(i, j + 1, k, imy), fmy)
            tmp = atomicadd(fw(i, j + 1, k, imz), fmz)
            tmp = atomicadd(fw(i, j + 1, k, irhoE), frhoE)
            ! fw(i, j + 1, k, imx) = fw(i, j + 1, k, imx) + fmx
            ! fw(i, j + 1, k, imy) = fw(i, j + 1, k, imy) + fmy
            ! fw(i, j + 1, k, imz) = fw(i, j + 1, k, imz) + fmz
            ! fw(i, j + 1, k, irhoE) = fw(i, j + 1, k, irhoE) + frhoE
        end if  
        !
        !         Viscous fluxes in the i-direction.
        !
        !loop k 2 to kl, j 2 to jl, i 1 to il
        if (k <= kl .and. j <= jl .and. i <= il .and. k>=2 .and. j>=2 .and. i>=1) then
            ! Set the value of the porosity. If not zero, it is set
            ! to average the eddy-viscosity and to take the factor
            ! rFilv into account.

            por = half * rFilv
            if (porI(i, j, k) == noFlux) por = zero

            ! Compute the laminar and (if present) the eddy viscosities
            ! multiplied the porosity. Compute the factor in front of
            ! the gradients of the speed of sound squared for the heat
            ! flux.

            mul = por * (rlv(i, j, k) + rlv(i + 1, j, k))
            mue = por * (rev(i, j, k) + rev(i + 1, j, k))
            mut = mul + mue

            gm1 = half * (gamma(i, j, k) + gamma(i + 1, j, k)) - one
            factLamHeat = one / (prandtl * gm1)
            factTurbHeat = one / (prandtlTurb * gm1)

            heatCoef = mul * factLamHeat + mue * factTurbHeat

            ! Compute the gradients at the face by averaging the four
            ! nodal values.

            u_x = fourth * (ux(i, j - 1, k - 1) + ux(i, j, k - 1) &
                            + ux(i, j - 1, k) + ux(i, j, k))
            u_y = fourth * (uy(i, j - 1, k - 1) + uy(i, j, k - 1) &
                            + uy(i, j - 1, k) + uy(i, j, k))
            u_z = fourth * (uz(i, j - 1, k - 1) + uz(i, j, k - 1) &
                            + uz(i, j - 1, k) + uz(i, j, k))

            v_x = fourth * (vx(i, j - 1, k - 1) + vx(i, j, k - 1) &
                            + vx(i, j - 1, k) + vx(i, j, k))
            v_y = fourth * (vy(i, j - 1, k - 1) + vy(i, j, k - 1) &
                            + vy(i, j - 1, k) + vy(i, j, k))
            v_z = fourth * (vz(i, j - 1, k - 1) + vz(i, j, k - 1) &
                            + vz(i, j - 1, k) + vz(i, j, k))

            w_x = fourth * (wx(i, j - 1, k - 1) + wx(i, j, k - 1) &
                            + wx(i, j - 1, k) + wx(i, j, k))
            w_y = fourth * (wy(i, j - 1, k - 1) + wy(i, j, k - 1) &
                            + wy(i, j - 1, k) + wy(i, j, k))
            w_z = fourth * (wz(i, j - 1, k - 1) + wz(i, j, k - 1) &
                            + wz(i, j - 1, k) + wz(i, j, k))

            q_x = fourth * (qx(i, j - 1, k - 1) + qx(i, j, k - 1) &
                            + qx(i, j - 1, k) + qx(i, j, k))
            q_y = fourth * (qy(i, j - 1, k - 1) + qy(i, j, k - 1) &
                            + qy(i, j - 1, k) + qy(i, j, k))
            q_z = fourth * (qz(i, j - 1, k - 1) + qz(i, j, k - 1) &
                            + qz(i, j - 1, k) + qz(i, j, k))

            ! The gradients in the normal direction are corrected, such
            ! that no averaging takes places here.
            ! First determine the vector in the direction from the
            ! cell center i to cell center i+1.

            ssx = eighth * (x(i + 1, j - 1, k - 1, 1) - x(i - 1, j - 1, k - 1, 1) &
                            + x(i + 1, j - 1, k, 1) - x(i - 1, j - 1, k, 1) &
                            + x(i + 1, j, k - 1, 1) - x(i - 1, j, k - 1, 1) &
                            + x(i + 1, j, k, 1) - x(i - 1, j, k, 1))
            ssy = eighth * (x(i + 1, j - 1, k - 1, 2) - x(i - 1, j - 1, k - 1, 2) &
                            + x(i + 1, j - 1, k, 2) - x(i - 1, j - 1, k, 2) &
                            + x(i + 1, j, k - 1, 2) - x(i - 1, j, k - 1, 2) &
                            + x(i + 1, j, k, 2) - x(i - 1, j, k, 2))
            ssz = eighth * (x(i + 1, j - 1, k - 1, 3) - x(i - 1, j - 1, k - 1, 3) &
                            + x(i + 1, j - 1, k, 3) - x(i - 1, j - 1, k, 3) &
                            + x(i + 1, j, k - 1, 3) - x(i - 1, j, k - 1, 3) &
                            + x(i + 1, j, k, 3) - x(i - 1, j, k, 3))

            ! Determine the length of this vector and create the
            ! unit normal.

            snrm = one / sqrt(ssx * ssx + ssy * ssy + ssz * ssz)
            ssx = snrm * ssx
            ssy = snrm * ssy
            ssz = snrm * ssz

            ! Correct the gradients.

            corr = u_x * ssx + u_y * ssy + u_z * ssz &
                   - (w(i + 1, j, k, ivx) - w(i, j, k, ivx)) * snrm
            u_x = u_x - corr * ssx
            u_y = u_y - corr * ssy
            u_z = u_z - corr * ssz

            corr = v_x * ssx + v_y * ssy + v_z * ssz &
                   - (w(i + 1, j, k, ivy) - w(i, j, k, ivy)) * snrm
            v_x = v_x - corr * ssx
            v_y = v_y - corr * ssy
            v_z = v_z - corr * ssz

            corr = w_x * ssx + w_y * ssy + w_z * ssz &
                   - (w(i + 1, j, k, ivz) - w(i, j, k, ivz)) * snrm
            w_x = w_x - corr * ssx
            w_y = w_y - corr * ssy
            w_z = w_z - corr * ssz

            corr = q_x * ssx + q_y * ssy + q_z * ssz &
                   + (aa(i + 1, j, k) - aa(i, j, k)) * snrm
            q_x = q_x - corr * ssx
            q_y = q_y - corr * ssy
            q_z = q_z - corr * ssz

            ! Compute the stress tensor and the heat flux vector.
            ! We remove the viscosity from the stress tensor (tau)
            ! to define tauS since we still need to separate between
            ! laminar and turbulent stress for QCR.
            ! Therefore, laminar tau = mue*tauS, turbulent
            ! tau = mue*tauS, and total tau = mut*tauS.

            fracDiv = twoThird * (u_x + v_y + w_z)

            tauxxS = two * u_x - fracDiv
            tauyyS = two * v_y - fracDiv
            tauzzS = two * w_z - fracDiv

            tauxyS = u_y + v_x
            tauxzS = u_z + w_x
            tauyzS = v_z + w_y

            q_x = heatCoef * q_x
            q_y = heatCoef * q_y
            q_z = heatCoef * q_z

            ! Add QCR corrections if necessary
            if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x * u_x + u_y * u_y + u_z * u_z + &
                           v_x * v_x + v_y * v_y + v_z * v_z + &
                           w_x * w_x + w_y * w_y + w_z * w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue * Ccr1 / den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact * (Wxy * tauxyS + Wxz * tauxzS) * two
                eyy = fact * (Wyx * tauxyS + Wyz * tauyzS) * two
                ezz = fact * (Wzx * tauxzS + Wzy * tauyzS) * two

                exy = fact * (Wxy * tauyyS + Wxz * tauyzS + &
                              Wyx * tauxxS + Wyz * tauxzS)
                exz = fact * (Wxy * tauyzS + Wxz * tauzzS + &
                              Wzx * tauxxS + Wzy * tauxyS)
                eyz = fact * (Wyx * tauxzS + Wyz * tauzzS + &
                              Wzx * tauxyS + Wzy * tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut * tauxxS - exx
                tauyy = mut * tauyyS - eyy
                tauzz = mut * tauzzS - ezz
                tauxy = mut * tauxyS - exy
                tauxz = mut * tauxzS - exz
                tauyz = mut * tauyzS - eyz

            else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut * tauxxS
                tauyy = mut * tauyyS
                tauzz = mut * tauzzS
                tauxy = mut * tauxyS
                tauxz = mut * tauxzS
                tauyz = mut * tauyzS

            end if

            ! Compute the average velocities for the face. Remember that
            ! the velocities are stored and not the momentum.

            ubar = half * (w(i, j, k, ivx) + w(i + 1, j, k, ivx))
            vbar = half * (w(i, j, k, ivy) + w(i + 1, j, k, ivy))
            wbar = half * (w(i, j, k, ivz) + w(i + 1, j, k, ivz))

            ! Compute the viscous fluxes for this i-face.

            fmx = tauxx * si(i, j, k, 1) + tauxy * si(i, j, k, 2) &
                  + tauxz * si(i, j, k, 3)
            fmy = tauxy * si(i, j, k, 1) + tauyy * si(i, j, k, 2) &
                  + tauyz * si(i, j, k, 3)
            fmz = tauxz * si(i, j, k, 1) + tauyz * si(i, j, k, 2) &
                  + tauzz * si(i, j, k, 3)
            frhoE = (ubar * tauxx + vbar * tauxy + wbar * tauxz) * si(i, j, k, 1) &
                    + (ubar * tauxy + vbar * tauyy + wbar * tauyz) * si(i, j, k, 2) &
                    + (ubar * tauxz + vbar * tauyz + wbar * tauzz) * si(i, j, k, 3) &
                    - q_x * si(i, j, k, 1) - q_y * si(i, j, k, 2) - q_z * si(i, j, k, 3)

            ! Update the residuals of cell i and i+1.
            tmp = atomicsub(fw(i, j, k, imx), fmx)
            tmp = atomicsub(fw(i, j, k, imy), fmy)
            tmp = atomicsub(fw(i, j, k, imz), fmz)
            tmp = atomicsub(fw(i, j, k, irhoE), frhoE)
            ! fw(i, j, k, imx) = fw(i, j, k, imx) - fmx
            ! fw(i, j, k, imy) = fw(i, j, k, imy) - fmy
            ! fw(i, j, k, imz) = fw(i, j, k, imz) - fmz
            ! fw(i, j, k, irhoE) = fw(i, j, k, irhoE) - frhoE
            
            tmp = atomicadd(fw(i + 1, j, k, imx), fmx)
            tmp = atomicadd(fw(i + 1, j, k, imy), fmy)
            tmp = atomicadd(fw(i + 1, j, k, imz), fmz)
            tmp = atomicadd(fw(i + 1, j, k, irhoE), frhoE)
            ! fw(i + 1, j, k, imx) = fw(i + 1, j, k, imx) + fmx
            ! fw(i + 1, j, k, imy) = fw(i + 1, j, k, imy) + fmy
            ! fw(i + 1, j, k, imz) = fw(i + 1, j, k, imz) + fmz
            ! fw(i + 1, j, k, irhoE) = fw(i + 1, j, k, irhoE) + frhoE

        end if
    end subroutine viscousFlux_v1

     ! Miles
    attributes(global) subroutine sumDwandFw_v1

      use constants,       only: zero
      use cudaFlowVarRefState, only: nwf, nt1, nt2
      use precision,       only: intType, realType

      implicit none

      integer(kind=intType) :: i, j, k, l, nTurb
      real(kind=realType)   :: rBlank

      i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + 1  ! starting at 2
      j = (blockIdx%y - 1) * blockDim%y + threadIdx%y + 1
      k = (blockIdx%z - 1) * blockDim%z + threadIdx%z + 1
      !loop 2 to il for cell centers
      if (i>=2 .AND. i <= il .AND. j>=2 .AND. j <= jl .AND. k>=2 .AND. k <= kl) then

        nTurb = nt2 - nt1 + 1

        do l = 1, nwf
          rblank = max(real(iblank(i, j, k), kind=realType), zero)
          dw(i, j, k, l) = (dw(i, j, k, l) + fw(i, j, k, l)) * rBlank
        end do

      end if

    end subroutine sumDwandFw_v1

    ! ========================================================================================================
    ! NEWER WORK THAT IS CURRENTLY BROKEN
    ! ========================================================================================================
    !all non flux subroutines 
    attributes(global) subroutine computeSpeedOfSoundSquared_v2
        use constants, only: two, third,irho,itu1
        implicit none
        integer :: i,j,k,dom,sps

        ! Hard-coded for now
        dom = 1
        sps = 1
        
        ! --- GPU thread indices ---
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z
        
        if (i <= cudaDoms(dom,sps)%ie .AND. j <= cudaDoms(dom,sps)%je .AND. k <= cudaDoms(dom,sps)%ke) then
            !TODO: this will not work if TKE is present check 
            !computeSpeedOfSoundSquared in blockette.F90 need
            !check to correct for K  
            ! cudaDoms(dom,sps)%aa(i,j,k) = cudaDoms(dom,sps)%gamma(i,j,k)*cudaDoms(dom,sps)%p(i,j,k) &
            !                                             / cudaDoms(dom,sps)%w(i, j, k, irho
            cudaDoms(dom,sps)%aa(i,j,k) = (cudaDoms(dom,sps)%gamma(i,j,k)*cudaDoms(dom,sps)%p(i,j,k)) &
                                                          / cudaDoms(dom,sps)%w(i, j, k, irho)

        end if

    end subroutine computeSpeedOfSoundSquared_v2

    !alex
    attributes(global) subroutine allNodalGradients_v2
        use constants, only: zero,one,fourth,ivx,ivy,ivz
        use precision, only: realType
        implicit none

        real(kind=realType) :: a2, oVol, uBar,vBar,wBar,sx,sy,sz
        real(kind=realType) :: tmp

        integer(kind=intType) :: i, j, k,dom,sps
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z

        ! Hard-coded for now
        dom = 1
        sps = 1

        ! First part. Contribution in the k-direction.
        ! The contribution is scattered to both the left and right node
        ! in k-direction.
        
        !loop k 1 to cudaDoms(dom,sps)%ke, j 1 to cudaDoms(dom,sps)%jl and i 1 to cudaDoms(dom,sps)%il
        if (i <= cudaDoms(dom,sps)%il .AND. j <= cudaDoms(dom,sps)%jl .AND. k <= cudaDoms(dom,sps)%ke) then
            ! Compute 8 times the average normal for this part of
            ! the control volume. The factor 8 is taken care of later
            ! on when the division by the volume takes place.
            sx = cudaDoms(dom,sps)%sK(i, j, k - 1, 1) + cudaDoms(dom,sps)%sK(i + 1, j, k - 1, 1) &
                + cudaDoms(dom,sps)%sK(i, j + 1, k - 1, 1) + cudaDoms(dom,sps)%sK(i + 1, j + 1, k - 1, 1) &
                + cudaDoms(dom,sps)%sK(i, j, k, 1) + cudaDoms(dom,sps)%sK(i + 1, j, k, 1) &
                + cudaDoms(dom,sps)%sK(i, j + 1, k, 1) + cudaDoms(dom,sps)%sK(i + 1, j + 1, k, 1)
            sy = cudaDoms(dom,sps)%sK(i, j, k - 1, 2) + cudaDoms(dom,sps)%sK(i + 1, j, k - 1, 2) &
                + cudaDoms(dom,sps)%sK(i, j + 1, k - 1, 2) + cudaDoms(dom,sps)%sK(i + 1, j + 1, k - 1, 2) &
                + cudaDoms(dom,sps)%sK(i, j, k, 2) + cudaDoms(dom,sps)%sK(i + 1, j, k, 2) &
                + cudaDoms(dom,sps)%sK(i, j + 1, k, 2) + cudaDoms(dom,sps)%sK(i + 1, j + 1, k, 2)
            sz = cudaDoms(dom,sps)%sK(i, j, k - 1, 3) + cudaDoms(dom,sps)%sK(i + 1, j, k - 1, 3) &
                + cudaDoms(dom,sps)%sK(i, j + 1, k - 1, 3) + cudaDoms(dom,sps)%sK(i + 1, j + 1, k - 1, 3) &
                + cudaDoms(dom,sps)%sK(i, j, k, 3) + cudaDoms(dom,sps)%sK(i + 1, j, k, 3) &
                + cudaDoms(dom,sps)%sK(i, j + 1, k, 3) + cudaDoms(dom,sps)%sK(i + 1, j + 1, k, 3)
            ! Compute the average velocities and speed of sound squared
            ! for this integration point. Node that these variables are
            ! stored in w(ivx), w(ivy), w(ivz) and p.

            ubar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivx) + cudaDoms(dom,sps)%w(i + 1, j, k, ivx) &
                            + cudaDoms(dom,sps)%w(i, j + 1, k, ivx) + cudaDoms(dom,sps)%w(i + 1, j + 1, k, ivx))
            vbar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivy) + cudaDoms(dom,sps)%w(i + 1, j, k, ivy) &
                            + cudaDoms(dom,sps)%w(i, j + 1, k, ivy) + cudaDoms(dom,sps)%w(i + 1, j + 1, k, ivy))
            wbar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivz) + cudaDoms(dom,sps)%w(i + 1, j, k, ivz) &
                            + cudaDoms(dom,sps)%w(i, j + 1, k, ivz) + cudaDoms(dom,sps)%w(i + 1, j + 1, k, ivz))

            a2 = fourth * (cudaDoms(dom,sps)%aa(i, j, k) + cudaDoms(dom,sps)%aa(i + 1, j, k) + &
                         cudaDoms(dom,sps)%aa(i, j + 1, k) + cudaDoms(dom,sps)%aa(i + 1, j + 1, k))

            ! print *,"k:", i,j,k, a2
            
            ! Add the contributions to the surface integral to the node
            ! j-1 and substract it from the node j. For the heat flux it
            ! is reversed, because the negative of the gradient of the
            ! speed of sound must be computed.
            if (k > 1) then

                tmp = atomicadd(cudaDoms(dom,sps)%ux(i, j, k - 1), ubar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%uy(i, j, k - 1), ubar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%uz(i, j, k - 1), ubar * sz)
                
                tmp = atomicadd(cudaDoms(dom,sps)%vx(i, j, k - 1), vbar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%vy(i, j, k - 1), vbar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%vz(i, j, k - 1), vbar * sz)

                tmp = atomicadd(cudaDoms(dom,sps)%wx(i, j, k - 1), wbar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%wy(i, j, k - 1), wbar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%wz(i, j, k - 1), wbar * sz)
                
                tmp = atomicsub(cudaDoms(dom,sps)%qx(i, j, k - 1), a2 * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%qy(i, j, k - 1), a2 * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%qz(i, j, k - 1), a2 * sz)

                ! cudaDoms(dom,sps)%ux(i, j, k - 1) = cudaDoms(dom,sps)%ux(i, j, k - 1) + ubar * sx
                ! cudaDoms(dom,sps)%uy(i, j, k - 1) = cudaDoms(dom,sps)%uy(i, j, k - 1) + ubar * sy
                ! cudaDoms(dom,sps)%uz(i, j, k - 1) = cudaDoms(dom,sps)%uz(i, j, k - 1) + ubar * sz

                ! cudaDoms(dom,sps)%vx(i, j, k - 1) = cudaDoms(dom,sps)%vx(i, j, k - 1) + vbar * sx
                ! cudaDoms(dom,sps)%vy(i, j, k - 1) = cudaDoms(dom,sps)%vy(i, j, k - 1) + vbar * sy
                ! cudaDoms(dom,sps)%vz(i, j, k - 1) = cudaDoms(dom,sps)%vz(i, j, k - 1) + vbar * sz

                ! cudaDoms(dom,sps)%wx(i, j, k - 1) = cudaDoms(dom,sps)%wx(i, j, k - 1) + wbar * sx
                ! cudaDoms(dom,sps)%wy(i, j, k - 1) = cudaDoms(dom,sps)%wy(i, j, k - 1) + wbar * sy
                ! cudaDoms(dom,sps)%wz(i, j, k - 1) = cudaDoms(dom,sps)%wz(i, j, k - 1) + wbar * sz

                ! cudaDoms(dom,sps)%qx(i, j, k - 1) = cudaDoms(dom,sps)%qx(i, j, k - 1) - a2 * sx
                ! cudaDoms(dom,sps)%qy(i, j, k - 1) = cudaDoms(dom,sps)%qy(i, j, k - 1) - a2 * sy
                ! cudaDoms(dom,sps)%qz(i, j, k - 1) = cudaDoms(dom,sps)%qz(i, j, k - 1) - a2 * sz

            end if 

            if (k < cudaDoms(dom,sps)%ke) then

                tmp = atomicsub(cudaDoms(dom,sps)%ux(i, j, k), ubar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%uy(i, j, k), ubar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%uz(i, j, k), ubar * sz)

                tmp = atomicsub(cudaDoms(dom,sps)%vx(i, j, k), vbar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%vy(i, j, k), vbar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%vz(i, j, k), vbar * sz)

                tmp = atomicsub(cudaDoms(dom,sps)%wx(i, j, k), wbar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%wy(i, j, k), wbar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%wz(i, j, k), wbar * sz)
                
                tmp = atomicadd(cudaDoms(dom,sps)%qx(i, j, k), a2 * sx) 
                tmp = atomicadd(cudaDoms(dom,sps)%qy(i, j, k), a2 * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%qz(i, j, k), a2 * sz)

                ! cudaDoms(dom,sps)%ux(i, j, k) = cudaDoms(dom,sps)%ux(i, j, k) - ubar * sx
                ! cudaDoms(dom,sps)%uy(i, j, k) = cudaDoms(dom,sps)%uy(i, j, k) - ubar * sy
                ! cudaDoms(dom,sps)%uz(i, j, k) = cudaDoms(dom,sps)%uz(i, j, k) - ubar * sz

                ! cudaDoms(dom,sps)%vx(i, j, k) = cudaDoms(dom,sps)%vx(i, j, k) - vbar * sx
                ! cudaDoms(dom,sps)%vy(i, j, k) = cudaDoms(dom,sps)%vy(i, j, k) - vbar * sy
                ! cudaDoms(dom,sps)%vz(i, j, k) = cudaDoms(dom,sps)%vz(i, j, k) - vbar * sz

                ! cudaDoms(dom,sps)%wx(i, j, k) = cudaDoms(dom,sps)%wx(i, j, k) - wbar * sx
                ! cudaDoms(dom,sps)%wy(i, j, k) = cudaDoms(dom,sps)%wy(i, j, k) - wbar * sy
                ! cudaDoms(dom,sps)%wz(i, j, k) = cudaDoms(dom,sps)%wz(i, j, k) - wbar * sz

                ! cudaDoms(dom,sps)%qx(i, j, k) = cudaDoms(dom,sps)%qx(i, j, k) + a2 * sx
                ! cudaDoms(dom,sps)%qy(i, j, k) = cudaDoms(dom,sps)%qy(i, j, k) + a2 * sy
                ! cudaDoms(dom,sps)%qz(i, j, k) = cudaDoms(dom,sps)%qz(i, j, k) + a2 * sz

            end if 
        end if


        ! Second part. Contribution in the j-direction.
        ! The contribution is scattered to both the left and right node
        ! in j-direction.
        !loop k 1 to cudaDoms(dom,sps)%kl j 1 to cudaDoms(dom,sps)%je and i 1 to cudaDoms(dom,sps)%il
        if (i <= cudaDoms(dom,sps)%il .AND. j <= cudaDoms(dom,sps)%je .AND. k <= cudaDoms(dom,sps)%kl) then
            ! Compute 8 times the average normal for this part of
            ! the control volume. The factor 8 is taken care of later
            ! on when the division by the volume takes place.
            sx = cudaDoms(dom,sps)%sJ(i, j - 1, k, 1) + cudaDoms(dom,sps)%sJ(i + 1, j - 1, k, 1) &
                    + cudaDoms(dom,sps)%sJ(i, j - 1, k + 1, 1) + cudaDoms(dom,sps)%sJ(i + 1, j - 1, k + 1, 1) &
                    + cudaDoms(dom,sps)%sJ(i, j, k, 1) + cudaDoms(dom,sps)%sJ(i + 1, j, k, 1) &
                    + cudaDoms(dom,sps)%sJ(i, j, k + 1, 1) + cudaDoms(dom,sps)%sJ(i + 1, j, k + 1, 1)
            sy = cudaDoms(dom,sps)%sJ(i, j - 1, k, 2) + cudaDoms(dom,sps)%sJ(i + 1, j - 1, k, 2) &
                    + cudaDoms(dom,sps)%sJ(i, j - 1, k + 1, 2) + cudaDoms(dom,sps)%sJ(i + 1, j - 1, k + 1, 2) &
                    + cudaDoms(dom,sps)%sJ(i, j, k, 2) + cudaDoms(dom,sps)%sJ(i + 1, j, k, 2) &
                    + cudaDoms(dom,sps)%sJ(i, j, k + 1, 2) + cudaDoms(dom,sps)%sJ(i + 1, j, k + 1, 2)
            sz = cudaDoms(dom,sps)%sJ(i, j - 1, k, 3) + cudaDoms(dom,sps)%sJ(i + 1, j - 1, k, 3) &
                    + cudaDoms(dom,sps)%sJ(i, j - 1, k + 1, 3) + cudaDoms(dom,sps)%sJ(i + 1, j - 1, k + 1, 3) &
                    + cudaDoms(dom,sps)%sJ(i, j, k, 3) + cudaDoms(dom,sps)%sJ(i + 1, j, k, 3) &
                    + cudaDoms(dom,sps)%sJ(i, j, k + 1, 3) + cudaDoms(dom,sps)%sJ(i + 1, j, k + 1, 3)
            ! Compute the average velocities and speed of sound squared
            ! for this integration point. Node that these variables are
            ! stored in w(ivx), w(ivy), w(ivz) and p.
            ubar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivx) + cudaDoms(dom,sps)%w(i + 1, j, k, ivx) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivx) + cudaDoms(dom,sps)%w(i + 1, j, k + 1, ivx))
            vbar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivy) + cudaDoms(dom,sps)%w(i + 1, j, k, ivy) &
                                + cudaDoms(dom,sps)%w(i, j, k + 1, ivy) + cudaDoms(dom,sps)%w(i + 1, j, k + 1, ivy))
            wbar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivz) + cudaDoms(dom,sps)%w(i + 1, j, k, ivz) &
                                + cudaDoms(dom,sps)%w(i, j, k + 1, ivz) + cudaDoms(dom,sps)%w(i + 1, j, k + 1, ivz))

            a2 = fourth * (cudaDoms(dom,sps)%aa(i, j, k) + cudaDoms(dom,sps)%aa(i + 1, j, k) + &
                            cudaDoms(dom,sps)%aa(i, j, k + 1) + cudaDoms(dom,sps)%aa(i + 1, j, k + 1))
            ! print *,"j:", i,j,k, a2
            ! Add the contributions to the surface integral to the node
            ! j-1 and substract it from the node j. For the heat flux it
            ! is reversed, because the negative of the gradient of the
            ! speed of sound must be computed.
    
            if (j > 1) then

                tmp = atomicadd(cudaDoms(dom,sps)%ux(i, j - 1, k), ubar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%uy(i, j - 1, k), ubar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%uz(i, j - 1, k), ubar * sz)

                tmp = atomicadd(cudaDoms(dom,sps)%vx(i, j - 1, k), vbar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%vy(i, j - 1, k), vbar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%vz(i, j - 1, k), vbar * sz)

                tmp = atomicadd(cudaDoms(dom,sps)%wx(i, j - 1, k), wbar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%wy(i, j - 1, k), wbar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%wz(i, j - 1, k), wbar * sz)
                
                tmp = atomicsub(cudaDoms(dom,sps)%qx(i, j - 1, k), a2 * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%qy(i, j - 1, k), a2 * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%qz(i, j - 1, k), a2 * sz)

                ! cudaDoms(dom,sps)%ux(i, j - 1, k) = cudaDoms(dom,sps)%ux(i, j - 1, k) + ubar * sx
                ! cudaDoms(dom,sps)%uy(i, j - 1, k) = cudaDoms(dom,sps)%uy(i, j - 1, k) + ubar * sy
                ! cudaDoms(dom,sps)%uz(i, j - 1, k) = cudaDoms(dom,sps)%uz(i, j - 1, k) + ubar * sz

                ! cudaDoms(dom,sps)%vx(i, j - 1, k) = cudaDoms(dom,sps)%vx(i, j - 1, k) + vbar * sx
                ! cudaDoms(dom,sps)%vy(i, j - 1, k) = cudaDoms(dom,sps)%vy(i, j - 1, k) + vbar * sy
                ! cudaDoms(dom,sps)%vz(i, j - 1, k) = cudaDoms(dom,sps)%vz(i, j - 1, k) + vbar * sz

                ! cudaDoms(dom,sps)%wx(i, j - 1, k) = cudaDoms(dom,sps)%wx(i, j - 1, k) + wbar * sx
                ! cudaDoms(dom,sps)%wy(i, j - 1, k) = cudaDoms(dom,sps)%wy(i, j - 1, k) + wbar * sy
                ! cudaDoms(dom,sps)%wz(i, j - 1, k) = cudaDoms(dom,sps)%wz(i, j - 1, k) + wbar * sz

                ! cudaDoms(dom,sps)%qx(i, j - 1, k) = cudaDoms(dom,sps)%qx(i, j - 1, k) - a2 * sx
                ! cudaDoms(dom,sps)%qy(i, j - 1, k) = cudaDoms(dom,sps)%qy(i, j - 1, k) - a2 * sy
                ! cudaDoms(dom,sps)%qz(i, j - 1, k) = cudaDoms(dom,sps)%qz(i, j - 1, k) - a2 * sz

            end if

            if (j < cudaDoms(dom,sps)%je) then

                tmp = atomicsub(cudaDoms(dom,sps)%ux(i, j, k), ubar * sx) 
                tmp = atomicsub(cudaDoms(dom,sps)%uy(i, j, k), ubar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%uz(i, j, k), ubar * sz)

                tmp = atomicsub(cudaDoms(dom,sps)%vx(i, j, k), vbar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%vy(i, j, k), vbar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%vz(i, j, k), vbar * sz)

                tmp = atomicsub(cudaDoms(dom,sps)%wx(i, j, k), wbar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%wy(i, j, k), wbar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%wz(i, j, k), wbar * sz)
                
                tmp = atomicadd(cudaDoms(dom,sps)%qx(i, j, k), a2 * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%qy(i, j, k), a2 * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%qz(i, j, k), a2 * sz)

                ! cudaDoms(dom,sps)%ux(i, j, k) = cudaDoms(dom,sps)%ux(i, j, k) - ubar * sx
                ! cudaDoms(dom,sps)%uy(i, j, k) = cudaDoms(dom,sps)%uy(i, j, k) - ubar * sy
                ! cudaDoms(dom,sps)%uz(i, j, k) = cudaDoms(dom,sps)%uz(i, j, k) - ubar * sz

                ! cudaDoms(dom,sps)%vx(i, j, k) = cudaDoms(dom,sps)%vx(i, j, k) - vbar * sx
                ! cudaDoms(dom,sps)%vy(i, j, k) = cudaDoms(dom,sps)%vy(i, j, k) - vbar * sy
                ! cudaDoms(dom,sps)%vz(i, j, k) = cudaDoms(dom,sps)%vz(i, j, k) - vbar * sz

                ! cudaDoms(dom,sps)%wx(i, j, k) = cudaDoms(dom,sps)%wx(i, j, k) - wbar * sx
                ! cudaDoms(dom,sps)%wy(i, j, k) = cudaDoms(dom,sps)%wy(i, j, k) - wbar * sy
                ! cudaDoms(dom,sps)%wz(i, j, k) = cudaDoms(dom,sps)%wz(i, j, k) - wbar * sz

                ! cudaDoms(dom,sps)%qx(i, j, k) = cudaDoms(dom,sps)%qx(i, j, k) + a2 * sx
                ! cudaDoms(dom,sps)%qy(i, j, k) = cudaDoms(dom,sps)%qy(i, j, k) + a2 * sy
                ! cudaDoms(dom,sps)%qz(i, j, k) = cudaDoms(dom,sps)%qz(i, j, k) + a2 * sz

            end if 
        end if
        
        ! Third part. Contribution in the i-direction.
        ! The contribution is scattered to both the left and right node
        ! in i-direction.
        
        ! loop k 1 to cudaDoms(dom,sps)%kl j 1 to cudaDoms(dom,sps)%jl and i 1 to cudaDoms(dom,sps)%ie
        if (i <= cudaDoms(dom,sps)%ie .AND. j <= cudaDoms(dom,sps)%jl .AND. k <= cudaDoms(dom,sps)%kl) then
            ! Compute 8 times the average normal for this part of
            ! the control volume. The factor 8 is taken care of later
            ! on when the division by the volume takes place.
            sx = cudaDoms(dom,sps)%sI(i - 1, j, k, 1) + cudaDoms(dom,sps)%sI(i - 1, j + 1, k, 1) &
                    + cudaDoms(dom,sps)%sI(i - 1, j, k + 1, 1) + cudaDoms(dom,sps)%sI(i - 1, j + 1, k + 1, 1) &
                    + cudaDoms(dom,sps)%sI(i, j, k, 1) + cudaDoms(dom,sps)%sI(i, j + 1, k, 1) &
                    + cudaDoms(dom,sps)%sI(i, j, k + 1, 1) + cudaDoms(dom,sps)%sI(i, j + 1, k + 1, 1)
            sy = cudaDoms(dom,sps)%sI(i - 1, j, k, 2) + cudaDoms(dom,sps)%sI(i - 1, j + 1, k, 2) &
                    + cudaDoms(dom,sps)%sI(i - 1, j, k + 1, 2) + cudaDoms(dom,sps)%sI(i - 1, j + 1, k + 1, 2) &
                    + cudaDoms(dom,sps)%sI(i, j, k, 2) + cudaDoms(dom,sps)%sI(i, j + 1, k, 2) &
                    + cudaDoms(dom,sps)%sI(i, j, k + 1, 2) + cudaDoms(dom,sps)%sI(i, j + 1, k + 1, 2)
            sz = cudaDoms(dom,sps)%sI(i - 1, j, k, 3) + cudaDoms(dom,sps)%sI(i - 1, j + 1, k, 3) &
                    + cudaDoms(dom,sps)%sI(i - 1, j, k + 1, 3) + cudaDoms(dom,sps)%sI(i - 1, j + 1, k + 1, 3) &
                    + cudaDoms(dom,sps)%sI(i, j, k, 3) + cudaDoms(dom,sps)%sI(i, j + 1, k, 3) &
                    + cudaDoms(dom,sps)%sI(i, j, k + 1, 3) + cudaDoms(dom,sps)%sI(i, j + 1, k + 1, 3)
            ! Compute the average velocities and speed of sound squared
            ! for this integration point. Node that these variables are
            ! stored in w(ivx), w(ivy), w(ivz) and p.
            ubar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivx) + cudaDoms(dom,sps)%w(i, j + 1, k, ivx) &
                                + cudaDoms(dom,sps)%w(i, j, k + 1, ivx) + cudaDoms(dom,sps)%w(i, j + 1, k + 1, ivx))
            vbar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivy) + cudaDoms(dom,sps)%w(i, j + 1, k, ivy) &
                                + cudaDoms(dom,sps)%w(i, j, k + 1, ivy) + cudaDoms(dom,sps)%w(i, j + 1, k + 1, ivy))
            wbar = fourth * (cudaDoms(dom,sps)%w(i, j, k, ivz) + cudaDoms(dom,sps)%w(i, j + 1, k, ivz) &
                                + cudaDoms(dom,sps)%w(i, j, k + 1, ivz) + cudaDoms(dom,sps)%w(i, j + 1, k + 1, ivz))

            a2 = fourth * (cudaDoms(dom,sps)%aa(i, j, k) + cudaDoms(dom,sps)%aa(i, j + 1, k) + &
                        cudaDoms(dom,sps)%aa(i, j, k + 1) + cudaDoms(dom,sps)%aa(i, j + 1, k + 1))
            ! print *,"i:", i,j,k, a2
            ! Add the contributions to the surface integral to the node
            ! j-1 and substract it from the node j. For the heat flux it
            ! is reversed, because the negative of the gradient of the
            ! speed of sound must be computed.

            if (i > 1) then

                tmp = atomicadd(cudaDoms(dom,sps)%ux(i - 1, j, k), ubar * sx) 
                tmp = atomicadd(cudaDoms(dom,sps)%uy(i - 1, j, k), ubar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%uz(i - 1, j, k), ubar * sz)

                tmp = atomicadd(cudaDoms(dom,sps)%vx(i - 1, j, k), vbar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%vy(i - 1, j, k), vbar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%vz(i - 1, j, k), vbar * sz)

                tmp = atomicadd(cudaDoms(dom,sps)%wx(i - 1, j, k), wbar * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%wy(i - 1, j, k), wbar * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%wz(i - 1, j, k), wbar * sz)
                
                tmp = atomicsub(cudaDoms(dom,sps)%qx(i - 1, j, k), a2 * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%qy(i - 1, j, k), a2 * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%qz(i - 1, j, k), a2 * sz)

                ! cudaDoms(dom,sps)%ux(i - 1, j, k) = cudaDoms(dom,sps)%ux(i - 1, j, k) + ubar * sx
                ! cudaDoms(dom,sps)%uy(i - 1, j, k) = cudaDoms(dom,sps)%uy(i - 1, j, k) + ubar * sy
                ! cudaDoms(dom,sps)%uz(i - 1, j, k) = cudaDoms(dom,sps)%uz(i - 1, j, k) + ubar * sz

                ! cudaDoms(dom,sps)%vx(i - 1, j, k) = cudaDoms(dom,sps)%vx(i - 1, j, k) + vbar * sx
                ! cudaDoms(dom,sps)%vy(i - 1, j, k) = cudaDoms(dom,sps)%vy(i - 1, j, k) + vbar * sy
                ! cudaDoms(dom,sps)%vz(i - 1, j, k) = cudaDoms(dom,sps)%vz(i - 1, j, k) + vbar * sz

                ! cudaDoms(dom,sps)%wx(i - 1, j, k) = cudaDoms(dom,sps)%wx(i - 1, j, k) + wbar * sx
                ! cudaDoms(dom,sps)%wy(i - 1, j, k) = cudaDoms(dom,sps)%wy(i - 1, j, k) + wbar * sy
                ! cudaDoms(dom,sps)%wz(i - 1, j, k) = cudaDoms(dom,sps)%wz(i - 1, j, k) + wbar * sz

                ! cudaDoms(dom,sps)%qx(i - 1, j, k) = cudaDoms(dom,sps)%qx(i - 1, j, k) - a2 * sx
                ! cudaDoms(dom,sps)%qy(i - 1, j, k) = cudaDoms(dom,sps)%qy(i - 1, j, k) - a2 * sy
                ! cudaDoms(dom,sps)%qz(i - 1, j, k) = cudaDoms(dom,sps)%qz(i - 1, j, k) - a2 * sz

            end if

            if (i < cudaDoms(dom,sps)%ie) then

                tmp = atomicsub(cudaDoms(dom,sps)%ux(i, j, k), ubar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%uy(i, j, k), ubar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%uz(i, j, k), ubar * sz)

                tmp = atomicsub(cudaDoms(dom,sps)%vx(i, j, k), vbar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%vy(i, j, k), vbar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%vz(i, j, k), vbar * sz)

                tmp = atomicsub(cudaDoms(dom,sps)%wx(i, j, k), wbar * sx)
                tmp = atomicsub(cudaDoms(dom,sps)%wy(i, j, k), wbar * sy)
                tmp = atomicsub(cudaDoms(dom,sps)%wz(i, j, k), wbar * sz)
                
                tmp = atomicadd(cudaDoms(dom,sps)%qx(i, j, k), a2 * sx)
                tmp = atomicadd(cudaDoms(dom,sps)%qy(i, j, k), a2 * sy)
                tmp = atomicadd(cudaDoms(dom,sps)%qz(i, j, k), a2 * sz)

                ! cudaDoms(dom,sps)%ux(i, j, k) = cudaDoms(dom,sps)%ux(i, j, k) - ubar * sx
                ! cudaDoms(dom,sps)%uy(i, j, k) = cudaDoms(dom,sps)%uy(i, j, k) - ubar * sy
                ! cudaDoms(dom,sps)%uz(i, j, k) = cudaDoms(dom,sps)%uz(i, j, k) - ubar * sz

                ! cudaDoms(dom,sps)%vx(i, j, k) = cudaDoms(dom,sps)%vx(i, j, k) - vbar * sx
                ! cudaDoms(dom,sps)%vy(i, j, k) = cudaDoms(dom,sps)%vy(i, j, k) - vbar * sy
                ! cudaDoms(dom,sps)%vz(i, j, k) = cudaDoms(dom,sps)%vz(i, j, k) - vbar * sz

                ! cudaDoms(dom,sps)%wx(i, j, k) = cudaDoms(dom,sps)%wx(i, j, k) - wbar * sx
                ! cudaDoms(dom,sps)%wy(i, j, k) = cudaDoms(dom,sps)%wy(i, j, k) - wbar * sy
                ! cudaDoms(dom,sps)%wz(i, j, k) = cudaDoms(dom,sps)%wz(i, j, k) - wbar * sz

                ! cudaDoms(dom,sps)%qx(i, j, k) = cudaDoms(dom,sps)%qx(i, j, k) + a2 * sx
                ! cudaDoms(dom,sps)%qy(i, j, k) = cudaDoms(dom,sps)%qy(i, j, k) + a2 * sy
                ! cudaDoms(dom,sps)%qz(i, j, k) = cudaDoms(dom,sps)%qz(i, j, k) + a2 * sz

            end if
        end if 

        ! NOTE: DIFFERENT FROM CPU VERSION
        ! need to get rid of this to avoid race condition
        !each cell only scales one node 
        !however need to wait for each nodes accumulation to be done
        !if thread on edge of boundary need to wait for another block to be done
        !therefore split this part into a separate kernel

        ! ! Divide by 8 times the volume to obtain the correct gradients.
        ! if (i <= cudaDoms(dom,sps)%il .AND. j <= cudaDoms(dom,sps)%jl .AND. k <= cudaDoms(dom,sps)%kl) then
        !     ! Compute the inverse of 8 times the volume for this node.
        !     oVol = one / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i, j, k + 1) &
        !                           + cudaDoms(dom,sps)%vol(i + 1, j, k) + cudaDoms(dom,sps)%vol(i + 1, j, k + 1) &
        !                           + cudaDoms(dom,sps)%vol(i, j + 1, k) + cudaDoms(dom,sps)%vol(i, j + 1, k + 1) &
        !                           + cudaDoms(dom,sps)%vol(i + 1, j + 1, k) + cudaDoms(dom,sps)%vol(i + 1, j + 1, k + 1))
        !     ! Compute the correct velocity gradients and "unit" heat
        !     ! fluxes. The velocity gradients are stored in cudaDoms(dom,sps)%ux, etc.
        !     cudaDoms(dom,sps)%ux(i, j, k) = cudaDoms(dom,sps)%ux(i, j, k) * oVol
        !     cudaDoms(dom,sps)%uy(i, j, k) = cudaDoms(dom,sps)%uy(i, j, k) * oVol
        !     cudaDoms(dom,sps)%uz(i, j, k) = cudaDoms(dom,sps)%uz(i, j, k) * oVol

        !     cudaDoms(dom,sps)%vx(i, j, k) = cudaDoms(dom,sps)%vx(i, j, k) * oVol
        !     cudaDoms(dom,sps)%vy(i, j, k) = cudaDoms(dom,sps)%vy(i, j, k) * oVol
        !     cudaDoms(dom,sps)%vz(i, j, k) = cudaDoms(dom,sps)%vz(i, j, k) * oVol

        !     cudaDoms(dom,sps)%wx(i, j, k) = cudaDoms(dom,sps)%wx(i, j, k) * oVol
        !     cudaDoms(dom,sps)%wy(i, j, k) = cudaDoms(dom,sps)%wy(i, j, k) * oVol
        !     cudaDoms(dom,sps)%wz(i, j, k) = cudaDoms(dom,sps)%wz(i, j, k) * oVol

        !     cudaDoms(dom,sps)%qx(i, j, k) = cudaDoms(dom,sps)%qx(i, j, k) * oVol
        !     cudaDoms(dom,sps)%qy(i, j, k) = cudaDoms(dom,sps)%qy(i, j, k) * oVol
        !     cudaDoms(dom,sps)%qz(i, j, k) = cudaDoms(dom,sps)%qz(i, j, k) * oVol
        ! end if 
    end subroutine allNodalGradients_v2

    attributes(global) subroutine scaleNodalGradients_v2
        use constants, only: zero,fourth,ivx,ivy,ivz,one
        use precision, only: realType
        implicit none

        real(kind=realType) :: a2, oVol, uBar,vBar,wBar,sx,sy,sz

        integer(kind=intType) :: i, j, k,dom,sps

        ! Hard-coded for now
        dom = 1
        sps = 1

        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z

        if (i <= cudaDoms(dom,sps)%il .AND. j <= cudaDoms(dom,sps)%jl .AND. k <= cudaDoms(dom,sps)%kl) then
            ! Compute the inverse of 8 times the volume for this node.
            oVol = one / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i, j, k + 1) &
                                  + cudaDoms(dom,sps)%vol(i + 1, j, k) + cudaDoms(dom,sps)%vol(i + 1, j, k + 1) &
                                  + cudaDoms(dom,sps)%vol(i, j + 1, k) + cudaDoms(dom,sps)%vol(i, j + 1, k + 1) &
                                  + cudaDoms(dom,sps)%vol(i + 1, j + 1, k) + cudaDoms(dom,sps)%vol(i + 1, j + 1, k + 1))
            ! Compute the correct velocity gradients and "unit" heat
            ! fluxes. The velocity gradients are stored in cudaDoms(dom,sps)%ux, etc.
            cudaDoms(dom,sps)%ux(i, j, k) = cudaDoms(dom,sps)%ux(i, j, k) * oVol
            cudaDoms(dom,sps)%uy(i, j, k) = cudaDoms(dom,sps)%uy(i, j, k) * oVol
            cudaDoms(dom,sps)%uz(i, j, k) = cudaDoms(dom,sps)%uz(i, j, k) * oVol

            cudaDoms(dom,sps)%vx(i, j, k) = cudaDoms(dom,sps)%vx(i, j, k) * oVol
            cudaDoms(dom,sps)%vy(i, j, k) = cudaDoms(dom,sps)%vy(i, j, k) * oVol
            cudaDoms(dom,sps)%vz(i, j, k) = cudaDoms(dom,sps)%vz(i, j, k) * oVol

            cudaDoms(dom,sps)%wx(i, j, k) = cudaDoms(dom,sps)%wx(i, j, k) * oVol
            cudaDoms(dom,sps)%wy(i, j, k) = cudaDoms(dom,sps)%wy(i, j, k) * oVol
            cudaDoms(dom,sps)%wz(i, j, k) = cudaDoms(dom,sps)%wz(i, j, k) * oVol

            cudaDoms(dom,sps)%qx(i, j, k) = cudaDoms(dom,sps)%qx(i, j, k) * oVol
            cudaDoms(dom,sps)%qy(i, j, k) = cudaDoms(dom,sps)%qy(i, j, k) * oVol
            cudaDoms(dom,sps)%qz(i, j, k) = cudaDoms(dom,sps)%qz(i, j, k) * oVol
        end if 

    end subroutine scaleNodalGradients_v2

    !sabet
    attributes(global) subroutine metrics
        use constants, only: half
        implicit none

        integer(kind=intType) :: i, j, k, l, m, n,dom,sps
        real(kind=realType), dimension(3) :: v1, v2
        real(kind=realType) :: fact

        dom = 1
        sps = 1
        if (cudaDoms(dom,sps)%rightHanded) then
            fact = half
        else
            fact = -half
        end if

        ! Projected areas of cell faces in the i direction.    

        i = (blockIdx%x-1)*blockDim%x + threadIdx%x - 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z

        n = k - 1
        m = j - 1   
        ! Divide by 8 times the volume to obtain the correct gradients.
        if (i <= cudaDoms(dom,sps)%ie .and. j <= cudaDoms(dom,sps)%je .and. k <= cudaDoms(dom,sps)%ke) then

            ! Determine the two diagonal vectors of the face.

            v1(1) = cudaDoms(dom,sps)%x(i, j, n, 1) - cudaDoms(dom,sps)%x(i, m, k, 1)
            v1(2) = cudaDoms(dom,sps)%x(i, j, n, 2) - cudaDoms(dom,sps)%x(i, m, k, 2)
            v1(3) = cudaDoms(dom,sps)%x(i, j, n, 3) - cudaDoms(dom,sps)%x(i, m, k, 3)

            v2(1) = cudaDoms(dom,sps)%x(i, j, k, 1) - cudaDoms(dom,sps)%x(i, m, n, 1)
            v2(2) = cudaDoms(dom,sps)%x(i, j, k, 2) - cudaDoms(dom,sps)%x(i, m, n, 2)
            v2(3) = cudaDoms(dom,sps)%x(i, j, k, 3) - cudaDoms(dom,sps)%x(i, m, n, 3)

            ! The face normal, which is the cross product of the two
            ! diagonal vectors times fact; remember that fact is
            ! either -0.5 or 0.5.

            cudaDoms(dom,sps)%sI(i, j, k, 1) = fact * (v1(2) * v2(3) - v1(3) * v2(2))
            cudaDoms(dom,sps)%sI(i, j, k, 2) = fact * (v1(3) * v2(1) - v1(1) * v2(3))
            cudaDoms(dom,sps)%sI(i, j, k, 3) = fact * (v1(1) * v2(2) - v1(2) * v2(1))

        end if

        ! Projected areas of cell faces in the j direction.    

        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y - 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z

        n = k - 1
        l = i - 1   

        if (i <= cudaDoms(dom,sps)%ie .and. j <= cudaDoms(dom,sps)%je .and. k <= cudaDoms(dom,sps)%ke) then

            ! Determine the two diagonal vectors of the face.

            v1(1) = cudaDoms(dom,sps)%x(i, j, n, 1) - cudaDoms(dom,sps)%x(l, j, k, 1)
            v1(2) = cudaDoms(dom,sps)%x(i, j, n, 2) - cudaDoms(dom,sps)%x(l, j, k, 2)
            v1(3) = cudaDoms(dom,sps)%x(i, j, n, 3) - cudaDoms(dom,sps)%x(l, j, k, 3)

            v2(1) = cudaDoms(dom,sps)%x(l, j, n, 1) - cudaDoms(dom,sps)%x(i, j, k, 1)
            v2(2) = cudaDoms(dom,sps)%x(l, j, n, 2) - cudaDoms(dom,sps)%x(i, j, k, 2)
            v2(3) = cudaDoms(dom,sps)%x(l, j, n, 3) - cudaDoms(dom,sps)%x(i, j, k, 3)

            ! The face normal, which is the cross product of the two
            ! diagonal vectors times fact; remember that fact is
            ! either -0.5 or 0.5.

            cudaDoms(dom,sps)%sJ(i, j, k, 1) = fact * (v1(2) * v2(3) - v1(3) * v2(2))
            cudaDoms(dom,sps)%sJ(i, j, k, 2) = fact * (v1(3) * v2(1) - v1(1) * v2(3))
            cudaDoms(dom,sps)%sJ(i, j, k, 3) = fact * (v1(1) * v2(2) - v1(2) * v2(1))

        end if

        ! Projected areas of cell faces in the k direction.    

        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z - 1

        m = j - 1
        l = i - 1   

        if (i <= cudaDoms(dom,sps)%ie .and. j <= cudaDoms(dom,sps)%je .and. k <= cudaDoms(dom,sps)%ke) then

            ! Determine the two diagonal vectors of the face.

            v1(1) = cudaDoms(dom,sps)%x(i, j, k, 1) - cudaDoms(dom,sps)%x(l, m, k, 1)
            v1(2) = cudaDoms(dom,sps)%x(i, j, k, 2) - cudaDoms(dom,sps)%x(l, m, k, 2)
            v1(3) = cudaDoms(dom,sps)%x(i, j, k, 3) - cudaDoms(dom,sps)%x(l, m, k, 3)

            v2(1) = cudaDoms(dom,sps)%x(l, j, k, 1) - cudaDoms(dom,sps)%x(i, m, k, 1)
            v2(2) = cudaDoms(dom,sps)%x(l, j, k, 2) - cudaDoms(dom,sps)%x(i, m, k, 2)
            v2(3) = cudaDoms(dom,sps)%x(l, j, k, 3) - cudaDoms(dom,sps)%x(i, m, k, 3)

            ! The face normal, which is the cross product of the two
            ! diagonal vectors times fact; remember that fact is
            ! either -0.5 or 0.5.

            cudaDoms(dom,sps)%sK(i, j, k, 1) = fact * (v1(2) * v2(3) - v1(3) * v2(2))
            cudaDoms(dom,sps)%sK(i, j, k, 2) = fact * (v1(3) * v2(1) - v1(1) * v2(3))
            cudaDoms(dom,sps)%sK(i, j, k, 3) = fact * (v1(1) * v2(2) - v1(2) * v2(1))

        end if

    end subroutine metrics

    !sabet
    attributes(global) subroutine timeStep(updateDtl)

        use constants, only: eps, irho, ivx, ivy, ivz, irhoE, one, half, zero
        use cudaFlowVarRefState, only: pInfCorr, rhoInf, gammaInf, viscous
        use cudaInputDiscretization, only: adis, acousticScaleFactor

        implicit none

        logical,value,intent(in) :: updateDtl

        ! Local parameters.
        real(kind=realType), parameter :: b = 2.0_realType

        ! Variables for spectral Radius
        real(kind=realType) :: plim, rlim, clim2
        real(kind=realType) :: cc2, qsi, qsj, qsk, sx, sy, sz, rmu
        real(kind=realType) :: ri, rj, rk, rij, rjk, rki
        real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
        real(kind=realType) :: sFace, tmp, uux, uuy, uuz
        logical :: doScaling, updateDt
        integer(kind=intType) :: i, j, k, dom, sps

        dom = 1 
        sps = 1
        !no nonger optional argument must be passed in
        updateDt = updateDtl


        ! Set the value of plim. To be fully consistent this must have
        ! the dimension of a pressure. Therefore a fraction of pInfCorr
        ! is used. Idem for rlim; compute clim2 as well.

        plim = 0.001_realType * pInfCorr
        rlim = 0.001_realType * rhoInf
        clim2 = 0.000001_realType * gammaInf * pInfCorr / rhoInf
        doScaling = .True.

        ! Initialize sFace to zero. This value will be used if the
        ! block is not moving.

        sFace = zero
        !
        !           Inviscid contribution, depending on the preconditioner.
        !           Compute the cell centered values of the spectral radii.
        !
        ! Note:  DON'T change the ranges for i. It will mess up dtl.
        ! we don't copy the spectral-radii and dtl to keep the
        ! code simple, therefore this loop needs the full single
        ! halo range.

        i = (blockIdx%x-1)*blockDim%x + threadIdx%x
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z
        !loop 1 to cudaDoms(dom,sps)%ie
        if (i <= cudaDoms(dom,sps)%ie .and. j <= cudaDoms(dom,sps)%je .and. k <= cudaDoms(dom,sps)%ke) then

            ! Compute the velocities and speed of sound squared.

            uux = cudaDoms(dom,sps)%w(i, j, k, ivx)
            uuy = cudaDoms(dom,sps)%w(i, j, k, ivy)
            uuz = cudaDoms(dom,sps)%w(i, j, k, ivz)
            cc2 = cudaDoms(dom,sps)%gamma(i, j, k) * cudaDoms(dom,sps)%P(i, j, k) / cudaDoms(dom,sps)%w(i, j, k, irho)
            cc2 = max(cc2, clim2)

            ! Set the dot product of the grid velocity and the
            ! normal in i-direction for a moving face. To avoid
            ! a number of multiplications by 0.5 simply the sum
            ! is taken.

            sFace = cudaDoms(dom,sps)%sFaceI(i - 1, j, k) + cudaDoms(dom,sps)%sFaceI(i, j, k)

            ! Spectral radius in i-direction.

            sx = cudaDoms(dom,sps)%sI(i - 1, j, k, 1) + cudaDoms(dom,sps)%sI(i, j, k, 1)
            sy = cudaDoms(dom,sps)%sI(i - 1, j, k, 2) + cudaDoms(dom,sps)%sI(i, j, k, 2)
            sz = cudaDoms(dom,sps)%sI(i - 1, j, k, 3) + cudaDoms(dom,sps)%sI(i, j, k, 3)

            qsi = uux * sx + uuy * sy + uuz * sz - sFace

            ri = half * (abs(qsi) &
                            + acousticScaleFactor * sqrt(cc2 * (sx**2 + sy**2 + sz**2)))

            ! The grid velocity in j-direction.
            sFace = cudaDoms(dom,sps)%sFaceJ(i, j - 1, k) + cudaDoms(dom,sps)%sFaceJ(i, j, k)

            ! Spectral radius in j-direction.

            sx = cudaDoms(dom,sps)%sJ(i, j - 1, k, 1) + cudaDoms(dom,sps)%sJ(i, j, k, 1)
            sy = cudaDoms(dom,sps)%sJ(i, j - 1, k, 2) + cudaDoms(dom,sps)%sJ(i, j, k, 2)
            sz = cudaDoms(dom,sps)%sJ(i, j - 1, k, 3) + cudaDoms(dom,sps)%sJ(i, j, k, 3)

            qsj = uux * sx + uuy * sy + uuz * sz - sFace

            rj = half * (abs(qsj) &
                            + acousticScaleFactor * sqrt(cc2 * (sx**2 + sy**2 + sz**2)))

            ! The grid velocity in k-direction.
            sFace = cudaDoms(dom,sps)%sFaceK(i, j, k - 1) + cudaDoms(dom,sps)%sFaceK(i, j, k)

            ! Spectral radius in k-direction.

            sx = cudaDoms(dom,sps)%sK(i, j, k - 1, 1) + cudaDoms(dom,sps)%sK(i, j, k, 1)
            sy = cudaDoms(dom,sps)%sK(i, j, k - 1, 2) + cudaDoms(dom,sps)%sK(i, j, k, 2)
            sz = cudaDoms(dom,sps)%sK(i, j, k - 1, 3) + cudaDoms(dom,sps)%sK(i, j, k, 3)

            qsk = uux * sx + uuy * sy + uuz * sz - sFace

            rk = half * (abs(qsk) &
                            + acousticScaleFactor * sqrt(cc2 * (sx**2 + sy**2 + sz**2)))

            ! Store in tdl if required
            if (updateDt) then
                cudaDoms(dom,sps)%dtl(i, j, k) = ri + rj + rk
            end if

            ! Avoid division by zero by clipping radi, radJ and
            ! radK.

            ri = max(ri, eps)
            rj = max(rj, eps)
            rk = max(rk, eps)

            ! Compute the scaling in the three coordinate
            ! directions.

            rij = (ri / rj)**adis
            rjk = (rj / rk)**adis
            rki = (rk / ri)**adis

            ! Create the scaled versions of the aspect ratios.
            ! Note that the multiplication is done with radi, radJ
            ! and radK, such that the influence of the clipping
            ! is negligible.

            cudaDoms(dom,sps)%radi(i, j, k) = ri * (one + one / rij + rki)
            cudaDoms(dom,sps)%radj(i, j, k) = rj * (one + one / rjk + rij)
            cudaDoms(dom,sps)%radk(i, j, k) = rk * (one + one / rki + rjk)
        end if

        ! The rest is only necessary if the timeStep needs to be computed
        if (updateDt) then

            viscousTerm: if (viscous) then

                ! Loop over the owned cell centers.
                !loop starts at 2
                i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
                j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
                k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
                !loop at 2 to cudaDoms(dom,sps)%il
                if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then

                    ! Compute the effective viscosity coefficient. The
                    ! factor 0.5 is a combination of two things. In the
                    ! standard central discretization of a second
                    ! derivative there is a factor 2 multiplying the
                    ! central node. However in the code below not the
                    ! average but the sum of the left and the right face
                    ! is taken and squared. This leads to a factor 4.
                    ! Combining both effects leads to 0.5. Furthermore,
                    ! it is divided by the volume and density to obtain
                    ! the correct dimensions and multiplied by the
                    ! non-dimensional factor factVis.

                    rmu = cudaDoms(dom,sps)%rlv(i, j, k)
                    rmu = rmu + cudaDoms(dom,sps)%rev(i, j, k)
                    rmu = half * rmu / (cudaDoms(dom,sps)%w(i, j, k, irho) * cudaDoms(dom,sps)%vol(i, j, k))

                    ! Add the viscous contribution in i-direction to the
                    ! (inverse) of the time step.

                    sx = cudaDoms(dom,sps)%sI(i, j, k, 1) + cudaDoms(dom,sps)%sI(i - 1, j, k, 1)
                    sy = cudaDoms(dom,sps)%sI(i, j, k, 2) + cudaDoms(dom,sps)%sI(i - 1, j, k, 2)
                    sz = cudaDoms(dom,sps)%sI(i, j, k, 3) + cudaDoms(dom,sps)%sI(i - 1, j, k, 3)

                    vsi = rmu * (sx * sx + sy * sy + sz * sz)
                    cudaDoms(dom,sps)%dtl(i, j, k) = cudaDoms(dom,sps)%dtl(i, j, k) + vsi

                    ! Add the viscous contribution in j-direction to the
                    ! (inverse) of the time step.

                    sx = cudaDoms(dom,sps)%sJ(i, j, k, 1) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 1)
                    sy = cudaDoms(dom,sps)%sJ(i, j, k, 2) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 2)
                    sz = cudaDoms(dom,sps)%sJ(i, j, k, 3) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 3)

                    vsj = rmu * (sx * sx + sy * sy + sz * sz)
                    cudaDoms(dom,sps)%dtl(i, j, k) = cudaDoms(dom,sps)%dtl(i, j, k) + vsj

                    ! Add the viscous contribution in k-direction to the
                    ! (inverse) of the time step.

                    sx = cudaDoms(dom,sps)%sK(i, j, k, 1) + cudaDoms(dom,sps)%sK(i, j, k - 1, 1)
                    sy = cudaDoms(dom,sps)%sK(i, j, k, 2) + cudaDoms(dom,sps)%sK(i, j, k - 1, 2)
                    sz = cudaDoms(dom,sps)%sK(i, j, k, 3) + cudaDoms(dom,sps)%sK(i, j, k - 1, 3)

                    vsk = rmu * (sx * sx + sy * sy + sz * sz)
                    cudaDoms(dom,sps)%dtl(i, j, k) = cudaDoms(dom,sps)%dtl(i, j, k) + vsk

                end if
            end if viscousTerm
        end if
        ! Skipping time spectral

        ! Currently the inverse of dt/vol is stored in dtl. Invert
        ! this value such that the time step per unit cfl number is
        ! stored and correct in cases of high gradients.

        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
        !loop 2 to cudaDoms(dom,sps)%il (cell centers)
        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then
            dpi = abs(cudaDoms(dom,sps)%P(i + 1, j, k) - two * cudaDoms(dom,sps)%P(i, j, k) + cudaDoms(dom,sps)%P(i - 1, j, k)) &
                    / (cudaDoms(dom,sps)%P(i + 1, j, k) + two * cudaDoms(dom,sps)%P(i, j, k) + cudaDoms(dom,sps)%P(i - 1, j, k) + plim)
            dpj = abs(cudaDoms(dom,sps)%P(i, j + 1, k) - two * cudaDoms(dom,sps)%P(i, j, k) + cudaDoms(dom,sps)%P(i, j - 1, k)) &
                    / (cudaDoms(dom,sps)%P(i, j + 1, k) + two * cudaDoms(dom,sps)%P(i, j, k) + cudaDoms(dom,sps)%P(i, j - 1, k) + plim)
            dpk = abs(cudaDoms(dom,sps)%P(i, j, k + 1) - two * cudaDoms(dom,sps)%P(i, j, k) + cudaDoms(dom,sps)%P(i, j, k - 1)) &
                    / (cudaDoms(dom,sps)%P(i, j, k + 1) + two * cudaDoms(dom,sps)%P(i, j, k) + cudaDoms(dom,sps)%P(i, j, k - 1) + plim)
            rfl = one / (one + b * (dpi + dpj + dpk))

            cudaDoms(dom,sps)%dtl(i, j, k) = rfl / cudaDoms(dom,sps)%dtl(i, j, k)
        end if
    
    end subroutine timeStep
    
    ! Miles
    attributes(global) subroutine sumDwandFw_v2

      use constants,       only: zero
      use cudaFlowVarRefState, only: nwf, nt1, nt2
      use precision,       only: intType, realType

      implicit none

      integer(kind=intType) :: i, j, k, l, nTurb,dom ,sps
      real(kind=realType)   :: rBlank

      dom = 1
      sps = 1

      i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + 1  ! starting at 2
      j = (blockIdx%y - 1) * blockDim%y + threadIdx%y + 1
      k = (blockIdx%z - 1) * blockDim%z + threadIdx%z + 1
      !loop 2 to cudaDoms(dom,sps)%il for cell centers
      if (i>=2 .AND. i <= cudaDoms(dom,sps)%il .AND. j>=2 .AND. j <= cudaDoms(dom,sps)%jl .AND. k>=2 &
                                                                        .AND. k <= cudaDoms(dom,sps)%kl) then

        nTurb = nt2 - nt1 + 1

        do l = 1, nwf
          rblank = max(real(cudaDoms(dom,sps)%iblank(i, j, k), kind=realType), zero)
          cudaDoms(dom,sps)%dw(i, j, k, l) = (cudaDoms(dom,sps)%dw(i, j, k, l) + cudaDoms(dom,sps)%fw(i, j, k, l)) * rBlank
        end do

      end if

    end subroutine sumDwandFw_v2

    ! Galen
    attributes(global) subroutine resScale()
        use constants
        use cudaFlowVarRefState, only: nwf, nt1, nt2
        use cudaInputIteration, only: turbResScale
        implicit none

        ! Local Variables
        integer(kind=intType) :: i, j, k, ii, nTurb, dom, sps,l
        real(kind=realType) :: ovol

        dom = 1 
        sps = 1
        ! Get iteration indices from GPU
        ! These loops start at '2'
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x+1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y+1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z+1

        ! Divide through by the reference volume
        nTurb = nt2 - nt1 + 1
        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then
            oVol = one / cudaDoms(dom,sps)%volRef(i, j, k)
            do l = 1,nwf
                cudaDoms(dom,sps)%dw(i, j, k, l) = cudaDoms(dom,sps)%dw(i, j, k, l) * ovol
            end do 
            do l = nt1,nt2
                cudaDoms(dom,sps)%dw(i, j, k, l) = cudaDoms(dom,sps)%dw(i, j, k, l) * ovol * turbResScale(l)
            end do 
        end if
    end subroutine resScale

    attributes(global) subroutine inviscidCentralFlux
        ! ---------------------------------------------
        !               Inviscid central flux
        ! ---------------------------------------------
        use precision, only: realType, intType
        use constants, only: zero, one, two, third, fourth, eighth, ivx, ivy, ivz, irhoE, irho, itu1, imx, imy, imz
        use blockPointers, only: blockIsMoving, nBkGlobal
        use cudaFlowVarRefState, only: timeRef
        use cgnsGrid, only: cgnsDoms
        use cudaInputPhysics, only: equationMode
        implicit none
        !real(kind=realType), intent(in),optional :: wwx,wwy,wwz
        ! Variables for inviscid central flux
        real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
        real(kind=realType) :: pa, vnp, vnm, fs, sFace
        integer(kind=intType) :: i, j, k, dom, sps
        real(kind=realtype) :: tmp  
        ! real(kind=realType) :: wwx, wwy, wwz, rvol

        dom = 1
        sps = 1
        !thread indices start at 1
        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x 
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y 
        k = (blockIdx%z - 1) * blockDim%z + threadIdx%z 

        !ifaces flux
        !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i1 to cudaDoms(dom,sps)%il
        if (k <= cudaDoms(dom,sps)%kl .and. j <= cudaDoms(dom,sps)%jl .and. i <= cudaDoms(dom,sps)%il .and. j>=2 .and. k>=2) then
            ! Set the dot product of the grid velocity and the
            ! normal in i-direction for a moving face.
            sFace = cudaDoms(dom,sps)%sFaceI(i, j, k)
            ! Compute the normal velocities of the left and right state.

            vnp = cudaDoms(dom,sps)%w(i + 1, j, k, ivx) * cudaDoms(dom,sps)%sI(i, j, k, 1) &
                    + cudaDoms(dom,sps)%w(i + 1, j, k, ivy) * cudaDoms(dom,sps)%sI(i, j, k, 2) &
                    + cudaDoms(dom,sps)%w(i + 1, j, k, ivz) * cudaDoms(dom,sps)%sI(i, j, k, 3)
            vnm = cudaDoms(dom,sps)%w(i, j, k, ivx) * cudaDoms(dom,sps)%sI(i, j, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j, k, ivy) * cudaDoms(dom,sps)%sI(i, j, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j, k, ivz) * cudaDoms(dom,sps)%sI(i, j, k, 3)
            ! Set the values of the porosities for this face.
            ! porVel defines the porosity w.r.t. velocity;
            ! porFlux defines the porosity w.r.t. the entire flux.
            ! The latter is only zero for a discontinuous block
            ! boundary that must be treated conservatively.
            ! The default value of porFlux is 0.5, such that the
            ! correct central flux is scattered to both cells.
            ! In case of a boundFlux the normal velocity is set
            ! to sFace.

            porVel = one
            porFlux = half
            if (cudaDoms(dom,sps)%porI(i, j, k) == noFlux) porFlux = zero
            if (cudaDoms(dom,sps)%porI(i, j, k) == boundFlux) then
                porVel = zero
                vnp = sFace
                vnm = sFace
            end if

            ! Incorporate porFlux in porVel.

            porVel = porVel * porFlux

            ! Compute the normal velocities relative to the grid for
            ! the face as well as the mass fluxes.

            qsp = (vnp - sFace) * porVel
            qsm = (vnm - sFace) * porVel

            rqsp = qsp * cudaDoms(dom,sps)%w(i + 1, j, k, irho)
            rqsm = qsm * cudaDoms(dom,sps)%w(i, j, k, irho)

            ! Compute the sum of the pressure multiplied by porFlux.
            ! For the default value of porFlux, 0.5, this leads to
            ! the average pressure.

            pa = porFlux * (cudaDoms(dom,sps)%P(i + 1, j, k) + cudaDoms(dom,sps)%P(i, j, k))

            ! Compute the fluxes and scatter them to the cells
            ! i,j,k and i+1,j,k. Store the density flux in the
            ! mass flow of the appropriate sliding mesh interface.

            fs = rqsp + rqsm
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i + 1, j, k, irho), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, irho), fs)
            ! cudaDoms(dom,sps)%dw(i + 1, j, k, irho) = cudaDoms(dom,sps)%dw(i + 1, j, k, irho) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, irho) = cudaDoms(dom,sps)%dw(i, j, k, irho) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i + 1, j, k, ivx) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivx) &
                + pa * cudaDoms(dom,sps)%sI(i, j, k, 1)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i + 1, j, k, imx), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imx), fs)
            ! cudaDoms(dom,sps)%dw(i + 1, j, k, imx) = cudaDoms(dom,sps)%dw(i + 1, j, k, imx) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imx) = cudaDoms(dom,sps)%dw(i, j, k, imx) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i + 1, j, k, ivy) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivy) &
                + pa * cudaDoms(dom,sps)%sI(i, j, k, 2)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i + 1, j, k, imy), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imy), fs)
            ! cudaDoms(dom,sps)%dw(i + 1, j, k, imy) = cudaDoms(dom,sps)%dw(i + 1, j, k, imy) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imy) = cudaDoms(dom,sps)%dw(i, j, k, imy) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i + 1, j, k, ivz) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivz) &
                + pa * cudaDoms(dom,sps)%sI(i, j, k, 3)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i + 1, j, k, imz), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imz), fs)
            ! cudaDoms(dom,sps)%dw(i + 1, j, k, imz) = cudaDoms(dom,sps)%dw(i + 1, j, k, imz) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imz) = cudaDoms(dom,sps)%dw(i, j, k, imz) + fs

            fs = qsp * cudaDoms(dom,sps)%w(i + 1, j, k, irhoE) + qsm * cudaDoms(dom,sps)%w(i, j, k, irhoE) &
                + porFlux * (vnp * cudaDoms(dom,sps)%P(i + 1, j, k) + vnm * cudaDoms(dom,sps)%P(i, j, k))
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i + 1, j, k, irhoE), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, irhoE), fs)
            ! cudaDoms(dom,sps)%dw(i + 1, j, k, irhoE) = cudaDoms(dom,sps)%dw(i + 1, j, k, irhoE) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, irhoE) = cudaDoms(dom,sps)%dw(i, j, k, irhoE) + fs
        end if
        
        !loop k 2 to cudaDoms(dom,sps)%kl, j 1 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        !j faces flux
        if (k <= cudaDoms(dom,sps)%kl .and. j <= cudaDoms(dom,sps)%jl .and. i <= cudaDoms(dom,sps)%il .and. i>=2 .and. k>=2) then
            ! Set the dot product of the grid velocity and the
            ! normal in j-direction for a moving face.

            sFace = cudaDoms(dom,sps)%sFaceJ(i, j, k)
            ! Compute the normal velocities of the left and right state.

            vnp = cudaDoms(dom,sps)%w(i, j + 1, k, ivx) * cudaDoms(dom,sps)%sJ(i, j, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivy) * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivz) * cudaDoms(dom,sps)%sJ(i, j, k, 3)
            vnm = cudaDoms(dom,sps)%w(i, j, k, ivx) * cudaDoms(dom,sps)%sJ(i, j, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j, k, ivy) * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j, k, ivz) * cudaDoms(dom,sps)%sJ(i, j, k, 3)

            ! Set the values of the porosities for this face.
            ! porVel defines the porosity w.r.t. velocity;
            ! porFlux defines the porosity w.r.t. the entire flux.
            ! The latter is only zero for a discontinuous block
            ! boundary that must be treated conservatively.
            ! The default value of porFlux is 0.5, such that the
            ! correct central flux is scattered to both cells.
            ! In case of a boundFlux the normal velocity is set
            ! to sFace.

            porVel = one
            porFlux = half
            if (cudaDoms(dom,sps)%porJ(i, j, k) == noFlux) porFlux = zero
            if (cudaDoms(dom,sps)%porJ(i, j, k) == boundFlux) then
                porVel = zero
                vnp = sFace
                vnm = sFace
            end if

            ! Incorporate porFlux in porVel.

            porVel = porVel * porFlux

            ! Compute the normal velocities for the face as well as the
            ! mass fluxes.

            qsp = (vnp - sFace) * porVel
            qsm = (vnm - sFace) * porVel

            rqsp = qsp * cudaDoms(dom,sps)%w(i, j + 1, k, irho)
            rqsm = qsm * cudaDoms(dom,sps)%w(i, j, k, irho)

            ! Compute the sum of the pressure multiplied by porFlux.
            ! For the default value of porFlux, 0.5, this leads to
            ! the average pressure.

            pa = porFlux * (cudaDoms(dom,sps)%P(i, j + 1, k) + cudaDoms(dom,sps)%P(i, j, k))

            ! Compute the fluxes and scatter them to the cells
            ! i,j,k and i,j+1,k. Store the density flux in the
            ! mass flow of the appropriate sliding mesh interface.

            fs = rqsp + rqsm
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j + 1, k, irho), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, irho), fs)
            ! cudaDoms(dom,sps)%dw(i, j + 1, k, irho) = cudaDoms(dom,sps)%dw(i, j + 1, k, irho) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, irho) = cudaDoms(dom,sps)%dw(i, j, k, irho) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i, j + 1, k, ivx) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivx) &
                + pa * cudaDoms(dom,sps)%sJ(i, j, k, 1)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j + 1, k, imx), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imx), fs)
            ! cudaDoms(dom,sps)%dw(i, j + 1, k, imx) = cudaDoms(dom,sps)%dw(i, j + 1, k, imx) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imx) = cudaDoms(dom,sps)%dw(i, j, k, imx) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i, j + 1, k, ivy) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivy) &
                + pa * cudaDoms(dom,sps)%sJ(i, j, k, 2)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j + 1, k, imy), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imy), fs)
            ! cudaDoms(dom,sps)%dw(i, j + 1, k, imy) = cudaDoms(dom,sps)%dw(i, j + 1, k, imy) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imy) = cudaDoms(dom,sps)%dw(i, j, k, imy) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i, j + 1, k, ivz) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivz) &
                + pa * cudaDoms(dom,sps)%sJ(i, j, k, 3)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j + 1, k, imz), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imz), fs)
            ! cudaDoms(dom,sps)%dw(i, j + 1, k, imz) = cudaDoms(dom,sps)%dw(i, j + 1, k, imz) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imz) = cudaDoms(dom,sps)%dw(i, j, k, imz) + fs

            fs = qsp * cudaDoms(dom,sps)%w(i, j + 1, k, irhoE) + qsm * cudaDoms(dom,sps)%w(i, j, k, irhoE) &
                + porFlux * (vnp * cudaDoms(dom,sps)%P(i, j + 1, k) + vnm * cudaDoms(dom,sps)%P(i, j, k))
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j + 1, k, irhoE), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, irhoE), fs)
            ! cudaDoms(dom,sps)%dw(i, j + 1, k, irhoE) = cudaDoms(dom,sps)%dw(i, j + 1, k, irhoE) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, irhoE) = cudaDoms(dom,sps)%dw(i, j, k, irhoE) + fs
        end if

        !k face flux
        !loop k 1 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (k <= cudaDoms(dom,sps)%kl .and. j <= cudaDoms(dom,sps)%jl .and. i <= cudaDoms(dom,sps)%il .and. i>=2 .and. j>=2) then
            ! Set the dot product of the grid velocity and the
                    ! normal in k-direction for a moving face.

            sFace = cudaDoms(dom,sps)%sFaceK(i, j, k)

            ! Compute the normal velocities of the left and right state.

            vnp = cudaDoms(dom,sps)%w(i, j, k + 1, ivx) * cudaDoms(dom,sps)%sK(i, j, k, 1) &
                  + cudaDoms(dom,sps)%w(i, j, k + 1, ivy) * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                  + cudaDoms(dom,sps)%w(i, j, k + 1, ivz) * cudaDoms(dom,sps)%sK(i, j, k, 3)
            vnm = cudaDoms(dom,sps)%w(i, j, k, ivx) * cudaDoms(dom,sps)%sK(i, j, k, 1) &
                  + cudaDoms(dom,sps)%w(i, j, k, ivy) * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                  + cudaDoms(dom,sps)%w(i, j, k, ivz) * cudaDoms(dom,sps)%sK(i, j, k, 3)

            ! Set the values of the porosities for this face.
            ! porVel defines the porosity w.r.t. velocity;
            ! porFlux defines the porosity w.r.t. the entire flux.
            ! The latter is only zero for a discontinuous block
            ! block boundary that must be treated conservatively.
            ! The default value of porFlux is 0.5, such that the
            ! correct central flux is scattered to both cells.
            ! In case of a boundFlux the normal velocity is set
            ! to sFace.

            porVel = one
            porFlux = half

            if (cudaDoms(dom,sps)%porK(i, j, k) == noFlux) porFlux = zero
            if (cudaDoms(dom,sps)%porK(i, j, k) == boundFlux) then
                porVel = zero
                vnp = sFace
                vnm = sFace
            end if

            ! Incorporate porFlux in porVel.

            porVel = porVel * porFlux

            ! Compute the normal velocities for the face as well as the
            ! mass fluxes.

            qsp = (vnp - sFace) * porVel
            qsm = (vnm - sFace) * porVel

            rqsp = qsp * cudaDoms(dom,sps)%w(i, j, k + 1, irho)
            rqsm = qsm * cudaDoms(dom,sps)%w(i, j, k, irho)

            ! Compute the sum of the pressure multiplied by porFlux.
            ! For the default value of porFlux, 0.5, this leads to
            ! the average pressure.

            pa = porFlux * (cudaDoms(dom,sps)%P(i, j, k + 1) + cudaDoms(dom,sps)%P(i, j, k))

            ! Compute the fluxes and scatter them to the cells
            ! i,j,k and i,j,k+1. Store the density flux in the
            ! mass flow of the appropriate sliding mesh interface.

            fs = rqsp + rqsm
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j, k + 1, irho), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, irho), fs)
            ! cudaDoms(dom,sps)%dw(i, j, k + 1, irho) = cudaDoms(dom,sps)%dw(i, j, k + 1, irho) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, irho) = cudaDoms(dom,sps)%dw(i, j, k, irho) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i, j, k + 1, ivx) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivx) &
                 + pa * cudaDoms(dom,sps)%sK(i, j, k, 1)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j, k + 1, imx), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imx), fs)
            ! cudaDoms(dom,sps)%dw(i, j, k + 1, imx) = cudaDoms(dom,sps)%dw(i, j, k + 1, imx) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imx) = cudaDoms(dom,sps)%dw(i, j, k, imx) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i, j, k + 1, ivy) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivy) &
                 + pa * cudaDoms(dom,sps)%sK(i, j, k, 2)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j, k + 1, imy), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imy), fs)
            ! cudaDoms(dom,sps)%dw(i, j, k + 1, imy) = cudaDoms(dom,sps)%dw(i, j, k + 1, imy) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imy) = cudaDoms(dom,sps)%dw(i, j, k, imy) + fs

            fs = rqsp * cudaDoms(dom,sps)%w(i, j, k + 1, ivz) + rqsm * cudaDoms(dom,sps)%w(i, j, k, ivz) &
                 + pa * cudaDoms(dom,sps)%sK(i, j, k, 3)
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j, k + 1, imz), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, imz), fs)
            ! cudaDoms(dom,sps)%dw(i, j, k + 1, imz) = cudaDoms(dom,sps)%dw(i, j, k + 1, imz) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, imz) = cudaDoms(dom,sps)%dw(i, j, k, imz) + fs

            fs = qsp * cudaDoms(dom,sps)%w(i, j, k + 1, irhoE) + qsm * cudaDoms(dom,sps)%w(i, j, k, irhoE) &
                 + porFlux * (vnp * cudaDoms(dom,sps)%P(i, j, k + 1) + vnm * cudaDoms(dom,sps)%P(i, j, k))
            tmp = atomicsub(cudaDoms(dom,sps)%dw(i, j, k + 1, irhoE), fs)
            tmp = atomicadd(cudaDoms(dom,sps)%dw(i, j, k, irhoE), fs)
            ! cudaDoms(dom,sps)%dw(i, j, k + 1, irhoE) = cudaDoms(dom,sps)%dw(i, j, k + 1, irhoE) - fs
            ! cudaDoms(dom,sps)%dw(i, j, k, irhoE) = cudaDoms(dom,sps)%dw(i, j, k, irhoE) + fs
        end if 
        
        !we left out the rotationof the block


    end subroutine inviscidCentralFlux

    attributes(global) subroutine computeSS
        use constants,           only: zero, two, half,EulerEquations
        use precision,           only: intType, realType
        use cudaInputDiscretization, only: vis2, vis4
        use cudaInputIteration,      only: useDissContinuation, dissContMagnitude, dissContMidpoint, dissContSharpness
        use cudaInputPhysics,        only: equations
        use cudaIteration,           only: rFil, totalR0, totalR
        use cudaFlowVarRefState,     only: gammaInf, pInfCorr, rhoInf, nwf

        implicit none
    
        ! Variables for inviscid diss flux scalar
        real(kind=realType), parameter :: dssMax = 0.25_realType
        real(kind=realType)            :: sslim, rhoi
        real(kind=realType)            :: sfil, fis2, fis4
        real(kind=realType)            :: ppor, rrad, dis2, dis4, fs
        real(kind=realType)            :: ddw1, ddw2, ddw3, ddw4, ddw5
        real(kind=realType)            :: tmp
        integer(kind=intType)          :: l
        integer(kind=intType) :: i, j, k , dom, sps 

        dom = 1
        sps = 1
        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1  ! starting at 0
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1  ! for entropy part
        k = (blockIdx%z - 1) * blockDim%z + threadIdx%z - 1


        ! Determine the variables used to compute the switch.
      ! For the inviscid case this is the pressure; for the viscous
      ! case it is the entropy.
    
      if(equations == EulerEquations) then
    
        ! Inviscid case. Pressure switch is based on the pressure.
        ! Also set the value of sslim. To be fully consistent this
        ! must have the dimension of pressure and it is therefore
        ! set to a fraction of the free stream value.
    
        sslim = 0.001_realType * pInfCorr
    
        ! Copy the pressure in ss. Only need the entries used in the
        ! discretization, i.e. not including the corner halo's, but we'll
        ! just copy all anyway.
        
        !TODO this was a copy as ss=P 
        !on gpu each thread needs to copy its own value
        !need to check we copied enough values this was temporary
        cudaDoms(dom,sps)%ss(i,j,k) = cudaDoms(dom,sps)%P(i,j,k)
        !===============================================================
    
      else if (equations == NSEquations .or. equations == RANSEquations) then
    
        ! Viscous case. Pressure switch is based on the entropy.
        ! Also set the value of sslim. To be fully consistent this
        ! must have the dimension of entropy and it is therefore
        ! set to a fraction of the free stream value.
    
        sslim = 0.001_realType * pInfCorr / (rhoInf**gammaInf)
    
        ! Store the entropy in ss. See above.
        !loop k 0 to cudaDoms(dom,sps)%kb, j 0 to cudaDoms(dom,sps)%jb, i 0 to cudaDoms(dom,sps)%ib
        if (((i >=               0) .AND. (i <= cudaDoms(dom,sps)%ib)) .AND. &
            ((j >=               0) .AND. (j <= cudaDoms(dom,sps)%jb)) .AND. &
            ((k >=               0) .AND. (k <= cudaDoms(dom,sps)%kb))) then 
    
          cudaDoms(dom,sps)%ss(i, j, k) = cudaDoms(dom,sps)%P(i, j, k) &
                                            / (cudaDoms(dom,sps)%w(i, j, k, irho)**cudaDoms(dom,sps)%gamma(i, j, k))
    
        end if 
      end if
        
    end subroutine computeSS


    attributes(global) subroutine inviscidDissFluxScalar
    
      use constants,           only: zero, two, half,EulerEquations
      use precision,           only: intType, realType
      use cudaInputDiscretization, only: vis2, vis4
      use cudaInputIteration,      only: useDissContinuation, dissContMagnitude, dissContMidpoint, dissContSharpness
      use cudaInputPhysics,        only: equations
      use cudaIteration,           only: rFil, totalR0, totalR
      use cudaFlowVarRefState,     only: gammaInf, pInfCorr, rhoInf, nwf
    
      implicit none
    
      ! Variables for inviscid diss flux scalar
      real(kind=realType), parameter :: dssMax = 0.25_realType
      real(kind=realType)            :: sslim, rhoi
      real(kind=realType)            :: sfil, fis2, fis4
      real(kind=realType)            :: ppor, rrad, dis2, dis4, fs
      real(kind=realType)            :: ddw1, ddw2, ddw3, ddw4, ddw5
      real(kind=realType)            :: tmp
      integer(kind=intType)          :: l
      integer(kind=intType) :: i, j, k, dom, sps 
      dom = 1
      sps = 1
      i = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1  ! starting at 0
      j = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1  ! for entropy part
      k = (blockIdx%z - 1) * blockDim%z + threadIdx%z - 1
    
    
      
        if(equations == EulerEquations) then
            sslim = 0.001_realType * pInfCorr
        else if (equations == NSEquations .or. equations == RANSEquations) then
            sslim = 0.001_realType * pInfCorr / (rhoInf**gammaInf)
        end if
    
      ! Compute the pressure sensor for each cell, in each direction:
      !loop k 1 to cudaDoms(dom,sps)%ke, j 1 to cudaDoms(dom,sps)%je, i 1 to cudaDoms(dom,sps)%ie
      if (((i >=               1) .AND. (i <= cudaDoms(dom,sps)%ie)) .AND. &
          ((j >=               1) .AND. (j <= cudaDoms(dom,sps)%je)) .AND. &
          ((k >=               1) .AND. (k <= cudaDoms(dom,sps)%ke))) then 
    
        cudaDoms(dom,sps)%dss(i, j, k, 1) = abs((cudaDoms(dom,sps)%ss(i + 1, j, k) - two * cudaDoms(dom,sps)%ss(i, j, k) &
                            + cudaDoms(dom,sps)%ss(i - 1, j, k)) &
            / (cudaDoms(dom,sps)%ss(i + 1, j, k) + two * cudaDoms(dom,sps)%ss(i, j, k) + cudaDoms(dom,sps)%ss(i - 1, j, k) + sslim))
        cudaDoms(dom,sps)%dss(i, j, k, 2) = abs((cudaDoms(dom,sps)%ss(i, j + 1, k) - two * cudaDoms(dom,sps)%ss(i, j, k) &
        + cudaDoms(dom,sps)%ss(i, j - 1, k)) &
            / (cudaDoms(dom,sps)%ss(i, j + 1, k) + two * cudaDoms(dom,sps)%ss(i, j, k) + cudaDoms(dom,sps)%ss(i, j - 1, k) + sslim))
        cudaDoms(dom,sps)%dss(i, j, k, 3) = abs((cudaDoms(dom,sps)%ss(i, j, k + 1) - two * cudaDoms(dom,sps)%ss(i, j, k) &
        + cudaDoms(dom,sps)%ss(i, j, k - 1)) &
            / (cudaDoms(dom,sps)%ss(i, j, k + 1) + two * cudaDoms(dom,sps)%ss(i, j, k) + cudaDoms(dom,sps)%ss(i, j, k - 1) + sslim))
    
      end if
    
    
      ! Set the dissipation constants for the scheme.
      ! rFil and sFil are fractions used by the Runge-Kutta solver to compute residuals at intermediate steps.
      ! For the blockette code, rFil is always one, so sFil==0, fis2==vis2, and fis4==vis4.
    
      ! The sigmoid function used for dissipation-based continuation is described in Eq. 28 and Eq. 29 from the paper:
      ! "Improving the Performance of a Compressible RANS Solver for Low and High Mach Number Flows" (Seraj2022c).
      ! The options documentation also has information on the parameters in this formulation.
      if (useDissContinuation) then
          if (totalR == zero .or. totalR0 == zero) then
              fis2 = rFil * (vis2 + dissContMagnitude / (1 + exp(-dissContSharpness * dissContMidpoint)))
          else
              fis2 = rFil * (vis2 + dissContMagnitude / &
                             (1 + exp(-dissContSharpness * (log10(totalR / totalR0) + dissContMidpoint))))
          end if
      else
          fis2 = rFil * vis2
      end if
      fis4 = rFil * vis4
      sfil = one - rFil
    
      ! Initialize the dissipative residual to a certain times,
      ! possibly zero, the previously stored value. Owned cells
      ! only, because the halo values do not matter.
     
     if (((i >= 1) .AND. (i <= cudaDoms(dom,sps)%ie)) .AND. &
          ((j >= 1) .AND. (j <= cudaDoms(dom,sps)%je)) .AND. &
          ((k >= 1) .AND. (k <= cudaDoms(dom,sps)%ke))) then 
        do l = 1, nwf
            cudaDoms(dom,sps)%fw(i,j,k,l)  = sfil * cudaDoms(dom,sps)%fw(i,j,k,l)
        end do
     end if
      !
      !       Dissipative fluxes in the i-direction.
      !
      !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 1 to cudaDoms(dom,sps)%il
      if (((i >= 1) .AND. (i <= cudaDoms(dom,sps)%il)) .AND. &
          ((j >= 2) .AND. (j <= cudaDoms(dom,sps)%jl)) .AND. &
          ((k >= 2) .AND. (k <= cudaDoms(dom,sps)%kl))) then 
    
        ! Compute the dissipation coefficients for this face.
    
        ppor = zero
        if (cudaDoms(dom,sps)%porI(i, j, k) == normalFlux) ppor = half
        rrad = ppor * (cudaDoms(dom,sps)%radi(i, j, k) + cudaDoms(dom,sps)%radi(i + 1, j, k))
    
        dis2 = fis2 * rrad * min(dssMax, max(cudaDoms(dom,sps)%dss(i, j, k, 1), cudaDoms(dom,sps)%dss(i + 1, j, k, 1)))
        dis4 = dim(fis4 * rrad, dis2)
    
        ! Compute and scatter the dissipative flux.
        ! Density. Store it in the mass flow of the
        ! appropriate sliding mesh interface.
    
        ddw1 = cudaDoms(dom,sps)%w(i + 1, j, k, irho) - cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw1 &
             - dis4 * (cudaDoms(dom,sps)%w(i + 2, j, k, irho) - cudaDoms(dom,sps)%w(i - 1, j, k, irho) - three * ddw1)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, irho), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(    i, j, k, irho), fs)

        ! cudaDoms(dom,sps)%fw(i + 1, j, k, irho) = cudaDoms(dom,sps)%fw(i + 1, j, k, irho) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, irho) = cudaDoms(dom,sps)%fw(i, j, k, irho) - fs
    
        ! X-momentum.
    
        ddw2 = cudaDoms(dom,sps)%w(i + 1, j, k, ivx) * cudaDoms(dom,sps)%w(i + 1, j, k, irho) &
        - cudaDoms(dom,sps)%w(i, j, k, ivx) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw2 &
             - dis4 * (cudaDoms(dom,sps)%w(i + 2, j, k, ivx) * cudaDoms(dom,sps)%w(i + 2, j, k, irho) - &
                       cudaDoms(dom,sps)%w(i - 1, j, k, ivx) * cudaDoms(dom,sps)%w(i - 1, j, k, irho) - three * ddw2)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, imx), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(    i, j, k, imx), fs)

        ! cudaDoms(dom,sps)%fw(i + 1, j, k, imx) = cudaDoms(dom,sps)%fw(i + 1, j, k, imx) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imx) = cudaDoms(dom,sps)%fw(i, j, k, imx) - fs
    
        ! Y-momentum.
    
        ddw3 = cudaDoms(dom,sps)%w(i + 1, j, k, ivy) * cudaDoms(dom,sps)%w(i + 1, j, k, irho) &
        - cudaDoms(dom,sps)%w(i, j, k, ivy) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw3 &
             - dis4 * (cudaDoms(dom,sps)%w(i + 2, j, k, ivy) * cudaDoms(dom,sps)%w(i + 2, j, k, irho) - &
                       cudaDoms(dom,sps)%w(i - 1, j, k, ivy) * cudaDoms(dom,sps)%w(i - 1, j, k, irho) - three * ddw3)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, imy), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(    i, j, k, imy), fs)

        ! cudaDoms(dom,sps)%fw(i + 1, j, k, imy) = cudaDoms(dom,sps)%fw(i + 1, j, k, imy) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imy) = cudaDoms(dom,sps)%fw(i, j, k, imy) - fs
    
        ! Z-momentum.
    
        ddw4 = cudaDoms(dom,sps)%w(i + 1, j, k, ivz) * cudaDoms(dom,sps)%w(i + 1, j, k, irho) &
        - cudaDoms(dom,sps)%w(i, j, k, ivz) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw4 &
             - dis4 * (cudaDoms(dom,sps)%w(i + 2, j, k, ivz) * cudaDoms(dom,sps)%w(i + 2, j, k, irho) - &
                       cudaDoms(dom,sps)%w(i - 1, j, k, ivz) * cudaDoms(dom,sps)%w(i - 1, j, k, irho) - three * ddw4)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, imz), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(    i, j, k, imz), fs)
    
        ! cudaDoms(dom,sps)%fw(i + 1, j, k, imz) = cudaDoms(dom,sps)%fw(i + 1, j, k, imz) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imz) = cudaDoms(dom,sps)%fw(i, j, k, imz) - fs

        ! Energy.
    
        ddw5 = (cudaDoms(dom,sps)%w(i + 1, j, k, irhoE) + cudaDoms(dom,sps)%P(i + 1, j, K)) &
                                - (cudaDoms(dom,sps)%w(i, j, k, irhoE) + cudaDoms(dom,sps)%P(i, j, k))
        fs = dis2 * ddw5 &
             - dis4 * ((cudaDoms(dom,sps)%w(i + 2, j, k, irhoE) + cudaDoms(dom,sps)%P(i + 2, j, k)) - &
                       (cudaDoms(dom,sps)%w(i - 1, j, k, irhoE) + cudaDoms(dom,sps)%P(i - 1, j, k)) - three * ddw5)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, irhoE), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(    i, j, k, irhoE), fs)

        ! cudaDoms(dom,sps)%fw(i + 1, j, k, irhoE) = cudaDoms(dom,sps)%fw(i + 1, j, k, irhoE) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, irhoE) = cudaDoms(dom,sps)%fw(i, j, k, irhoE) - fs
    
      end if
    
      !
      !       Dissipative fluxes in the j-direction.
      !
      !loop k 2 to cudaDoms(dom,sps)%kl, j 1 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
      if (((i >= 2) .AND. (i <= cudaDoms(dom,sps)%il)) .AND. &
          ((j >= 1) .AND. (j <= cudaDoms(dom,sps)%jl)) .AND. &
          ((k >= 2) .AND. (k <= cudaDoms(dom,sps)%kl))) then 
    
        ! Compute the dissipation coefficients for this face.
    
        ppor = zero
        if (cudaDoms(dom,sps)%porJ(i, j, k) == normalFlux) ppor = half
        rrad = ppor * (cudaDoms(dom,sps)%radj(i, j, k) + cudaDoms(dom,sps)%radj(i, j + 1, k))
    
        dis2 = fis2 * rrad * min(dssMax, max(cudaDoms(dom,sps)%dss(i, j, k, 2), cudaDoms(dom,sps)%dss(i, j + 1, k, 2)))
        dis4 = dim(fis4 * rrad, dis2)
    
        ! Compute and scatter the dissipative flux.
        ! Density. Store it in the mass flow of the
        ! appropriate sliding mesh interface.
    
        ddw1 = cudaDoms(dom,sps)%w(i, j + 1, k, irho) - cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw1 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j + 2, k, irho) - cudaDoms(dom,sps)%w(i, j - 1, k, irho) - three * ddw1)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, irho), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i,     j, k, irho), fs)

        ! cudaDoms(dom,sps)%fw(i, j + 1, k, irho) = cudaDoms(dom,sps)%fw(i, j + 1, k, irho) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, irho) = cudaDoms(dom,sps)%fw(i, j, k, irho) - fs
    
        ! X-momentum.
    
        ddw2 = cudaDoms(dom,sps)%w(i, j + 1, k, ivx) * cudaDoms(dom,sps)%w(i, j + 1, k, irho) &
                                            - cudaDoms(dom,sps)%w(i, j, k, ivx) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw2 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j + 2, k, ivx) * cudaDoms(dom,sps)%w(i, j + 2, k, irho) - &
                       cudaDoms(dom,sps)%w(i, j - 1, k, ivx) * cudaDoms(dom,sps)%w(i, j - 1, k, irho) - three * ddw2)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, imx), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i,     j, k, imx), fs)

        ! cudaDoms(dom,sps)%fw(i, j + 1, k, imx) = cudaDoms(dom,sps)%fw(i, j + 1, k, imx) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imx) = cudaDoms(dom,sps)%fw(i, j, k, imx) - fs
    
        ! Y-momentum.
    
        ddw3 = cudaDoms(dom,sps)%w(i, j + 1, k, ivy) * cudaDoms(dom,sps)%w(i, j + 1, k, irho) &
                                            - cudaDoms(dom,sps)%w(i, j, k, ivy) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw3 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j + 2, k, ivy) * cudaDoms(dom,sps)%w(i, j + 2, k, irho) - &
                       cudaDoms(dom,sps)%w(i, j - 1, k, ivy) * cudaDoms(dom,sps)%w(i, j - 1, k, irho) - three * ddw3)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, imy), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i,     j, k, imy), fs)

        ! cudaDoms(dom,sps)%fw(i, j + 1, k, imy) = cudaDoms(dom,sps)%fw(i, j + 1, k, imy) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imy) = cudaDoms(dom,sps)%fw(i, j, k, imy) - fs
    
        ! Z-momentum.
    
        ddw4 = cudaDoms(dom,sps)%w(i, j + 1, k, ivz) * cudaDoms(dom,sps)%w(i, j + 1, k, irho) &
                                                    - cudaDoms(dom,sps)%w(i, j, k, ivz) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw4 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j + 2, k, ivz) * cudaDoms(dom,sps)%w(i, j + 2, k, irho) - &
                       cudaDoms(dom,sps)%w(i, j - 1, k, ivz) * cudaDoms(dom,sps)%w(i, j - 1, k, irho) - three * ddw4)

        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, imz), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i,     j, k, imz), fs)

        ! cudaDoms(dom,sps)%fw(i, j + 1, k, imz) = cudaDoms(dom,sps)%fw(i, j + 1, k, imz) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imz) = cudaDoms(dom,sps)%fw(i, j, k, imz) - fs
    
        ! Energy.
    
        ddw5 = (cudaDoms(dom,sps)%w(i, j + 1, k, irhoE) + cudaDoms(dom,sps)%P(i, j + 1, k)) &
                                                    - (cudaDoms(dom,sps)%w(i, j, k, irhoE) + cudaDoms(dom,sps)%P(i, j, k))
        fs = dis2 * ddw5 &
             - dis4 * ((cudaDoms(dom,sps)%w(i, j + 2, k, irhoE) + cudaDoms(dom,sps)%P(i, j + 2, k)) - &
                       (cudaDoms(dom,sps)%w(i, j - 1, k, irhoE) + cudaDoms(dom,sps)%P(i, j - 1, k)) - three * ddw5)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, irhoE), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i,     j, k, irhoE), fs)

        ! cudaDoms(dom,sps)%fw(i, j + 1, k, irhoE) = cudaDoms(dom,sps)%fw(i, j + 1, k, irhoE) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, irhoE) = cudaDoms(dom,sps)%fw(i, j, k, irhoE) - fs
    
      end if
      !
      !       Dissipative fluxes in the k-direction.
      !
      !loop k 1 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
      if (((i >= 2) .AND. (i <= cudaDoms(dom,sps)%il)) .AND. &
          ((j >= 2) .AND. (j <= cudaDoms(dom,sps)%jl)) .AND. &
          ((k >= 1) .AND. (k <= cudaDoms(dom,sps)%kl))) then 
    
        ! Compute the dissipation coefficients for this face.
    
        ppor = zero
        if (cudaDoms(dom,sps)%porK(i, j, k) == normalFlux) ppor = half
        rrad = ppor * (cudaDoms(dom,sps)%radk(i, j, k) + cudaDoms(dom,sps)%radk(i, j, k + 1))
    
        dis2 = fis2 * rrad * min(dssMax, max(cudaDoms(dom,sps)%dss(i, j, k, 3), cudaDoms(dom,sps)%dss(i, j, k + 1, 3)))
        dis4 = dim(fis4 * rrad, dis2)
    
        ! Compute and scatter the dissipative flux.
        ! Density. Store it in the mass flow of the
        ! appropriate sliding mesh interface.
    
        ddw1 = cudaDoms(dom,sps)%w(i, j, k + 1, irho) - cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw1 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j, k + 2, irho) - cudaDoms(dom,sps)%w(i, j, k - 1, irho) - three * ddw1)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, irho), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j,     k, irho), fs)

        ! cudaDoms(dom,sps)%fw(i, j, k + 1, irho) = cudaDoms(dom,sps)%fw(i, j, k + 1, irho) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, irho) = cudaDoms(dom,sps)%fw(i, j, k, irho) - fs
    
        ! X-momentum.
    
        ddw2 = cudaDoms(dom,sps)%w(i, j, k + 1, ivx) * cudaDoms(dom,sps)%w(i, j, k + 1, irho) &
                                                    - cudaDoms(dom,sps)%w(i, j, k, ivx) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw2 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j, k + 2, ivx) * cudaDoms(dom,sps)%w(i, j, k + 2, irho) - &
                       cudaDoms(dom,sps)%w(i, j, k - 1, ivx) * cudaDoms(dom,sps)%w(i, j, k - 1, irho) - three * ddw2)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, imx), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j,     k, imx), fs)

        ! cudaDoms(dom,sps)%fw(i, j, k + 1, imx) = cudaDoms(dom,sps)%fw(i, j, k + 1, imx) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imx) = cudaDoms(dom,sps)%fw(i, j, k, imx) - fs
    
        ! Y-momentum.
    
        ddw3 = cudaDoms(dom,sps)%w(i, j, k + 1, ivy) * cudaDoms(dom,sps)%w(i, j, k + 1, irho) &
                                                    - cudaDoms(dom,sps)%w(i, j, k, ivy) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw3 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j, k + 2, ivy) * cudaDoms(dom,sps)%w(i, j, k + 2, irho) - &
                       cudaDoms(dom,sps)%w(i, j, k - 1, ivy) * cudaDoms(dom,sps)%w(i, j, k - 1, irho) - three * ddw3)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, imy), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j,     k, imy), fs)

        ! cudaDoms(dom,sps)%fw(i, j, k + 1, imy) = cudaDoms(dom,sps)%fw(i, j, k + 1, imy) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imy) = cudaDoms(dom,sps)%fw(i, j, k, imy) - fs
    
        ! Z-momentum.
    
        ddw4 = cudaDoms(dom,sps)%w(i, j, k + 1, ivz) * cudaDoms(dom,sps)%w(i, j, k + 1, irho) &
                                                    - cudaDoms(dom,sps)%w(i, j, k, ivz) * cudaDoms(dom,sps)%w(i, j, k, irho)
        fs = dis2 * ddw4 &
             - dis4 * (cudaDoms(dom,sps)%w(i, j, k + 2, ivz) * cudaDoms(dom,sps)%w(i, j, k + 2, irho) - &
                       cudaDoms(dom,sps)%w(i, j, k - 1, ivz) * cudaDoms(dom,sps)%w(i, j, k - 1, irho) - three * ddw4)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, imz), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j,     k, imz), fs)

        ! cudaDoms(dom,sps)%fw(i, j, k + 1, imz) = cudaDoms(dom,sps)%fw(i, j, k + 1, imz) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, imz) = cudaDoms(dom,sps)%fw(i, j, k, imz) - fs
    
        ! Energy.
    
        ddw5 = (cudaDoms(dom,sps)%w(i, j, k + 1, irhoE) + cudaDoms(dom,sps)%P(i, j, k + 1)) &
                                                    - (cudaDoms(dom,sps)%w(i, j, k, irhoE) + cudaDoms(dom,sps)%P(i, j, k))
        fs = dis2 * ddw5 &
             - dis4 * ((cudaDoms(dom,sps)%w(i, j, k + 2, irhoE) + cudaDoms(dom,sps)%P(i, j, k + 2)) - &
                       (cudaDoms(dom,sps)%w(i, j, k - 1, irhoE) + cudaDoms(dom,sps)%P(i, j, k - 1)) - three * ddw5)
    
        tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, irhoE), fs)
        tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j,     k, irhoE), fs)

        ! cudaDoms(dom,sps)%fw(i, j, k + 1, irhoE) = cudaDoms(dom,sps)%fw(i, j, k + 1, irhoE) + fs
        ! cudaDoms(dom,sps)%fw(i, j, k, irhoE) = cudaDoms(dom,sps)%fw(i, j, k, irhoE) - fs
    
      end if
    
    end subroutine inviscidDissFluxScalar

    ! !alex
    attributes(global) subroutine viscousFlux
        use precision, only: realType, intType
        use constants, only: half, zero,one,two, third,fourth,eighth,ivx,ivy,ivz,irhoE,irho,itu1,imx,imy,imz,noFlux
        use cudaInputPhysics, only: useQCR, prandtl, prandtlturb
        use cudaFlowvarRefState, only: eddyModel
        use cudaIteration, only: rFil
        implicit none

        ! Variables for viscous flux
        real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
        real(kind=realType) :: gm1, factLamHeat, factTurbHeat
        real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
        real(kind=realType) :: q_x, q_y, q_z
        real(kind=realType) :: corr, ssx, ssy, ssz, fracDiv, snrm
        real(kind=realType) :: tauxx, tauyy, tauzz
        real(kind=realType) :: tauxy, tauxz, tauyz
        real(kind=realType) :: tauxxS, tauyyS, tauzzS
        real(kind=realType) :: tauxyS, tauxzS, tauyzS
        real(kind=realType) :: ubar, vbar, wbar
        real(kind=realType) :: exx, eyy, ezz
        real(kind=realType) :: exy, exz, eyz
        real(kind=realType) :: Wxx, Wyy, Wzz
        real(kind=realType) :: Wxy, Wxz, Wyz, Wyx, Wzx, Wzy
        real(kind=realType) :: den, Ccr1
        real(kind=realType) :: fmx, fmy, fmz, frhoE, fact
        integer(kind=intType) :: i, j, k, io, jo, ko, dom, sps
        real(kind=realType), parameter :: xminn = 1.e-10_realType
        real(kind=realType), parameter :: twoThird = two * third
        real(kind=realType) :: tmp

        dom = 1
        sps = 1
        !thread indices start at 1
        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x 
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y 
        k = (blockIdx%z - 1) * blockDim%z + threadIdx%z 

        ! Set QCR parameters
        Ccr1 = 0.3_realType
        rFilv = rFil

        ! The diagonals of the vorticity tensor components are always zero
        Wxx = zero
        Wyy = zero
        Wzz = zero

        !
        !         viscous fluxes in the k-direction.
        !
        mue = zero
        !loop k 1 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (k <= cudaDoms(dom,sps)%kl .and. j <= cudaDoms(dom,sps)%jl &
                .and. i <= cudaDoms(dom,sps)%il .and. j>=2 .and. i>=2 .and. k>=1) then
            ! Set the value of the porosity. If not zero, it is set
            ! to average the eddy-viscosity and to take the factor
            ! rFilv into account.
            por = half * rFilv
            if (cudaDoms(dom,sps)%porK(i, j, k) == noFlux) por = zero

            ! Compute the laminar and (if present) the eddy viscosities
            ! multiplied by the porosity. Compute the factor in front of
            ! the gradients of the speed of sound squared for the heat
            ! flux.
            mul = por * (cudaDoms(dom,sps)%rlv(i, j, k) + cudaDoms(dom,sps)%rlv(i, j, k + 1))
            mue = por * (cudaDoms(dom,sps)%rev(i, j, k) + cudaDoms(dom,sps)%rev(i, j, k + 1))
            mut = mul + mue

            gm1 = half * (cudaDoms(dom,sps)%gamma(i, j, k) + cudaDoms(dom,sps)%gamma(i, j, k + 1)) - one
            factLamHeat = one / (prandtl * gm1)
            factTurbHeat = one / (prandtlTurb * gm1)

            heatCoef = mul * factLamHeat + mue * factTurbHeat
            ! Compute the gradients at the face by averaging the four
                    ! nodal values.

            u_x = fourth * (cudaDoms(dom,sps)%ux(i - 1, j - 1, k) + cudaDoms(dom,sps)%ux(i, j - 1, k) &
                        + cudaDoms(dom,sps)%ux(i - 1, j, k) + cudaDoms(dom,sps)%ux(i, j, k))
            u_y = fourth * (cudaDoms(dom,sps)%uy(i - 1, j - 1, k) + cudaDoms(dom,sps)%uy(i, j - 1, k) &
                        + cudaDoms(dom,sps)%uy(i - 1, j, k) + cudaDoms(dom,sps)%uy(i, j, k))
            u_z = fourth * (cudaDoms(dom,sps)%uz(i - 1, j - 1, k) + cudaDoms(dom,sps)%uz(i, j - 1, k) &
                        + cudaDoms(dom,sps)%uz(i - 1, j, k) + cudaDoms(dom,sps)%uz(i, j, k))

            v_x = fourth * (cudaDoms(dom,sps)%vx(i - 1, j - 1, k) + cudaDoms(dom,sps)%vx(i, j - 1, k) &
                        + cudaDoms(dom,sps)%vx(i - 1, j, k) + cudaDoms(dom,sps)%vx(i, j, k))
            v_y = fourth * (cudaDoms(dom,sps)%vy(i - 1, j - 1, k) + cudaDoms(dom,sps)%vy(i, j - 1, k) &
                        + cudaDoms(dom,sps)%vy(i - 1, j, k) + cudaDoms(dom,sps)%vy(i, j, k))
            v_z = fourth * (cudaDoms(dom,sps)%vz(i - 1, j - 1, k) + cudaDoms(dom,sps)%vz(i, j - 1, k) &
                        + cudaDoms(dom,sps)%vz(i - 1, j, k) + cudaDoms(dom,sps)%vz(i, j, k))

            w_x = fourth * (cudaDoms(dom,sps)%wx(i - 1, j - 1, k) + cudaDoms(dom,sps)%wx(i, j - 1, k) &
                        + cudaDoms(dom,sps)%wx(i - 1, j, k) + cudaDoms(dom,sps)%wx(i, j, k))
            w_y = fourth * (cudaDoms(dom,sps)%wy(i - 1, j - 1, k) + cudaDoms(dom,sps)%wy(i, j - 1, k) &
                        + cudaDoms(dom,sps)%wy(i - 1, j, k) + cudaDoms(dom,sps)%wy(i, j, k))
            w_z = fourth * (cudaDoms(dom,sps)%wz(i - 1, j - 1, k) + cudaDoms(dom,sps)%wz(i, j - 1, k) &
                        + cudaDoms(dom,sps)%wz(i - 1, j, k) + cudaDoms(dom,sps)%wz(i, j, k))

            q_x = fourth * (cudaDoms(dom,sps)%qx(i - 1, j - 1, k) + cudaDoms(dom,sps)%qx(i, j - 1, k) &
                        + cudaDoms(dom,sps)%qx(i - 1, j, k) + cudaDoms(dom,sps)%qx(i, j, k))
            q_y = fourth * (cudaDoms(dom,sps)%qy(i - 1, j - 1, k) + cudaDoms(dom,sps)%qy(i, j - 1, k) &
                        + cudaDoms(dom,sps)%qy(i - 1, j, k) + cudaDoms(dom,sps)%qy(i, j, k))
            q_z = fourth * (cudaDoms(dom,sps)%qz(i - 1, j - 1, k) + cudaDoms(dom,sps)%qz(i, j - 1, k) &
                        + cudaDoms(dom,sps)%qz(i - 1, j, k) + cudaDoms(dom,sps)%qz(i, j, k))
            ! The gradients in the normal direction are corrected, such
            ! that no averaging takes places here.
            ! First determine the vector in the direction from the
            ! cell center k to cell center k+1.

            ssx = eighth * (cudaDoms(dom,sps)%x(i - 1, j - 1, k + 1, 1) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 1) &
                        + cudaDoms(dom,sps)%x(i - 1, j, k + 1, 1) - cudaDoms(dom,sps)%x(i - 1, j, k - 1, 1) &
                        + cudaDoms(dom,sps)%x(i, j - 1, k + 1, 1) - cudaDoms(dom,sps)%x(i, j - 1, k - 1, 1) &
                        + cudaDoms(dom,sps)%x(i, j, k + 1, 1) - cudaDoms(dom,sps)%x(i, j, k - 1, 1))
            ssy = eighth * (cudaDoms(dom,sps)%x(i - 1, j - 1, k + 1, 2) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 2) &
                        + cudaDoms(dom,sps)%x(i - 1, j, k + 1, 2) - cudaDoms(dom,sps)%x(i - 1, j, k - 1, 2) &
                        + cudaDoms(dom,sps)%x(i, j - 1, k + 1, 2) - cudaDoms(dom,sps)%x(i, j - 1, k - 1, 2) &
                        + cudaDoms(dom,sps)%x(i, j, k + 1, 2) - cudaDoms(dom,sps)%x(i, j, k - 1, 2))
            ssz = eighth * (cudaDoms(dom,sps)%x(i - 1, j - 1, k + 1, 3) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 3) &
                        + cudaDoms(dom,sps)%x(i - 1, j, k + 1, 3) - cudaDoms(dom,sps)%x(i - 1, j, k - 1, 3) &
                        + cudaDoms(dom,sps)%x(i, j - 1, k + 1, 3) - cudaDoms(dom,sps)%x(i, j - 1, k - 1, 3) &
                        + cudaDoms(dom,sps)%x(i, j, k + 1, 3) - cudaDoms(dom,sps)%x(i, j, k - 1, 3))

            ! Determine the length of this vector and create the
            ! unit normal.

            snrm = one / sqrt(ssx * ssx + ssy * ssy + ssz * ssz)
            ssx = snrm * ssx
            ssy = snrm * ssy
            ssz = snrm * ssz
    
            ! Correct the gradients.

            corr = u_x * ssx + u_y * ssy + u_z * ssz &
                    - (cudaDoms(dom,sps)%w(i, j, k + 1, ivx) - cudaDoms(dom,sps)%w(i, j, k, ivx)) * snrm
            u_x = u_x - corr * ssx
            u_y = u_y - corr * ssy
            u_z = u_z - corr * ssz

            corr = v_x * ssx + v_y * ssy + v_z * ssz &
                    - (cudaDoms(dom,sps)%w(i, j, k + 1, ivy) - cudaDoms(dom,sps)%w(i, j, k, ivy)) * snrm
            v_x = v_x - corr * ssx
            v_y = v_y - corr * ssy
            v_z = v_z - corr * ssz

            corr = w_x * ssx + w_y * ssy + w_z * ssz &
                    - (cudaDoms(dom,sps)%w(i, j, k + 1, ivz) - cudaDoms(dom,sps)%w(i, j, k, ivz)) * snrm
            w_x = w_x - corr * ssx
            w_y = w_y - corr * ssy
            w_z = w_z - corr * ssz

            corr = q_x * ssx + q_y * ssy + q_z * ssz &
                    + (cudaDoms(dom,sps)%aa(i, j, k + 1) - cudaDoms(dom,sps)%aa(i, j, k)) * snrm
            q_x = q_x - corr * ssx
            q_y = q_y - corr * ssy
            q_z = q_z - corr * ssz
            
            ! Compute the stress tensor and the heat flux vector.
            ! We remove the viscosity from the stress tensor (tau)
            ! to define tauS since we still need to separate between
            ! laminar and turbulent stress for QCR.
            ! Therefore, laminar tau = mue*tauS, turbulent
            ! tau = mue*tauS, and total tau = mut*tauS.

            fracDiv = twoThird * (u_x + v_y + w_z)
            tauxxS = two * u_x - fracDiv
            tauyyS = two * v_y - fracDiv
            tauzzS = two * w_z - fracDiv

            tauxyS = u_y + v_x
            tauxzS = u_z + w_x
            tauyzS = v_z + w_y

            q_x = heatCoef * q_x
            q_y = heatCoef * q_y
            q_z = heatCoef * q_z
            
            ! Add QCR corrections if necessary
            if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x * u_x + u_y * u_y + u_z * u_z + &
                           v_x * v_x + v_y * v_y + v_z * v_z + &
                           w_x * w_x + w_y * w_y + w_z * w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue * Ccr1 / den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact * (Wxy * tauxyS + Wxz * tauxzS) * two
                eyy = fact * (Wyx * tauxyS + Wyz * tauyzS) * two
                ezz = fact * (Wzx * tauxzS + Wzy * tauyzS) * two

                exy = fact * (Wxy * tauyyS + Wxz * tauyzS + &
                              Wyx * tauxxS + Wyz * tauxzS)
                exz = fact * (Wxy * tauyzS + Wxz * tauzzS + &
                              Wzx * tauxxS + Wzy * tauxyS)
                eyz = fact * (Wyx * tauxzS + Wyz * tauzzS + &
                              Wzx * tauxyS + Wzy * tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut * tauxxS - exx
                tauyy = mut * tauyyS - eyy
                tauzz = mut * tauzzS - ezz
                tauxy = mut * tauxyS - exy
                tauxz = mut * tauxzS - exz
                tauyz = mut * tauyzS - eyz

            else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut * tauxxS
                tauyy = mut * tauyyS
                tauzz = mut * tauzzS
                tauxy = mut * tauxyS
                tauxz = mut * tauxzS
                tauyz = mut * tauyzS

            end if

            ! Compute the average velocities for the face. Remember that
            ! the velocities are stored and not the momentum.
            ubar = half * (cudaDoms(dom,sps)%w(i, j, k, ivx) + cudaDoms(dom,sps)%w(i, j, k + 1, ivx))
            vbar = half * (cudaDoms(dom,sps)%w(i, j, k, ivy) + cudaDoms(dom,sps)%w(i, j, k + 1, ivy))
            wbar = half * (cudaDoms(dom,sps)%w(i, j, k, ivz) + cudaDoms(dom,sps)%w(i, j, k + 1, ivz))

            ! Compute the viscous fluxes for this k-face.
            fmx = tauxx * cudaDoms(dom,sps)%sK(i, j, k, 1) + tauxy * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                    + tauxz * cudaDoms(dom,sps)%sK(i, j, k, 3)
            fmy = tauxy * cudaDoms(dom,sps)%sK(i, j, k, 1) + tauyy * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                    + tauyz * cudaDoms(dom,sps)%sK(i, j, k, 3)
            fmz = tauxz * cudaDoms(dom,sps)%sK(i, j, k, 1) + tauyz * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                    + tauzz * cudaDoms(dom,sps)%sK(i, j, k, 3)
            frhoE = (ubar * tauxx + vbar * tauxy + wbar * tauxz) * cudaDoms(dom,sps)%sK(i, j, k, 1)
            frhoE = frhoE + (ubar * tauxy + vbar * tauyy + wbar * tauyz) * cudaDoms(dom,sps)%sK(i, j, k, 2)
            frhoE = frhoE + (ubar * tauxz + vbar * tauyz + wbar * tauzz) * cudaDoms(dom,sps)%sK(i, j, k, 3)
            frhoE = frhoE - q_x * cudaDoms(dom,sps)%sK(i, j, k, 1) - q_y * cudaDoms(dom,sps)%sK(i, j, k, 2) - q_z * cudaDoms(dom,sps)%sK(i, j, k, 3)

            ! Update the residuals of cell k and k+1.
            
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imx), fmx)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imy), fmy)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imz), fmz)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, irhoE), frhoE)
            ! cudaDoms(dom,sps)%fw(i, j, k, imx) = cudaDoms(dom,sps)%fw(i, j, k, imx) - fmx
            ! cudaDoms(dom,sps)%fw(i, j, k, imy) = cudaDoms(dom,sps)%fw(i, j, k, imy) - fmy
            ! cudaDoms(dom,sps)%fw(i, j, k, imz) = cudaDoms(dom,sps)%fw(i, j, k, imz) - fmz
            ! cudaDoms(dom,sps)%fw(i, j, k, irhoE) = cudaDoms(dom,sps)%fw(i, j, k, irhoE) - frhoE

            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, imx), fmx)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, imy), fmy)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, imz), fmz)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j, k + 1, irhoE), frhoE)
            ! cudaDoms(dom,sps)%fw(i, j, k + 1, imx) = cudaDoms(dom,sps)%fw(i, j, k + 1, imx) + fmx
            ! cudaDoms(dom,sps)%fw(i, j, k + 1, imy) = cudaDoms(dom,sps)%fw(i, j, k + 1, imy) + fmy
            ! cudaDoms(dom,sps)%fw(i, j, k + 1, imz) = cudaDoms(dom,sps)%fw(i, j, k + 1, imz) + fmz
            ! cudaDoms(dom,sps)%fw(i, j, k + 1, irhoE) = cudaDoms(dom,sps)%fw(i, j, k + 1, irhoE) + frhoE

            ! Temporarily store the shear stress and heat flux, even
            ! if we won't need it. This can still vectorize
        end if

        !
        !         Viscous fluxes in the j-direction.
        !
        !loop k 2 to cudaDoms(dom,sps)%kl, j 1 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (k <= cudaDoms(dom,sps)%kl .and. j <= cudaDoms(dom,sps)%jl .and. i <= cudaDoms(dom,sps)%il .and. k>=2 .and. i>=2 .and. j>=1) then
            ! Set the value of the porosity. If not zero, it is set
                    ! to average the eddy-viscosity and to take the factor
                    ! rFilv into account.

            por = half * rFilv
            if (cudaDoms(dom,sps)%porJ(i, j, k) == noFlux) por = zero

            ! Compute the laminar and (if present) the eddy viscosities
            ! multiplied by the porosity. Compute the factor in front of
            ! the gradients of the speed of sound squared for the heat
            ! flux.

            mul = por * (cudaDoms(dom,sps)%rlv(i, j, k) + cudaDoms(dom,sps)%rlv(i, j + 1, k))
            mue = por * (cudaDoms(dom,sps)%rev(i, j, k) + cudaDoms(dom,sps)%rev(i, j + 1, k))
            mut = mul + mue

            gm1 = half * (cudaDoms(dom,sps)%gamma(i, j, k) + cudaDoms(dom,sps)%gamma(i, j + 1, k)) - one
            factLamHeat = one / (prandtl * gm1)
            factTurbHeat = one / (prandtlTurb * gm1)

            heatCoef = mul * factLamHeat + mue * factTurbHeat

            ! Compute the gradients at the face by averaging the four
            ! nodal values.

            u_x = fourth * (cudaDoms(dom,sps)%ux(i - 1, j, k - 1) + cudaDoms(dom,sps)%ux(i, j, k - 1) &
                            + cudaDoms(dom,sps)%ux(i - 1, j, k) + cudaDoms(dom,sps)%ux(i, j, k))
            u_y = fourth * (cudaDoms(dom,sps)%uy(i - 1, j, k - 1) + cudaDoms(dom,sps)%uy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%uy(i - 1, j, k) + cudaDoms(dom,sps)%uy(i, j, k))
            u_z = fourth * (cudaDoms(dom,sps)%uz(i - 1, j, k - 1) + cudaDoms(dom,sps)%uz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%uz(i - 1, j, k) + cudaDoms(dom,sps)%uz(i, j, k))

            v_x = fourth * (cudaDoms(dom,sps)%vx(i - 1, j, k - 1) + cudaDoms(dom,sps)%vx(i, j, k - 1) &
                            + cudaDoms(dom,sps)%vx(i - 1, j, k) + cudaDoms(dom,sps)%vx(i, j, k))
            v_y = fourth * (cudaDoms(dom,sps)%vy(i - 1, j, k - 1) + cudaDoms(dom,sps)%vy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%vy(i - 1, j, k) + cudaDoms(dom,sps)%vy(i, j, k))
            v_z = fourth * (cudaDoms(dom,sps)%vz(i - 1, j, k - 1) + cudaDoms(dom,sps)%vz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%vz(i - 1, j, k) + cudaDoms(dom,sps)%vz(i, j, k))

            w_x = fourth * (cudaDoms(dom,sps)%wx(i - 1, j, k - 1) + cudaDoms(dom,sps)%wx(i, j, k - 1) &
                            + cudaDoms(dom,sps)%wx(i - 1, j, k) + cudaDoms(dom,sps)%wx(i, j, k))
            w_y = fourth * (cudaDoms(dom,sps)%wy(i - 1, j, k - 1) + cudaDoms(dom,sps)%wy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%wy(i - 1, j, k) + cudaDoms(dom,sps)%wy(i, j, k))
            w_z = fourth * (cudaDoms(dom,sps)%wz(i - 1, j, k - 1) + cudaDoms(dom,sps)%wz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%wz(i - 1, j, k) + cudaDoms(dom,sps)%wz(i, j, k))

            q_x = fourth * (cudaDoms(dom,sps)%qx(i - 1, j, k - 1) + cudaDoms(dom,sps)%qx(i, j, k - 1) &
                            + cudaDoms(dom,sps)%qx(i - 1, j, k) + cudaDoms(dom,sps)%qx(i, j, k))
            q_y = fourth * (cudaDoms(dom,sps)%qy(i - 1, j, k - 1) + cudaDoms(dom,sps)%qy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%qy(i - 1, j, k) + cudaDoms(dom,sps)%qy(i, j, k))
            q_z = fourth * (cudaDoms(dom,sps)%qz(i - 1, j, k - 1) + cudaDoms(dom,sps)%qz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%qz(i - 1, j, k) + cudaDoms(dom,sps)%qz(i, j, k))

            ! The gradients in the normal direction are corrected, such
            ! that no averaging takes places here.
            ! First determine the vector in the direction from the
            ! cell center j to cell center j+1.

            ssx = eighth * (cudaDoms(dom,sps)%x(i - 1, j + 1, k - 1, 1) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 1) &
                            + cudaDoms(dom,sps)%x(i - 1, j + 1, k, 1) - cudaDoms(dom,sps)%x(i - 1, j - 1, k, 1) &
                            + cudaDoms(dom,sps)%x(i, j + 1, k - 1, 1) - cudaDoms(dom,sps)%x(i, j - 1, k - 1, 1) &
                            + cudaDoms(dom,sps)%x(i, j + 1, k, 1) - cudaDoms(dom,sps)%x(i, j - 1, k, 1))
            ssy = eighth * (cudaDoms(dom,sps)%x(i - 1, j + 1, k - 1, 2) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 2) &
                            + cudaDoms(dom,sps)%x(i - 1, j + 1, k, 2) - cudaDoms(dom,sps)%x(i - 1, j - 1, k, 2) &
                            + cudaDoms(dom,sps)%x(i, j + 1, k - 1, 2) - cudaDoms(dom,sps)%x(i, j - 1, k - 1, 2) &
                            + cudaDoms(dom,sps)%x(i, j + 1, k, 2) - cudaDoms(dom,sps)%x(i, j - 1, k, 2))
            ssz = eighth * (cudaDoms(dom,sps)%x(i - 1, j + 1, k - 1, 3) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 3) &
                            + cudaDoms(dom,sps)%x(i - 1, j + 1, k, 3) - cudaDoms(dom,sps)%x(i - 1, j - 1, k, 3) &
                            + cudaDoms(dom,sps)%x(i, j + 1, k - 1, 3) - cudaDoms(dom,sps)%x(i, j - 1, k - 1, 3) &
                            + cudaDoms(dom,sps)%x(i, j + 1, k, 3) - cudaDoms(dom,sps)%x(i, j - 1, k, 3))

            ! Determine the length of this vector and create the
            ! unit normal.

            snrm = one / sqrt(ssx * ssx + ssy * ssy + ssz * ssz)
            ssx = snrm * ssx
            ssy = snrm * ssy
            ssz = snrm * ssz

            ! Correct the gradients.

            corr = u_x * ssx + u_y * ssy + u_z * ssz &
                    - (cudaDoms(dom,sps)%w(i, j + 1, k, ivx) - cudaDoms(dom,sps)%w(i, j, k, ivx)) * snrm
            u_x = u_x - corr * ssx
            u_y = u_y - corr * ssy
            u_z = u_z - corr * ssz

            corr = v_x * ssx + v_y * ssy + v_z * ssz &
                    - (cudaDoms(dom,sps)%w(i, j + 1, k, ivy) - cudaDoms(dom,sps)%w(i, j, k, ivy)) * snrm
            v_x = v_x - corr * ssx
            v_y = v_y - corr * ssy
            v_z = v_z - corr * ssz

            corr = w_x * ssx + w_y * ssy + w_z * ssz &
                    - (cudaDoms(dom,sps)%w(i, j + 1, k, ivz) - cudaDoms(dom,sps)%w(i, j, k, ivz)) * snrm
            w_x = w_x - corr * ssx
            w_y = w_y - corr * ssy
            w_z = w_z - corr * ssz

            corr = q_x * ssx + q_y * ssy + q_z * ssz &
                    + (cudaDoms(dom,sps)%aa(i, j + 1, k) - cudaDoms(dom,sps)%aa(i, j, k)) * snrm
            q_x = q_x - corr * ssx
            q_y = q_y - corr * ssy
            q_z = q_z - corr * ssz

            ! Compute the stress tensor and the heat flux vector.
            ! We remove the viscosity from the stress tensor (tau)
            ! to define tauS since we still need to separate between
            ! laminar and turbulent stress for QCR.
            ! Therefore, laminar tau = mue*tauS, turbulent
            ! tau = mue*tauS, and total tau = mut*tauS.

            fracDiv = twoThird * (u_x + v_y + w_z)

            tauxxS = two * u_x - fracDiv
            tauyyS = two * v_y - fracDiv
            tauzzS = two * w_z - fracDiv

            tauxyS = u_y + v_x
            tauxzS = u_z + w_x
            tauyzS = v_z + w_y

            q_x = heatCoef * q_x
            q_y = heatCoef * q_y
            q_z = heatCoef * q_z

            ! Add QCR corrections if necessary
            if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x * u_x + u_y * u_y + u_z * u_z + &
                            v_x * v_x + v_y * v_y + v_z * v_z + &
                            w_x * w_x + w_y * w_y + w_z * w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue * Ccr1 / den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact * (Wxy * tauxyS + Wxz * tauxzS) * two
                eyy = fact * (Wyx * tauxyS + Wyz * tauyzS) * two
                ezz = fact * (Wzx * tauxzS + Wzy * tauyzS) * two

                exy = fact * (Wxy * tauyyS + Wxz * tauyzS + &
                                Wyx * tauxxS + Wyz * tauxzS)
                exz = fact * (Wxy * tauyzS + Wxz * tauzzS + &
                                Wzx * tauxxS + Wzy * tauxyS)
                eyz = fact * (Wyx * tauxzS + Wyz * tauzzS + &
                                Wzx * tauxyS + Wzy * tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut * tauxxS - exx
                tauyy = mut * tauyyS - eyy
                tauzz = mut * tauzzS - ezz
                tauxy = mut * tauxyS - exy
                tauxz = mut * tauxzS - exz
                tauyz = mut * tauyzS - eyz

            else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut * tauxxS
                tauyy = mut * tauyyS
                tauzz = mut * tauzzS
                tauxy = mut * tauxyS
                tauxz = mut * tauxzS
                tauyz = mut * tauyzS

            end if

            ! Compute the average velocities for the face. Remember that
            ! the velocities are stored and not the momentum.

            ubar = half * (cudaDoms(dom,sps)%w(i, j, k, ivx) + cudaDoms(dom,sps)%w(i, j + 1, k, ivx))
            vbar = half * (cudaDoms(dom,sps)%w(i, j, k, ivy) + cudaDoms(dom,sps)%w(i, j + 1, k, ivy))
            wbar = half * (cudaDoms(dom,sps)%w(i, j, k, ivz) + cudaDoms(dom,sps)%w(i, j + 1, k, ivz))

            ! Compute the viscous fluxes for this j-face.

            fmx = tauxx * cudaDoms(dom,sps)%sJ(i, j, k, 1) + tauxy * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    + tauxz * cudaDoms(dom,sps)%sJ(i, j, k, 3)
            fmy = tauxy * cudaDoms(dom,sps)%sJ(i, j, k, 1) + tauyy * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    + tauyz * cudaDoms(dom,sps)%sJ(i, j, k, 3)
            fmz = tauxz * cudaDoms(dom,sps)%sJ(i, j, k, 1) + tauyz * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    + tauzz * cudaDoms(dom,sps)%sJ(i, j, k, 3)
            frhoE = (ubar * tauxx + vbar * tauxy + wbar * tauxz) * cudaDoms(dom,sps)%sJ(i, j, k, 1) &
                    + (ubar * tauxy + vbar * tauyy + wbar * tauyz) * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    + (ubar * tauxz + vbar * tauyz + wbar * tauzz) * cudaDoms(dom,sps)%sJ(i, j, k, 3) &
                    - q_x * cudaDoms(dom,sps)%sJ(i, j, k, 1) - q_y * cudaDoms(dom,sps)%sJ(i, j, k, 2) - q_z * cudaDoms(dom,sps)%sJ(i, j, k, 3)

            ! Update the residuals of cell j and j+1.

            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imx), fmx)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imy), fmy)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imz), fmz)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, irhoE), frhoE)
            ! cudaDoms(dom,sps)%fw(i, j, k, imx) = cudaDoms(dom,sps)%fw(i, j, k, imx) - fmx
            ! cudaDoms(dom,sps)%fw(i, j, k, imy) = cudaDoms(dom,sps)%fw(i, j, k, imy) - fmy
            ! cudaDoms(dom,sps)%fw(i, j, k, imz) = cudaDoms(dom,sps)%fw(i, j, k, imz) - fmz
            ! cudaDoms(dom,sps)%fw(i, j, k, irhoE) = cudaDoms(dom,sps)%fw(i, j, k, irhoE) - frhoE

            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, imx), fmx)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, imy), fmy)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, imz), fmz)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i, j + 1, k, irhoE), frhoE)
            ! cudaDoms(dom,sps)%fw(i, j + 1, k, imx) = cudaDoms(dom,sps)%fw(i, j + 1, k, imx) + fmx
            ! cudaDoms(dom,sps)%fw(i, j + 1, k, imy) = cudaDoms(dom,sps)%fw(i, j + 1, k, imy) + fmy
            ! cudaDoms(dom,sps)%fw(i, j + 1, k, imz) = cudaDoms(dom,sps)%fw(i, j + 1, k, imz) + fmz
            ! cudaDoms(dom,sps)%fw(i, j + 1, k, irhoE) = cudaDoms(dom,sps)%fw(i, j + 1, k, irhoE) + frhoE
        end if  
        !
        !         Viscous fluxes in the i-direction.
        !
        !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 1 to cudaDoms(dom,sps)%il
        if (k <= cudaDoms(dom,sps)%kl .and. j <= cudaDoms(dom,sps)%jl .and. i <= cudaDoms(dom,sps)%il .and. k>=2 .and. j>=2 .and. i>=1) then
            ! Set the value of the porosity. If not zero, it is set
            ! to average the eddy-viscosity and to take the factor
            ! rFilv into account.

            por = half * rFilv
            if (cudaDoms(dom,sps)%porI(i, j, k) == noFlux) por = zero

            ! Compute the laminar and (if present) the eddy viscosities
            ! multiplied the porosity. Compute the factor in front of
            ! the gradients of the speed of sound squared for the heat
            ! flux.

            mul = por * (cudaDoms(dom,sps)%rlv(i, j, k) + cudaDoms(dom,sps)%rlv(i + 1, j, k))
            mue = por * (cudaDoms(dom,sps)%rev(i, j, k) + cudaDoms(dom,sps)%rev(i + 1, j, k))
            mut = mul + mue

            gm1 = half * (cudaDoms(dom,sps)%gamma(i, j, k) + cudaDoms(dom,sps)%gamma(i + 1, j, k)) - one
            factLamHeat = one / (prandtl * gm1)
            factTurbHeat = one / (prandtlTurb * gm1)

            heatCoef = mul * factLamHeat + mue * factTurbHeat

            ! Compute the gradients at the face by averaging the four
            ! nodal values.

            u_x = fourth * (cudaDoms(dom,sps)%ux(i, j - 1, k - 1) + cudaDoms(dom,sps)%ux(i, j, k - 1) &
                            + cudaDoms(dom,sps)%ux(i, j - 1, k) + cudaDoms(dom,sps)%ux(i, j, k))
            u_y = fourth * (cudaDoms(dom,sps)%uy(i, j - 1, k - 1) + cudaDoms(dom,sps)%uy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%uy(i, j - 1, k) + cudaDoms(dom,sps)%uy(i, j, k))
            u_z = fourth * (cudaDoms(dom,sps)%uz(i, j - 1, k - 1) + cudaDoms(dom,sps)%uz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%uz(i, j - 1, k) + cudaDoms(dom,sps)%uz(i, j, k))

            v_x = fourth * (cudaDoms(dom,sps)%vx(i, j - 1, k - 1) + cudaDoms(dom,sps)%vx(i, j, k - 1) &
                            + cudaDoms(dom,sps)%vx(i, j - 1, k) + cudaDoms(dom,sps)%vx(i, j, k))
            v_y = fourth * (cudaDoms(dom,sps)%vy(i, j - 1, k - 1) + cudaDoms(dom,sps)%vy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%vy(i, j - 1, k) + cudaDoms(dom,sps)%vy(i, j, k))
            v_z = fourth * (cudaDoms(dom,sps)%vz(i, j - 1, k - 1) + cudaDoms(dom,sps)%vz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%vz(i, j - 1, k) + cudaDoms(dom,sps)%vz(i, j, k))

            w_x = fourth * (cudaDoms(dom,sps)%wx(i, j - 1, k - 1) + cudaDoms(dom,sps)%wx(i, j, k - 1) &
                            + cudaDoms(dom,sps)%wx(i, j - 1, k) + cudaDoms(dom,sps)%wx(i, j, k))
            w_y = fourth * (cudaDoms(dom,sps)%wy(i, j - 1, k - 1) + cudaDoms(dom,sps)%wy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%wy(i, j - 1, k) + cudaDoms(dom,sps)%wy(i, j, k))
            w_z = fourth * (cudaDoms(dom,sps)%wz(i, j - 1, k - 1) + cudaDoms(dom,sps)%wz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%wz(i, j - 1, k) + cudaDoms(dom,sps)%wz(i, j, k))

            q_x = fourth * (cudaDoms(dom,sps)%qx(i, j - 1, k - 1) + cudaDoms(dom,sps)%qx(i, j, k - 1) &
                            + cudaDoms(dom,sps)%qx(i, j - 1, k) + cudaDoms(dom,sps)%qx(i, j, k))
            q_y = fourth * (cudaDoms(dom,sps)%qy(i, j - 1, k - 1) + cudaDoms(dom,sps)%qy(i, j, k - 1) &
                            + cudaDoms(dom,sps)%qy(i, j - 1, k) + cudaDoms(dom,sps)%qy(i, j, k))
            q_z = fourth * (cudaDoms(dom,sps)%qz(i, j - 1, k - 1) + cudaDoms(dom,sps)%qz(i, j, k - 1) &
                            + cudaDoms(dom,sps)%qz(i, j - 1, k) + cudaDoms(dom,sps)%qz(i, j, k))

            ! The gradients in the normal direction are corrected, such
            ! that no averaging takes places here.
            ! First determine the vector in the direction from the
            ! cell center i to cell center i+1.

            ssx = eighth * (cudaDoms(dom,sps)%x(i + 1, j - 1, k - 1, 1) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 1) &
                            + cudaDoms(dom,sps)%x(i + 1, j - 1, k, 1) - cudaDoms(dom,sps)%x(i - 1, j - 1, k, 1) &
                            + cudaDoms(dom,sps)%x(i + 1, j, k - 1, 1) - cudaDoms(dom,sps)%x(i - 1, j, k - 1, 1) &
                            + cudaDoms(dom,sps)%x(i + 1, j, k, 1) - cudaDoms(dom,sps)%x(i - 1, j, k, 1))
            ssy = eighth * (cudaDoms(dom,sps)%x(i + 1, j - 1, k - 1, 2) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 2) &
                            + cudaDoms(dom,sps)%x(i + 1, j - 1, k, 2) - cudaDoms(dom,sps)%x(i - 1, j - 1, k, 2) &
                            + cudaDoms(dom,sps)%x(i + 1, j, k - 1, 2) - cudaDoms(dom,sps)%x(i - 1, j, k - 1, 2) &
                            + cudaDoms(dom,sps)%x(i + 1, j, k, 2) - cudaDoms(dom,sps)%x(i - 1, j, k, 2))
            ssz = eighth * (cudaDoms(dom,sps)%x(i + 1, j - 1, k - 1, 3) - cudaDoms(dom,sps)%x(i - 1, j - 1, k - 1, 3) &
                            + cudaDoms(dom,sps)%x(i + 1, j - 1, k, 3) - cudaDoms(dom,sps)%x(i - 1, j - 1, k, 3) &
                            + cudaDoms(dom,sps)%x(i + 1, j, k - 1, 3) - cudaDoms(dom,sps)%x(i - 1, j, k - 1, 3) &
                            + cudaDoms(dom,sps)%x(i + 1, j, k, 3) - cudaDoms(dom,sps)%x(i - 1, j, k, 3))

            ! Determine the length of this vector and create the
            ! unit normal.

            snrm = one / sqrt(ssx * ssx + ssy * ssy + ssz * ssz)
            ssx = snrm * ssx
            ssy = snrm * ssy
            ssz = snrm * ssz

            ! Correct the gradients.

            corr = u_x * ssx + u_y * ssy + u_z * ssz &
                   - (cudaDoms(dom,sps)%w(i + 1, j, k, ivx) - cudaDoms(dom,sps)%w(i, j, k, ivx)) * snrm
            u_x = u_x - corr * ssx
            u_y = u_y - corr * ssy
            u_z = u_z - corr * ssz

            corr = v_x * ssx + v_y * ssy + v_z * ssz &
                   - (cudaDoms(dom,sps)%w(i + 1, j, k, ivy) - cudaDoms(dom,sps)%w(i, j, k, ivy)) * snrm
            v_x = v_x - corr * ssx
            v_y = v_y - corr * ssy
            v_z = v_z - corr * ssz

            corr = w_x * ssx + w_y * ssy + w_z * ssz &
                   - (cudaDoms(dom,sps)%w(i + 1, j, k, ivz) - cudaDoms(dom,sps)%w(i, j, k, ivz)) * snrm
            w_x = w_x - corr * ssx
            w_y = w_y - corr * ssy
            w_z = w_z - corr * ssz

            corr = q_x * ssx + q_y * ssy + q_z * ssz &
                   + (cudaDoms(dom,sps)%aa(i + 1, j, k) - cudaDoms(dom,sps)%aa(i, j, k)) * snrm
            q_x = q_x - corr * ssx
            q_y = q_y - corr * ssy
            q_z = q_z - corr * ssz

            ! Compute the stress tensor and the heat flux vector.
            ! We remove the viscosity from the stress tensor (tau)
            ! to define tauS since we still need to separate between
            ! laminar and turbulent stress for QCR.
            ! Therefore, laminar tau = mue*tauS, turbulent
            ! tau = mue*tauS, and total tau = mut*tauS.

            fracDiv = twoThird * (u_x + v_y + w_z)

            tauxxS = two * u_x - fracDiv
            tauyyS = two * v_y - fracDiv
            tauzzS = two * w_z - fracDiv

            tauxyS = u_y + v_x
            tauxzS = u_z + w_x
            tauyzS = v_z + w_y

            q_x = heatCoef * q_x
            q_y = heatCoef * q_y
            q_z = heatCoef * q_z

            ! Add QCR corrections if necessary
            if (useQCR) then

                ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                !
                ! tau_ij,QCR = tau_ij - e_ij
                !
                ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                !
                ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                !
                ! We are computing O_ik as follows:
                !
                ! O_ik = 2*W_ik/den
                !
                ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                ! Compute denominator
                den = sqrt(u_x * u_x + u_y * u_y + u_z * u_z + &
                           v_x * v_x + v_y * v_y + v_z * v_z + &
                           w_x * w_x + w_y * w_y + w_z * w_z)

                ! Denominator should be limited to avoid division by zero in regions with
                ! no gradients
                den = max(den, xminn)

                ! Compute factor that will multiply all tensor components.
                ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                ! components as well.
                fact = mue * Ccr1 / den

                ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                ! The diagonals of the vorticity tensor components are always zero
                Wxy = u_y - v_x
                Wxz = u_z - w_x
                Wyz = v_z - w_y
                Wyx = -Wxy
                Wzx = -Wxz
                Wzy = -Wyz

                ! Compute the extra terms of the Boussinesq relation
                exx = fact * (Wxy * tauxyS + Wxz * tauxzS) * two
                eyy = fact * (Wyx * tauxyS + Wyz * tauyzS) * two
                ezz = fact * (Wzx * tauxzS + Wzy * tauyzS) * two

                exy = fact * (Wxy * tauyyS + Wxz * tauyzS + &
                              Wyx * tauxxS + Wyz * tauxzS)
                exz = fact * (Wxy * tauyzS + Wxz * tauzzS + &
                              Wzx * tauxxS + Wzy * tauxyS)
                eyz = fact * (Wyx * tauxzS + Wyz * tauzzS + &
                              Wzx * tauxyS + Wzy * tauyyS)

                ! Apply the total viscosity to the stress tensor and add extra terms
                tauxx = mut * tauxxS - exx
                tauyy = mut * tauyyS - eyy
                tauzz = mut * tauzzS - ezz
                tauxy = mut * tauxyS - exy
                tauxz = mut * tauxzS - exz
                tauyz = mut * tauyzS - eyz

            else

                ! Just apply the total viscosity to the stress tensor
                tauxx = mut * tauxxS
                tauyy = mut * tauyyS
                tauzz = mut * tauzzS
                tauxy = mut * tauxyS
                tauxz = mut * tauxzS
                tauyz = mut * tauyzS

            end if

            ! Compute the average velocities for the face. Remember that
            ! the velocities are stored and not the momentum.

            ubar = half * (cudaDoms(dom,sps)%w(i, j, k, ivx) + cudaDoms(dom,sps)%w(i + 1, j, k, ivx))
            vbar = half * (cudaDoms(dom,sps)%w(i, j, k, ivy) + cudaDoms(dom,sps)%w(i + 1, j, k, ivy))
            wbar = half * (cudaDoms(dom,sps)%w(i, j, k, ivz) + cudaDoms(dom,sps)%w(i + 1, j, k, ivz))

            ! Compute the viscous fluxes for this i-face.

            fmx = tauxx * cudaDoms(dom,sps)%sI(i, j, k, 1) + tauxy * cudaDoms(dom,sps)%sI(i, j, k, 2) &
                  + tauxz * cudaDoms(dom,sps)%sI(i, j, k, 3)
            fmy = tauxy * cudaDoms(dom,sps)%sI(i, j, k, 1) + tauyy * cudaDoms(dom,sps)%sI(i, j, k, 2) &
                  + tauyz * cudaDoms(dom,sps)%sI(i, j, k, 3)
            fmz = tauxz * cudaDoms(dom,sps)%sI(i, j, k, 1) + tauyz * cudaDoms(dom,sps)%sI(i, j, k, 2) &
                  + tauzz * cudaDoms(dom,sps)%sI(i, j, k, 3)
            frhoE = (ubar * tauxx + vbar * tauxy + wbar * tauxz) * cudaDoms(dom,sps)%sI(i, j, k, 1) &
                    + (ubar * tauxy + vbar * tauyy + wbar * tauyz) * cudaDoms(dom,sps)%sI(i, j, k, 2) &
                    + (ubar * tauxz + vbar * tauyz + wbar * tauzz) * cudaDoms(dom,sps)%sI(i, j, k, 3) &
                    - q_x * cudaDoms(dom,sps)%sI(i, j, k, 1) - q_y * cudaDoms(dom,sps)%sI(i, j, k, 2) - q_z * cudaDoms(dom,sps)%sI(i, j, k, 3)

            ! Update the residuals of cell i and i+1.
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imx), fmx)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imy), fmy)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, imz), fmz)
            tmp = atomicsub(cudaDoms(dom,sps)%fw(i, j, k, irhoE), frhoE)
            ! cudaDoms(dom,sps)%fw(i, j, k, imx) = cudaDoms(dom,sps)%fw(i, j, k, imx) - fmx
            ! cudaDoms(dom,sps)%fw(i, j, k, imy) = cudaDoms(dom,sps)%fw(i, j, k, imy) - fmy
            ! cudaDoms(dom,sps)%fw(i, j, k, imz) = cudaDoms(dom,sps)%fw(i, j, k, imz) - fmz
            ! cudaDoms(dom,sps)%fw(i, j, k, irhoE) = cudaDoms(dom,sps)%fw(i, j, k, irhoE) - frhoE
            
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, imx), fmx)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, imy), fmy)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, imz), fmz)
            tmp = atomicadd(cudaDoms(dom,sps)%fw(i + 1, j, k, irhoE), frhoE)
            ! cudaDoms(dom,sps)%fw(i + 1, j, k, imx) = cudaDoms(dom,sps)%fw(i + 1, j, k, imx) + fmx
            ! cudaDoms(dom,sps)%fw(i + 1, j, k, imy) = cudaDoms(dom,sps)%fw(i + 1, j, k, imy) + fmy
            ! cudaDoms(dom,sps)%fw(i + 1, j, k, imz) = cudaDoms(dom,sps)%fw(i + 1, j, k, imz) + fmz
            ! cudaDoms(dom,sps)%fw(i + 1, j, k, irhoE) = cudaDoms(dom,sps)%fw(i + 1, j, k, irhoE) + frhoE

        end if
    end subroutine viscousFlux

    ! (hannah)
    attributes(global) subroutine saSource
        use constants, only: fourth, two, one, zero, sixth
        use paramTurb, only: rsaCv1, rsaCw3, rsaK, rsaCb3, rsaCw2, rsaCb1, rsaCw1, rsact3,rsact4
        use cudaInputPhysics, only: useft2SA, turbProd, equations
        use cudaInputDiscretization, only: approxSA
        !use sa, only: cv13, kar2Inv, cw36, cb3Inv
        use cudaFlowvarRefState, only: timeRef

        implicit none

        ! variables for SA source
        real(kind=realType) :: fv1, fv2, ft2
        real(kind=realType) :: sst, nu, dist2Inv, chi, chi2, chi3
        real(kind=realType) :: rr, gg, gg6, termFw, fwSa, term1, term2
        real(kind=realType) :: dfv1, dfv2, dft2, drr, dgg, dfw, sqrtProd
        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
        real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
        real(kind=realType) :: vortx, vorty, vortz
        real(kind=realType) :: strainMag2, prod
        real(kind=realType), parameter :: xminn = 1.e-10_realType
        real(kind=realType), parameter :: f23 = two * third
        integer(kind=intType) :: i, j, k, dom, sps
        real(kind=realType) :: term1Fact

        !this is new no longer read these from sa module
        !since they are just computed here anyways
        real(kind=realType) :: cv13, kar2Inv, cw36, cb3Inv


        dom = 1 
        sps = 1


        ! Set model constants
        cv13 = rsaCv1**3
        kar2Inv = one / (rsaK**2)
        cw36 = rsaCw3**6
        cb3Inv = one / rsaCb3

        ! set the approximate multiplier here
        term1Fact = one
        if (approxSA) term1Fact = zero
        
        ! skipping omega calculations because we are ignoring sections and rotation

        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1

        !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then

            ! Compute the gradient of u in the cell center. Use is made
            ! of the fact that the surrounding normals sum up to zero,
            ! such that the cell i,j,k does not give a contribution.
            ! The gradient is scaled by the factor 2*vol.

            uux = cudaDoms(dom,sps)%w(i + 1, j, k, ivx) * cudaDoms(dom,sps)%sI(i, j, k, 1) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivx) * cudaDoms(dom,sps)%sI(i - 1, j, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivx) * cudaDoms(dom,sps)%sJ(i, j, k, 1) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivx) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivx) * cudaDoms(dom,sps)%sK(i, j, k, 1) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivx) * cudaDoms(dom,sps)%sK(i, j, k - 1, 1)
            uuy = cudaDoms(dom,sps)%w(i + 1, j, k, ivx) * cudaDoms(dom,sps)%sI(i, j, k, 2) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivx) * cudaDoms(dom,sps)%sI(i - 1, j, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivx) * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivx) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivx) * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivx) * cudaDoms(dom,sps)%sK(i, j, k - 1, 2)
            uuz = cudaDoms(dom,sps)%w(i + 1, j, k, ivx) * cudaDoms(dom,sps)%sI(i, j, k, 3) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivx) * cudaDoms(dom,sps)%sI(i - 1, j, k, 3) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivx) * cudaDoms(dom,sps)%sJ(i, j, k, 3) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivx) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 3) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivx) * cudaDoms(dom,sps)%sK(i, j, k, 3) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivx) * cudaDoms(dom,sps)%sK(i, j, k - 1, 3)

            ! Idem for the gradient of v.

            vvx = cudaDoms(dom,sps)%w(i + 1, j, k, ivy) * cudaDoms(dom,sps)%sI(i, j, k, 1) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivy) * cudaDoms(dom,sps)%sI(i - 1, j, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivy) * cudaDoms(dom,sps)%sJ(i, j, k, 1) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivy) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivy) * cudaDoms(dom,sps)%sK(i, j, k, 1) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivy) * cudaDoms(dom,sps)%sK(i, j, k - 1, 1)
            vvy = cudaDoms(dom,sps)%w(i + 1, j, k, ivy) * cudaDoms(dom,sps)%sI(i, j, k, 2) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivy) * cudaDoms(dom,sps)%sI(i - 1, j, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivy) * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivy) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivy) * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivy) * cudaDoms(dom,sps)%sK(i, j, k - 1, 2)
            vvz = cudaDoms(dom,sps)%w(i + 1, j, k, ivy) * cudaDoms(dom,sps)%sI(i, j, k, 3) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivy) * cudaDoms(dom,sps)%sI(i - 1, j, k, 3) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivy) * cudaDoms(dom,sps)%sJ(i, j, k, 3) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivy) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 3) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivy) * cudaDoms(dom,sps)%sK(i, j, k, 3) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivy) * cudaDoms(dom,sps)%sK(i, j, k - 1, 3)

            ! And for the gradient of w.

            wwx = cudaDoms(dom,sps)%w(i + 1, j, k, ivz) * cudaDoms(dom,sps)%sI(i, j, k, 1) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivz) * cudaDoms(dom,sps)%sI(i - 1, j, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivz) * cudaDoms(dom,sps)%sJ(i, j, k, 1) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivz) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 1) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivz) * cudaDoms(dom,sps)%sK(i, j, k, 1) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivz) * cudaDoms(dom,sps)%sK(i, j, k - 1, 1)
            wwy = cudaDoms(dom,sps)%w(i + 1, j, k, ivz) * cudaDoms(dom,sps)%sI(i, j, k, 2) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivz) * cudaDoms(dom,sps)%sI(i - 1, j, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivz) * cudaDoms(dom,sps)%sJ(i, j, k, 2) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivz) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 2) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivz) * cudaDoms(dom,sps)%sK(i, j, k, 2) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivz) * cudaDoms(dom,sps)%sK(i, j, k - 1, 2)
            wwz = cudaDoms(dom,sps)%w(i + 1, j, k, ivz) * cudaDoms(dom,sps)%sI(i, j, k, 3) &
            - cudaDoms(dom,sps)%w(i - 1, j, k, ivz) * cudaDoms(dom,sps)%sI(i - 1, j, k, 3) &
                    + cudaDoms(dom,sps)%w(i, j + 1, k, ivz) * cudaDoms(dom,sps)%sJ(i, j, k, 3) &
                    - cudaDoms(dom,sps)%w(i, j - 1, k, ivz) * cudaDoms(dom,sps)%sJ(i, j - 1, k, 3) &
                    + cudaDoms(dom,sps)%w(i, j, k + 1, ivz) * cudaDoms(dom,sps)%sK(i, j, k, 3) &
                    - cudaDoms(dom,sps)%w(i, j, k - 1, ivz) * cudaDoms(dom,sps)%sK(i, j, k - 1, 3)

            ! Compute the components of the stress tensor.
            ! The combination of the current scaling of the velocity
            ! gradients (2*vol) and the definition of the stress tensor,
            ! leads to the factor 1/(4*vol).

            fact = fourth / cudaDoms(dom,sps)%vol(i, j, k)

            ! -- Calcs for strain --
            sxx = two * fact * uux
            syy = two * fact * vvy
            szz = two * fact * wwz

            sxy = fact * (uuy + vvx)
            sxz = fact * (uuz + wwx)
            syz = fact * (vvz + wwy)

            ! Compute 2/3 * divergence of velocity squared

            div2 = f23 * (sxx + syy + szz)**2

            ! Compute strain production term

            strainMag2 = two * (sxy**2 + sxz**2 + syz**2) &
                            + sxx**2 + syy**2 + szz**2

            ! -- Calcs for vorticity --

            ! Compute the three components of the vorticity vector.
            ! Substract the part coming from the rotating frame.

            ! took omega out of this calculation
            vortx = two * fact * (wwy - vvz)
            vorty = two * fact * (uuz - wwx)
            vortz = two * fact * (vvx - uuy)

            if (turbProd == strain) then
                sqrtProd = sqrt(max(two * strainMag2 - div2, eps))
            else
                sqrtProd = sqrt(vortx**2 + vorty**2 + vortz**2)
            end if

            ! Compute the laminar kinematic viscosity, the inverse of
            ! wall distance squared, the ratio chi (ratio of nuTilde
            ! and nu) and the functions fv1 and fv2. The latter corrects
            ! the production term near a viscous wall.

            nu = cudaDoms(dom,sps)%rlv(i, j, k) / cudaDoms(dom,sps)%w(i, j, k, irho)
            dist2Inv = one / (cudaDoms(dom,sps)%d2wall(i, j, k)**2)
            chi = cudaDoms(dom,sps)%w(i, j, k, itu1) / nu
            chi2 = chi * chi
            chi3 = chi * chi2
            fv1 = chi3 / (chi3 + cv13)
            fv2 = one - chi / (one + chi * fv1)

            ! The function ft2, which is designed to keep a laminar
            ! solution laminar. When running in fully turbulent mode
            ! this function should be set to 0.0.

            ft2 = zero
            if (useft2SA) then
                ft2 = rsaCt3 * exp(-rsaCt4 * chi2)
            end if

            ! Correct the production term to account for the influence
            ! of the wall.

            sst = sqrtProd + cudaDoms(dom,sps)%w(i, j, k, itu1) * fv2 * kar2Inv * dist2Inv

            ! removed useRotationSA

            ! Make sure that this term remains positive
            ! (the function fv2 is negative between chi = 1 and 18.4,
            ! which can cause sst to go negative, which is undesirable).

            sst = max(sst, xminn)

            ! Compute the function cudaDoms(dom,sps)%fw. The argument rr is cut off at 10
            ! to avoid numerical problems. This is ok, because the
            ! asymptotical value of cudaDoms(dom,sps)%fw is then already reached.

            rr = cudaDoms(dom,sps)%w(i, j, k, itu1) * kar2Inv * dist2Inv / sst
            rr = min(rr, 10.0_realType)
            gg = rr + rsaCw2 * (rr**6 - rr)
            gg6 = gg**6
            termFw = ((one + cw36) / (gg6 + cw36))**sixth
            fwSa = gg * termFw

            ! Compute the source term; some terms are saved for the
            ! linearization. The source term is stored in dvt.

            term1 = rsaCb1 * (one - ft2) * sqrtProd * term1Fact
            term2 = dist2Inv * (kar2Inv * rsaCb1 * ((one - ft2) * fv2 + ft2) &
                                - rsaCw1 * fwSa)

            cudaDoms(dom,sps)%dw(i, j, k, itu1) = cudaDoms(dom,sps)%dw(i, j, k, itu1) &
            + (term1 + term2 * cudaDoms(dom,sps)%w(i, j, k, itu1)) * cudaDoms(dom,sps)%w(i, j, k, itu1)
        
        end if

    end subroutine saSource

    ! (HANNAH)
    attributes(global) subroutine saAdvection
        use constants, only: zero, half
        use cudaInputDiscretization, only: orderTurb
        use cudaiteration, only: groundlevel
        use cudaTurbMod, only: secondOrd
        implicit none

        ! Variables for sa Advection
        real(kind=realType) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk, qs
        real(kind=realType) :: voli, xa, ya, za
        integer(kind=intType), parameter :: nAdv = 1
        integer(kind=intType) :: offset, i, j, k, ii, jj, dom, sps

        dom = 1 
        sps = 1

        ! Determine whether or not a second order discretization for the
        ! advective terms must be used.
        secondOrd = .false.
        if (groundLevel == 1_intType .and. &
            orderTurb == secondOrder) secondOrd = .true.

        offset = itu1 - 1
        !+1 loop sart at 2
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
        !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then

            ! Compute the grid velocity if present.
            ! It is taken as the average of k and k-1,

            voli = half / cudaDoms(dom,sps)%vol(i, j, k)
            qs = (cudaDoms(dom,sps)%sFaceK(i, j, k) + cudaDoms(dom,sps)%sFaceK(i, j, k - 1)) * voli

            ! Compute the normal velocity, where the normal direction
            ! is taken as the average of faces k and k-1.

            xa = (cudaDoms(dom,sps)%sK(i, j, k, 1) + cudaDoms(dom,sps)%sK(i, j, k - 1, 1)) * voli
            ya = (cudaDoms(dom,sps)%sK(i, j, k, 2) + cudaDoms(dom,sps)%sK(i, j, k - 1, 2)) * voli
            za = (cudaDoms(dom,sps)%sK(i, j, k, 3) + cudaDoms(dom,sps)%sK(i, j, k - 1, 3)) * voli

            uu = xa * cudaDoms(dom,sps)%w(i, j, k, ivx) + ya * cudaDoms(dom,sps)%w(i, j, k, ivy) &
                                        + za * cudaDoms(dom,sps)%w(i, j, k, ivz) - qs

            ! Determine the situation we are having here, i.e. positive
            ! or negative normal velocity.

            velKdir: if (uu > zero) then

                ! Velocity has a component in positive k-direction.
                ! Loop over the number of advection equations.

                do ii = 1, nAdv

                    ! Set the value of jj such that it corresponds to the
                    ! turbulent entry in w.

                    jj = ii + offset

                    ! Check whether a first or a second order discretization
                    ! must be used.

                    if (secondOrd) then

                        ! Second order; store the three differences for the
                        ! discretization of the derivative in k-direction.

                        dwtm1 = cudaDoms(dom,sps)%w(i, j, k - 1, jj) - cudaDoms(dom,sps)%w(i, j, k - 2, jj)
                        dwt = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i, j, k - 1, jj)
                        dwtp1 = cudaDoms(dom,sps)%w(i, j, k + 1, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)

                        ! Construct the derivative in this cell center. This
                        ! is the first order upwind derivative with two
                        ! nonlinear corrections.

                        dwtk = dwt

                        if (dwt * dwtp1 > zero) then
                            if (abs(dwt) < abs(dwtp1)) then
                                dwtk = dwtk + half * dwt
                            else
                                dwtk = dwtk + half * dwtp1
                            end if
                        end if

                        if (dwt * dwtm1 > zero) then
                            if (abs(dwt) < abs(dwtm1)) then
                                dwtk = dwtk - half * dwt
                            else
                                dwtk = dwtk - half * dwtm1
                            end if
                        end if

                    else

                        ! 1st order upwind scheme.

                        dwtk = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i, j, k - 1, jj)

                    end if

                    ! Update the residual. The convective term must be
                    ! substracted, because it appears on the other side of
                    ! the equation as the source and viscous terms.

                    cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) = cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) - uu * dwtk
                end do

            else velKdir

                ! Velocity has a component in negative k-direction.
                ! Loop over the number of advection equations
                do ii = 1, nAdv

                    ! Set the value of jj such that it corresponds to the
                    ! turbulent entry in w.

                    jj = ii + offset

                    ! Check whether a first or a second order discretization
                    ! must be used.

                    if (secondOrd) then

                        ! Store the three differences for the discretization of
                        ! the derivative in k-direction.

                        dwtm1 = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i, j, k - 1, jj)
                        dwt = cudaDoms(dom,sps)%w(i, j, k + 1, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)
                        dwtp1 = cudaDoms(dom,sps)%w(i, j, k + 2, jj) - cudaDoms(dom,sps)%w(i, j, k + 1, jj)

                        ! Construct the derivative in this cell center. This is
                        ! the first order upwind derivative with two nonlinear
                        ! corrections.

                        dwtk = dwt

                        if (dwt * dwtp1 > zero) then
                            if (abs(dwt) < abs(dwtp1)) then
                                dwtk = dwtk - half * dwt
                            else
                                dwtk = dwtk - half * dwtp1
                            end if
                        end if

                        if (dwt * dwtm1 > zero) then
                            if (abs(dwt) < abs(dwtm1)) then
                                dwtk = dwtk + half * dwt
                            else
                                dwtk = dwtk + half * dwtm1
                            end if
                        end if

                    else

                        ! 1st order upwind scheme.

                        dwtk = cudaDoms(dom,sps)%w(i, j, k + 1, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)

                    end if

                    ! Update the residual. The convective term must be
                    ! substracted, because it appears on the other side
                    ! of the equation as the source and viscous terms.

                    cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) = cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) - uu * dwtk
                end do
            end if velKdir
        end if

        !
        !       Upwind discretization of the convective term in j (eta)
        !       direction. Either the 1st order upwind or the second order
        !       fully upwind interpolation scheme, kappa = -1, is used in
        !       combination with the minmod limiter.
        !       The possible grid velocity must be taken into account.
        !
        
        !these start at 2
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
        !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then

            ! Compute the grid velocity if present.
            ! It is taken as the average of j and j-1,

            voli = half / cudaDoms(dom,sps)%vol(i, j, k)
            qs = (cudaDoms(dom,sps)%sFaceJ(i, j, k) + cudaDoms(dom,sps)%sFaceJ(i, j - 1, k)) * voli

            ! Compute the normal velocity, where the normal direction
            ! is taken as the average of faces j and j-1.

            xa = (cudaDoms(dom,sps)%sJ(i, j, k, 1) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 1)) * voli
            ya = (cudaDoms(dom,sps)%sJ(i, j, k, 2) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 2)) * voli
            za = (cudaDoms(dom,sps)%sJ(i, j, k, 3) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 3)) * voli

            uu = xa * cudaDoms(dom,sps)%w(i, j, k, ivx) + ya * cudaDoms(dom,sps)%w(i, j, k, ivy) + za * cudaDoms(dom,sps)%w(i, j, k, ivz) - qs

            ! Determine the situation we are having here, i.e. positive
            ! or negative normal velocity.

            velJdir: if (uu > zero) then

                ! Velocity has a component in positive j-direction.
                ! Loop over the number of advection equations.
                do ii = 1, nAdv

                    ! Set the value of jj such that it corresponds to the
                    ! turbulent entry in w.

                    jj = ii + offset

                    ! Check whether a first or a second order discretization
                    ! must be used.

                    if (secondOrd) then

                        ! Second order; store the three differences for the
                        ! discretization of the derivative in j-direction.

                        dwtm1 = cudaDoms(dom,sps)%w(i, j - 1, k, jj) - cudaDoms(dom,sps)%w(i, j - 2, k, jj)
                        dwt = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i, j - 1, k, jj)
                        dwtp1 = cudaDoms(dom,sps)%w(i, j + 1, k, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)

                        ! Construct the derivative in this cell center. This is
                        ! the first order upwind derivative with two nonlinear
                        ! corrections.

                        dwtj = dwt

                        if (dwt * dwtp1 > zero) then
                            if (abs(dwt) < abs(dwtp1)) then
                                dwtj = dwtj + half * dwt
                            else
                                dwtj = dwtj + half * dwtp1
                            end if
                        end if

                        if (dwt * dwtm1 > zero) then
                            if (abs(dwt) < abs(dwtm1)) then
                                dwtj = dwtj - half * dwt
                            else
                                dwtj = dwtj - half * dwtm1
                            end if
                        end if

                    else

                        ! 1st order upwind scheme.

                        dwtj = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i, j - 1, k, jj)

                    end if

                    ! Update the residual. The convective term must be
                    ! substracted, because it appears on the other side of
                    ! the equation as the source and viscous terms.

                    cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) = cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) - uu * dwtj
                end do

            else velJdir

                ! Velocity has a component in negative j-direction.
                ! Loop over the number of advection equations.
                do ii = 1, nAdv

                    ! Set the value of jj such that it corresponds to the
                    ! turbulent entry in w.

                    jj = ii + offset

                    ! Check whether a first or a second order discretization
                    ! must be used.

                    if (secondOrd) then

                        ! Store the three differences for the discretization of
                        ! the derivative in j-direction.

                        dwtm1 = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i, j - 1, k, jj)
                        dwt = cudaDoms(dom,sps)%w(i, j + 1, k, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)
                        dwtp1 = cudaDoms(dom,sps)%w(i, j + 2, k, jj) - cudaDoms(dom,sps)%w(i, j + 1, k, jj)

                        ! Construct the derivative in this cell center. This is
                        ! the first order upwind derivative with two nonlinear
                        ! corrections.

                        dwtj = dwt

                        if (dwt * dwtp1 > zero) then
                            if (abs(dwt) < abs(dwtp1)) then
                                dwtj = dwtj - half * dwt
                            else
                                dwtj = dwtj - half * dwtp1
                            end if
                        end if

                        if (dwt * dwtm1 > zero) then
                            if (abs(dwt) < abs(dwtm1)) then
                                dwtj = dwtj + half * dwt
                            else
                                dwtj = dwtj + half * dwtm1
                            end if
                        end if

                    else

                        ! 1st order upwind scheme.

                        dwtj = cudaDoms(dom,sps)%w(i, j + 1, k, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)

                    end if

                    ! Update the residual. The convective term must be
                    ! substracted, because it appears on the other side
                    ! of the equation as the source and viscous terms.

                    cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) = cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) - uu * dwtj
                end do
            end if velJdir

        end if
        !
        !       Upwind discretization of the convective term in i (xi)
        !       direction. Either the 1st order upwind or the second order
        !       fully upwind interpolation scheme, kappa = -1, is used in
        !       combination with the minmod limiter.
        !       The possible grid velocity must be taken into account.
        !
        qs = zero

        !loop starts at 2
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
        !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then
        
            ! Compute the grid velocity if present.
            ! It is taken as the average of i and i-1,

            voli = half / cudaDoms(dom,sps)%vol(i, j, k)
            qs = (cudaDoms(dom,sps)%sFaceI(i, j, k) + cudaDoms(dom,sps)%sFaceI(i - 1, j, k)) * voli

            ! Compute the normal velocity, where the normal direction
            ! is taken as the average of faces i and i-1.

            xa = (cudaDoms(dom,sps)%sI(i, j, k, 1) + cudaDoms(dom,sps)%sI(i - 1, j, k, 1)) * voli
            ya = (cudaDoms(dom,sps)%sI(i, j, k, 2) + cudaDoms(dom,sps)%sI(i - 1, j, k, 2)) * voli
            za = (cudaDoms(dom,sps)%sI(i, j, k, 3) + cudaDoms(dom,sps)%sI(i - 1, j, k, 3)) * voli

            uu = xa * cudaDoms(dom,sps)%w(i, j, k, ivx) + ya * cudaDoms(dom,sps)%w(i, j, k, ivy) + za * cudaDoms(dom,sps)%w(i, j, k, ivz) - qs

            ! Determine the situation we are having here, i.e. positive
            ! or negative normal velocity.

            velIdir: if (uu > zero) then

                ! Velocity has a component in positive i-direction.
                ! Loop over the number of advection equations.
                do ii = 1, nAdv

                    ! Set the value of jj such that it corresponds to the
                    ! turbulent entry in w.

                    jj = ii + offset

                    ! Check whether a first or a second order discretization
                    ! must be used.

                    if (secondOrd) then

                        ! Second order; store the three differences for the
                        ! discretization of the derivative in i-direction.

                        dwtm1 = cudaDoms(dom,sps)%w(i - 1, j, k, jj) - cudaDoms(dom,sps)%w(i - 2, j, k, jj)
                        dwt = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i - 1, j, k, jj)
                        dwtp1 = cudaDoms(dom,sps)%w(i + 1, j, k, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)

                        ! Construct the derivative in this cell center. This is
                        ! the first order upwind derivative with two nonlinear
                        ! corrections.

                        dwti = dwt

                        if (dwt * dwtp1 > zero) then
                            if (abs(dwt) < abs(dwtp1)) then
                                dwti = dwti + half * dwt
                            else
                                dwti = dwti + half * dwtp1
                            end if
                        end if

                        if (dwt * dwtm1 > zero) then
                            if (abs(dwt) < abs(dwtm1)) then
                                dwti = dwti - half * dwt
                            else
                                dwti = dwti - half * dwtm1
                            end if
                        end if

                    else

                        ! 1st order upwind scheme.

                        dwti = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i - 1, j, k, jj)

                    end if

                    ! Update the residual. The convective term must be
                    ! substracted, because it appears on the other side of
                    ! the equation as the source and viscous terms.

                    cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) = cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) - uu * dwti
                end do

            else velIdir

                ! Velocity has a component in negative i-direction.
                ! Loop over the number of advection equations.
                do ii = 1, nAdv

                    ! Set the value of jj such that it corresponds to the
                    ! turbulent entry in w.

                    jj = ii + offset

                    ! Check whether a first or a second order discretization
                    ! must be used.

                    if (secondOrd) then

                        ! Second order; store the three differences for the
                        ! discretization of the derivative in i-direction.

                        dwtm1 = cudaDoms(dom,sps)%w(i, j, k, jj) - cudaDoms(dom,sps)%w(i - 1, j, k, jj)
                        dwt = cudaDoms(dom,sps)%w(i + 1, j, k, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)
                        dwtp1 = cudaDoms(dom,sps)%w(i + 2, j, k, jj) - cudaDoms(dom,sps)%w(i + 1, j, k, jj)

                        ! Construct the derivative in this cell center. This is
                        ! the first order upwind derivative with two nonlinear
                        ! corrections.

                        dwti = dwt

                        if (dwt * dwtp1 > zero) then
                            if (abs(dwt) < abs(dwtp1)) then
                                dwti = dwti - half * dwt
                            else
                                dwti = dwti - half * dwtp1
                            end if
                        end if

                        if (dwt * dwtm1 > zero) then
                            if (abs(dwt) < abs(dwtm1)) then
                                dwti = dwti + half * dwt
                            else
                                dwti = dwti + half * dwtm1
                            end if
                        end if

                    else

                        ! 1st order upwind scheme.

                        dwti = cudaDoms(dom,sps)%w(i + 1, j, k, jj) - cudaDoms(dom,sps)%w(i, j, k, jj)

                    end if

                    ! Update the residual. The convective term must be
                    ! substracted, because it appears on the other side
                    ! of the equation as the source and viscous terms.

                    cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) = cudaDoms(dom,sps)%dw(i, j, k, itu1 + ii - 1) - uu * dwti

                    ! Update the central jacobian. First the term which is
                    ! always present, i.e. -uu.
                end do

            end if velIdir

        end if

    end subroutine saAdvection

    attributes(global) subroutine saViscous
        use constants, only: one, itu1, irho
        ! use sa, only: kar2Inv, cw36, cb3Inv
        use paramTurb, only: rsaCv1, rsaCw3, rsaCb3, rsaK, rsacb2
        implicit none

        ! Variables for sa Viscous
        real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
        real(kind=realType) :: xa, ya, za, ttm, ttp, cnud, cam, cap
        real(kind=realType) :: nutm, nutp, num, nup, cdm, cdp
        real(kind=realType) :: c1m, c1p, c10, b1, c1, d1, qs, nu
        integer(Kind=intType) :: i, j, k, dom, sps

        !this is new no longer read these from sa module
        !since they are just computed here anyways
        real(kind=realType) :: cv13, kar2Inv,cw36, cb3Inv


        dom = 1 
        sps = 1

        ! Set model constants
        cv13 = rsaCv1**3
        kar2Inv = one / (rsaK**2)
        cw36 = rsaCw3**6
        cb3Inv = one / rsaCb3
       
        !loops start at 2
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1

        !loop k 2 to cudaDoms(dom,sps)%kl, j 2 to cudaDoms(dom,sps)%jl, i 2 to cudaDoms(dom,sps)%il
        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then

            !
            !       Viscous terms in k-direction.
            !

            ! Compute the metrics in zeta-direction, i.e. along the
            ! line k = constant.

            voli = one / cudaDoms(dom,sps)%vol(i, j, k)
            volmi = two / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i, j, k - 1))
            volpi = two / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i, j, k + 1))

            xm = cudaDoms(dom,sps)%sK(i, j, k - 1, 1) * volmi
            ym = cudaDoms(dom,sps)%sK(i, j, k - 1, 2) * volmi
            zm = cudaDoms(dom,sps)%sK(i, j, k - 1, 3) * volmi
            xp = cudaDoms(dom,sps)%sK(i, j, k, 1) * volpi
            yp = cudaDoms(dom,sps)%sK(i, j, k, 2) * volpi
            zp = cudaDoms(dom,sps)%sK(i, j, k, 3) * volpi

            xa = half * (cudaDoms(dom,sps)%sK(i, j, k, 1) + cudaDoms(dom,sps)%sK(i, j, k - 1, 1)) * voli
            ya = half * (cudaDoms(dom,sps)%sK(i, j, k, 2) + cudaDoms(dom,sps)%sK(i, j, k - 1, 2)) * voli
            za = half * (cudaDoms(dom,sps)%sK(i, j, k, 3) + cudaDoms(dom,sps)%sK(i, j, k - 1, 3)) * voli
            ttm = xm * xa + ym * ya + zm * za
            ttp = xp * xa + yp * ya + zp * za

            ! Computation of the viscous terms in zeta-direction; note
            ! that cross-derivatives are neglected, i.e. the mesh is
            ! assumed to be orthogonal.
            ! Furthermore, the grad(nu)**2 has been rewritten as
            ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
            ! The second derivative in zeta-direction is constructed as
            ! the central difference of the first order derivatives, i.e.
            ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
            ! In this way the metric can be taken into account.

            ! Compute the diffusion coefficients multiplying the nodes
            ! k+1, k and k-1 in the second derivative. Make sure that
            ! these coefficients are nonnegative.

            cnud = -rsaCb2 * cudaDoms(dom,sps)%w(i, j, k, itu1) * cb3Inv
            cam = ttm * cnud
            cap = ttp * cnud

            nutm = half * (cudaDoms(dom,sps)%w(i, j, k - 1, itu1) + cudaDoms(dom,sps)%w(i, j, k, itu1))
            nutp = half * (cudaDoms(dom,sps)%w(i, j, k + 1, itu1) + cudaDoms(dom,sps)%w(i, j, k, itu1))
            nu = cudaDoms(dom,sps)%rlv(i, j, k) / cudaDoms(dom,sps)%w(i, j, k, irho)
            num = half * (cudaDoms(dom,sps)%rlv(i, j, k - 1) / cudaDoms(dom,sps)%w(i, j, k - 1, irho) + nu)
            nup = half * (cudaDoms(dom,sps)%rlv(i, j, k + 1) / cudaDoms(dom,sps)%w(i, j, k + 1, irho) + nu)
            cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
            cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

            c1m = max(cdm + cam, zero)
            c1p = max(cdp + cap, zero)
            c10 = c1m + c1p

            ! Update the residual for this cell and store the possible
            ! coefficients for the matrix in b1, c1 and d1.

            cudaDoms(dom,sps)%dw(i, j, k, itu1) = cudaDoms(dom,sps)%dw(i, j, k, itu1) + c1m * cudaDoms(dom,sps)%w(i, j, k - 1, itu1) &
                                - c10 * cudaDoms(dom,sps)%w(i, j, k, itu1) + c1p * cudaDoms(dom,sps)%w(i, j, k + 1, itu1)

            !
            !       Viscous terms in j-direction.
            !

            ! Compute the metrics in eta-direction, i.e. along the
            ! line j = constant.

            voli = one / cudaDoms(dom,sps)%vol(i, j, k)
            volmi = two / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i, j - 1, k))
            volpi = two / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i, j + 1, k))

            xm = cudaDoms(dom,sps)%sJ(i, j - 1, k, 1) * volmi
            ym = cudaDoms(dom,sps)%sJ(i, j - 1, k, 2) * volmi
            zm = cudaDoms(dom,sps)%sJ(i, j - 1, k, 3) * volmi
            xp = cudaDoms(dom,sps)%sJ(i, j, k, 1) * volpi
            yp = cudaDoms(dom,sps)%sJ(i, j, k, 2) * volpi
            zp = cudaDoms(dom,sps)%sJ(i, j, k, 3) * volpi

            xa = half * (cudaDoms(dom,sps)%sJ(i, j, k, 1) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 1)) * voli
            ya = half * (cudaDoms(dom,sps)%sJ(i, j, k, 2) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 2)) * voli
            za = half * (cudaDoms(dom,sps)%sJ(i, j, k, 3) + cudaDoms(dom,sps)%sJ(i, j - 1, k, 3)) * voli
            ttm = xm * xa + ym * ya + zm * za
            ttp = xp * xa + yp * ya + zp * za

            ! Computation of the viscous terms in eta-direction; note
            ! that cross-derivatives are neglected, i.e. the mesh is
            ! assumed to be orthogonal.
            ! Furthermore, the grad(nu)**2 has been rewritten as
            ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
            ! The second derivative in eta-direction is constructed as
            ! the central difference of the first order derivatives, i.e.
            ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
            ! In this way the metric can be taken into account.

            ! Compute the diffusion coefficients multiplying the nodes
            ! j+1, j and j-1 in the second derivative. Make sure that
            ! these coefficients are nonnegative.

            cnud = -rsaCb2 * cudaDoms(dom,sps)%w(i, j, k, itu1) * cb3Inv
            cam = ttm * cnud
            cap = ttp * cnud

            nutm = half * (cudaDoms(dom,sps)%w(i, j - 1, k, itu1) + cudaDoms(dom,sps)%w(i, j, k, itu1))
            nutp = half * (cudaDoms(dom,sps)%w(i, j + 1, k, itu1) + cudaDoms(dom,sps)%w(i, j, k, itu1))
            nu = cudaDoms(dom,sps)%rlv(i, j, k) / cudaDoms(dom,sps)%w(i, j, k, irho)
            num = half * (cudaDoms(dom,sps)%rlv(i, j - 1, k) / cudaDoms(dom,sps)%w(i, j - 1, k, irho) + nu)
            nup = half * (cudaDoms(dom,sps)%rlv(i, j + 1, k) / cudaDoms(dom,sps)%w(i, j + 1, k, irho) + nu)
            cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
            cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

            c1m = max(cdm + cam, zero)
            c1p = max(cdp + cap, zero)
            c10 = c1m + c1p

            ! Update the residual for this cell and store the possible
            ! coefficients for the matrix in b1, c1 and d1.

            cudaDoms(dom,sps)%dw(i, j, k, itu1) = cudaDoms(dom,sps)%dw(i, j, k, itu1) + c1m * cudaDoms(dom,sps)%w(i, j - 1, k, itu1) &
                                - c10 * cudaDoms(dom,sps)%w(i, j, k, itu1) + c1p * cudaDoms(dom,sps)%w(i, j + 1, k, itu1)

            !
            !       Viscous terms in i-direction.
            !
            
            ! Compute the metrics in xi-direction, i.e. along the
            ! line i = constant.

            voli = one / cudaDoms(dom,sps)%vol(i, j, k)
            volmi = two / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i - 1, j, k))
            volpi = two / (cudaDoms(dom,sps)%vol(i, j, k) + cudaDoms(dom,sps)%vol(i + 1, j, k))

            xm = cudaDoms(dom,sps)%sI(i - 1, j, k, 1) * volmi
            ym = cudaDoms(dom,sps)%sI(i - 1, j, k, 2) * volmi
            zm = cudaDoms(dom,sps)%sI(i - 1, j, k, 3) * volmi
            xp = cudaDoms(dom,sps)%sI(i, j, k, 1) * volpi
            yp = cudaDoms(dom,sps)%sI(i, j, k, 2) * volpi
            zp = cudaDoms(dom,sps)%sI(i, j, k, 3) * volpi

            xa = half * (cudaDoms(dom,sps)%sI(i, j, k, 1) + cudaDoms(dom,sps)%sI(i - 1, j, k, 1)) * voli
            ya = half * (cudaDoms(dom,sps)%sI(i, j, k, 2) + cudaDoms(dom,sps)%sI(i - 1, j, k, 2)) * voli
            za = half * (cudaDoms(dom,sps)%sI(i, j, k, 3) + cudaDoms(dom,sps)%sI(i - 1, j, k, 3)) * voli
            ttm = xm * xa + ym * ya + zm * za
            ttp = xp * xa + yp * ya + zp * za

            ! Computation of the viscous terms in xi-direction; note
            ! that cross-derivatives are neglected, i.e. the mesh is
            ! assumed to be orthogonal.
            ! Furthermore, the grad(nu)**2 has been rewritten as
            ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
            ! The second derivative in xi-direction is constructed as
            ! the central difference of the first order derivatives, i.e.
            ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
            ! In this way the metric can be taken into account.

            ! Compute the diffusion coefficients multiplying the nodes
            ! i+1, i and i-1 in the second derivative. Make sure that
            ! these coefficients are nonnegative.

            cnud = -rsaCb2 * cudaDoms(dom,sps)%w(i, j, k, itu1) * cb3Inv
            cam = ttm * cnud
            cap = ttp * cnud

            nutm = half * (cudaDoms(dom,sps)%w(i - 1, j, k, itu1) + cudaDoms(dom,sps)%w(i, j, k, itu1))
            nutp = half * (cudaDoms(dom,sps)%w(i + 1, j, k, itu1) + cudaDoms(dom,sps)%w(i, j, k, itu1))
            nu = cudaDoms(dom,sps)%rlv(i, j, k) / cudaDoms(dom,sps)%w(i, j, k, irho)
            num = half * (cudaDoms(dom,sps)%rlv(i - 1, j, k) / cudaDoms(dom,sps)%w(i - 1, j, k, irho) + nu)
            nup = half * (cudaDoms(dom,sps)%rlv(i + 1, j, k) / cudaDoms(dom,sps)%w(i + 1, j, k, irho) + nu)
            cdm = (num + (one + rsaCb2) * nutm) * ttm * cb3Inv
            cdp = (nup + (one + rsaCb2) * nutp) * ttp * cb3Inv

            c1m = max(cdm + cam, zero)
            c1p = max(cdp + cap, zero)
            c10 = c1m + c1p

            ! Update the residual for this cell and store the possible
            ! coefficients for the matrix in b1, c1 and d1.

            cudaDoms(dom,sps)%dw(i, j, k, itu1) = cudaDoms(dom,sps)%dw(i, j, k, itu1) + c1m * cudaDoms(dom,sps)%w(i - 1, j, k, itu1) &
                                - c10 * cudaDoms(dom,sps)%w(i, j, k, itu1) + c1p * cudaDoms(dom,sps)%w(i + 1, j, k, itu1)
        end if
        !TODO finish copying i and j directions
    end subroutine saViscous


    attributes(global) subroutine saResScale
        !
        !  Multiply the residual by the volume and store this in cudaDoms(dom,sps)%dw; this
        ! * is done for monitoring reasons only. The multiplication with the
        ! * volume is present to be consistent with the flow residuals; also
        !  the negative value is taken, again to be consistent with the
        ! * flow equations. Also multiply by iblank so that no updates occur
        !  in holes or the overset boundary.
        use constants, only: zero
        implicit none

        ! Local variables
        integer(kind=intType) :: i, j, k, dom, sps
        real(kind=realType) :: rblank

        dom = 1 
        sps = 1
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
        j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1

        if (i <= cudaDoms(dom,sps)%il .and. j <= cudaDoms(dom,sps)%jl .and. k <= cudaDoms(dom,sps)%kl) then
            rblank = max(real(cudaDoms(dom,sps)%iblank(i, j, k), realType), zero)
            cudaDoms(dom,sps)%dw(i, j, k, itu1) = -cudaDoms(dom,sps)%volRef(i, j, k) * cudaDoms(dom,sps)%dw(i, j, k, itu1) * rblank
        end if    
    end subroutine saResScale


    subroutine calculateCudaResidual(updateIntermed, varStart, varEnd)
      use constants, only: zero, spalartAllmaras
      use cudafor,   only: dim3
      use precision, only: intType
      use cudaFlowVarRefState, only: cudaCopyFlowVarRefState
      use cudaInputDiscretization, only: cudaCopyInputDiscretization
      use cudaInputIteration,      only: cudaCopyInputIteration
      use cudaInputPhysics,        only: cudaCopyInputPhysics
      use cudaParamTurb,           only: cudaCopyParamTurb
      use cudaIteration,           only: cudaCopyIteration
      use cudaTurbMod,             only: cudaCopyTurbMod
      use inputPhysics, only: equationMode, equations, turbModel
      use inputDiscretization, only: spaceDiscr
      use flowvarRefState, only: viscous
      use blockPointers, only: nDom, bib=>ib, bjb=>jb, bkb=>kb
      use inputTimeSpectral, only: nTimeIntervalsSpectral
      use cudaBlock, only: copyCudaBlock
      use utils, only: setPointers



      implicit none

      logical, intent(in)               :: updateIntermed
      integer(kind=intType), intent(in) :: varStart, varEnd
      type(dim3)                        :: grid_size, block_size
      integer(kind=intType) :: istat
      real(kind=realType) :: start, finish, startInv, finishInv, startVisc,finishVisc,startSA,finishSA
      integer(kind = intType) :: dom, sps
      integer(kind = intType) :: ibmax, jbmax, kbmax
    
      ! copy constants
      call cudaCopyFlowVarRefState
      call cudaCopyInputDiscretization
      call cudaCopyInputIteration
      call cudaCopyInputPhysics
      call cudaCopyParamTurb
      call cudaCopyIteration
      call cudaCopyTurbMod

      !copy data from cpu to gpu
      call copyCudaBlock
      call copyData

        ibmax = 0
        jbmax = 0
        kbmax = 0
    ! zeroing
      do dom = 1,nDom
        do sps = 1, ntimeIntervalsSpectral
            call setPointers(dom,1,sps)
            ibmax = max(ibmax,bib)
            jbmax = max(jbmax,bjb)
            kbmax = max(kbmax,bkb)
            h_cudaDoms(dom,sps)%fw                           = zero
            h_cudaDoms(dom,sps)%dw(:, :, :, varStart:varEnd) = zero
        end do
      end do 
     
      istat = cudaDeviceSynchronize()

      ! metrics
      print *,ibmax,jbmax,kbmax
      block_size = dim3(8, 8, 8)
      grid_size  = dim3(ceiling(real(ibmax+2) / block_size%x), ceiling(real(jbmax+2) / block_size%y), ceiling(real(kbmax+2) / block_size%z))
      

      call CPU_TIME(start)
    !   call metrics<<<grid_size, block_size>>>
    !   istat = cudaDeviceSynchronize()

      ! call SA routines
      if (equations == RANSEquations .and. turbModel == spalartAllmaras) then
        call CPU_TIME(startSA)

        call saSource<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()

        call saAdvection<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()

        call saViscous<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()

        call saResScale<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()
        call CPU_TIME(finishSA)
                print  '("CUDA SA = ",E22.16," seconds.")', finishSA-startSA


      end if
      
      ! call time step routine
        call CPU_TIME(startInv)
      call timeStep<<<grid_size, block_size>>>(updateIntermed)
      istat = cudaDeviceSynchronize()
      
      ! inviscid central flux
      call computeSS<<<grid_size, block_size>>>
      istat = cudaDeviceSynchronize()

      call inviscidCentralFlux<<<grid_size, block_size>>>
      istat = cudaDeviceSynchronize()

      ! inviscid diss flux scalar
      call inviscidDissFluxScalar<<<grid_size, block_size>>>
      istat = cudaDeviceSynchronize()
      call CPU_TIME(finishInv)

        print  '("CUDA Inviscid Time = ",E22.16," seconds.")', finishInv-startInv

      ! viscousFlux
      if (viscous) then
      call CPU_TIME(startVisc)
        call computeSpeedOfSoundSquared_v1<<<grid_size, block_size>>>
        call computeSpeedOfSoundSquared_v2<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()

        call allNodalGradients_v1<<<grid_size, block_size>>>
        call allNodalGradients_v2<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()

        call scaleNodalGradients_v1<<<grid_size, block_size>>>
        call scaleNodalGradients_v2<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()
        block_size = dim3(4, 4, 2)
        grid_size  = dim3(ceiling(real(ibmax+2) / block_size%x), ceiling(real(jbmax+2) / block_size%y), ceiling(real(kbmax+2) / block_size%z))
        call viscousFlux<<<grid_size, block_size>>>
        istat = cudaDeviceSynchronize()
        call CPU_TIME(finishVisc)
        print  '("CUDA Visc = ",E22.16," seconds.")', finishVisc-startVisc

      end if 
      

    !   sumDwAndFw
    !   call sumDwandFw_v1<<<grid_size, block_size>>>
      call sumDwandFw_v2<<<grid_size, block_size>>>
      istat = cudaDeviceSynchronize()
      call CPU_TIME(finish)


        print  '("CUDA Time = ",E22.16," seconds.")', finish-start
      ! TODO: update residual?

    end subroutine calculateCudaResidual

    subroutine cudaResAPI(res,ndimw)

        use constants
        use blockPointers, only: bvolRef=>volRef,bie=>ie, bje=>je, bke=>ke, &
                                    bil=>il, bjl=>jl, bkl=>kl, &
                                
                                    bux=>ux, buy=>uy, buz=>uz, &
                                    bvx=>vx, bvy=>vy, bvz=>vz, &
                                    bwx=>wx, bwy=>wy, bwz=>wz, &
                                    bqx=>qx, bqy=>qy, bqz=>qz, &
                                    bfw=>fw, bw=>w, &
                                    bib=>ib, bjb=>jb, bkb=>kb, &
                                    bsI=>sI, bsJ=>sJ, bsK=>sK, &
                                    bvol=>vol, baa=>aa, bp=>p, bgamma=>gamma

        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowvarrefstate, only: nw,nwf,nt1,nt2
        use utils, only: setPointers
        use flowutils, only: ballNodalGradients=>allNodalGradients, bcomputeSpeedOfSoundSquared=>computeSpeedOfSoundSquared

        implicit none
        integer(kind=intType), intent(in) :: ndimw
        real(kind=realType), dimension(ndimw), intent(inout) :: res(ndimw)
        
        ! Local Variables
        integer(kind=intType) :: nn, i, j, k, l, counter, sps, dom
        real(kind=realType) :: ovv
        real(kind=realType), dimension(:,:,:,:), allocatable :: h_dw, h_dw2
        ! These are v2 versions of the variables
        real(kind=realType), dimension(:, :, :), allocatable :: ux_v2, uy_v2, uz_v2
        real(kind=realType), dimension(:, :, :), allocatable :: vx_v2, vy_v2, vz_v2
        real(kind=realType), dimension(:, :, :), allocatable :: wx_v2, wy_v2, wz_v2
        real(kind=realType), dimension(:, :, :), allocatable :: qx_v2, qy_v2, qz_v2

        !versions from old working version (v1)
        real(kind=realType), dimension(:, :, :), allocatable :: ux_v1, uy_v1, uz_v1
        real(kind=realType), dimension(:, :, :), allocatable :: vx_v1, vy_v1, vz_v1
        real(kind=realType), dimension(:, :, :), allocatable :: wx_v1, wy_v1, wz_v1
        real(kind=realType), dimension(:, :, :), allocatable :: qx_v1, qy_v1, qz_v1

        real(kind=realType), dimension(:,:,:,:), allocatable :: h_w
        real(kind=realType), dimension(:,:,:,:), allocatable :: h_sI, h_sJ, h_sK
        real(kind=realType), dimension(:,:,:), allocatable :: h_vol, h_aa,h_p,h_gamma, h_aa_v1,h_vol_v1
        integer(kind=intType) :: h_ie,h_je,h_ke

        dom = 1 
        sps = 1
        call setPointers(1, 1, 1)
        call calculateCudaResidual(.True.,1,nw)

        

        allocate(ux_v1(1:bil,1:bjl,1:bkl))
        allocate(uy_v1(1:bil,1:bjl,1:bkl))
        allocate(uz_v1(1:bil,1:bjl,1:bkl))
        allocate(vx_v1(1:bil,1:bjl,1:bkl))
        allocate(vy_v1(1:bil,1:bjl,1:bkl))
        allocate(vz_v1(1:bil,1:bjl,1:bkl))
        allocate(wx_v1(1:bil,1:bjl,1:bkl))
        allocate(wy_v1(1:bil,1:bjl,1:bkl))
        allocate(wz_v1(1:bil,1:bjl,1:bkl))
        allocate(qx_v1(1:bil,1:bjl,1:bkl))
        allocate(qy_v1(1:bil,1:bjl,1:bkl))
        allocate(qz_v1(1:bil,1:bjl,1:bkl))

        allocate(ux_v2(1:bil,1:bjl,1:bkl))
        allocate(uy_v2(1:bil,1:bjl,1:bkl))
        allocate(uz_v2(1:bil,1:bjl,1:bkl))
        allocate(vx_v2(1:bil,1:bjl,1:bkl))
        allocate(vy_v2(1:bil,1:bjl,1:bkl))
        allocate(vz_v2(1:bil,1:bjl,1:bkl))
        allocate(wx_v2(1:bil,1:bjl,1:bkl))
        allocate(wy_v2(1:bil,1:bjl,1:bkl))
        allocate(wz_v2(1:bil,1:bjl,1:bkl))
        allocate(qx_v2(1:bil,1:bjl,1:bkl))
        allocate(qy_v2(1:bil,1:bjl,1:bkl))
        allocate(qz_v2(1:bil,1:bjl,1:bkl))

        allocate(h_sI(0:bie,1:bje,1:bke,3))
        allocate(h_sJ(1:bie,0:bje,1:bke,3))
        allocate(h_sK(1:bie,1:bje,0:bke,3))
        allocate(h_vol(1:bie,1:bje,1:bke))

        allocate(h_vol_v1(1:bie,1:bje,1:bke))

        allocate(h_aa(1:bie,1:bje,1:bke))
        allocate(h_aa_v1(1:bie,1:bje,1:bke))

        allocate(h_p(0:bib,0:bjb,0:bkb))
        allocate(h_gamma(0:bib,0:bjb,0:bkb))
        

        ! allocate memory for cudaDoms(dom,sps)%dw
        allocate(h_dw(1:bie,1:bje,1:bke,1:nw))
        allocate(h_dw2(1:bie,1:bje,1:bke,1:nw))

        allocate(h_w(0:bib,0:bjb,0:bkb,1:nw))

        
        !copy from gpu to cpu
        h_cudaDoms(1,1) = cudaDoms(1,1)
        h_dw(1:bie,1:bje,1:bke,1:nw) = h_cudaDoms(dom,sps)%dw(1:bie,1:bje,1:bke,1:nw)
        ! h_dw2 = dw

        ! h_dw = h_dw + h_dw2
        print *, bie,bje,bke
        h_ie = h_cudaDoms(1,1)%ie
        h_je =  h_cudaDoms(1,1)%je
        h_ke = h_cudaDoms(1,1)%ke
        print *, h_ie,h_je,h_ke
        baa = zero
        call bcomputeSpeedOfSoundSquared
        call ballNodalGradients


        ux_v2 = h_cudaDoms(dom,sps)%ux
        uy_v2 = h_cudaDoms(dom,sps)%uy
        uz_v2 = h_cudaDoms(dom,sps)%uz
        
        ux_v1 = ux
        uy_v1 = uy
        uz_v1 = uz
        
        print *,ux_v1(:,:,1)
        print *,"=================="
        print *, ux_v2(:,:,1)
        print *,"=================="
        print *, ux_v1(:,:,1)-ux_v2(:,:,1)

        print  '("Max DIFF UX = ",E22.16," seconds.")', maxval(abs(ux_v2-ux_v1))
        print  '("Max DIFF UY = ",E22.16," seconds.")', maxval(abs(uy_v2-uy_v1))
        print  '("Max DIFF UZ = ",E22.16," seconds.")', maxval(abs(uz_v2-uz_v1))

        vx_v2 = h_cudaDoms(dom,sps)%vx
        vy_v2 = h_cudaDoms(dom,sps)%vy
        vz_v2 = h_cudaDoms(dom,sps)%vz
        vx_v1 = vx
        vy_v1 = vy
        vz_v1 = vz
        print  '("Max DIFF VX = ",E22.16," seconds.")', maxval(abs(vx_v2-vx_v1))
        print  '("Max DIFF Vy = ",E22.16," seconds.")', maxval(abs(vy_v2-vy_v1))
        print  '("Max DIFF vz = ",E22.16," seconds.")', maxval(abs(vz_v2-vz_v1))

        wx_v2 = h_cudaDoms(dom,sps)%wx
        wy_v2 = h_cudaDoms(dom,sps)%wy
        wz_v2 = h_cudaDoms(dom,sps)%wz
        wx_v1 = wx
        wy_v1 = wy
        wz_v1 = wz
        print  '("Max DIFF wX = ",E22.16," seconds.")', maxval(abs(wx_v2-wx_v1))
        print  '("Max DIFF wy = ",E22.16," seconds.")', maxval(abs(wy_v2-wy_v1))
        print  '("Max DIFF wz = ",E22.16," seconds.")', maxval(abs(wz_v2-wz_v1))

        qx_v2 = h_cudaDoms(dom,sps)%qx
        qy_v2 = h_cudaDoms(dom,sps)%qy
        qz_v2 = h_cudaDoms(dom,sps)%qz
        qx_v1 = qx
        qy_v1 = qy
        qz_v1 = qz
        print  '("Max DIFF qX = ",E22.16," seconds.")', maxval(abs(qx_v2-qx_v1))
        print  '("Max DIFF qy = ",E22.16," seconds.")', maxval(abs(qy_v2-qy_v1))
        print  '("Max DIFF qz = ",E22.16," seconds.")', maxval(abs(qz_v2-qz_v1))


        h_w = h_cudaDoms(dom,sps)%w
        print  '("Max DIFF w = ",E22.16," seconds.")', maxval(abs(h_w-bw))


        ! Surface normals are good
        h_sI = h_cudaDoms(dom,sps)%sI
        h_sJ = h_cudaDoms(dom,sps)%sJ
        h_sK = h_cudaDoms(dom,sps)%sK
        print  '("Max DIFF sI = ",E22.16," seconds.")', maxval(abs(h_SI-bsI))
        print  '("Max DIFF sJ = ",E22.16," seconds.")', maxval(abs(h_SJ-bsJ))
        print  '("Max DIFF sK = ",E22.16," seconds.")', maxval(abs(h_SK-bsK))

        h_vol = h_cudaDoms(dom,sps)%vol
        print  '("Max DIFF vol = ",E22.16," seconds.")', maxval(abs(h_vol-bvol))

        ! print *, "volumes:"
        ! print *,"========"
        ! print *, h_vol
        ! print *, "========"
        ! print *, bvol
        ! print *,"========"
        ! print *, size(h_vol), size(bvol)
        ! counter = 0
        
        ! Speed of sound 
        h_aa = h_cudaDoms(dom,sps)%aa
        h_aa_v1 = aa
        print  '("Max DIFF aa = ",E22.16," seconds.")', maxval(abs(h_aa-h_aa_v1))
        print  '("Max DIFF aa cpu vs gpu= ",E22.16," seconds.")', maxval(abs(h_aa_v1-baa))

        print *, maxval(h_aa), maxval(baa),minval(h_aa),minval(baa)

        h_p = h_cudaDoms(dom,sps)%p
        print  '("Max DIFF p = ",E22.16," seconds.")', maxval(abs(h_p-bp))

        h_gamma = h_cudaDoms(dom,sps)%gamma
        print  '("Max DIFF gamma = ",E22.16," seconds.")', maxval(abs(h_gamma-bgamma))
        counter = 0
        !copy cudaDoms(dom,sps)%dw to res
        do k = 2,bkl
            do j = 2, bjl
                do i = 2,bil
                    ovv = one / bvolRef(i, j, k)
                    do l = 1, nw
                        counter = counter + 1
                        res(counter) = h_dw(i, j, k, l) * ovv
                    end do
                end do
            end do
        end do

        deallocate(h_dw)
    end subroutine cudaResAPI
    
end module cudaResidual


