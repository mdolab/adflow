! @File    :   cudaBCPointers.f90
! @Desc    :   This module contains the BC data type and provides subroutines for copying data between the CPU and GPU

module cudaBCPointers
    use constants, only: intType, realType
    use cudaFlowVarRefState, only: nw
    implicit none
    save

    real(kind=realType), dimension(:, :, :), pointer, device :: ww0, ww1, ww2, ww3
    real(kind=realType), dimension(:, :), pointer, device :: pp0, pp1, pp2, pp3
    real(kind=realType), dimension(:, :), pointer, device :: rlv0, rlv1, rlv2, rlv3
    real(kind=realType), dimension(:, :), pointer, device :: rev0, rev1, rev2, rev3
    real(kind=realType), dimension(:, :), pointer, device :: gamma0, gamma1, gamma2, gamma3
    
    real(kind=realType), dimension(:, :, :), pointer, device :: ssi, ssj, ssk
    real(kind=realType), dimension(:, :, :), pointer, device :: ss, xx
    real(kind=realType), dimension(:, :), pointer, device :: dd2wall, sFace
    integer(kind=intType), dimension(:, :), pointer, device :: gcp

    integer(kind=intType), pointer, device :: iStart, iEnd, iSize
    integer(kind=intType), pointer, device :: jStart, jEnd, jSize


    type cudaBCType
        ! Data structure similar to BCType
        real(kind=realType), dimension(:, :, :), pointer, device :: ww0, ww1, ww2, ww3
        real(kind=realType), dimension(:, :), pointer, device :: pp0, pp1, pp2, pp3
        real(kind=realType), dimension(:, :), pointer, device :: rlv0, rlv1, rlv2, rlv3
        real(kind=realType), dimension(:, :), pointer, device :: rev0, rev1, rev2, rev3
        real(kind=realType), dimension(:, :), pointer, device :: gamma0, gamma1, gamma2, gamma3
        real(kind=realType), dimension(:, :, :), pointer, device :: ssi, ssj, ssk
        real(kind=realType), dimension(:, :, :), pointer, device :: ss, xx
        real(kind=realType), dimension(:, :), pointer, device :: dd2wall, sFace
        integer(kind=intType), dimension(:, :), pointer, device :: gcp
    
        integer(kind=intType), pointer, device :: iStart, iEnd, iSize
        integer(kind=intType), pointer, device :: jStart, jEnd, jSize

    ! Commenting out differentiated code for now

    ! #ifndef USE_TAPENADE
    !     real(kind=realType), dimension(:, :, :), pointer, device :: ww0d, ww1d, ww2d, ww3d
    !     real(kind=realType), dimension(:, :), pointer, device :: pp0d, pp1d, pp2d, pp3d
    !     real(kind=realType), dimension(:, :), pointer, device :: rlv0d, rlv1d, rlv2d, rlv3d
    !     real(kind=realType), dimension(:, :), pointer, device :: rev0d, rev1d, rev2d, rev3d
    !     real(kind=realType), dimension(:, :), pointer, device :: gamma0d, gamma1d, gamma2d, gamma3d
    !     real(kind=realType), dimension(:, :, :), pointer, device :: ssid, ssjd, sskd, xxd
    !     real(kind=realType), dimension(:, :, :), pointer, device :: ssd
    !     real(kind=realType), dimension(:, :), pointer, device :: dd2walld, sFaced
    ! #endif
    
    end type cudaBCType

    ! Device (GPU) and host (CPU) domains
    type(cudaBCType), device, allocatable,dimension(:,:) :: d_cudaBCDoms(:,:)
    type(cudaBCType), allocatable,dimension(:,:) :: h_cudaBCDoms(:,:)

    contains 
    
    subroutine copyCudaBC
        use blockPointers, only: nDom
        use BCPointers, only: iStart, iEnd, iSize, jStart, jEnd, jSize, &
            bww0=>ww0, bww1=>ww1, bww2=>ww2, bww3=>ww3, &
            bpp0=>pp0, bpp1=>pp1, bpp2=>pp2, bpp3=>pp3, &
            brlv0=>rlv0, brlv1=>rlv1, brlv2=>rlv2, brlv3=>rlv3, &
            brev0=>rev0, brev1=>rev1, brev2=>rev2, brev3=>rev3, &
            bgamma0=>gamma0, bgamma1=>gamma1, bgamma2=>gamma2, bgamma3=>gamma3, &
            bssi=>ssi, bssj=>ssj, bssk=>ssk, bss=>ss, bxx=>xx, bdd2wall=>dd2wall, &
            bsFace=>sFace, bgcp=>gcp

        use utils, only: setBCPointers
        
        implicit none
        
        integer(kind=intType) :: nn, sps    ! domain and time spectral interval indices
        integer(kind=intType) :: nTimeIntervalsSpectral=1

        ! Allocate device 
        allocate(d_cudaBCDoms(nDom, nTimeIntervalsSpectral))
        ! Allocate host
        allocate(h_cudaBCDoms(nDom, nTimeIntervalsSpectral))
        
        do nn = 1, nDom ! loop over domains
            do sps = 1, nTimeIntervalsSpectral ! loop over time spectral intervals
                
                call setBCPointers(nn, .True.)
                
                ! Indices
                h_cudaBCDoms(nn, sps)%iStart = iStart
                h_cudaBCDoms(nn, sps)%iEnd = iEnd
                h_cudaBCDoms(nn, sps)%iSize = iSize
                h_cudaBCDoms(nn, sps)%jStart = jStart
                h_cudaBCDoms(nn, sps)%jEnd = jEnd
                h_cudaBCDoms(nn, sps)%jSize = jSize

                ! Allocate BC data
                allocate(h_cudaBCDoms(nn,sps)%ww0(iSize, jSize, nw)) ! 2nd halo flowfield states
                allocate(h_cudaBCDoms(nn,sps)%ww1(iSize, jSize, nw)) ! 1st halo flowfield states
                allocate(h_cudaBCDoms(nn,sps)%ww2(iSize, jSize, nw)) ! 1st interior flowfield states
                allocate(h_cudaBCDoms(nn,sps)%ww3(iSize, jSize, nw)) ! 2nd interior flowfield states
                allocate(h_cudaBCDoms(nn,sps)%pp0(iSize, jSize)) ! 2nd halo pressure states
                allocate(h_cudaBCDoms(nn,sps)%pp1(iSize, jSize)) ! 1st halo pressure states
                allocate(h_cudaBCDoms(nn,sps)%pp2(iSize, jSize)) ! 1st interior pressure states
                allocate(h_cudaBCDoms(nn,sps)%pp3(iSize, jSize)) ! 2nd interior pressure states
                allocate(h_cudaBCDoms(nn,sps)%rlv0(iSize, jSize)) ! 2nd halo laminar viscosity states
                allocate(h_cudaBCDoms(nn,sps)%rlv1(iSize, jSize)) ! 1st halo laminar viscosity states
                allocate(h_cudaBCDoms(nn,sps)%rlv2(iSize, jSize)) ! 1st interior laminar viscosity states
                allocate(h_cudaBCDoms(nn,sps)%rlv3(iSize, jSize)) ! 2nd interior laminar viscosity states
                allocate(h_cudaBCDoms(nn,sps)%rev0(iSize, jSize)) ! 2nd halo eddy viscosity states
                allocate(h_cudaBCDoms(nn,sps)%rev1(iSize, jSize)) ! 1st halo eddy viscosity states
                allocate(h_cudaBCDoms(nn,sps)%rev2(iSize, jSize)) ! 1st interior eddy viscosity states
                allocate(h_cudaBCDoms(nn,sps)%rev3(iSize, jSize)) ! 2nd interior eddy viscosity states
                allocate(h_cudaBCDoms(nn,sps)%gamma0(iSize, jSize)) ! 2nd halo gamma states
                allocate(h_cudaBCDoms(nn,sps)%gamma1(iSize, jSize)) ! 1st halo gamma states
                allocate(h_cudaBCDoms(nn,sps)%gamma2(iSize, jSize)) ! 1st interior gamma states
                allocate(h_cudaBCDoms(nn,sps)%gamma3(iSize, jSize)) ! 2nd interior gamma states
                allocate(h_cudaBCDoms(nn,sps)%ssi(iSize, jSize, 3)) ! i face normal
                allocate(h_cudaBCDoms(nn,sps)%ssj(iSize, jSize, 3)) ! j face normal
                allocate(h_cudaBCDoms(nn,sps)%ssk(iSize, jSize, 3)) ! k face normal
                allocate(h_cudaBCDoms(nn,sps)%ss(iSize, jSize, 3)) ! mesh velocities (see block.F90) this var naming sucks
                allocate(h_cudaBCDoms(nn,sps)%xx(iSize, jSize, 3)) ! coordinates
                allocate(h_cudaBCDoms(nn,sps)%dd2wall(iSize, jSize)) ! distance squared to wall
                allocate(h_cudaBCDoms(nn,sps)%sFace(iSize, jSize)) ! dot product of face velocity & normal 
                allocate(h_cudaBCDoms(nn,sps)%gcp(iSize, jSize)) ! global cell index

                ! Set the BC data
                h_cudaBCDoms(nn,sps)%ww0 = bww0
                h_cudaBCDoms(nn,sps)%ww1 = bww1
                h_cudaBCDoms(nn,sps)%ww2 = bww2
                h_cudaBCDoms(nn,sps)%ww3 = bww3
                h_cudaBCDoms(nn,sps)%pp0 = bpp0
                h_cudaBCDoms(nn,sps)%pp1 = bpp1
                h_cudaBCDoms(nn,sps)%pp2 = bpp2
                h_cudaBCDoms(nn,sps)%pp3 = bpp3
                h_cudaBCDoms(nn,sps)%rlv0 = brlv0
                h_cudaBCDoms(nn,sps)%rlv1 = brlv1
                h_cudaBCDoms(nn,sps)%rlv2 = brlv2
                h_cudaBCDoms(nn,sps)%rlv3 = brlv3
                h_cudaBCDoms(nn,sps)%rev0 = brev0
                h_cudaBCDoms(nn,sps)%rev1 = brev1
                h_cudaBCDoms(nn,sps)%rev2 = brev2
                h_cudaBCDoms(nn,sps)%rev3 = brev3
                h_cudaBCDoms(nn,sps)%gamma0 = bgamma0
                h_cudaBCDoms(nn,sps)%gamma1 = bgamma1
                h_cudaBCDoms(nn,sps)%gamma2 = bgamma2
                h_cudaBCDoms(nn,sps)%gamma3 = bgamma3
                h_cudaBCDoms(nn,sps)%ssi = bssi
                h_cudaBCDoms(nn,sps)%ssj = bssj
                h_cudaBCDoms(nn,sps)%ssk = bssk
                h_cudaBCDoms(nn,sps)%ss = bss
                h_cudaBCDoms(nn,sps)%xx = bxx
                h_cudaBCDoms(nn,sps)%dd2wall = bdd2wall
                h_cudaBCDoms(nn,sps)%sFace = bsFace
                h_cudaBCDoms(nn,sps)%gcp = bgcp
                
            end do
        end do

        ! --- Copy data from host back to device---
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                d_cudaBCDoms(nn, sps) = h_cudaBCDoms(nn, sps)
            end do
        end do

    end subroutine copyCudaBC
end module cudaBCPointers