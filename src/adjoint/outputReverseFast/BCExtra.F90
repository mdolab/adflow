subroutine applyAllBC_block_fast_b(secondHalo)

  ! Apply BC's for a single block
  
  use blockPointers
  use flowVarRefState
  use inputDiscretization
  use inputTimeSpectral
  use iteration
  use bcTypes
  use bcroutines_fast_b
  implicit none

    ! Subroutine arguments.
    logical, intent(in) :: secondHalo

    ! Local variables.
    logical :: correctForK
    integer(kind=intType) :: nn
    !
    ! Determine whether or not the total energy must be corrected
    ! for the presence of the turbulent kinetic energy.
    if( kPresent ) then
       if((currentLevel <= groundLevel) .or. turbCoupled) then
          correctForK = .true.
       else
          correctForK = .false.
       endif
    else
       correctForK = .false.
    endif

    ! Apply all the boundary conditions. The order is important!  Only
    ! some of them have been AD'ed
    ! ------------------------------------
    !  Symmetry Boundary Condition 
    ! ------------------------------------
    DO nn=1,nbocos
       IF (bctype(nn) .EQ. symm) THEN
          CALL SETBCPOINTERS(nn, .false.)
          CALL BCSYMM(nn, secondhalo)
          END IF
    END DO

    ! ------------------------------------
    !  Adibatic Wall Boundary Condition 
    ! ------------------------------------
    DO nn=1,nviscbocos
       IF (bctype(nn) .EQ. nswalladiabatic) THEN
          CALL SETBCPOINTERS(nn, .false.)
          CALL BCNSWALLADIABATIC(nn, secondhalo, correctfork)
       END IF
    END DO

    ! ------------------------------------
    !  Farfield Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == farField) then
          call setBCPointers(nn, .False.)
          call bcFarField(nn, secondHalo, correctForK)
       end if
    end do

    ! ------------------------------------
    !  Farfield Boundary Condition 
    ! ------------------------------------
    do nn=1,nBocos
       if (bcType(nn) == farField) then
          call setBCPointers_fast_b(nn, .False.)
          call bcFarField_fast_b(nn, secondHalo, correctForK)
       end if
    end do

    do nn=1,nBocos
       if (bcType(nn) == eulerWall) then
          call setBCPointers_fast_b(nn, .False.)
          call bcEulerWall_fast_b(nn, secondHalo, correctForK)
       end if
    end do


    ! And now do the reverse pass
    ! ------------------------------------
    !  Adibatic Wall Boundary Condition 
    ! ------------------------------------
    DO nn=1,nviscbocos
       IF (bctype(nn) .EQ. nswalladiabatic) THEN
          CALL SETBCPOINTERS_fast_b(nn, .false.)
          CALL BCNSWALLADIABATIC_FAST_B(nn, secondhalo, correctfork)
       END IF
    END DO
   
    ! ------------------------------------
    !  Symmetry Boundary Condition 
    ! ------------------------------------
    DO nn=1,nbocos
       IF (bctype(nn) .EQ. symm) THEN
          CALL SETBCPOINTERS_fast_b(nn, .false.)
          CALL BCSYMM_fast_b(nn, secondhalo)
          END IF
    END DO


    ! ! ------------------------------------
    ! !  Euler Wall Boundary Condition 
    ! ! ------------------------------------
    ! do nn=nBocos,1,-1
    !    if (bcType(nn) == EulerWall) then
    !       call setBCPointers_fast_b(nn, .True.)
    !       call bcEulerWall_fast_b(nn, secondHalo, correctForK)
    !    end if
    ! end do

    ! ! ------------------------------------
    ! !  Farfield Boundary Condition 
    ! ! ------------------------------------
    ! if (precond == Turkel .or. precond == ChoiMerkle) then 
    !    call returnFail("applyAllBC", &
    !         "Farfield Turkel and Coid/Merkle preconditioners not implemented")
    ! end if
    ! do nn=nBocos,1,-1
    !    if (bcType(nn) == farField) then
    !       call setBCPointers_fast_b(nn, .False.)
    !       call bcFarField_fast_b(nn, secondHalo, correctForK)
    !    end if
    ! end do

    ! ! ! ------------------------------------
    ! ! !  Isotermal Wall Boundary Condition 
    ! ! ! ------------------------------------
    ! ! do nn=nViscBocos,1,-1
    ! !    if (bcType(nn) == NSWallIsoThermal) then 
    ! !       call setBCPointers_fast_b(nn, .False.)
    ! !       call bcNSWallIsothermal_fast_b(nn, secondHalo, correctForK)
    ! !    end if
    ! ! end do

    ! ! ------------------------------------
    ! !  adibatic Wall Boundary Condition 
    ! ! ------------------------------------
    ! do nn=nViscBocos,1,-1
    !    if (bcType(nn) == NSWallAdiabatic) then 
    !       call setBCPointers_fast_b(nn, .False.)
    !       call bcNSWallAdiabatic_fast_b(nn, secondHalo, correctForK)
    !    end if
    ! end do

    ! ! ------------------------------------
    ! !  Symmetry Boundary Condition 
    ! ! ------------------------------------
    ! do nn=nBocos,1,-1
    !    if (bcType(nn) == symm) then 
    !       call setBCPointers_fast_b(nn, .False.)
    !       call bcSymm_fast_b(nn, secondHalo)
    !    end if
    ! end do
  end subroutine applyAllBC_block_fast_b


 subroutine setBCPointers_fast_b(nn, spatialPointers)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * setBCPointers sets the pointers needed for the boundary        *
    !      * condition treatment on a general face, such that the boundary  *
    !      * routines are only implemented once instead of 6 times.         *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use flowVarRefState
    use bcroutines_fast_b
    implicit none

    ! Subroutine arguments.
    integer(kind=intType), intent(in) :: nn
    logical, intent(in) :: spatialPointers

    ! Determine the sizes of each face and point to just the range we
    ! need on each face. 
    iStart = BCData(nn)%icBeg
    iEnd   = BCData(nn)%icEnd
    jStart = BCData(nn)%jcBeg
    jEnd   = BCData(nn)%jcEnd

    ! Set the size of the subface
    isize = iEnd-iStart + 1
    jsize = jEnd-jStart + 1

    ! Determine the face id on which the subface is located and set
    ! the pointers accordinly.

    select case (BCFaceID(nn))

       !===============================================================
    case (iMin)

       ww3 => w(3, 1:, 1:, :)
       ww2 => w(2, 1:, 1:, :)
       ww1 => w(1, 1:, 1:, :)
       ww0 => w(0, 1:, 1:, :)

       pp3 => p(3, 1:, 1:)
       pp2 => p(2, 1:, 1:)
       pp1 => p(1, 1:, 1:)
       pp0 => p(0, 1:, 1:)

       rlv3 => rlv(3, 1:, 1:)
       rlv2 => rlv(2, 1:, 1:)
       rlv1 => rlv(1, 1:, 1:)
       rlv0 => rlv(0, 1:, 1:)

       rev3 => rev(3, 1:, 1:)
       rev2 => rev(2, 1:, 1:)
       rev1 => rev(1, 1:, 1:)
       rev0 => rev(0, 1:, 1:)

       gamma3 => gamma(3, 1:, 1:)
       gamma2 => gamma(2, 1:, 1:)
       gamma1 => gamma(1, 1:, 1:)
       gamma0 => gamma(0, 1:, 1:)

       !===============================================================

    case (iMax)

       ww3 => w(nx, 1:, 1:, :)
       ww2 => w(il, 1:, 1:, :)
       ww1 => w(ie, 1:, 1:, :)
       ww0 => w(ib, 1:, 1:, :)

       pp3 => p(nx, 1:, 1:)
       pp2 => p(il, 1:, 1:)
       pp1 => p(ie, 1:, 1:)
       pp0 => p(ib, 1:, 1:)

       rlv3 => rlv(nx, 1:, 1:)
       rlv2 => rlv(il, 1:, 1:)
       rlv1 => rlv(ie, 1:, 1:)
       rlv0 => rlv(ib, 1:, 1:)

       rev3 => rev(nx, 1:, 1:)
       rev2 => rev(il, 1:, 1:)
       rev1 => rev(ie, 1:, 1:)
       rev0 => rev(ib, 1:, 1:)

       gamma3 => gamma(nx, 1:, 1:)
       gamma2 => gamma(il, 1:, 1:)
       gamma1 => gamma(ie, 1:, 1:)
       gamma0 => gamma(ib, 1:, 1:)

       !===============================================================

    case (jMin)

       ww3 => w(1:, 3, 1:, :)
       ww2 => w(1:, 2, 1:, :)
       ww1 => w(1:, 1, 1:, :)
       ww0 => w(1:, 0, 1:, :)

       pp3 => p(1:, 3, 1:)
       pp2 => p(1:, 2, 1:)
       pp1 => p(1:, 1, 1:)
       pp0 => p(1:, 0, 1:)

       rlv3 => rlv(1:, 3, 1:)
       rlv2 => rlv(1:, 2, 1:)
       rlv1 => rlv(1:, 1, 1:)
       rlv0 => rlv(1:, 0, 1:)

       rev3 => rev(1:, 3, 1:)
       rev2 => rev(1:, 2, 1:)
       rev1 => rev(1:, 1, 1:)
       rev0 => rev(1:, 0, 1:)

       gamma3 => gamma(1:, 3, 1:)
       gamma2 => gamma(1:, 2, 1:)
       gamma1 => gamma(1:, 1, 1:)
       gamma0 => gamma(1:, 0, 1:)

       !===============================================================

    case (jMax)

       ww3 => w(1:, ny, 1:, :)
       ww2 => w(1:, jl, 1:, :)
       ww1 => w(1:, je, 1:, :)
       ww0 => w(1:, jb, 1:, :)

       pp3 => p(1:, ny, 1:)
       pp2 => p(1:, jl, 1:)
       pp1 => p(1:, je, 1:)
       pp0 => p(1:, jb, 1:)

       rlv3 => rlv(1:, ny, 1:)
       rlv2 => rlv(1:, jl, 1:)
       rlv1 => rlv(1:, je, 1:)
       rlv0 => rlv(1:, jb, 1:)

       rev3 => rev(1:, ny, 1:)
       rev2 => rev(1:, jl, 1:)
       rev1 => rev(1:, je, 1:)
       rev0 => rev(1:, jb, 1:)

       gamma3 => gamma(1:, ny, 1:)
       gamma2 => gamma(1:, jl, 1:)
       gamma1 => gamma(1:, je, 1:)
       gamma0 => gamma(1:, jb, 1:)

       !===============================================================

    case (kMin)

       ww3 => w(1:, 1:, 3, :)
       ww2 => w(1:, 1:, 2, :)
       ww1 => w(1:, 1:, 1, :)
       ww0 => w(1:, 1:, 0, :)

       pp3 => p(1:, 1:, 3)
       pp2 => p(1:, 1:, 2)
       pp1 => p(1:, 1:, 1)
       pp0 => p(1:, 1:, 0)

       rlv3 => rlv(1:, 1:, 3)
       rlv2 => rlv(1:, 1:, 2)
       rlv1 => rlv(1:, 1:, 1)
       rlv0 => rlv(1:, 1:, 0)

       rev3 => rev(1:, 1:, 3)
       rev2 => rev(1:, 1:, 2)
       rev1 => rev(1:, 1:, 1)
       rev0 => rev(1:, 1:, 0)

       gamma3 => gamma(1:, 1:, 3)
       gamma2 => gamma(1:, 1:, 2)
       gamma1 => gamma(1:, 1:, 1)
       gamma0 => gamma(1:, 1:, 0)

       !===============================================================

    case (kMax)

       ww3 => w(1:, 1:, nz, :)
       ww2 => w(1:, 1:, kl, :)
       ww1 => w(1:, 1:, ke, :)
       ww0 => w(1:, 1:, kb, :)

       pp3 => p(1:, 1:, nz)
       pp2 => p(1:, 1:, kl)
       pp1 => p(1:, 1:, ke)
       pp0 => p(1:, 1:, kb)

       rlv3 => rlv(1:, 1:, nz)
       rlv2 => rlv(1:, 1:, kl)
       rlv1 => rlv(1:, 1:, ke)
       rlv0 => rlv(1:, 1:, kb)

       rev3 => rev(1:, 1:, nz)
       rev2 => rev(1:, 1:, kl)
       rev1 => rev(1:, 1:, ke)
       rev0 => rev(1:, 1:, kb)

       gamma3 => gamma(1:, 1:, nz)
       gamma2 => gamma(1:, 1:, kl)
       gamma1 => gamma(1:, 1:, ke)
       gamma0 => gamma(1:, 1:, kb)

    end select

    if (spatialPointers) then 
       select case (BCFaceID(nn))
       case (iMin)
          xx => x(1,:,:,:)
          ssi => si(1,:,:,:)
          ssj => sj(2,:,:,:)
          ssk => sk(2,:,:,:)
          ss  => s (2,:,:,:)
          dd2Wall => d2Wall(2,:,:)
       case (iMax)
          xx => x(il,:,:,:)
          ssi => si(il,:,:,:)
          ssj => sj(il,:,:,:)
          ssk => sk(il,:,:,:)
          ss  =>  s(il,:,:,:)
          dd2Wall => d2Wall(il,:,:)
       case (jMin)
          xx => x(:,1,:,:)
          ssi => sj(:,1,:,:)
          ssj => si(:,2,:,:)
          ssk => sk(:,2,:,:)
          ss   => s(:,2,:,:)
          dd2Wall => d2Wall(:,2,:)
       case (jMax)
          xx => x(:,jl,:,:)
          ssi => sj(:,jl,:,:)
          ssj => si(:,jl,:,:)
          ssk => sk(:,jl,:,:)
          ss  =>  s(:,jl,:,:)
          dd2Wall => d2Wall(:,jl,:)
       case (kMin)
          xx => x(:,:,1,:)
          ssi => sk(:,:,1,:)
          ssj => si(:,:,2,:)
          ssk => sj(:,:,2,:)
          ss  =>  s(:,:,2,:)
          dd2Wall => d2Wall(:,:,2)
       case (kMax)
          xx => x(:,:,kl,:)
          ssi => sk(:,:,kl,:)
          ssj => si(:,:,kl,:)
          ssk => sj(:,:,kl,:)
          ss  =>  s(:,:,kl,:)
          dd2Wall => d2Wall(:,:,kl)
       end select
    end if

    ! And now all the derivative values
    select case (BCFaceID(nn))

       !===============================================================
    case (iMin)

       ww3d => wd(3, 1:, 1:, :)
       ww2d => wd(2, 1:, 1:, :)
       ww1d => wd(1, 1:, 1:, :)
       ww0d => wd(0, 1:, 1:, :)

       pp3d => pd(3, 1:, 1:)
       pp2d => pd(2, 1:, 1:)
       pp1d => pd(1, 1:, 1:)
       pp0d => pd(0, 1:, 1:)

       rlv3d => rlvd(3, 1:, 1:)
       rlv2d => rlvd(2, 1:, 1:)
       rlv1d => rlvd(1, 1:, 1:)
       rlv0d => rlvd(0, 1:, 1:)

       rev3d => revd(3, 1:, 1:)
       rev2d => revd(2, 1:, 1:)
       rev1d => revd(1, 1:, 1:)
       rev0d => revd(0, 1:, 1:)

       !gamma3d => gammad(3, 1:, 1:)
       !gamma2d => gammad(2, 1:, 1:)
       !gamma1d => gammad(1, 1:, 1:)
       !gamma0d => gammad(0, 1:, 1:)

       !===============================================================

    case (iMax)

       ww3d => wd(nx, 1:, 1:, :)
       ww2d => wd(il, 1:, 1:, :)
       ww1d => wd(ie, 1:, 1:, :)
       ww0d => wd(ib, 1:, 1:, :)

       pp3d => pd(nx, 1:, 1:)
       pp2d => pd(il, 1:, 1:)
       pp1d => pd(ie, 1:, 1:)
       pp0d => pd(ib, 1:, 1:)

       rlv3d => rlvd(nx, 1:, 1:)
       rlv2d => rlvd(il, 1:, 1:)
       rlv1d => rlvd(ie, 1:, 1:)
       rlv0d => rlvd(ib, 1:, 1:)

       rev3d => revd(nx, 1:, 1:)
       rev2d => revd(il, 1:, 1:)
       rev1d => revd(ie, 1:, 1:)
       rev0d => revd(ib, 1:, 1:)

       !gamma3d => gammad(nx, 1:, 1:)
       !gamma2d => gammad(il, 1:, 1:)
       !gamma1d => gammad(ie, 1:, 1:)
       !gamma0d => gammad(ib, 1:, 1:)

       !===============================================================

    case (jMin)

       ww3d => wd(1:, 3, 1:, :)
       ww2d => wd(1:, 2, 1:, :)
       ww1d => wd(1:, 1, 1:, :)
       ww0d => wd(1:, 0, 1:, :)

       pp3d => pd(1:, 3, 1:)
       pp2d => pd(1:, 2, 1:)
       pp1d => pd(1:, 1, 1:)
       pp0d => pd(1:, 0, 1:)

       rlv3d => rlvd(1:, 3, 1:)
       rlv2d => rlvd(1:, 2, 1:)
       rlv1d => rlvd(1:, 1, 1:)
       rlv0d => rlvd(1:, 0, 1:)

       rev3d => revd(1:, 3, 1:)
       rev2d => revd(1:, 2, 1:)
       rev1d => revd(1:, 1, 1:)
       rev0d => revd(1:, 0, 1:)

       !gamma3d => gammad(1:, 3, 1:)
       !gamma2d => gammad(1:, 2, 1:)
       !gamma1d => gammad(1:, 1, 1:)
       !gamma0d => gammad(1:, 0, 1:)

       !===============================================================

    case (jMax)

       ww3d => wd(1:, ny, 1:, :)
       ww2d => wd(1:, jl, 1:, :)
       ww1d => wd(1:, je, 1:, :)
       ww0d => wd(1:, jb, 1:, :)

       pp3d => pd(1:, ny, 1:)
       pp2d => pd(1:, jl, 1:)
       pp1d => pd(1:, je, 1:)
       pp0d => pd(1:, jb, 1:)

       rlv3d => rlvd(1:, ny, 1:)
       rlv2d => rlvd(1:, jl, 1:)
       rlv1d => rlvd(1:, je, 1:)
       rlv0d => rlvd(1:, jb, 1:)

       rev3d => revd(1:, ny, 1:)
       rev2d => revd(1:, jl, 1:)
       rev1d => revd(1:, je, 1:)
       rev0d => revd(1:, jb, 1:)

       !gamma3d => gammad(1:, ny, 1:)
       !gamma2d => gammad(1:, jl, 1:)
       !gamma1d => gammad(1:, je, 1:)
       !gamma0d => gammad(1:, jb, 1:)

       !===============================================================

    case (kMin)

       ww3d => wd(1:, 1:, 3, :)
       ww2d => wd(1:, 1:, 2, :)
       ww1d => wd(1:, 1:, 1, :)
       ww0d => wd(1:, 1:, 0, :)

       pp3d => pd(1:, 1:, 3)
       pp2d => pd(1:, 1:, 2)
       pp1d => pd(1:, 1:, 1)
       pp0d => pd(1:, 1:, 0)

       rlv3d => rlvd(1:, 1:, 3)
       rlv2d => rlvd(1:, 1:, 2)
       rlv1d => rlvd(1:, 1:, 1)
       rlv0d => rlvd(1:, 1:, 0)

       rev3d => revd(1:, 1:, 3)
       rev2d => revd(1:, 1:, 2)
       rev1d => revd(1:, 1:, 1)
       rev0d => revd(1:, 1:, 0)

       !gamma3d => gammad(1:, 1:, 3)
       !gamma2d => gammad(1:, 1:, 2)
       !gamma1d => gammad(1:, 1:, 1)
       !gamma0d => gammad(1:, 1:, 0)

       !===============================================================

    case (kMax)

       ww3d => wd(1:, 1:, nz, :)
       ww2d => wd(1:, 1:, kl, :)
       ww1d => wd(1:, 1:, ke, :)
       ww0d => wd(1:, 1:, kb, :)

       pp3d => pd(1:, 1:, nz)
       pp2d => pd(1:, 1:, kl)
       pp1d => pd(1:, 1:, ke)
       pp0d => pd(1:, 1:, kb)

       rlv3d => rlvd(1:, 1:, nz)
       rlv2d => rlvd(1:, 1:, kl)
       rlv1d => rlvd(1:, 1:, ke)
       rlv0d => rlvd(1:, 1:, kb)

       rev3d => revd(1:, 1:, nz)
       rev2d => revd(1:, 1:, kl)
       rev1d => revd(1:, 1:, ke)
       rev0d => revd(1:, 1:, kb)

       !gamma3d => gammad(1:, 1:, nz)
       !gamma2d => gammad(1:, 1:, kl)
       !gamma1d => gammad(1:, 1:, ke)
       !gamma0d => gammad(1:, 1:, kb)

    end select

    ! if (spatialPointers) then 
    !    select case (BCFaceID(nn))
    !    case (iMin)
    !       xxd => xd(1,:,:,:)
    !       ssid => sid(1,:,:,:)
    !       ssjd => sjd(2,:,:,:)
    !       sskd => skd(2,:,:,:)
    !       ssd  => sd (2,:,:,:)
    !       dd2Wall => d2Walld(2,:,:)
    !    case (iMax)
    !       xxd => xd(il,:,:,:)
    !       ssid => sid(il,:,:,:)
    !       ssjd => sjd(il,:,:,:)
    !       sskd => skd(il,:,:,:)
    !       ssd  =>  sd(il,:,:,:)
    !       dd2Wall => d2Walld(il,:,:)
    !    case (jMin)
    !       xxd => xd(:,1,:,:)
    !       ssid => sjd(:,1,:,:)
    !       ssjd => sid(:,2,:,:)
    !       sskd => skd(:,2,:,:)
    !       ssd   => sd(:,2,:,:)
    !       dd2Wall => d2Walld(:,2,:)
    !    case (jMax)
    !       xxd => xd(:,jl,:,:)
    !       ssid => sjd(:,jl,:,:)
    !       ssjd => sid(:,jl,:,:)
    !       sskd => skd(:,jl,:,:)
    !       ssd  =>  sd(:,jl,:,:)
    !       dd2Wall => d2Walld(:,jl,:)
    !    case (kMin)
    !       xxd => xd(:,:,1,:)
    !       ssid => skd(:,:,1,:)
    !       ssjd => sid(:,:,2,:)
    !       sskd => sjd(:,:,2,:)
    !       ssd  =>  sd(:,:,2,:)
    !       dd2Wall => d2Walld(:,:,2)
    !    case (kMax)
    !       xxd => xd(:,:,kl,:)
    !       ssid => skd(:,:,kl,:)
    !       ssjd => sid(:,:,kl,:)
    !       sskd => sjd(:,:,kl,:)
    !       ssd  =>  sd(:,:,kl,:)
    !       dd2Walld => d2Walld(:,:,kl)
    !    end select
    ! end if
  end subroutine setBCPointers_fast_b

  

