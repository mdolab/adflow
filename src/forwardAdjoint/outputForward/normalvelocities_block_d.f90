!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of normalvelocities_block in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *(*bcdata.rface)
!   with respect to varying inputs: *sfacei *sfacej *sfacek *si
!                *sj *sk
!   plus diff mem management of: sfacei:in sfacej:in sfacek:in
!                si:in sj:in sk:in bcdata:in *bcdata.rface:in
!
!      ******************************************************************
!      *                                                                *
!      * file:          normalvelocities.f90                            *
!      * author:        edwin van der weide                             *
!      * starting date: 02-23-2004                                      *
!      * last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine normalvelocities_block_d(sps)
!
!      ******************************************************************
!      *                                                                *
!      * normalvelocitiesalllevels computes the normal grid             *
!      * velocities of some boundary faces of the moving blocks for     *
!      * spectral mode sps. all grid levels from ground level to the    *
!      * coarsest level are considered.                                 *
!      *                                                                *
!      ******************************************************************
!
  use bctypes
  use blockpointers
  use iteration
  use diffsizes
!  hint: isize1ofdrfbcdata should be the size of dimension 1 of array *bcdata
  implicit none
!
!      subroutine arguments.
!
  integer(kind=inttype), intent(in) :: sps
!
!      local variables.
!
  integer(kind=inttype) :: mm
  integer(kind=inttype) :: i, j
  real(kind=realtype) :: weight, mult
  real(kind=realtype) :: weightd
  real(kind=realtype), dimension(:, :), pointer :: sface
  real(kind=realtype), dimension(:, :), pointer :: sfaced
  real(kind=realtype), dimension(:, :, :), pointer :: ss
  real(kind=realtype), dimension(:, :, :), pointer :: ssd
  intrinsic associated
  intrinsic sqrt
  real(kind=realtype) :: arg1
  real(kind=realtype) :: arg1d
  integer :: ii1
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! check for a moving block. as it is possible that in a
! multidisicplinary environment additional grid velocities
! are set, the test should be done on addgridvelocities
! and not on blockismoving.
  if (addgridvelocities) then
    do ii1=1,isize1ofdrfbcdata
      bcdatad(ii1)%rface = 0.0_8
    end do
!
!            ************************************************************
!            *                                                          *
!            * determine the normal grid velocities of the boundaries.  *
!            * as these values are based on the unit normal. a division *
!            * by the length of the normal is needed.                   *
!            * furthermore the boundary unit normals are per definition *
!            * outward pointing, while on the imin, jmin and kmin       *
!            * boundaries the face normals are inward pointing. this    *
!            * is taken into account by the factor mult.                *
!            *                                                          *
!            ************************************************************
!
! loop over the boundary subfaces.
bocoloop:do mm=1,nbocos
! check whether rface is allocated.
      if (associated(bcdata(mm)%rface)) then
! determine the block face on which the subface is
! located and set some variables accordingly.
        select case  (bcfaceid(mm)) 
        case (imin) 
          mult = -one
          ssd => sid(1, :, :, :)
          ss => si(1, :, :, :)
          sfaced => sfaceid(1, :, :)
          sface => sfacei(1, :, :)
        case (imax) 
          mult = one
          ssd => sid(il, :, :, :)
          ss => si(il, :, :, :)
          sfaced => sfaceid(il, :, :)
          sface => sfacei(il, :, :)
        case (jmin) 
          mult = -one
          ssd => sjd(:, 1, :, :)
          ss => sj(:, 1, :, :)
          sfaced => sfacejd(:, 1, :)
          sface => sfacej(:, 1, :)
        case (jmax) 
          mult = one
          ssd => sjd(:, jl, :, :)
          ss => sj(:, jl, :, :)
          sfaced => sfacejd(:, jl, :)
          sface => sfacej(:, jl, :)
        case (kmin) 
          mult = -one
          ssd => skd(:, :, 1, :)
          ss => sk(:, :, 1, :)
          sfaced => sfacekd(:, :, 1)
          sface => sfacek(:, :, 1)
        case (kmax) 
          mult = one
          ssd => skd(:, :, kl, :)
          ss => sk(:, :, kl, :)
          sfaced => sfacekd(:, :, kl)
          sface => sfacek(:, :, kl)
        end select
! loop over the faces of the subface.
        do j=bcdata(mm)%jcbeg,bcdata(mm)%jcend
          do i=bcdata(mm)%icbeg,bcdata(mm)%icend
! compute the inverse of the length of the normal
! vector and possibly correct for inward pointing.
            arg1d = 2*ss(i, j, 1)*ssd(i, j, 1) + 2*ss(i, j, 2)*ssd(i, j&
&             , 2) + 2*ss(i, j, 3)*ssd(i, j, 3)
            arg1 = ss(i, j, 1)**2 + ss(i, j, 2)**2 + ss(i, j, 3)**2
            if (arg1 .eq. 0.0_8) then
              weightd = 0.0_8
            else
              weightd = arg1d/(2.0*sqrt(arg1))
            end if
            weight = sqrt(arg1)
            if (weight .gt. zero) then
              weightd = -(mult*weightd/weight**2)
              weight = mult/weight
            end if
! compute the normal velocity based on the outward
! pointing unit normal.
            bcdatad(mm)%rface(i, j) = weightd*sface(i, j) + weight*&
&             sfaced(i, j)
            bcdata(mm)%rface(i, j) = weight*sface(i, j)
          end do
        end do
      end if
    end do bocoloop
  else
! block is not moving. loop over the boundary faces and set
! the normal grid velocity to zero if allocated.
    do mm=1,nbocos
      if (associated(bcdata(mm)%rface)) then
        bcdatad(mm)%rface = 0.0_8
        bcdata(mm)%rface = zero
      end if
    end do
    do ii1=1,isize1ofdrfbcdata
      bcdatad(ii1)%rface = 0.0_8
    end do
  end if
end subroutine normalvelocities_block_d
