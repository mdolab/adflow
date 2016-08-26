!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of getcostfunction2 in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: funcvalues
!   with respect to varying inputs: gammainf pinf rhoinfdim pinfdim
!                pref machgrid lengthref machcoef dragdirection
!                liftdirection pointref moment sepsensoravg force
!                cavitation sepsensor
subroutine getcostfunction2_d(force, forced, moment, momentd, sepsensor&
& , sepsensord, sepsensoravg, sepsensoravgd, cavitation, cavitationd, &
& alpha, beta, liftindex)
! compute the value of the actual objective function based on the
! (summed) forces and moments and any other "extra" design
! variables. the index of the objective is determined by 'idv'. this
! function is intended to be ad'ed in reverse mode. 
  use constants
  use inputtimespectral
  use costfunctions
  use inputphysics
  use flowvarrefstate
  use inputtsstabderiv
  use utils_d, only : computetsderivatives, computetsderivatives_d, &
& computerootbendingmoment, computerootbendingmoment_d
  implicit none
! input 
  integer(kind=inttype), intent(in) :: liftindex
  real(kind=realtype), dimension(3, ntimeintervalsspectral), intent(in) &
& :: force, moment
  real(kind=realtype), dimension(3, ntimeintervalsspectral), intent(in) &
& :: forced, momentd
  real(kind=realtype), intent(in) :: sepsensor, cavitation, sepsensoravg&
& (3)
  real(kind=realtype), intent(in) :: sepsensord, cavitationd, &
& sepsensoravgd(3)
  real(kind=realtype), intent(in) :: alpha, beta
! working
  real(kind=realtype) :: fact, factmoment, scaledim, ovrnts
  real(kind=realtype) :: factd, factmomentd, scaledimd
  real(kind=realtype), dimension(3) :: cf, cm
  real(kind=realtype), dimension(3) :: cfd, cmd
  real(kind=realtype) :: elasticmomentx, elasticmomenty, elasticmomentz
  real(kind=realtype), dimension(ntimeintervalsspectral, 8) :: basecoef
  real(kind=realtype), dimension(8) :: coef0, dcdalpha, dcdalphadot, &
& dcdq, dcdqdot
  real(kind=realtype), dimension(8) :: coef0d, dcdalphad, dcdalphadotd
  real(kind=realtype) :: bendingmoment
  real(kind=realtype) :: bendingmomentd
  integer(kind=inttype) :: sps
! generate constants
  scaledimd = (prefd*pinf-pref*pinfd)/pinf**2
  scaledim = pref/pinf
  factd = -(two*surfaceref*lref**2*((gammainfd*pinf+gammainf*pinfd)*&
&   machcoef**2*scaledim+gammainf*pinf*(2*machcoef*machcoefd*scaledim+&
&   machcoef**2*scaledimd))/(gammainf*pinf*machcoef**2*surfaceref*lref**&
&   2*scaledim)**2)
  fact = two/(gammainf*pinf*machcoef**2*surfaceref*lref**2*scaledim)
  factmomentd = (factd*lengthref*lref-fact*lref*lengthrefd)/(lengthref*&
&   lref)**2
  factmoment = fact/(lengthref*lref)
  ovrnts = one/ntimeintervalsspectral
! pre-compute ts stability info if required:
  if (tsstability) then
    call computetsderivatives_d(force, forced, moment, momentd, &
&                         liftindex, coef0, coef0d, dcdalpha, dcdalphad&
&                         , dcdalphadot, dcdalphadotd, dcdq, dcdqdot)
  else
    dcdalphadotd = 0.0_8
    coef0d = 0.0_8
    dcdalphad = 0.0_8
  end if
  funcvalues = zero
  funcvaluesd = 0.0_8
! now we just compute each cost function:
  do sps=1,ntimeintervalsspectral
    funcvaluesd(costfuncforcex) = funcvaluesd(costfuncforcex) + ovrnts*&
&     forced(1, sps)
    funcvalues(costfuncforcex) = funcvalues(costfuncforcex) + ovrnts*&
&     force(1, sps)
    funcvaluesd(costfuncforcey) = funcvaluesd(costfuncforcey) + ovrnts*&
&     forced(2, sps)
    funcvalues(costfuncforcey) = funcvalues(costfuncforcey) + ovrnts*&
&     force(2, sps)
    funcvaluesd(costfuncforcez) = funcvaluesd(costfuncforcez) + ovrnts*&
&     forced(3, sps)
    funcvalues(costfuncforcez) = funcvalues(costfuncforcez) + ovrnts*&
&     force(3, sps)
    funcvaluesd(costfuncmomx) = funcvaluesd(costfuncmomx) + ovrnts*&
&     momentd(1, sps)
    funcvalues(costfuncmomx) = funcvalues(costfuncmomx) + ovrnts*moment(&
&     1, sps)
    funcvaluesd(costfuncmomy) = funcvaluesd(costfuncmomy) + ovrnts*&
&     momentd(2, sps)
    funcvalues(costfuncmomy) = funcvalues(costfuncmomy) + ovrnts*moment(&
&     2, sps)
    funcvaluesd(costfuncmomz) = funcvaluesd(costfuncmomz) + ovrnts*&
&     momentd(3, sps)
    funcvalues(costfuncmomz) = funcvalues(costfuncmomz) + ovrnts*moment(&
&     3, sps)
    funcvaluesd(costfuncsepsensor) = funcvaluesd(costfuncsepsensor) + &
&     ovrnts*sepsensord
    funcvalues(costfuncsepsensor) = funcvalues(costfuncsepsensor) + &
&     ovrnts*sepsensor
    funcvaluesd(costfunccavitation) = funcvaluesd(costfunccavitation) + &
&     ovrnts*cavitationd
    funcvalues(costfunccavitation) = funcvalues(costfunccavitation) + &
&     ovrnts*cavitation
    funcvaluesd(costfuncsepsensoravgx) = funcvaluesd(&
&     costfuncsepsensoravgx) + ovrnts*sepsensoravgd(1)
    funcvalues(costfuncsepsensoravgx) = funcvalues(costfuncsepsensoravgx&
&     ) + ovrnts*sepsensoravg(1)
    funcvaluesd(costfuncsepsensoravgy) = funcvaluesd(&
&     costfuncsepsensoravgy) + ovrnts*sepsensoravgd(2)
    funcvalues(costfuncsepsensoravgy) = funcvalues(costfuncsepsensoravgy&
&     ) + ovrnts*sepsensoravg(2)
    funcvaluesd(costfuncsepsensoravgz) = funcvaluesd(&
&     costfuncsepsensoravgz) + ovrnts*sepsensoravgd(3)
    funcvalues(costfuncsepsensoravgz) = funcvalues(costfuncsepsensoravgz&
&     ) + ovrnts*sepsensoravg(3)
! bending moment calc
    cmd = factmomentd*moment(:, sps) + factmoment*momentd(:, sps)
    cm = factmoment*moment(:, sps)
    cfd = factd*force(:, sps) + fact*forced(:, sps)
    cf = fact*force(:, sps)
    call computerootbendingmoment_d(cf, cfd, cm, cmd, liftindex, &
&                             bendingmoment, bendingmomentd)
    funcvaluesd(costfuncbendingcoef) = funcvaluesd(costfuncbendingcoef) &
&     + ovrnts*bendingmomentd
    funcvalues(costfuncbendingcoef) = funcvalues(costfuncbendingcoef) + &
&     ovrnts*bendingmoment
  end do
  funcvaluesd(costfuncforcexcoef) = funcvaluesd(costfuncforcex)*fact + &
&   funcvalues(costfuncforcex)*factd
  funcvalues(costfuncforcexcoef) = funcvalues(costfuncforcex)*fact
  funcvaluesd(costfuncforceycoef) = funcvaluesd(costfuncforcey)*fact + &
&   funcvalues(costfuncforcey)*factd
  funcvalues(costfuncforceycoef) = funcvalues(costfuncforcey)*fact
  funcvaluesd(costfuncforcezcoef) = funcvaluesd(costfuncforcez)*fact + &
&   funcvalues(costfuncforcez)*factd
  funcvalues(costfuncforcezcoef) = funcvalues(costfuncforcez)*fact
  funcvaluesd(costfuncmomxcoef) = funcvaluesd(costfuncmomx)*factmoment +&
&   funcvalues(costfuncmomx)*factmomentd
  funcvalues(costfuncmomxcoef) = funcvalues(costfuncmomx)*factmoment
  funcvaluesd(costfuncmomycoef) = funcvaluesd(costfuncmomy)*factmoment +&
&   funcvalues(costfuncmomy)*factmomentd
  funcvalues(costfuncmomycoef) = funcvalues(costfuncmomy)*factmoment
  funcvaluesd(costfuncmomzcoef) = funcvaluesd(costfuncmomz)*factmoment +&
&   funcvalues(costfuncmomz)*factmomentd
  funcvalues(costfuncmomzcoef) = funcvalues(costfuncmomz)*factmoment
  funcvaluesd(costfunclift) = funcvaluesd(costfuncforcex)*liftdirection(&
&   1) + funcvalues(costfuncforcex)*liftdirectiond(1) + funcvaluesd(&
&   costfuncforcey)*liftdirection(2) + funcvalues(costfuncforcey)*&
&   liftdirectiond(2) + funcvaluesd(costfuncforcez)*liftdirection(3) + &
&   funcvalues(costfuncforcez)*liftdirectiond(3)
  funcvalues(costfunclift) = funcvalues(costfuncforcex)*liftdirection(1)&
&   + funcvalues(costfuncforcey)*liftdirection(2) + funcvalues(&
&   costfuncforcez)*liftdirection(3)
  funcvaluesd(costfuncdrag) = funcvaluesd(costfuncforcex)*dragdirection(&
&   1) + funcvalues(costfuncforcex)*dragdirectiond(1) + funcvaluesd(&
&   costfuncforcey)*dragdirection(2) + funcvalues(costfuncforcey)*&
&   dragdirectiond(2) + funcvaluesd(costfuncforcez)*dragdirection(3) + &
&   funcvalues(costfuncforcez)*dragdirectiond(3)
  funcvalues(costfuncdrag) = funcvalues(costfuncforcex)*dragdirection(1)&
&   + funcvalues(costfuncforcey)*dragdirection(2) + funcvalues(&
&   costfuncforcez)*dragdirection(3)
  funcvaluesd(costfuncliftcoef) = funcvaluesd(costfunclift)*fact + &
&   funcvalues(costfunclift)*factd
  funcvalues(costfuncliftcoef) = funcvalues(costfunclift)*fact
  funcvaluesd(costfuncdragcoef) = funcvaluesd(costfuncdrag)*fact + &
&   funcvalues(costfuncdrag)*factd
  funcvalues(costfuncdragcoef) = funcvalues(costfuncdrag)*fact
! -------------------- time spectral objectives ------------------
  funcvaluesd(costfunccl0) = coef0d(1)
  funcvalues(costfunccl0) = coef0(1)
  funcvaluesd(costfunccd0) = coef0d(2)
  funcvalues(costfunccd0) = coef0(2)
  funcvaluesd(costfunccm0) = coef0d(8)
  funcvalues(costfunccm0) = coef0(8)
  funcvaluesd(costfuncclalpha) = dcdalphad(1)
  funcvalues(costfuncclalpha) = dcdalpha(1)
  funcvaluesd(costfunccdalpha) = dcdalphad(2)
  funcvalues(costfunccdalpha) = dcdalpha(2)
  funcvaluesd(costfunccmzalpha) = dcdalphad(8)
  funcvalues(costfunccmzalpha) = dcdalpha(8)
  funcvaluesd(costfuncclalphadot) = dcdalphadotd(1)
  funcvalues(costfuncclalphadot) = dcdalphadot(1)
  funcvaluesd(costfunccdalphadot) = dcdalphadotd(2)
  funcvalues(costfunccdalphadot) = dcdalphadot(2)
  funcvaluesd(costfunccmzalphadot) = dcdalphadotd(8)
  funcvalues(costfunccmzalphadot) = dcdalphadot(8)
  funcvaluesd(costfuncclq) = 0.0_8
  funcvalues(costfuncclq) = dcdq(1)
  funcvaluesd(costfunccdq) = 0.0_8
  funcvalues(costfunccdq) = dcdq(2)
  funcvaluesd(costfunccmzq) = 0.0_8
  funcvalues(costfunccmzq) = dcdq(8)
  funcvaluesd(costfuncclqdot) = 0.0_8
  funcvalues(costfuncclqdot) = dcdqdot(1)
  funcvaluesd(costfunccdqdot) = 0.0_8
  funcvalues(costfunccdqdot) = dcdqdot(2)
  funcvaluesd(costfunccmzqdot) = 0.0_8
  funcvalues(costfunccmzqdot) = dcdqdot(8)
end subroutine getcostfunction2_d
