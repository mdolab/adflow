   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of referencestate in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: gammainf pinf timeref rhoinf
   !                muref tref muinf uinf rgas pref
   !   with respect to varying inputs: pref mach tempfreestream veldirfreestream
   !                machcoef
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          referenceState.f90                              *
   !      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
   !      * Starting date: 05-29-2003                                      *
   !      * Last modified: 04-22-2006                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE REFERENCESTATE_D()
   !
   !      ******************************************************************
   !      *                                                                *
   !      * referenceState computes the reference state values in case     *
   !      * these have not been specified. A distinction is made between   *
   !      * internal and external flows. In case nothing has been          *
   !      * specified for the former a dimensional computation will be     *
   !      * made. For the latter the reference state is set to an          *
   !      * arbitrary state for an inviscid computation and computed for a *
   !      * viscous computation. Furthermore for internal flows an average *
   !      * velocity direction is computed from the boundary conditions,   *
   !      * which is used for initialization.                              *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BCTYPES
   USE BLOCK
   USE COMMUNICATION
   USE CONSTANTS
   USE FLOWVARREFSTATE
   USE INPUTMOTION
   USE INPUTPHYSICS
   USE INPUTTIMESPECTRAL
   USE ITERATION
   IMPLICIT NONE
   !
   !      Local variables.
   !
   INTEGER :: ierr
   INTEGER(kind=inttype) :: sps, nn, mm
   REAL(kind=realtype) :: gm1, ratio, tmp
   REAL(kind=realtype) :: mx, my, mz, re, v, tinfdim
   REAL(kind=realtype) :: mxd, myd, mzd, vd, tinfdimd
   REAL(kind=realtype), DIMENSION(3) :: dirloc, dirglob
   REAL(kind=realtype), DIMENSION(5) :: valloc, valglob
   TYPE(BCDATATYPE), DIMENSION(:), POINTER :: bcdata
   INTERFACE 
   SUBROUTINE VELMAGNANDDIRECTIONSUBFACE(vmag, dir, bcdata, mm)
   USE BLOCK
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: mm
   REAL(kind=realtype), INTENT(OUT) :: vmag
   REAL(kind=realtype), DIMENSION(3), INTENT(INOUT) :: dir
   TYPE(BCDATATYPE), DIMENSION(:), POINTER :: bcdata
   END SUBROUTINE VELMAGNANDDIRECTIONSUBFACE
   END INTERFACE
      INTRINSIC SQRT
   REAL(kind=realtype) :: arg1
   REAL(kind=realtype) :: arg1d
   REAL(kind=realtype) :: result1
   REAL(kind=realtype) :: result1d
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Initialize the dimensional free stream temperature and pressure.
   ! From these values the density and viscosity is computed. For
   ! external viscous and internal computation this is corrected
   ! later on.
   pinfdimd = prefd
   pinfdim = pref
   IF (pref .LE. zero) THEN
   pinfdim = 101325.0_realType
   pinfdimd = 0.0_8
   END IF
   tinfdimd = tempfreestreamd
   tinfdim = tempfreestream
   rhoinfdimd = (pinfdimd*rgasdim*tinfdim-pinfdim*rgasdim*tinfdimd)/(&
   &   rgasdim*tinfdim)**2
   rhoinfdim = pinfdim/(rgasdim*tinfdim)
   mudimd = musuthdim*((tsuthdim+ssuthdim)*1.5_realType*(tinfdim/tsuthdim&
   &   )**0.5*tinfdimd/((tinfdim+ssuthdim)*tsuthdim)-(tsuthdim+ssuthdim)*&
   &   tinfdimd*(tinfdim/tsuthdim)**1.5_realType/(tinfdim+ssuthdim)**2)
   mudim = musuthdim*((tsuthdim+ssuthdim)/(tinfdim+ssuthdim))*(tinfdim/&
   &   tsuthdim)**1.5_realType
   ! Check the flow type we are having here.
   IF (flowtype .EQ. internalflow) THEN
   gammainfd = 0.0_8
   trefd = 0.0_8
   rhorefd = 0.0_8
   ELSE
   ! External flow. Compute the value of gammaInf.
   CALL COMPUTEGAMMA_D(tempfreestream, tempfreestreamd, gammainf, &
   &                 gammainfd, 1)
   ! In case of a viscous problem, compute the
   ! dimensional free stream density and pressure.
   IF (equations .EQ. nsequations .OR. equations .EQ. ransequations) &
   &   THEN
   ! Compute the x, y, and z-components of the Mach number
   ! relative to the body; i.e. the mesh velocity must be
   ! taken into account here.
   mxd = machcoefd*veldirfreestream(1) + machcoef*veldirfreestreamd(1&
   &       )
   mx = machcoef*veldirfreestream(1)
   myd = machcoefd*veldirfreestream(2) + machcoef*veldirfreestreamd(2&
   &       )
   my = machcoef*veldirfreestream(2)
   mzd = machcoefd*veldirfreestream(3) + machcoef*veldirfreestreamd(3&
   &       )
   mz = machcoef*veldirfreestream(3)
   ! Reynolds number per meter, the viscosity using sutherland's
   ! law and the free stream velocity relative to the body.
   re = reynolds/reynoldslength
   mudimd = musuthdim*((tsuthdim+ssuthdim)*1.5*(tempfreestream/&
   &       tsuthdim)**0.5*tempfreestreamd/((tempfreestream+ssuthdim)*&
   &       tsuthdim)-(tsuthdim+ssuthdim)*tempfreestreamd*(tempfreestream/&
   &       tsuthdim)**1.5/(tempfreestream+ssuthdim)**2)
   mudim = musuthdim*((tsuthdim+ssuthdim)/(tempfreestream+ssuthdim))*&
   &       (tempfreestream/tsuthdim)**1.5
   arg1d = rgasdim*((mxd*mx+mx*mxd+myd*my+my*myd+mzd*mz+mz*mzd)*&
   &       gammainf*tempfreestream+(mx*mx+my*my+mz*mz)*(gammainfd*&
   &       tempfreestream+gammainf*tempfreestreamd))
   arg1 = (mx*mx+my*my+mz*mz)*gammainf*rgasdim*tempfreestream
   IF (arg1 .EQ. 0.0_8) THEN
   vd = 0.0_8
   ELSE
   vd = arg1d/(2.0*SQRT(arg1))
   END IF
   v = SQRT(arg1)
   ! Compute the free stream density and pressure.
   ! Set TInfDim to tempFreestream.
   rhoinfdimd = (re*mudimd*v-re*mudim*vd)/v**2
   rhoinfdim = re*mudim/v
   pinfdimd = rgasdim*(rhoinfdimd*tempfreestream+rhoinfdim*&
   &       tempfreestreamd)
   pinfdim = rhoinfdim*rgasdim*tempfreestream
   tinfdimd = tempfreestreamd
   tinfdim = tempfreestream
   END IF
   ! In case the reference pressure, density and temperature were
   ! not specified, set them to the infinity values.
   IF (pref .LE. zero) THEN
   prefd = pinfdimd
   pref = pinfdim
   END IF
   IF (rhoref .LE. zero) THEN
   rhorefd = rhoinfdimd
   rhoref = rhoinfdim
   ELSE
   rhorefd = 0.0_8
   END IF
   IF (tref .LE. zero) THEN
   trefd = tinfdimd
   tref = tinfdim
   ELSE
   trefd = 0.0_8
   END IF
   END IF
   ! Compute the value of muRef, such that the nonDimensional
   ! equations are identical to the dimensional ones.
   ! Note that in the non-dimensionalization of muRef there is
   ! a reference length. However this reference length is 1.0
   ! in this code, because the coordinates are converted to
   ! meters.
   IF (pref*rhoref .EQ. 0.0_8) THEN
   murefd = 0.0_8
   ELSE
   murefd = (prefd*rhoref+pref*rhorefd)/(2.0*SQRT(pref*rhoref))
   END IF
   muref = SQRT(pref*rhoref)
   ! Compute timeRef for a correct nonDimensionalization of the
   ! unsteady equations. Some story as for the reference viscosity
   ! concerning the reference length.
   IF (rhoref/pref .EQ. 0.0_8) THEN
   timerefd = 0.0_8
   ELSE
   timerefd = (rhorefd*pref-rhoref*prefd)/(pref**2*2.0*SQRT(rhoref/pref&
   &     ))
   END IF
   timeref = SQRT(rhoref/pref)
   ! Compute the nonDimensional pressure, density, velocity,
   ! viscosity and gas constant.
   pinfd = (pinfdimd*pref-pinfdim*prefd)/pref**2
   pinf = pinfdim/pref
   rhoinfd = (rhoinfdimd*rhoref-rhoinfdim*rhorefd)/rhoref**2
   rhoinf = rhoinfdim/rhoref
   arg1d = ((gammainfd*pinf+gammainf*pinfd)*rhoinf-gammainf*pinf*rhoinfd)&
   &   /rhoinf**2
   arg1 = gammainf*pinf/rhoinf
   IF (arg1 .EQ. 0.0_8) THEN
   result1d = 0.0_8
   ELSE
   result1d = arg1d/(2.0*SQRT(arg1))
   END IF
   result1 = SQRT(arg1)
   uinfd = machd*result1 + mach*result1d
   uinf = mach*result1
   rgasd = (rgasdim*(rhorefd*tref+rhoref*trefd)*pref-rgasdim*rhoref*tref*&
   &   prefd)/pref**2
   rgas = rgasdim*rhoref*tref/pref
   muinfd = (mudimd*muref-mudim*murefd)/muref**2
   muinf = mudim/muref
      CONTAINS
   !=================================================================
   !===============================================================
   FUNCTION MAXVALUESUBFACE(var)
   IMPLICIT NONE
   !
   !        Function type
   !
   REAL(kind=realtype) :: maxvaluesubface
   !
   !        Function argument.
   !
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: var
   !
   !        Local variables.
   !
   INTEGER(kind=inttype) :: i, j
   INTRINSIC ASSOCIATED
   INTRINSIC MAX
   !
   !        ****************************************************************
   !        *                                                              *
   !        * Begin execution                                              *
   !        *                                                              *
   !        ****************************************************************
   !
   ! Initialize the function to -1 and return immediately if
   ! var is not associated with data.
   maxvaluesubface = -one
   IF (.NOT.ASSOCIATED(var)) THEN
   RETURN
   ELSE
   ! Loop over the owned faces of the subface. As the cell range
   ! may contain halo values, the nodal range is used.
   DO j=bcdata(mm)%jnbeg+1,bcdata(mm)%jnend
   DO i=bcdata(mm)%inbeg+1,bcdata(mm)%inend
   IF (maxvaluesubface .LT. var(i, j)) THEN
   maxvaluesubface = var(i, j)
   ELSE
   maxvaluesubface = maxvaluesubface
   END IF
   END DO
   END DO
   END IF
   END FUNCTION MAXVALUESUBFACE
   END SUBROUTINE REFERENCESTATE_D
