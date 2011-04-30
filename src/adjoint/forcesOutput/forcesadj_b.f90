   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
   !  
   !  Differentiation of forcesadj in reverse (adjoint) mode:
   !   gradient, with respect to input variables: moment padj refpoint
   !                pts force normadj
   !   of linear combination of output variables: moment refpoint
   !                pts force
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          forcesAdj.f90                                   *
   !      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
   !      * Starting date: 08-17-2008                                      *
   !      * Last modified: 08-17-2008                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE FORCESADJ_B(padj, padjb, pts, ptsb, normadj, normadjb, &
   &  refpoint, refpointb, force, forceb, moment, momentb, fact, ibeg, iend&
   &  , jbeg, jend, inode, jnode)
   USE flowvarrefstate
   IMPLICIT NONE
   REAL(KIND=REALTYPE), INTENT(IN) :: fact
   REAL(KIND=REALTYPE) :: force(3), forceb(3)
   INTEGER(KIND=INTTYPE), INTENT(IN) :: ibeg
   INTEGER(KIND=INTTYPE), INTENT(IN) :: iend
   INTEGER(KIND=INTTYPE), INTENT(IN) :: inode
   INTEGER(KIND=INTTYPE), INTENT(IN) :: jbeg
   INTEGER(KIND=INTTYPE), INTENT(IN) :: jend
   INTEGER(KIND=INTTYPE), INTENT(IN) :: jnode
   REAL(KIND=REALTYPE) :: moment(3), momentb(3)
   REAL(KIND=REALTYPE), DIMENSION(3, 2, 2), INTENT(IN) :: normadj
   REAL(KIND=REALTYPE) :: normadjb(3, 2, 2)
   REAL(KIND=REALTYPE), DIMENSION(3, 2, 2), INTENT(IN) :: padj
   REAL(KIND=REALTYPE) :: padjb(3, 2, 2)
   REAL(KIND=REALTYPE), DIMENSION(3, 3, 3), INTENT(IN) :: pts
   REAL(KIND=REALTYPE) :: ptsb(3, 3, 3)
   REAL(KIND=REALTYPE), DIMENSION(3), INTENT(IN) :: refpoint
   REAL(KIND=REALTYPE) :: refpointb(3)
   INTEGER :: branch
   INTEGER(KIND=INTTYPE) :: i, j
   REAL(KIND=REALTYPE) :: tauxx, tauyy, tauzz
   REAL(KIND=REALTYPE) :: tauxy, tauxz, tauyz
   REAL(KIND=REALTYPE) :: addforce(3), addforceb(3), pp, ppb, r(3), rb(3)&
   &  , scaledim, tempb(3), xc(3), xcb(3)
   ! Subroutine Arguments
   ! Local Variables
   scaledim = pref
   ! Force is the contribution of each of 4 cells
   DO j=1,2
   DO i=1,2
   IF (.NOT.(inode + i - 2 .LT. ibeg .OR. inode + i - 1 .GT. iend &
   &          .OR. jnode + j - 2 .LT. jbeg .OR. jnode + j - 1 .GT. jend)) &
   &      THEN
   CALL PUSHREAL8(pp)
   ! Calculate the Pressure
   pp = half*(padj(1, i, j)+padj(2, i, j)) - pinf
   pp = fourth*fact*scaledim*pp
   CALL PUSHREAL8ARRAY(addforce, 3)
   ! Incremental Force to Add
   addforce = pp*normadj(:, i, j)
   ! Add increment to total Force for this node
   ! Cell Center, xc
   xc(:) = fourth*(pts(:, i, j)+pts(:, i+1, j)+pts(:, i, j+1)+pts(:&
   &          , i+1, j+1))
   CALL PUSHREAL8ARRAY(r, 3)
   ! Vector from center to refPoint
   r(:) = xc(:) - refpoint(:)
   ! Moment is F x r ( F cross r)
   !          ! Viscous Force: Below is the code one would use to compute
   !            ! the viscous forces. This code has NOT been verified, and
   !            ! will not run. tau MUST be passed to this function, which
   !            ! will have been computed by a routine similar to
   !            ! viscousFlux.  The outer driving routines, will remain the
   !            ! same (since this is just an extra w-dependance). 
   !            if  (viscousSubface) then
   !               tauXx = tau(i,j,1)
   !               tauYy = tau(i,j,2)
   !               tauZz = tau(i,j,3)
   !               tauXy = tau(i,j,4)
   !               tauXz = tau(i,j,5)
   !               tauYz = tau(i,j,6)
   !               ! Compute the viscous force on the face. A minus sign
   !               ! is now present, due to the definition of this force.
   !               addForce(1)= -fact*(tauXx*normAdj(1,i,j) + tauXy*normAdj(2,i,j) &
   !                    +        tauXz*normAdj(3,i,j))
   !               addForce(2) = -fact*(tauXy*normAdj(i,j,1) + tauYy*normAdj(2,i,j) &
   !                    +        tauYz*normAdj(3,i,j))
   !               addForce(3) = -fact*(tauXz*normAdj(1,i,j) + tauYz*normAdj(2,i,j) &
   !                    +        tauZz*normAdj(3,i,j))
   !               ! Increment the Force
   !               force(:) = force(:) + addForce(:)
   !               ! Increment the Moment
   !               moment(1) = moment(1) + r(2)*addForce(3) - r(3)*addForce(2)
   !               moment(2) = moment(2) + r(3)*addForce(1) - r(1)*addForce(3)
   !               moment(3) = moment(3) + r(1)*addForce(2) - r(2)*addForce(1)
   !            end if
   CALL PUSHINTEGER4(2)
   ELSE
   CALL PUSHINTEGER4(1)
   END IF
   END DO
   END DO
   padjb(1:3, 1:2, 1:2) = 0.0
   normadjb(1:3, 1:2, 1:2) = 0.0
   rb(1:3) = 0.0
   xcb(1:3) = 0.0
   DO j=2,1,-1
   DO i=2,1,-1
   CALL POPINTEGER4(branch)
   IF (.NOT.branch .LT. 2) THEN
   addforceb(1:3) = 0.0
   rb(1) = rb(1) + addforce(2)*momentb(3)
   addforceb(2) = r(1)*momentb(3)
   rb(2) = rb(2) - addforce(1)*momentb(3)
   addforceb(1) = r(3)*momentb(2) - r(2)*momentb(3)
   rb(3) = rb(3) + addforce(1)*momentb(2)
   rb(1) = rb(1) - addforce(3)*momentb(2)
   addforceb(3) = addforceb(3) + r(2)*momentb(1) - r(1)*momentb(2)
   rb(2) = rb(2) + addforce(3)*momentb(1)
   rb(3) = rb(3) - addforce(2)*momentb(1)
   addforceb(2) = addforceb(2) - r(3)*momentb(1)
   CALL POPREAL8ARRAY(r, 3)
   xcb(:) = xcb(:) + rb(:)
   refpointb(:) = refpointb(:) - rb(:)
   rb(:) = 0.0
   tempb = fourth*xcb(:)
   ptsb(:, i, j) = ptsb(:, i, j) + tempb
   ptsb(:, i+1, j) = ptsb(:, i+1, j) + tempb
   ptsb(:, i, j+1) = ptsb(:, i, j+1) + tempb
   ptsb(:, i+1, j+1) = ptsb(:, i+1, j+1) + tempb
   xcb(:) = 0.0
   addforceb(:) = addforceb(:) + forceb(:)
   CALL POPREAL8ARRAY(addforce, 3)
   ppb = SUM(normadj(:, i, j)*addforceb)
   normadjb(:, i, j) = normadjb(:, i, j) + pp*addforceb
   ppb = fourth*fact*scaledim*ppb
   CALL POPREAL8(pp)
   padjb(1, i, j) = padjb(1, i, j) + half*ppb
   padjb(2, i, j) = padjb(2, i, j) + half*ppb
   END IF
   END DO
   END DO
   momentb(:) = 0.0
   forceb(:) = 0.0
   END SUBROUTINE FORCESADJ_B
