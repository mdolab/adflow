      SUBROUTINE DELQ3DM_solid (IL, JL, KL, I1, I2, J1, J2, K1, K2, XYZ0, S0,&
           DFACEI, DFACEJ, DFACEK, XYZ)

!     ******************************************************************
!     *   DELQ3DM performs stage 1 of the WARPQ3DM 3-space surface     *
!     *   grid perturbation in a form which is reusable by WARP-BLK    *
!     *   It returns face perturbations rather than perturbed face     *
!     *   coordinates.  The three cases of a block face are handled    *
!     *   here by three similar code sections.  Special handling of    *
!     *   fixed corners is avoided to keep the bulk down.              *
!     *                                                                *
!     *   11/29/95  D.Saunders  Adaptation of DELQ3D for specialized   *
!     *                         WARP-BLK and WARPQ3DM used by          *
!     *                         FLO107-MB.                             *
!     *   04/04/96      "       DELQ3DM does only stage 1 now.         *
!     *   12/11/08  C.A.Mader   Converted to *.f90                     *
!     *                                                                *
!     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
!     ******************************************************************

      !IMPLICIT REAL*8 (A-H,O-Z) ! Take out when all compilers have a switch
      use precision
      implicit none
!     Arguments.

      INTEGER(kind=intType)   IL, JL, KL                  ! I Grid array dimensions.
      INTEGER(kind=intType)   I1, I2, J1, J2, K1, K2      ! I Define active face,
                                            !   one pair being equal.
      real(kind=realType) :: XYZ0(3,0:IL+1,0:JL+1,0:KL+1)! I Base face coordinates in
                                            !   appropriate places
      real(kind=realType) :: S0(3,0:IL+1,0:JL+1,0:KL+1)  ! I Base normalized arc-lengths
      real(kind=realType) :: DFACEI(3,JL,KL)             ! O Reqd. face perturbations:
      real(kind=realType) :: DFACEJ(3,IL,KL)             !   DFACEI(1:3,1:JL,1:KL) =
      real(kind=realType) :: DFACEK(3,IL,JL)             !   dX, dY, dZ for an I face, etc.
      real(kind=realType) :: XYZ(3,0:IL+1,0:JL+1,0:KL+1) ! I Grid coordinates: new edges
                                            !   of a face input; unchanged
                                            !   on output
!     Local constants.

      REAL(kind=realType)::       EPS, ONE

      PARAMETER (EPS = 1.E-14, ONE = 1.E+0) ! EPS safeguards a divide by zero -
                                            ! presumably only if result is zero.
!     Local variables.

      INTEGER(kind=intType)   I, J, K
      REAL(kind=realType):: DELI, DELJ, DELK, WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2


!     Execution.
!     ----------


      IF (I1 .EQ. I2) THEN

!        I plane case:
!        -------------

         I = I1

!        Set up the corner perturbations:

         DO K = 1, KL, KL - 1
            DO J = 1, JL, JL - 1
               DFACEI(1,J,K) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DFACEI(2,J,K) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DFACEI(3,J,K) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
            END DO
         END DO

!        Set up intermediate edge perturbations corresponding to the
!        final corners but otherwise derived from the original edges.

         DO J = 1, JL, JL - 1
            DO K = 2, KL - 1
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2
               DFACEI(1,J,K) = WTK1 * DFACEI(1,J, 1) +&
                    WTK2 * DFACEI(1,J,KL)
               DFACEI(2,J,K) = WTK1 * DFACEI(2,J, 1) +&
                    WTK2 * DFACEI(2,J,KL)
               DFACEI(3,J,K) = WTK1 * DFACEI(3,J, 1) +&
                    WTK2 * DFACEI(3,J,KL)
            END DO
         END DO

         DO K = 1, KL, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               DFACEI(1,J,K) = WTJ1 * DFACEI(1, 1,K) +&
                    WTJ2 * DFACEI(1,JL,K)
               DFACEI(2,J,K) = WTJ1 * DFACEI(2, 1,K) +&
                    WTJ2 * DFACEI(2,JL,K)
               DFACEI(3,J,K) = WTJ1 * DFACEI(3, 1,K) +&
                    WTJ2 * DFACEI(3,JL,K)
            END DO
         END DO

!        Interpolate the intermediate perturbations of interior points.
!        The contributions from each pair of edges are not independent.

         DO K = 2, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI(1,1,K) + WTJ2 * DFACEI(1,JL,K)
               DELK = WTK1 * DFACEI(1,J,1) + WTK2 * DFACEI(1,J,KL)

               DFACEI(1,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /&
                    MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI(2,1,K) + WTJ2 * DFACEI(2,JL,K)
               DELK = WTK1 * DFACEI(2,J,1) + WTK2 * DFACEI(2,J,KL)

               DFACEI(2,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /&
                    MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI(3,1,K) + WTJ2 * DFACEI(3,JL,K)
               DELK = WTK1 * DFACEI(3,J,1) + WTK2 * DFACEI(3,J,KL)

               DFACEI(3,J,K) = (ABS (DELJ) * DELJ + ABS (DELK) * DELK) /&
                    MAX (ABS (DELJ) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (J1 .EQ. J2) THEN

!        J plane case:
!        -------------

         J = J1

!        Corner perturbations:

         DO K = 1, KL, KL - 1
            DO I = 1, IL, IL - 1
               DFACEJ(1,I,K) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DFACEJ(2,I,K) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DFACEJ(3,I,K) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
            END DO
         END DO

!        Intermediate edge perturbations:

         DO I = 1, IL, IL - 1
            DO K = 2, KL - 1
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2
               DFACEJ(1,I,K) = WTK1 * DFACEJ(1,I, 1) +&
                              WTK2 * DFACEJ(1,I,KL)
               DFACEJ(2,I,K) = WTK1 * DFACEJ(2,I, 1) +&
                              WTK2 * DFACEJ(2,I,KL)
               DFACEJ(3,I,K) = WTK1 * DFACEJ(3,I, 1) +&
                              WTK2 * DFACEJ(3,I,KL)
            END DO
         END DO

         DO K = 1, KL, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               DFACEJ(1,I,K) = WTI1 * DFACEJ(1, 1,K) +&
                              WTI2 * DFACEJ(1,IL,K)
               DFACEJ(2,I,K) = WTI1 * DFACEJ(2, 1,K) +&
                              WTI2 * DFACEJ(2,IL,K)
               DFACEJ(3,I,K) = WTI1 * DFACEJ(3, 1,K) +&
                              WTI2 * DFACEJ(3,IL,K)
            END DO
         END DO

!        Intermediate perturbations of interior points:

         DO K = 2, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ(1,1,K) + WTI2 * DFACEJ(1,IL,K)
               DELK = WTK1 * DFACEJ(1,I,1) + WTK2 * DFACEJ(1,I,KL)

               DFACEJ(1,I,K) = (ABS (DELI) * DELI + ABS (DELK) * DELK) /&
                               MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ(2,1,K) + WTI2 * DFACEJ(2,IL,K)
               DELK = WTK1 * DFACEJ(2,I,1) + WTK2 * DFACEJ(2,I,KL)

               DFACEJ(2,I,K) = (ABS (DELI) * DELI + ABS (DELK) * DELK) /&
                               MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ(3,1,K) + WTI2 * DFACEJ(3,IL,K)
               DELK = WTK1 * DFACEJ(3,I,1) + WTK2 * DFACEJ(3,I,KL)

               DFACEJ(3,I,K) = (ABS (DELI) * DELI + ABS (DELK) * DELK) /&
                               MAX (ABS (DELI) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (K1 .EQ. K2) THEN

!        K plane case:
!        -------------

         K = K1

!        Corner perturbations:

         DO J = 1, JL, JL - 1
            DO I = 1, IL, IL - 1
               DFACEK(1,I,J) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DFACEK(2,I,J) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DFACEK(3,I,J) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
            END DO
         END DO

!        Intermediate edge perturbations:

         DO I = 1, IL, IL - 1
            DO J = 2, JL - 1
               WTJ2 = S0 (2,I,J,K)
               WTJ1 = ONE - WTJ2
               DFACEK(1,I,J) = WTJ1 * DFACEK(1,I, 1) +&
                              WTJ2 * DFACEK(1,I,JL)
               DFACEK(2,I,J) = WTJ1 * DFACEK(2,I, 1) +&
                              WTJ2 * DFACEK(2,I,JL)
               DFACEK(3,I,J) = WTJ1 * DFACEK(3,I, 1) +&
                              WTJ2 * DFACEK(3,I,JL)
            END DO
         END DO

         DO J = 1, JL, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               DFACEK(1,I,J) = WTI1 * DFACEK(1, 1,J) +&
                              WTI2 * DFACEK(1,IL,J)
               DFACEK(2,I,J) = WTI1 * DFACEK(2, 1,J) +&
                              WTI2 * DFACEK(2,IL,J)
               DFACEK(3,I,J) = WTI1 * DFACEK(3, 1,J) +&
                              WTI2 * DFACEK(3,IL,J)
            END DO
         END DO

!        Intermediate perturbations of interior points:

         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK(1,1,J) + WTI2 * DFACEK(1,IL,J)
               DELJ = WTJ1 * DFACEK(1,I,1) + WTJ2 * DFACEK(1,I,JL)

               DFACEK(1,I,J) = (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /&
                               MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK(2,1,J) + WTI2 * DFACEK(2,IL,J)
               DELJ = WTJ1 * DFACEK(2,I,1) + WTJ2 * DFACEK(2,I,JL)

               DFACEK(2,I,J) = (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /&
                               MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK(3,1,J) + WTI2 * DFACEK(3,IL,J)
               DELJ = WTJ1 * DFACEK(3,I,1) + WTJ2 * DFACEK(3,I,JL)

               DFACEK(3,I,J) = (ABS (DELI) * DELI + ABS (DELJ) * DELJ) /&
                               MAX (ABS (DELI) + ABS (DELJ), EPS)
            END DO
         END DO

      END IF

    END SUBROUTINE DELQ3DM_SOLID
    SUBROUTINE PARAM3DM_solid (IL, JL, KL, XYZ, S)

!     ******************************************************************
!     *   PARAM3DM parameterizes the volume of one block of a multi-   *
!     *   block grid structure by setting up the normalized arc-length *
!     *   increments in all three index directions.                    *
!     *                                                                *
!     *   11/29/95  D.Saunders  Adaptation of PARAMXYZ for specialized *
!     *                         WARP-BLK used by FLO107-MB.            *
!     *   06/19/96      "       Allow for degenerate edges.            *
!     *   12/11/08  C.A.Mader   Converted to *.f90                     *
!     *                                                                *
!     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
!     ******************************************************************

      !IMPLICIT REAL*8 (A-H,O-Z) ! Take out when all compilers have a switch
      use precision

      implicit none
!     Arguments.

      INTEGER(kind=intType)::   IL, JL, KL                  ! I Grid array dimensions.
      real(kind=realType):: XYZ(3,0:IL+1,0:JL+1,0:KL+1) ! I Grid coordinates
      real(kind=realType):: S(3,0:IL+1,0:JL+1,0:KL+1)   ! O Normalized arc-lengths:
                                            !   S(1,1,J,K) = 0.,
                                            !   S(2,I,1,K) = 0.,
                                            !   S(3,I,J,1) = 0.,
                                            !   S(1,IL,J,K) = 1.,etc.
!     Local constants.

      REAL      ONE, ZERO
      PARAMETER (ONE = 1.E+0, ZERO = 0.E+0)

!     Local variables.

      INTEGER(kind=intType)::   I, J, K

!     Local functions.

      REAL      DELI, DELJ, DELK

      DELI(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I-1,J,K)) ** 2 +&
           (XYZ(2,I,J,K) - XYZ(2,I-1,J,K)) ** 2 +&
           (XYZ(3,I,J,K) - XYZ(3,I-1,J,K)) ** 2)

      DELJ(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I,J-1,K)) ** 2 +&
           (XYZ(2,I,J,K) - XYZ(2,I,J-1,K)) ** 2 +&
           (XYZ(3,I,J,K) - XYZ(3,I,J-1,K)) ** 2)

      DELK(I,J,K) = SQRT ((XYZ(1,I,J,K) - XYZ(1,I,J,K-1)) ** 2 +&
           (XYZ(2,I,J,K) - XYZ(2,I,J,K-1)) ** 2 +&
           (XYZ(3,I,J,K) - XYZ(3,I,J,K-1)) ** 2)

!     Execution.
!     ----------

!     Zero the three low-end faces (or edges if one plane is specified).

      DO K = 1, KL
         DO J = 1, JL
            S(1,1,J,K) = ZERO
         END DO

         DO I = 1, IL
            S(2,I,1,K) = ZERO
         END DO
      END DO

      DO J = 1, JL
         DO I = 1, IL
            S(3,I,J,1) = ZERO
         END DO
      END DO

!     Set up the low-end edge lines because they are missed by the
!     following loops over most of the low-end faces:

      DO I = 2, IL
         S(1,I,1,1) = S(1,I-1,1,1) + DELI(I,1,1)
      END DO

      DO J = 2, JL
         S(2,1,J,1) = S(2,1,J-1,1) + DELJ(1,J,1)
      END DO

      DO K = 2, KL
         S(3,1,1,K) = S(3,1,1,K-1) + DELK(1,1,K)
      END DO

!     Set up the rest of the low-end face lines because they are
!     missed by the the main loop over most of the volume.

      DO K = 2, KL
         DO J = 2, JL
            S(2,1,J,K) = S(2,1,J-1,K) + DELJ(1,J,K)
            S(3,1,J,K) = S(3,1,J,K-1) + DELK(1,J,K)
         END DO
         DO I = 2, IL
            S(1,I,1,K) = S(1,I-1,1,K) + DELI(I,1,K)
            S(3,I,1,K) = S(3,I,1,K-1) + DELK(I,1,K)
         END DO
      END DO

      DO J = 2, JL
         DO I = 2, IL
            S(1,I,J,1) = S(1,I-1,J,1) + DELI(I,J,1)
            S(2,I,J,1) = S(2,I,J-1,1) + DELJ(I,J,1)
         END DO
      END DO

!     Traverse the block just once for all lines except those within
!     the low-end faces.

      DO K = 2, KL
         DO J = 2, JL
            DO I = 2, IL
               S(1,I,J,K) = S(1,I-1,J,K) + DELI(I,J,K)
               S(2,I,J,K) = S(2,I,J-1,K) + DELJ(I,J,K)
               S(3,I,J,K) = S(3,I,J,K-1) + DELK(I,J,K)
            END DO
         END DO
      END DO

!     Normalizing requires another pass through the volume.
!     Handle lines of zero length first by inserting uniform
!     distributions.  Then the standard normalization can be
!     applied safely everywhere.

      DO K = 1, KL

!        Zero-length lines in the I direction?

         DO J = 1, JL
            IF (S(1,IL,J,K) .EQ. ZERO) THEN
               DO I = 2, IL
                  S(1,I,J,K) = I - 1
               END DO
            END IF
         END DO

!        Zero-length lines in the J direction?

         DO I = 1, IL
            IF (S(2,I,JL,K) .EQ. ZERO) THEN
               DO J = 2, JL
                  S(2,I,J,K) = J - 1
               END DO
            END IF
         END DO
      END DO

!     Zero-length lines in the K direction?

      DO J = 1, JL
         DO I = 1, IL
            IF (S(3,I,J,KL) .EQ. ZERO) THEN
               DO K = 2, KL
                  S(3,I,J,K) = K - 1
               END DO
            END IF
         END DO
      END DO

!     Normalize:

      DO K = 1, KL
         DO J = 1, JL
            DO I = 1, IL
               S(1,I,J,K) = S(1,I,J,K) / S(1,IL,J,K)
               S(2,I,J,K) = S(2,I,J,K) / S(2,I,JL,K)
               S(3,I,J,K) = S(3,I,J,K) / S(3,I,J,KL)
            END DO
         END DO
      END DO

!     Finally, precise 1s for the three high-end faces:

      DO K = 1, KL
         DO J = 1, JL
            S(1,IL,J,K) = ONE
         END DO

         DO I = 1, IL
            S(2,I,JL,K) = ONE
         END DO
      END DO

      DO J = 1, JL
         DO I = 1, IL
            S(3,I,J,KL) = ONE
         END DO
      END DO

    END SUBROUTINE PARAM3DM_SOLID
      SUBROUTINE WARPBLK_solid (IFACEPTB, IEDGEPTB, NCALL, IG, IGO, &
           IL, JL, KL, XYZ0, S0, XYZ)

!     ******************************************************************
!     *   WARP-BLK completes the perturbation of one block of a multi- *
!     *   block grid structure given a base grid and all fixed and all *
!     *   explicitly-perturbed faces in-place for the desired block.   *
!     *                                                                *
!     *   Ancillary routines:   PARAM3DM, WARPQ3DM, DELQ3DM            *
!     *                                                                *
!     *   11/29/95  D.Saunders  Specialized adaptation of WARP3D and   *
!     *                         ancillary routines for FLO107-MB.      *
!     *   01/26/96  DAS/JJR     Edges affected implicitly by corner    *
!     *                         motion had been overlooked.            *
!     *   04/04/96     "        Algorithm is now 3-stage/1-pass, not   *
!     *                         2-stage/2-pass as it was originally.   *
!     *   05/24/96     "        Juan Alonso debugged it for us.        *
!     *   06/??/96    JJR       Implicit edge motion is specified at   *
!     *                         the higher level now via IEDGEPTB(*).  *
!     *   06/19/96    DAS       Handled degenerate edges in PARAM3DM.  *
!     *   12/11/08  C.A.Mader   Converted to *.f90                     *
!     *                                                                *
!     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
!     ******************************************************************

      !IMPLICIT REAL*8 (A-H,O-Z) ! Take out when all compilers have a switch
      use precision

      implicit none
!     Arguments:

      INTEGER(kind=intType)::IFACEPTB(6)    ! I  Controls for faces 1 - 6:
                                            !    0 means no perturbation;
                                            !    1 means implicit perturbation;
                                            !    2 means explicit perturbation
      INTEGER(kind=intType)::IEDGEPTB(12)   ! I  Controls for edges 1 - 12:
                                            !    0 means no perturbation;
                                            !    1 means implicit perturbation;
                                            !    2 means explicit perturbation
      INTEGER(kind=intType)::NCALL          ! I  First call if NCALL = -3
      INTEGER(kind=intType)::IG, IGO        ! I  Design point i.d.
      INTEGER(kind=intType)::IL, JL, KL     ! I  Block dimensions
      REAL(kind=realType)::XYZ0(3,0:IL+1,0:JL+1,0:KL+1)! I  Base block including halo
      REAL(kind=realType)::S0(3,0:IL+1,0:JL+1,0:KL+1)  !I/O Base normalized arc-lengths;
                                            !    halo allows use of XYZ pointers
    !   REAL(kind=realType)::DFACEI(3,JL,KL,2,4)         ! S  For face perturbations; e.g.,
!       REAL(kind=realType)::DFACEJ(3,IL,KL,2,4)         !    DFACEI(1:3,JL,KL,1,M) =
!       REAL(kind=realType)::DFACEK(3,IL,JL,2,4)         !    dX, dY, dZ on the I=1 face for
!                                             !    stage M, etc.; M = 1,2a 2b,3.
!                                             !    each can be (3,MDM,MDM,2,4) if
!                                             !    MDM=MAX(IL,JL,KL) (all blocks)
      REAL(kind=realType)::XYZ(3,0:IL+1,0:JL+1,0:KL+1) !I/O Perturbed block

!     Local constants.

      REAL(kind=realType):: EPS, ONE, ZERO

      PARAMETER (EPS = 1.E-14, ONE = 1.E+0,& ! EPS safeguards a divide by zero -
           ZERO = 0.E+0)              ! presumably only if result is zero.

!     Local variables.

      INTEGER(kind=inttype):: I, I1, I2, J, J1, J2, K, K1, K2, L, LC, LE, M
      integer(kind=inttype):: faceindex
      REAL(kind=realType)::     DEL1, DEL2, DELI, DELJ, DELK,&
           DELIJ, DELIK, DELJI, DELJK, DELKI, DELKJ,&
           WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2
      real(kind=realType):: DX(8), DY(8), DZ(8)         ! Corner perturbations
      real(kind=realType),allocatable,dimension(:,:,:,:,:)::dFaceI,dFaceJ,dFaceK

      ALLOCATE(DFACEI(3,JL,KL,2,4),DFACEJ(3,IL,KL,2,4),&
           DFACEK(3,IL,JL,2,4))

!     Execution.

!     Parameterize the block volume on the first call:

      IF ( (NCALL .LE. -3) .AND. (IG /= IGO) ) THEN
!         write(*,*) 'parameterizing'
         CALL PARAM3DM_solid (IL, JL, KL, XYZ0, S0)
      END IF

!     Calculate the DFACEI,DFACEJ,DFACEK

      DFACEI(1,1:JL,1:KL,1,1) = XYZ(1,1,1:JL,1:KL) - Xyz0(1,1,1:JL,1:KL)
      DFACEI(2,1:JL,1:KL,1,1) = XYZ(2,1,1:JL,1:KL) - Xyz0(2,1,1:JL,1:KL)
      DFACEI(3,1:JL,1:KL,1,1) = XYZ(3,1,1:JL,1:KL) - Xyz0(3,1,1:JL,1:KL)
      
      DFACEI(1,1:JL,1:KL,2,1) = XYZ(1,IL,1:JL,1:KL) - Xyz0(1,IL,1:JL,1:KL)
      DFACEI(2,1:JL,1:KL,2,1) = XYZ(2,IL,1:JL,1:KL) - Xyz0(2,IL,1:JL,1:KL)
      DFACEI(3,1:JL,1:KL,2,1) = XYZ(3,IL,1:JL,1:KL) - Xyz0(3,IL,1:JL,1:KL)

      DFACEJ(1,1:IL,1:KL,1,1) = XYZ(1,1:IL,1,1:KL) - Xyz0(1,1:IL,1,1:KL)
      DFACEJ(2,1:IL,1:KL,1,1) = XYZ(2,1:IL,1,1:KL) - Xyz0(2,1:IL,1,1:KL)
      DFACEJ(3,1:IL,1:KL,1,1) = XYZ(3,1:IL,1,1:KL) - Xyz0(3,1:IL,1,1:KL)

      DFACEJ(1,1:IL,1:KL,2,1) = XYZ(1,1:IL,JL,1:KL) - Xyz0(1,1:IL,JL,1:KL)
      DFACEJ(2,1:IL,1:KL,2,1) = XYZ(2,1:IL,JL,1:KL) - Xyz0(2,1:IL,JL,1:KL)
      DFACEJ(3,1:IL,1:KL,2,1) = XYZ(3,1:IL,JL,1:KL) - Xyz0(3,1:IL,JL,1:KL)

      DFACEK(1,1:IL,1:JL,1,1) = XYZ(1,1:IL,1:JL,1) - Xyz0(1,1:IL,1:JL,1)
      DFACEK(2,1:IL,1:JL,1,1) = XYZ(2,1:IL,1:JL,1) - Xyz0(2,1:IL,1:JL,1)
      DFACEK(3,1:IL,1:JL,1,1) = XYZ(3,1:IL,1:JL,1) - Xyz0(3,1:IL,1:JL,1)

      DFACEK(1,1:IL,1:JL,2,1) = XYZ(1,1:IL,1:JL,KL) - Xyz0(1,1:IL,1:JL,KL)
      DFACEK(2,1:IL,1:JL,2,1) = XYZ(2,1:IL,1:JL,KL) - Xyz0(2,1:IL,1:JL,KL)
      DFACEK(3,1:IL,1:JL,2,1) = XYZ(3,1:IL,1:JL,KL) - Xyz0(3,1:IL,1:JL,KL)

!     Calculate corner point motion:

      LC = 1
      DO K = 1, KL, KL - 1
         DO J = 1, JL, JL - 1
            DO I = 1, IL, IL - 1
               DX(LC) = XYZ(1,I,J,K) - XYZ0(1,I,J,K)
               DY(LC) = XYZ(2,I,J,K) - XYZ0(2,I,J,K)
               DZ(LC) = XYZ(3,I,J,K) - XYZ0(3,I,J,K)
               LC = LC + 1
            END DO
         END DO
      END DO
     
!     Perturb implicit block edges.

!     Edges in the I direction:

      LE = 1
      LC = 1

      DO K = 1, KL, KL - 1
         DO J = 1, JL, JL - 1
            IF (IEDGEPTB(LE) .EQ. 1) THEN
               DO I = 2, IL - 1
                  WTI2 = S0(1,I,J,K)
                  WTI1 = ONE - WTI2
                  XYZ(1,I,J,K)=WTI1*DX(LC)+WTI2*DX(LC+1) +XYZ0(1,I,J,K)
                  XYZ(2,I,J,K)=WTI1*DY(LC)+WTI2*DY(LC+1) +XYZ0(2,I,J,K)
                  XYZ(3,I,J,K)=WTI1*DZ(LC)+WTI2*DZ(LC+1) +XYZ0(3,I,J,K)
               END DO
            END IF
            LE = LE + 1
            LC = LC + 2
         END DO
      END DO

!     Edges in the J direction:

      LE = 5
      LC = 1

      DO K = 1, KL, KL - 1
         DO I = 1, IL, IL - 1
            IF (IEDGEPTB(LE) .EQ. 1) THEN
               DO J = 2, JL - 1
                  WTJ2 = S0(2,I,J,K)
                  WTJ1 = ONE - WTJ2
                  XYZ(1,I,J,K)=WTJ1*DX(LC)+WTJ2*DX(LC+2) +XYZ0(1,I,J,K)
                  XYZ(2,I,J,K)=WTJ1*DY(LC)+WTJ2*DY(LC+2) +XYZ0(2,I,J,K)
                  XYZ(3,I,J,K)=WTJ1*DZ(LC)+WTJ2*DZ(LC+2) +XYZ0(3,I,J,K)
               END DO
            END IF
            LE = LE + 1
            LC = LC + 1
         END DO
         LC = LC + 2
      END DO

!     Edges in the K direction:

      LE = 9 
      LC = 1

      DO J = 1, JL, JL - 1
         DO I = 1, IL, IL - 1
            IF (IEDGEPTB(LE) .EQ. 1 ) THEN
               DO K = 2, KL - 1
                  WTK2 = S0(3,I,J,K)
                  WTK1 = ONE - WTK2
                  XYZ(1,I,J,K)=WTK1*DX(LC)+WTK2*DX(LC+4) +XYZ0(1,I,J,K)
                  XYZ(2,I,J,K)=WTK1*DY(LC)+WTK2*DY(LC+4) +XYZ0(2,I,J,K)
                  XYZ(3,I,J,K)=WTK1*DZ(LC)+WTK2*DZ(LC+4) +XYZ0(3,I,J,K)
               END DO
            END IF
            LE = LE + 1
            LC = LC + 1
         END DO
      END DO

!     Perturb block faces implicitly affected by explicit edge changes:

      DO M = 1, 6

         IF (IFACEPTB(M) .EQ. 1) THEN

            I1 = 1   ! WARPQ3DM expects one pair of indices to be
            I2 = IL  ! equal to define the face
            J1 = 1
            J2 = JL
            K1 = 1
            K2 = KL

            IF (M .EQ. 1) THEN
               I2 = 1
               faceindex = 1
            ELSE IF (M .EQ. 2) THEN
               I1 = IL
               faceindex = 2
            ELSE IF (M .EQ. 3) THEN
               J2 = 1
               faceindex = 1
            ELSE IF (M .EQ. 4) THEN
               J1 = JL
               faceindex = 2
            ELSE IF (M .EQ. 5) THEN
               K2 = 1
               faceindex = 1
            ELSE !  (M .EQ. 6)
               K1 = KL
               faceindex = 2
            END IF

            !CALL WARPQ3DM (IL, JL, KL, I1, I2, J1, J2, K1, K2,&
            !     XYZ0, S0, DFACEI, DFACEJ, DFACEK, XYZ)
            CALL WARPQ3DM_solid (IL, JL, KL, I1, I2, J1, J2, K1, K2,&
                 XYZ0, S0, DFACEI(:,:,:,faceindex,1),&
                 DFACEJ(:,:,:,faceindex,1), &
                 DFACEK(:,:,:,faceindex,1), XYZ)
         END IF
      END DO


!     Perturb the volume grid as though all faces have changed.
!     The face perturbations are stored for all three stages so that
!     the volume points can be perturbed in a single pass.


!     Stage 1:  Corner motion (only). This stage is reusable by WARPQ3DM.
!     --------

!     I = 1 and IL intermediate faces.

      !CALL DELQ3DM (IL, JL, KL,  1,  1, 1, JL, 1, KL, XYZ0, S0,&
      !     DFACEI(1,1,1,1,1), DFACEJ, DFACEK, XYZ)
      CALL DELQ3DM_solid (IL, JL, KL,  1,  1, 1, JL, 1, KL, XYZ0, S0,&
           DFACEI(:,:,:,1,1), DFACEJ(:,:,:,1,1), DFACEK(:,:,:,1,1), XYZ)

      !CALL DELQ3DM (IL, JL, KL, IL, IL, 1, JL, 1, KL, XYZ0, S0,&
      !     DFACEI(1,1,1,2,1), DFACEJ, DFACEK, XYZ)
      CALL DELQ3DM_solid (IL, JL, KL, IL, IL, 1, JL, 1, KL, XYZ0, S0,&
           DFACEI(:,:,:,2,1), DFACEJ(:,:,:,2,1), DFACEK(:,:,:,2,1), XYZ)
      
!     J = 1 and JL intermediate faces.

      !CALL DELQ3DM (IL, JL, KL, 1, IL,  1,  1, 1, KL, XYZ0, S0,&
      !     DFACEI, DFACEJ(1,1,1,1,1), DFACEK, XYZ)
      CALL DELQ3DM_solid (IL, JL, KL, 1, IL,  1,  1, 1, KL, XYZ0, S0,&
           DFACEI(:,:,:,1,1), DFACEJ(:,:,:,1,1), DFACEK(:,:,:,1,1), XYZ)

      !CALL DELQ3DM (IL, JL, KL, 1, IL, JL, JL, 1, KL, XYZ0, S0,&
      !     DFACEI, DFACEJ(1,1,1,2,1), DFACEK, XYZ)
      CALL DELQ3DM_solid (IL, JL, KL, 1, IL, JL, JL, 1, KL, XYZ0, S0,&
           DFACEI(:,:,:,2,1), DFACEJ(:,:,:,2,1), DFACEK(:,:,:,2,1), XYZ)

!     K = 1 and KL intermediate faces.

      !CALL DELQ3DM (IL, JL, KL, 1, IL, 1, JL,  1,  1, XYZ0, S0,&
      !     DFACEI, DFACEJ, DFACEK(1,1,1,1,1), XYZ)
      CALL DELQ3DM_solid (IL, JL, KL, 1, IL, 1, JL,  1,  1, XYZ0, S0,&
           DFACEI(:,:,:,1,1), DFACEJ(:,:,:,1,1), DFACEK(:,:,:,1,1), XYZ)


      !CALL DELQ3DM (IL, JL, KL, 1, IL, 1, JL, KL, KL, XYZ0, S0,&
      !     DFACEI, DFACEJ, DFACEK(1,1,1,2,1), XYZ)
      CALL DELQ3DM_solid (IL, JL, KL, 1, IL, 1, JL, KL, KL, XYZ0, S0,&
           DFACEI(:,:,:,2,1), DFACEJ(:,:,:,2,1), DFACEK(:,:,:,2,1), XYZ)


!     Stage 2:  Handle edge motion from above interim edges to final edges.
!     --------

!     James's insight here: consider the 4 faces affected by the motion of
!     4 edges in a given (index) direction; furthermore, keep the two
!     directions of the resulting face perturbations separated for proper
!     weighted combination during the final interpolation into the interior.
!     The 3 index directions are then independent of each other (added).

!     I = 1 and IL faces:

      L = 1
      DO I = 1, IL, IL - 1

!        K = 1 and KL edge perturbations (J direction edges):

         DO J = 2, JL - 1
            DFACEI(1,J, 1,L,2) = XYZ(1,I,J, 1) - XYZ0(1,I,J, 1) -&
           DFACEI(1,J, 1,L,1)
            DFACEI(2,J, 1,L,2) = XYZ(2,I,J, 1) - XYZ0(2,I,J, 1) -&
           DFACEI(2,J, 1,L,1)
            DFACEI(3,J, 1,L,2) = XYZ(3,I,J, 1) - XYZ0(3,I,J, 1) -&
           DFACEI(3,J, 1,L,1)
            DFACEI(1,J,KL,L,2) = XYZ(1,I,J,KL) - XYZ0(1,I,J,KL) -&
           DFACEI(1,J,KL,L,1)
            DFACEI(2,J,KL,L,2) = XYZ(2,I,J,KL) - XYZ0(2,I,J,KL) -&
           DFACEI(2,J,KL,L,1)
            DFACEI(3,J,KL,L,2) = XYZ(3,I,J,KL) - XYZ0(3,I,J,KL) -&
           DFACEI(3,J,KL,L,1)
         END DO

!        J = 1 and JL edge perturbations (K direction edges):

         DO K = 2, KL - 1
            DFACEI(1, 1,K,L,3) = XYZ(1,I, 1,K) - XYZ0(1,I, 1,K) -&
           DFACEI(1, 1,K,L,1)
            DFACEI(2, 1,K,L,3) = XYZ(2,I, 1,K) - XYZ0(2,I, 1,K) -&
           DFACEI(2, 1,K,L,1)
            DFACEI(3, 1,K,L,3) = XYZ(3,I, 1,K) - XYZ0(3,I, 1,K) -&
           DFACEI(3, 1,K,L,1)
            DFACEI(1,JL,K,L,3) = XYZ(1,I,JL,K) - XYZ0(1,I,JL,K) -&
           DFACEI(1,JL,K,L,1)
            DFACEI(2,JL,K,L,3) = XYZ(2,I,JL,K) - XYZ0(2,I,JL,K) -&
           DFACEI(2,JL,K,L,1)
            DFACEI(3,JL,K,L,3) = XYZ(3,I,JL,K) - XYZ0(3,I,JL,K) -&
           DFACEI(3,JL,K,L,1)
         END DO

!        Interpolate stage 2 interior points for this I face, keeping the
!        two index directions separated.

         DO K = 2, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DFACEI(1,J,K,L,3) = WTJ1*DFACEI(1, 1,K,L,3) +&
                                  WTJ2*DFACEI(1,JL,K,L,3)
               DFACEI(2,J,K,L,3) = WTJ1*DFACEI(2, 1,K,L,3) +&
                                  WTJ2*DFACEI(2,JL,K,L,3)
               DFACEI(3,J,K,L,3) = WTJ1*DFACEI(3, 1,K,L,3) +&
                                  WTJ2*DFACEI(3,JL,K,L,3)

               DFACEI(1,J,K,L,2) = WTK1*DFACEI(1,J, 1,L,2) +&
                                  WTK2*DFACEI(1,J,KL,L,2)
               DFACEI(2,J,K,L,2) = WTK1*DFACEI(2,J, 1,L,2) +&
                                  WTK2*DFACEI(2,J,KL,L,2)
               DFACEI(3,J,K,L,2) = WTK1*DFACEI(3,J, 1,L,2) +&
                                  WTK2*DFACEI(3,J,KL,L,2)
            END DO
         END DO
         L = 2
      END DO

!     J = 1 and JL faces, stage 2:

      L = 1
      DO J = 1, JL, JL - 1

!        K = 1 and KL edge perturbations (I direction):

         DO I = 2, IL - 1
            DFACEJ(1,I, 1,L,2) = XYZ(1,I,J, 1) - XYZ0(1,I,J, 1) -&
           DFACEJ(1,I, 1,L,1)
            DFACEJ(2,I, 1,L,2) = XYZ(2,I,J, 1) - XYZ0(2,I,J, 1) -&
           DFACEJ(2,I, 1,L,1)
            DFACEJ(3,I, 1,L,2) = XYZ(3,I,J, 1) - XYZ0(3,I,J, 1) -&
           DFACEJ(3,I, 1,L,1)
            DFACEJ(1,I,KL,L,2) = XYZ(1,I,J,KL) - XYZ0(1,I,J,KL) -&
           DFACEJ(1,I,KL,L,1)
            DFACEJ(2,I,KL,L,2) = XYZ(2,I,J,KL) - XYZ0(2,I,J,KL) -&
           DFACEJ(2,I,KL,L,1)
            DFACEJ(3,I,KL,L,2) = XYZ(3,I,J,KL) - XYZ0(3,I,J,KL) -&
           DFACEJ(3,I,KL,L,1)
         END DO

!        I = 1 and IL edge perturbations (K direction):

         DO K = 2, KL - 1
            DFACEJ(1, 1,K,L,3) = XYZ(1, 1,J,K) - XYZ0(1, 1,J,K) -&
           DFACEJ(1, 1,K,L,1)
            DFACEJ(2, 1,K,L,3) = XYZ(2, 1,J,K) - XYZ0(2, 1,J,K) -&
           DFACEJ(2, 1,K,L,1)
            DFACEJ(3, 1,K,L,3) = XYZ(3, 1,J,K) - XYZ0(3, 1,J,K) -&
           DFACEJ(3, 1,K,L,1)
            DFACEJ(1,IL,K,L,3) = XYZ(1,IL,J,K) - XYZ0(1,IL,J,K) -&
           DFACEJ(1,IL,K,L,1)
            DFACEJ(2,IL,K,L,3) = XYZ(2,IL,J,K) - XYZ0(2,IL,J,K) -&
           DFACEJ(2,IL,K,L,1)
            DFACEJ(3,IL,K,L,3) = XYZ(3,IL,J,K) - XYZ0(3,IL,J,K) -&
           DFACEJ(3,IL,K,L,1)
         END DO

!        Interpolate stage 2 interior points for this J face.

         DO K = 2, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DFACEJ(1,I,K,L,3) = WTI1*DFACEJ(1, 1,K,L,3) +&
                                  WTI2*DFACEJ(1,IL,K,L,3)
               DFACEJ(2,I,K,L,3) = WTI1*DFACEJ(2, 1,K,L,3) +&
                                  WTI2*DFACEJ(2,IL,K,L,3)
               DFACEJ(3,I,K,L,3) = WTI1*DFACEJ(3, 1,K,L,3) +&
                                  WTI2*DFACEJ(3,IL,K,L,3)

               DFACEJ(1,I,K,L,2) = WTK1*DFACEJ(1,I, 1,L,2) +&
                                  WTK2*DFACEJ(1,I,KL,L,2)
               DFACEJ(2,I,K,L,2) = WTK1*DFACEJ(2,I, 1,L,2) +&
                                  WTK2*DFACEJ(2,I,KL,L,2)
               DFACEJ(3,I,K,L,2) = WTK1*DFACEJ(3,I, 1,L,2) +&
                                  WTK2*DFACEJ(3,I,KL,L,2)
            END DO
         END DO
         L = 2
      END DO

!     K = 1 and KL faces, stage 2:

      L = 1
      DO K = 1, KL, KL - 1

!        J = 1 and JL edge perturbations (I direction):

         DO I = 2, IL - 1
            DFACEK(1,I, 1,L,2) = XYZ(1,I, 1,K) - XYZ0(1,I, 1,K) -&
           DFACEK(1,I, 1,L,1)
            DFACEK(2,I, 1,L,2) = XYZ(2,I, 1,K) - XYZ0(2,I, 1,K) -&
           DFACEK(2,I, 1,L,1)
            DFACEK(3,I, 1,L,2) = XYZ(3,I, 1,K) - XYZ0(3,I, 1,K) -&
           DFACEK(3,I, 1,L,1)
            DFACEK(1,I,JL,L,2) = XYZ(1,I,JL,K) - XYZ0(1,I,JL,K) -&
           DFACEK(1,I,JL,L,1)
            DFACEK(2,I,JL,L,2) = XYZ(2,I,JL,K) - XYZ0(2,I,JL,K) -&
           DFACEK(2,I,JL,L,1)
            DFACEK(3,I,JL,L,2) = XYZ(3,I,JL,K) - XYZ0(3,I,JL,K) -&
           DFACEK(3,I,JL,L,1)
         END DO

!        I = 1 and IL edge perturbations (J direction):

         DO J = 2, JL - 1
            DFACEK(1, 1,J,L,3) = XYZ(1, 1,J,K) - XYZ0(1, 1,J,K) -&
           DFACEK(1, 1,J,L,1)
            DFACEK(2, 1,J,L,3) = XYZ(2, 1,J,K) - XYZ0(2, 1,J,K) -&
           DFACEK(2, 1,J,L,1)
            DFACEK(3, 1,J,L,3) = XYZ(3, 1,J,K) - XYZ0(3, 1,J,K) -&
           DFACEK(3, 1,J,L,1)
            DFACEK(1,IL,J,L,3) = XYZ(1,IL,J,K) - XYZ0(1,IL,J,K) -&
           DFACEK(1,IL,J,L,1)
            DFACEK(2,IL,J,L,3) = XYZ(2,IL,J,K) - XYZ0(2,IL,J,K) -&
           DFACEK(2,IL,J,L,1)
            DFACEK(3,IL,J,L,3) = XYZ(3,IL,J,K) - XYZ0(3,IL,J,K) -&
           DFACEK(3,IL,J,L,1)
         END DO

!        Interpolate stage 2 interior points for this K face.

         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2

               DFACEK(1,I,J,L,3) = WTI1*DFACEK(1, 1,J,L,3) +&
                                  WTI2*DFACEK(1,IL,J,L,3)
               DFACEK(2,I,J,L,3) = WTI1*DFACEK(2, 1,J,L,3) +&
                                  WTI2*DFACEK(2,IL,J,L,3)
               DFACEK(3,I,J,L,3) = WTI1*DFACEK(3, 1,J,L,3) +&
                                  WTI2*DFACEK(3,IL,J,L,3)

               DFACEK(1,I,J,L,2) = WTJ1*DFACEK(1,I, 1,L,2) +&
                                  WTJ2*DFACEK(1,I,JL,L,2)
               DFACEK(2,I,J,L,2) = WTJ1*DFACEK(2,I, 1,L,2) +&
                                  WTJ2*DFACEK(2,I,JL,L,2)
               DFACEK(3,I,J,L,2) = WTJ1*DFACEK(3,I, 1,L,2) +&
                                  WTJ2*DFACEK(3,I,JL,L,2)
            END DO
         END DO
         L = 2
      END DO


!     Stage 3:  Handle face motion from above interim faces to final faces.
!     --------

      DO K = 2, KL - 1
         DO J = 2, JL - 1
            DFACEI(1,J,K,1,4) = XYZ(1, 1,J,K) - XYZ0(1, 1,J,K) -&
              DFACEI(1,J,K,1,1) - DFACEI(1,J,K,1,2) - DFACEI(1,J,K,1,3)
            DFACEI(2,J,K,1,4) = XYZ(2, 1,J,K) - XYZ0(2, 1,J,K) -&
              DFACEI(2,J,K,1,1) - DFACEI(2,J,K,1,2) - DFACEI(2,J,K,1,3)
            DFACEI(3,J,K,1,4) = XYZ(3, 1,J,K) - XYZ0(3, 1,J,K) -&
              DFACEI(3,J,K,1,1) - DFACEI(3,J,K,1,2) - DFACEI(3,J,K,1,3)
            DFACEI(1,J,K,2,4) = XYZ(1,IL,J,K) - XYZ0(1,IL,J,K) -&
              DFACEI(1,J,K,2,1) - DFACEI(1,J,K,2,2) - DFACEI(1,J,K,2,3)
            DFACEI(2,J,K,2,4) = XYZ(2,IL,J,K) - XYZ0(2,IL,J,K) -&
              DFACEI(2,J,K,2,1) - DFACEI(2,J,K,2,2) - DFACEI(2,J,K,2,3)
            DFACEI(3,J,K,2,4) = XYZ(3,IL,J,K) - XYZ0(3,IL,J,K) -&
              DFACEI(3,J,K,2,1) - DFACEI(3,J,K,2,2) - DFACEI(3,J,K,2,3)
         END DO

         DO I = 2, IL - 1
            DFACEJ(1,I,K,1,4) = XYZ(1,I, 1,K) - XYZ0(1,I, 1,K) -&
              DFACEJ(1,I,K,1,1) - DFACEJ(1,I,K,1,2) - DFACEJ(1,I,K,1,3)
            DFACEJ(2,I,K,1,4) = XYZ(2,I, 1,K) - XYZ0(2,I, 1,K) -&
              DFACEJ(2,I,K,1,1) - DFACEJ(2,I,K,1,2) - DFACEJ(2,I,K,1,3)
            DFACEJ(3,I,K,1,4) = XYZ(3,I, 1,K) - XYZ0(3,I, 1,K) -&
              DFACEJ(3,I,K,1,1) - DFACEJ(3,I,K,1,2) - DFACEJ(3,I,K,1,3)
            DFACEJ(1,I,K,2,4) = XYZ(1,I,JL,K) - XYZ0(1,I,JL,K) -&
              DFACEJ(1,I,K,2,1) - DFACEJ(1,I,K,2,2) - DFACEJ(1,I,K,2,3)
            DFACEJ(2,I,K,2,4) = XYZ(2,I,JL,K) - XYZ0(2,I,JL,K) -&
              DFACEJ(2,I,K,2,1) - DFACEJ(2,I,K,2,2) - DFACEJ(2,I,K,2,3)
            DFACEJ(3,I,K,2,4) = XYZ(3,I,JL,K) - XYZ0(3,I,JL,K) -&
              DFACEJ(3,I,K,2,1) - DFACEJ(3,I,K,2,2) - DFACEJ(3,I,K,2,3)
         END DO
      END DO

      DO J = 2, JL - 1
         DO I = 2, IL - 1
            DFACEK(1,I,J,1,4) = XYZ(1,I,J, 1) - XYZ0(1,I,J, 1) -&
              DFACEK(1,I,J,1,1) - DFACEK(1,I,J,1,2) - DFACEK(1,I,J,1,3)
            DFACEK(2,I,J,1,4) = XYZ(2,I,J, 1) - XYZ0(2,I,J, 1) -&
              DFACEK(2,I,J,1,1) - DFACEK(2,I,J,1,2) - DFACEK(2,I,J,1,3)
            DFACEK(3,I,J,1,4) = XYZ(3,I,J, 1) - XYZ0(3,I,J, 1) -&
              DFACEK(3,I,J,1,1) - DFACEK(3,I,J,1,2) - DFACEK(3,I,J,1,3)
            DFACEK(1,I,J,2,4) = XYZ(1,I,J,KL) - XYZ0(1,I,J,KL) -&
              DFACEK(1,I,J,2,1) - DFACEK(1,I,J,2,2) - DFACEK(1,I,J,2,3)
            DFACEK(2,I,J,2,4) = XYZ(2,I,J,KL) - XYZ0(2,I,J,KL) -&
              DFACEK(2,I,J,2,1) - DFACEK(2,I,J,2,2) - DFACEK(2,I,J,2,3)
            DFACEK(3,I,J,2,4) = XYZ(3,I,J,KL) - XYZ0(3,I,J,KL) -&
              DFACEK(3,I,J,2,1) - DFACEK(3,I,J,2,2) - DFACEK(3,I,J,2,3)
         END DO
      END DO


!     Perturb the interior volume points.
!     All stages are performed at once via interpolation from the face
!     perturbations stored for each stage.  Note that the three stages
!     accumulate the contributions from the three subscript directions
!     with varying degrees of independence.

      DO K = 2, KL - 1
         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

!              Stage 1:

               DELI = WTI1*DFACEI (1,J,K,1,1) + WTI2*DFACEI (1,J,K,2,1)
               DELJ = WTJ1*DFACEJ (1,I,K,1,1) + WTJ2*DFACEJ (1,I,K,2,1)
               DELK = WTK1*DFACEK (1,I,J,1,1) + WTK2*DFACEK (1,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)&
                     / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

!              Stage 2:

               DELIJ= WTJ1*DFACEJ (1,I,K,1,2) + WTJ2*DFACEJ (1,I,K,2,2)
               DELIK= WTK1*DFACEK (1,I,J,1,2) + WTK2*DFACEK (1,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /&
                      MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (1,J,K,1,2) + WTI2*DFACEI (1,J,K,2,2)
               DELJK= WTK1*DFACEK (1,I,J,1,3) + WTK2*DFACEK (1,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /&
                      MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (1,J,K,1,3) + WTI2*DFACEI (1,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (1,I,K,1,3) + WTJ2*DFACEJ (1,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /&
                      MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

!              Stage 3:

               DELI = WTI1*DFACEI (1,J,K,1,4) + WTI2*DFACEI (1,J,K,2,4)
               DELJ = WTJ1*DFACEJ (1,I,K,1,4) + WTJ2*DFACEJ (1,I,K,2,4)
               DELK = WTK1*DFACEK (1,I,J,1,4) + WTK2*DFACEK (1,I,J,2,4)

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DEL1 +DEL2 +DELI+DELJ+DELK

!              Repeat all stages for Y perturbations:

               DELI = WTI1*DFACEI (2,J,K,1,1) + WTI2*DFACEI (2,J,K,2,1)
               DELJ = WTJ1*DFACEJ (2,I,K,1,1) + WTJ2*DFACEJ (2,I,K,2,1)
               DELK = WTK1*DFACEK (2,I,J,1,1) + WTK2*DFACEK (2,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)&
                     / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

               DELIJ= WTJ1*DFACEJ (2,I,K,1,2) + WTJ2*DFACEJ (2,I,K,2,2)
               DELIK= WTK1*DFACEK (2,I,J,1,2) + WTK2*DFACEK (2,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /&
                      MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (2,J,K,1,2) + WTI2*DFACEI (2,J,K,2,2)
               DELJK= WTK1*DFACEK (2,I,J,1,3) + WTK2*DFACEK (2,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /&
                      MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (2,J,K,1,3) + WTI2*DFACEI (2,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (2,I,K,1,3) + WTJ2*DFACEJ (2,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /&
                      MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

               DELI = WTI1*DFACEI (2,J,K,1,4) + WTI2*DFACEI (2,J,K,2,4)
               DELJ = WTJ1*DFACEJ (2,I,K,1,4) + WTJ2*DFACEJ (2,I,K,2,4)
               DELK = WTK1*DFACEK (2,I,J,1,4) + WTK2*DFACEK (2,I,J,2,4)

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DEL1 +DEL2 +DELI+DELJ+DELK

!              ... and for Z perturbations:

               DELI = WTI1*DFACEI (3,J,K,1,1) + WTI2*DFACEI (3,J,K,2,1)
               DELJ = WTJ1*DFACEJ (3,I,K,1,1) + WTJ2*DFACEJ (3,I,K,2,1)
               DELK = WTK1*DFACEK (3,I,J,1,1) + WTK2*DFACEK (3,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)&
                     / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

               DELIJ= WTJ1*DFACEJ (3,I,K,1,2) + WTJ2*DFACEJ (3,I,K,2,2)
               DELIK= WTK1*DFACEK (3,I,J,1,2) + WTK2*DFACEK (3,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /&
                      MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (3,J,K,1,2) + WTI2*DFACEI (3,J,K,2,2)
               DELJK= WTK1*DFACEK (3,I,J,1,3) + WTK2*DFACEK (3,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /&
                      MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (3,J,K,1,3) + WTI2*DFACEI (3,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (3,I,K,1,3) + WTJ2*DFACEJ (3,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /&
                      MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

               DELI = WTI1*DFACEI (3,J,K,1,4) + WTI2*DFACEI (3,J,K,2,4)
               DELJ = WTJ1*DFACEJ (3,I,K,1,4) + WTJ2*DFACEJ (3,I,K,2,4)
               DELK = WTK1*DFACEK (3,I,J,1,4) + WTK2*DFACEK (3,I,J,2,4)

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DEL1 +DEL2 +DELI+DELJ+DELK
            END DO
         END DO
      END DO

    END SUBROUTINE WARPBLK_SOLID

      SUBROUTINE WARPQ3DM_solid (IL, JL, KL, I1, I2, J1, J2, K1, K2, XYZ0, S0,&
           DFACEI, DFACEJ, DFACEK, XYZ)

!     ******************************************************************
!     *   WARPQ3DM perturbs the interior of one face of one block of a *
!     *   multiblock grid structure given perturbed edges of that face,*
!     *   which is indicated by one pair of equal index arguments.     *
!     *   (E.g.: I2 = I1 means a face in the J/K subspace.)            *
!     *                                                                *
!     *   The two-stage algorithm uses an intermediate perturbation to *
!     *   account for any corner motion, involving edges derived from  *
!     *   the original edges, then a second perturbation to account    *
!     *   for any differences between the intermediate edges and the   *
!     *   specified new edges.                                         *
!     *                                                                *
!     *   The perturbed edges should be input as edges of the desired  *
!     *   output face.  The original relative arc-length increments in *
!     *   each index direction should also be input.  See PARAM3DM for *
!     *   setting them up in preparation for multiple perturbations.   *
!     *                                                                *
!     *   11/29/95  D.Saunders  Adaptation of WARPQ3D for specialized  *
!     *                         WARP-BLK used by FLO107-MB.            *
!     *   04/04/96      "       DELQ3DM does stage 1 only now (all     *
!     *                         that WARP-BLK needs).                  *
!     *   12/11/08  C.A.Mader   Converted to *.f90                     *
!     *                                                                *
!     *   David Saunders/James Reuther, NASA Ames Research Center, CA. *
!     ******************************************************************

        !IMPLICIT REAL*8 (A-H,O-Z) ! Take out when all compilers have a switch
        use precision      

        implicit none
!     Arguments.

      INTEGER(kind=intType)::   IL, JL, KL                  ! I  Grid array dimensions.
      INTEGER(kind=intType)::   I1, I2, J1, J2, K1, K2      ! I  Define active face,
                                            !    one pair being equal.
      real(kind=realType):: XYZ0(3,0:IL+1,0:JL+1,0:KL+1)! I  Base face coordinates in
                                            !    appropriate places
      real(kind=realType):: S0(3,0:IL+1,0:JL+1,0:KL+1)  ! I  Base normalized arc-lengths
      real(kind=realType):: DFACEI(3,JL,KL)             ! S  For face perturbations; e.g.,
      real(kind=realType):: DFACEJ(3,IL,KL)             !    DFACEI(1:3,1:JL,1:KL) =
      real(kind=realType):: DFACEK(3,IL,JL)             !    dX, dY, dZ for an I face, etc.
      real(kind=realType):: XYZ(3,0:IL+1,0:JL+1,0:KL+1) !I/O Grid coordinates: new edges
                                            !    of a face in; full new face
                                            !    out
!     Local constants.

      real(kind=realType)::     ONE
      PARAMETER (ONE = 1.E+0)

!     Local variables.

      INTEGER(kind=intType)::   I, J, K
      real(kind=realType):: WTJ2,WTJ1,WTK2,WTK1,WTI2,WTI1,DELI,DELJ,DELK


!     Execution.
!     ----------

!     Stage 1:
!     Handle any corner motion by generating an intermediate face with
!     the final corners but otherwise derived from the original edges.
!     Actually, just set up the appropriate face perturbations.

      CALL DELQ3DM_solid (IL, JL, KL, I1, I2, J1, J2, K1, K2, XYZ0, S0,&
                   DFACEI, DFACEJ, DFACEK, XYZ)


!     Stage 2:
!     Set up the perturbations from the intermediate edges to the final
!     edges, then interpolate them into the interior points.

      IF (I1 .EQ. I2) THEN          ! I plane case:

         I = I1

!        J = 1 and JL edge perturbations:

         DO K = 2, KL - 1
            DFACEI(1, 1,K) = XYZ(1,I, 1,K)-XYZ0(1,I, 1,K)-DFACEI(1, 1,K)
            DFACEI(2, 1,K) = XYZ(2,I, 1,K)-XYZ0(2,I, 1,K)-DFACEI(2, 1,K)
            DFACEI(3, 1,K) = XYZ(3,I, 1,K)-XYZ0(3,I, 1,K)-DFACEI(3, 1,K)
            DFACEI(1,JL,K) = XYZ(1,I,JL,K)-XYZ0(1,I,JL,K)-DFACEI(1,JL,K)
            DFACEI(2,JL,K) = XYZ(2,I,JL,K)-XYZ0(2,I,JL,K)-DFACEI(2,JL,K)
            DFACEI(3,JL,K) = XYZ(3,I,JL,K)-XYZ0(3,I,JL,K)-DFACEI(3,JL,K)
         END DO

!        K = 1 and KL edge perturbations:

         DO J = 2, JL - 1
            DFACEI(1,J, 1) = XYZ(1,I,J, 1)-XYZ0(1,I,J, 1)-DFACEI(1,J, 1)
            DFACEI(2,J, 1) = XYZ(2,I,J, 1)-XYZ0(2,I,J, 1)-DFACEI(2,J, 1)
            DFACEI(3,J, 1) = XYZ(3,I,J, 1)-XYZ0(3,I,J, 1)-DFACEI(3,J, 1)
            DFACEI(1,J,KL) = XYZ(1,I,J,KL)-XYZ0(1,I,J,KL)-DFACEI(1,J,KL)
            DFACEI(2,J,KL) = XYZ(2,I,J,KL)-XYZ0(2,I,J,KL)-DFACEI(2,J,KL)
            DFACEI(3,J,KL) = XYZ(3,I,J,KL)-XYZ0(3,I,J,KL)-DFACEI(3,J,KL)
         END DO

!        Interior points: accumulate the (independent) contributions.

         DO K = 2, KL - 1
            DO J = 2, JL - 1
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI(1,1,K) + WTJ2 * DFACEI(1,JL,K)
               DELK = WTK1 * DFACEI(1,J,1) + WTK2 * DFACEI(1,J,KL)

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DFACEI(1,J,K) + DELJ +DELK

               DELJ = WTJ1 * DFACEI(2,1,K) + WTJ2 * DFACEI(2,JL,K)
               DELK = WTK1 * DFACEI(2,J,1) + WTK2 * DFACEI(2,J,KL)

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DFACEI(2,J,K) + DELJ +DELK

               DELJ = WTJ1 * DFACEI(3,1,K) + WTJ2 * DFACEI(3,JL,K)
               DELK = WTK1 * DFACEI(3,J,1) + WTK2 * DFACEI(3,J,KL)

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DFACEI(3,J,K) + DELJ +DELK
            END DO
         END DO

      ELSE IF (J1 .EQ. J2) THEN     ! J plane case:

         J = J1

!        I = 1 and IL edge perturbations:

         DO K = 2, KL - 1
            DFACEJ(1, 1,K) = XYZ(1, 1,J,K)-XYZ0(1, 1,J,K)-DFACEJ(1, 1,K)
            DFACEJ(2, 1,K) = XYZ(2, 1,J,K)-XYZ0(2, 1,J,K)-DFACEJ(2, 1,K)
            DFACEJ(3, 1,K) = XYZ(3, 1,J,K)-XYZ0(3, 1,J,K)-DFACEJ(3, 1,K)
            DFACEJ(1,IL,K) = XYZ(1,IL,J,K)-XYZ0(1,IL,J,K)-DFACEJ(1,IL,K)
            DFACEJ(2,IL,K) = XYZ(2,IL,J,K)-XYZ0(2,IL,J,K)-DFACEJ(2,IL,K)
            DFACEJ(3,IL,K) = XYZ(3,IL,J,K)-XYZ0(3,IL,J,K)-DFACEJ(3,IL,K)
         END DO

!        K = 1 and KL edge perturbations:

         DO I = 2, IL - 1
            DFACEJ(1,I, 1) = XYZ(1,I,J, 1)-XYZ0(1,I,J, 1)-DFACEJ(1,I, 1)
            DFACEJ(2,I, 1) = XYZ(2,I,J, 1)-XYZ0(2,I,J, 1)-DFACEJ(2,I, 1)
            DFACEJ(3,I, 1) = XYZ(3,I,J, 1)-XYZ0(3,I,J, 1)-DFACEJ(3,I, 1)
            DFACEJ(1,I,KL) = XYZ(1,I,J,KL)-XYZ0(1,I,J,KL)-DFACEJ(1,I,KL)
            DFACEJ(2,I,KL) = XYZ(2,I,J,KL)-XYZ0(2,I,J,KL)-DFACEJ(2,I,KL)
            DFACEJ(3,I,KL) = XYZ(3,I,J,KL)-XYZ0(3,I,J,KL)-DFACEJ(3,I,KL)
         END DO

!        Interior points: accumulate the (independent) contributions.

         DO K = 2, KL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTK2 = S0(3,I,J,K)
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ(1,1,K) + WTI2 * DFACEJ(1,IL,K)
               DELK = WTK1 * DFACEJ(1,I,1) + WTK2 * DFACEJ(1,I,KL)

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DFACEJ(1,I,K) + DELI +DELK

               DELI = WTI1 * DFACEJ(2,1,K) + WTI2 * DFACEJ(2,IL,K)
               DELK = WTK1 * DFACEJ(2,I,1) + WTK2 * DFACEJ(2,I,KL)

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DFACEJ(2,I,K) + DELI +DELK
               DELI = WTI1 * DFACEJ(3,1,K) + WTI2 * DFACEJ(3,IL,K)
               DELK = WTK1 * DFACEJ(3,I,1) + WTK2 * DFACEJ(3,I,KL)

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DFACEJ(3,I,K) + DELI +DELK
            END DO
         END DO

      ELSE IF (K1 .EQ. K2) THEN     ! K plane case:

         K = K1

!        I = 1 and IL edge perturbations:

         DO J = 2, JL - 1
            DFACEK(1, 1,J) = XYZ(1, 1,J,K)-XYZ0(1, 1,J,K)-DFACEK(1, 1,J)
            DFACEK(2, 1,J) = XYZ(2, 1,J,K)-XYZ0(2, 1,J,K)-DFACEK(2, 1,J)
            DFACEK(3, 1,J) = XYZ(3, 1,J,K)-XYZ0(3, 1,J,K)-DFACEK(3, 1,J)
            DFACEK(1,IL,J) = XYZ(1,IL,J,K)-XYZ0(1,IL,J,K)-DFACEK(1,IL,J)
            DFACEK(2,IL,J) = XYZ(2,IL,J,K)-XYZ0(2,IL,J,K)-DFACEK(2,IL,J)
            DFACEK(3,IL,J) = XYZ(3,IL,J,K)-XYZ0(3,IL,J,K)-DFACEK(3,IL,J)
         END DO

!        J = 1 and JL edge perturbations:

         DO I = 2, IL - 1
            DFACEK(1,I, 1) = XYZ(1,I, 1,K)-XYZ0(1,I, 1,K)-DFACEK(1,I, 1)
            DFACEK(2,I, 1) = XYZ(2,I, 1,K)-XYZ0(2,I, 1,K)-DFACEK(2,I, 1)
            DFACEK(3,I, 1) = XYZ(3,I, 1,K)-XYZ0(3,I, 1,K)-DFACEK(3,I, 1)
            DFACEK(1,I,JL) = XYZ(1,I,JL,K)-XYZ0(1,I,JL,K)-DFACEK(1,I,JL)
            DFACEK(2,I,JL) = XYZ(2,I,JL,K)-XYZ0(2,I,JL,K)-DFACEK(2,I,JL)
            DFACEK(3,I,JL) = XYZ(3,I,JL,K)-XYZ0(3,I,JL,K)-DFACEK(3,I,JL)
         END DO

!        Interior points: accumulate the (independent) contributions.

         DO J = 2, JL - 1
            DO I = 2, IL - 1
               WTI2 = S0(1,I,J,K)
               WTI1 = ONE - WTI2
               WTJ2 = S0(2,I,J,K)
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK(1,1,J) + WTI2 * DFACEK(1,IL,J)
               DELJ = WTJ1 * DFACEK(1,I,1) + WTJ2 * DFACEK(1,I,JL)

               XYZ(1,I,J,K) = XYZ0(1,I,J,K) + DFACEK(1,I,J) + DELI +DELJ

               DELI = WTI1 * DFACEK(2,1,J) + WTI2 * DFACEK(2,IL,J)
               DELJ = WTJ1 * DFACEK(2,I,1) + WTJ2 * DFACEK(2,I,JL)

               XYZ(2,I,J,K) = XYZ0(2,I,J,K) + DFACEK(2,I,J) + DELI +DELJ

               DELI = WTI1 * DFACEK(3,1,J) + WTI2 * DFACEK(3,IL,J)
               DELJ = WTJ1 * DFACEK(3,I,1) + WTJ2 * DFACEK(3,I,JL)

               XYZ(3,I,J,K) = XYZ0(3,I,J,K) + DFACEK(3,I,J) + DELI +DELJ
            END DO
         END DO

      END IF

    END SUBROUTINE WARPQ3DM_SOLID
