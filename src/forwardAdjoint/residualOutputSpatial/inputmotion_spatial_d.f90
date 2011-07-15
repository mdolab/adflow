   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !      ==================================================================
   MODULE INPUTMOTION_SPATIAL_D
   USE PRECISION
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Input parameters which are related to the rigid body motion of *
   !      * the entire mesh, i.e. translation and rotation.                *
   !      * These parameters can only be specified for an external flow    *
   !      * computation.                                                   *
   !      *                                                                *
   !      ******************************************************************
   !
   SAVE 
   ! rotPoint(3): Rotation point of the rigid body rotation.
   REAL(kind=realtype), DIMENSION(3) :: rotpoint
   ! degreePolXRot: Degree of the x-rotation polynomial.
   ! degreePolYRot: Degree of the y-rotation polynomial.
   ! degreePolZRot: Degree of the z-rotation polynomial.
   INTEGER(kind=inttype) :: degreepolxrot
   INTEGER(kind=inttype) :: degreepolyrot
   INTEGER(kind=inttype) :: degreepolzrot
   ! coefPolXRot(0:): coefficients of the x-rotation polynomial.
   ! coefPolYRot(0:): coefficients of the y-rotation polynomial.
   ! coefPolZRot(0:): coefficients of the z-rotation polynomial.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolxrot
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolyrot
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolzrot
   ! degreeFourXRot: Degree of the x-rotation fourier series.
   ! degreeFourYRot: Degree of the y-rotation fourier series.
   ! degreeFourZRot: Degree of the z-rotation fourier series.
   INTEGER(kind=inttype) :: degreefourxrot
   INTEGER(kind=inttype) :: degreefouryrot
   INTEGER(kind=inttype) :: degreefourzrot
   ! omegaFourXRot: Fourier frequency of the x-rotation; the
   !                   period of the motion is 2*pi/omega.
   ! omegaFourYRot: Fourier frequency of the y-rotation.
   ! omegaFourZRot: Fourier frequency of the z-rotation.
   REAL(kind=realtype) :: omegafourxrot, omegafourxrotb
   REAL(kind=realtype) :: omegafouryrot, omegafouryrotb
   REAL(kind=realtype) :: omegafourzrot, omegafourzrotb
   ! cosCoefFourXRot(0:): cosine coefficients of the
   !                      x-rotation fourier series.
   ! cosCoefFourYRot(0:): cosine coefficients of the
   !                      y-rotation fourier series.
   ! cosCoefFourZRot(0:): cosine coefficients of the
   !                      z-rotation fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffourxrot
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffouryrot
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffourzrot
   ! sinCoefFourXRot(1:): sine coefficients of the
   !                      x-rotation fourier series.
   ! sinCoefFourYRot(1:): sine coefficients of the
   !                      y-rotation fourier series.
   ! sinCoefFourZRot(1:): sine coefficients of the
   !                      z-rotation fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffourxrot
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffouryrot
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffourzrot
   ! degreePolAlpha: Degree of the Alpha polynomial.
   INTEGER(kind=inttype) :: degreepolalpha
   ! coefPolAlpha(0:): coefficients of the Alpha polynomial.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolalpha
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolalphab
   ! degreeFourAlpha: Degree of the Alpha fourier series.
   INTEGER(kind=inttype) :: degreefouralpha
   ! omegaFourAlpha: Fourier frequency of the Alpha; the
   !                   period of the motion is 2*pi/omega.
   REAL(kind=realtype) :: omegafouralpha, omegafouralphab
   ! cosCoefFourAlpha(0:): cosine coefficients of the
   !                      x-rotation fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffouralpha
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffouralphab
   ! sinCoefFourAlpha(1:): sine coefficients of the
   !                      Alpha fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffouralpha
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffouralphab
   ! degreePolXRot: Degree of the Beta polynomial.
   INTEGER(kind=inttype) :: degreepolbeta
   ! coefPolXRot(0:): coefficients of the Beta polynomial.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolbeta
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolbetab
   ! degreeFourBeta: Degree of the Beta fourier series.
   INTEGER(kind=inttype) :: degreefourbeta
   ! omegaFourBeta: Fourier frequency of the Beta; the
   !                   period of the motion is 2*pi/omega.
   REAL(kind=realtype) :: omegafourbeta, omegafourbetab
   ! cosCoefFourBeta(0:): cosine coefficients of the
   !                      Beta fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffourbeta
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffourbetab
   ! sinCoefFourBeta(1:): sine coefficients of the
   !                      Beta fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffourbeta
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffourbetab
   ! degreePolMach: Degree of the Mach polynomial.
   INTEGER(kind=inttype) :: degreepolmach
   ! coefPolMach(0:): coefficients of the Mach polynomial.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolmach
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coefpolmachb
   ! degreeFourMach: Degree of the Mach fourier series.
   INTEGER(kind=inttype) :: degreefourmach
   ! omegaFourMach: Fourier frequency of the Mach Number; the
   !                   period of the motion is 2*pi/omega.
   REAL(kind=realtype) :: omegafourmach, omegafourmachb
   ! cosCoefFourMach(0:): cosine coefficients of the
   !                      Mach Number fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffourmach
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: coscoeffourmachb
   ! sinCoefFourMach(1:): sine coefficients of the
   !                      Mach Number fourier series.
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffourmach
   REAL(kind=realtype), DIMENSION(:), ALLOCATABLE :: sincoeffourmachb
   ! gridMotionSpecified: Whether or not a rigid body motion of
   !                      the grid has been specified.
   LOGICAL :: gridmotionspecified
   END MODULE INPUTMOTION_SPATIAL_D
