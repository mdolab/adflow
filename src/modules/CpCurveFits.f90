!
!      ******************************************************************
!      *                                                                *
!      * File:          CpCurveFits.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-11-2003                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module CpCurveFits
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the curve fit data for Cp, or better      *
!      * Cp/R, as a function of the temperature. The temperature is     *
!      * assumed to be in Kelvin.                                       *
!      *                                                                *
!      ******************************************************************
!      
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type CpTempFitType, which   *
!      * stores the curve fit data for a certain temperature range.     *
!      *                                                                *
!      ******************************************************************
!
       type CpTempFitType

         ! nterm:     Number of terms in the polynomial expansion.
         ! exponents: The powers of the polynomial fit.
         ! constants: The constants in front of each term.

         integer(kind=intType) :: nTerm
         integer(kind=intType), dimension(:), pointer :: exponents
         real(kind=realType),   dimension(:), pointer :: constants

         ! Additional constant to compute the internal energy.

         real(kind=realType) :: eint0

         ! Values of the integrand of Cp/(R*T) at the lower (_1) and
         ! upper (_2) curve fit boundary. Needed to compute the total
         ! pressure.

         real(kind=realType) :: intCpOvRT_1, intCpOvRT_2

       end type CpTempFitType

       ! CpNParts:             Number of temperature ranges.
       ! CpTRange(0:CpNParts): The temperature ranges.
       ! CpEint(0:CpNParts):   Internal energies at the curve fit
       !                       boundaries.
       ! CpHint(0:CpNParts):   Internal enthalpies at the curve fit
       !                       boundaries.
       ! CpTempFit(CpNParts):  The actual curve fit data.

       integer(kind=intType) :: CpNParts

       real(kind=realType), dimension(:), allocatable :: CpTRange
       real(kind=realType), dimension(:), allocatable :: CpEint
       real(kind=realType), dimension(:), allocatable :: CpHint
       type(cpTempFitType), dimension(:), allocatable :: CpTempFit

       ! cv0:   The cv value at the the temperature CpTRange(0). If the
       !        temperature is lower than the lowest curve fit boundary
       !        this value is needed to extrapolate the energy (assuming
       !        constant cv).
       ! cvn:   Idem, but than at CpTRange(CpNParts). Used when the
       !        temperature is higher than the highest curve fit boundary.

       real(kind=realType) :: cv0, cvn

       end module CpCurveFits
