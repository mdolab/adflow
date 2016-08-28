module turbCurveFits
contains
  function curveUpRe(Re)
    !
    !       curveUpRe determines the value of the nonDimensional           
    !       tangential velocity (made nonDimensional with the skin         
    !       friction velocity) for the given Reynolds number.              
    !       This data has been curve fitted with cubic splines.            
    !
    use paramTurb
    implicit none
    !
    !      Function type.
    !
    real(kind=realType) :: curveUpRe
    !
    !      Function arguments.
    !
    real(kind=realType), intent(in) :: Re
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, nn, start
    real(kind=realType)   :: x, x2, x3, upRe

    ! Determine the situation we are dealing with.

    if(Re <= reT(0)) then

       ! Reynolds number is less than the smallest number in the curve
       ! fit. Use extrapolation.

       x     = sqrt(Re/reT(0))
       upRe = x*up0(1)

    else if(Re >= reT(nFit)) then

       ! Reynolds number is larger than the largest number in the curve
       ! fit. Set upRe to the largest value available.

       nn = nFit
       x  = reT(nn) - reT(nn-1)
       x2 = x*x
       x3 = x*x2

       upRe = up0(nn) + up1(nn)*x + up2(nn)*x2 + up3(nn)*x3

    else

       ! Reynolds number is in the range of the curve fits.
       ! First find the correct interval.

       ii    = nFit
       start = 1
       interval: do

          ! Next guess for the interval.

          nn = start + ii/2

          ! Determine the situation we are having here.

          if(Re > reT(nn)) then

             ! Reynoldls number is larger than the upper boundary of
             ! the current interval. Update the lower boundary.

             start = nn + 1
             ii    = ii - 1

          else if(Re >= reT(nn-1)) then

             ! This is the correct range. Exit the do-loop.

             exit

          endif

          ! Modify ii for the next branch to search.

          ii = ii/2

       enddo interval

       ! Compute upRe using the cubic polynomial for this interval.

       x  = Re - reT(nn-1)
       x2 = x*x
       x3 = x*x2

       upRe = up0(nn) + up1(nn)*x + up2(nn)*x2 + up3(nn)*x3

    endif

    ! And set the function value.

    curveUpRe = upRe

  end function curveUpRe
  !
  !      ==================================================================
  !
  subroutine curveTupYp(tup, yp, ntu1, ntu2)
    !
    !       CurveTupYp determines the value of the turbulent variables     
    !       ntu1 to ntu2 for the given yplus.                              
    !       This data has been curve fitted with cubic splines.            
    !
    use constants
    use inputPhysics
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: ntu1, ntu2
    real(kind=realType),   intent(in) :: yp

    real(kind=realType), dimension(ntu1:ntu2), intent(out) :: tup
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, nn, start, mm
    real(kind=realType)   :: x, x2, x3, epsWall, fWall

    ! Determine the situation we are dealing with.

    if(yp <= ypT(0)) then

       ! Yplus is less than the smallest number in the curve
       ! fit. The treatment is turbulence model dependent.

       select case(turbModel)

       case (spalartAllmaras, spalartAllmarasEdwards)

          ! Transport variable is zero on the wall. Use linear
          ! interpolation.

          x = yp/ypT(0)
          do mm=ntu1,ntu2
             tup(mm) = x*tup0(1,mm)
          enddo

          !=============================================================

       case (komegaWilcox, komegaModified, menterSST)

          ! Use the near wall expressions for k and omega.

          x = yp/ypT(0)
          do mm=ntu1,ntu2
             select case(mm)
             case (itu1)
                if( tuLogFit(mm) ) then
                   tup(mm) = exp(tup0(1,mm))*(x**3.23_realType)
                else
                   tup(mm) = tup0(1,mm)*(x**3.23_realType)
                endif

             case (itu2)
                if( tuLogFit(mm) ) then
                   tup(mm) = exp(tup0(1,mm))/(max(x,eps)**2)
                else
                   tup(mm) = tup0(1,mm)/(max(x,eps)**2)
                endif
             end select
          enddo

          !=============================================================

       case (ktau)

          ! Use the near wall expressions for k and tau.

          x = yp/ypT(0)
          do mm=ntu1,ntu2
             select case(mm)
             case (itu1)
                if( tuLogFit(mm) ) then
                   tup(mm) = exp(tup0(1,mm))*(x**3.23_realType)
                else
                   tup(mm) = tup0(1,mm)*(x**3.23_realType)
                endif

             case (itu2)
                if( tuLogFit(mm) ) then
                   tup(mm) = exp(tup0(1,mm))*x*x
                else
                   tup(mm) = tup0(1,mm)*x*x
                endif
             end select
          enddo

          !=============================================================

       case (v2f)

          ! Use the near wall expressions for k, epsilon, v2 and f.

          x = yp/ypT(0)
          do mm=ntu1,ntu2
             select case(mm)
             case (itu1)
                if( tuLogFit(mm) ) then
                   tup(mm) = exp(tup0(1,mm))*x**2
                else
                   tup(mm) = tup0(1,mm)*x**2
                endif

             case (itu2)  ! epsilon cannot be fitted logarithmically.
                if( tuLogFit(mm) ) then
                   call terminate(&
                        "curveTupYp", &
                        "Check curveFit, epsilon cannot be fitted with log")
                else
                   if(rvfN == 1) epsWall = 0.33_realType
                   if(rvfN == 6) epsWall = 0.27_realType
                   tup(mm) = epsWall + (tup0(1,mm)-epsWall)*x
                endif

             case (itu3)
                if( tuLogFit(mm) ) then
                   tup(mm) = exp(tup0(1,mm))*x**4
                else
                   tup(mm) = tup0(1,mm)*x**4
                endif

             case (itu4)
                if( tuLogFit(mm) ) then
                   if(rvfN == 1) &
                        call terminate(&
                        "curveTupYp", &
                        "Check curveFit, f cannot be fitted with log")
                   if(rvfN == 6) tup(mm) = exp(tup(mm))*x
                else
                   if(rvfN == 1) fWall =-0.0035_realType
                   if(rvfN == 6) fWall = zero
                   tup(mm) = fWall + (tup0(1,mm)-fWall)*x
                endif

             case (itu5)
                if( tuLogFit(mm) ) then
                   tup(mm) = exp(tup(mm))*x**4
                else
                   tup(mm) = tup0(1,mm)*x**4
                endif
             end select
          enddo

       end select

       !=================================================================

    else if(yp >= ypT(nFit)) then

       ! Yplus is larger than the largest number in the curve
       ! fit. Set tup to the largest value available.

       nn = nFit
       x  = ypT(nn) - ypT(nn-1)
       x2 = x*x
       x3 = x*x2

       do mm=ntu1,ntu2
          tup(mm) = tup0(nn,mm)    + tup1(nn,mm)*x &
               + tup2(nn,mm)*x2 + tup3(nn,mm)*x3
          if( tuLogFit(mm) ) tup(mm) = exp(tup(mm))
       enddo

       !=================================================================

    else

       ! y-plus is in the range of the curve fits.
       ! First find the correct interval.

       ii    = nFit
       start = 1
       interval: do

          ! Next guess for the interval.

          nn = start + ii/2

          ! Determine the situation we are having here.

          if(yp > ypT(nn)) then

             ! Yplus is larger than the upper boundary of
             ! the current interval. Update the lower boundary.

             start = nn + 1
             ii    = ii - 1

          else if(yp >= ypT(nn-1)) then

             ! This is the correct range. Exit the do-loop.

             exit

          endif

          ! Modify ii for the next branch to search.

          ii = ii/2

       enddo interval

       ! Compute tup using the cubic polynomial for this interval.

       x  = yp - ypT(nn-1)
       x2 = x*x
       x3 = x*x2

       do mm=ntu1,ntu2
          tup(mm) = tup0(nn,mm)    + tup1(nn,mm)*x &
               + tup2(nn,mm)*x2 + tup3(nn,mm)*x3
          if( tuLogFit(mm) ) tup(mm) = exp(tup(mm))
       enddo

    endif

  end subroutine curveTupYp

  subroutine initCurveFitDataKtau
    !
    !       initCurveFitDataKtau contains the curve fit constants for      
    !       the wall function data for the k-tau turbulence model.         
    !
    use flowVarRefState
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    ! integer :: ierr

    call terminate("initCurveFitDataKtau", &
         "Not implemented yet")

  end subroutine initCurveFitDataKtau

  subroutine initCurveFitDataKw
    !
    !       initCurveFitDataKw contains the curve fit constants for        
    !       the wall function data for the standard Wilcox k-omega model.  
    !
    use constants
    use flowVarRefState
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    ! Set the number of data points and allocate the memory for the
    ! arrays of the curve fits.

    nFit = 34

    allocate(ypT(0:nFit), reT(0:nFit),                     &
         up0(nFit),   up1(nFit), up2(nFit), up3(nFit), &
         tup0(nFit,nt1:nt2), tup1(nFit,nt1:nt2),       &
         tup2(nFit,nt1:nt2), tup3(nFit,nt1:nt2),       &
         tuLogFit(nt1:nt2), stat=ierr)
    if(ierr /= 0)                              &
         call terminate("initCurveFitDataKw", &
         "Memory allocation failure for curve fit &
         &coefficients")

    ! Set the values of the Reynolds numbers at interval boundaries.

    reT(0)  = 0.12529547e+00_realType
    reT(1)  = 0.44996057e+00_realType
    reT(2)  = 0.11581311e+01_realType
    reT(3)  = 0.25353238e+01_realType
    reT(4)  = 0.50446282e+01_realType
    reT(5)  = 0.94194631e+01_realType
    reT(6)  = 0.16766555e+02_realType
    reT(7)  = 0.28556753e+02_realType
    reT(8)  = 0.46274930e+02_realType
    reT(9)  = 0.71021000e+02_realType
    reT(10) = 0.10383163e+03_realType
    reT(11) = 0.14621738e+03_realType
    reT(12) = 0.20028019e+03_realType
    reT(13) = 0.26868298e+03_realType
    reT(14) = 0.35467049e+03_realType
    reT(15) = 0.46212508e+03_realType
    reT(16) = 0.59566097e+03_realType
    reT(17) = 0.76073076e+03_realType
    reT(18) = 0.96373333e+03_realType
    reT(19) = 0.12121761e+04_realType
    reT(20) = 0.15147917e+04_realType
    reT(21) = 0.18817196e+04_realType
    reT(22) = 0.23247121e+04_realType
    reT(23) = 0.28572322e+04_realType
    reT(24) = 0.34947840e+04_realType
    reT(25) = 0.42551444e+04_realType
    reT(26) = 0.51584529e+04_realType
    reT(27) = 0.62277581e+04_realType
    reT(28) = 0.74889831e+04_realType
    reT(29) = 0.89716314e+04_realType
    reT(30) = 0.10708764e+05_realType
    reT(31) = 0.12737815e+05_realType
    reT(32) = 0.15100490e+05_realType
    reT(33) = 0.17843939e+05_realType
    reT(34) = 0.21020534e+05_realType

    ! Set the values of the y+ values at interval boundaries.

    ypT(0)  = 0.35397100e+00_realType
    ypT(1)  = 0.67079200e+00_realType
    ypT(2)  = 0.10761700e+01_realType
    ypT(3)  = 0.15923000e+01_realType
    ypT(4)  = 0.22463000e+01_realType
    ypT(5)  = 0.30710600e+01_realType
    ypT(6)  = 0.41063800e+01_realType
    ypT(7)  = 0.54001200e+01_realType
    ypT(8)  = 0.70095900e+01_realType
    ypT(9)  = 0.90031400e+01_realType
    ypT(10) = 0.11461900e+02_realType
    ypT(11) = 0.14481700e+02_realType
    ypT(12) = 0.18175400e+02_realType
    ypT(13) = 0.22675200e+02_realType
    ypT(14) = 0.28135500e+02_realType
    ypT(15) = 0.34735800e+02_realType
    ypT(16) = 0.42683800e+02_realType
    ypT(17) = 0.52219300e+02_realType
    ypT(18) = 0.63617800e+02_realType
    ypT(19) = 0.77194900e+02_realType
    ypT(20) = 0.93310400e+02_realType
    ypT(21) = 0.11237300e+03_realType
    ypT(22) = 0.13484800e+03_realType
    ypT(23) = 0.16125700e+03_realType
    ypT(24) = 0.19218900e+03_realType
    ypT(25) = 0.22830600e+03_realType
    ypT(26) = 0.27034500e+03_realType
    ypT(27) = 0.31913000e+03_realType
    ypT(28) = 0.37557400e+03_realType
    ypT(29) = 0.44069100e+03_realType
    ypT(30) = 0.51559800e+03_realType
    ypT(31) = 0.60152700e+03_realType
    ypT(32) = 0.69982900e+03_realType
    ypT(33) = 0.81198500e+03_realType
    ypT(34) = 0.93960800e+03_realType

    ! Set the values of constants for the cubic fits of the
    ! non-dimensional tangential velocity.

    up0(1)  = 0.35397100e+00_realType
    up0(2)  = 0.67079000e+00_realType
    up0(3)  = 0.10761600e+01_realType
    up0(4)  = 0.15922400e+01_realType
    up0(5)  = 0.22457500e+01_realType
    up0(6)  = 0.30671700e+01_realType
    up0(7)  = 0.40830500e+01_realType
    up0(8)  = 0.52881700e+01_realType
    up0(9)  = 0.66016600e+01_realType
    up0(10) = 0.78884700e+01_realType
    up0(11) = 0.90588500e+01_realType
    up0(12) = 0.10096700e+02_realType
    up0(13) = 0.11019300e+02_realType
    up0(14) = 0.11849200e+02_realType
    up0(15) = 0.12605800e+02_realType
    up0(16) = 0.13304000e+02_realType
    up0(17) = 0.13955200e+02_realType
    up0(18) = 0.14568000e+02_realType
    up0(19) = 0.15148800e+02_realType
    up0(20) = 0.15702800e+02_realType
    up0(21) = 0.16233900e+02_realType
    up0(22) = 0.16745300e+02_realType
    up0(23) = 0.17239500e+02_realType
    up0(24) = 0.17718500e+02_realType
    up0(25) = 0.18184100e+02_realType
    up0(26) = 0.18637900e+02_realType
    up0(27) = 0.19081000e+02_realType
    up0(28) = 0.19514800e+02_realType
    up0(29) = 0.19940100e+02_realType
    up0(30) = 0.20358100e+02_realType
    up0(31) = 0.20769600e+02_realType
    up0(32) = 0.21175800e+02_realType
    up0(33) = 0.21577400e+02_realType
    up0(34) = 0.21975700e+02_realType

    up1(1)  = 0.12846958e+01_realType
    up1(2)  = 0.69922936e+00_realType
    up1(3)  = 0.44186548e+00_realType
    up1(4)  = 0.30093680e+00_realType
    up1(5)  = 0.21425046e+00_realType
    up1(6)  = 0.15674045e+00_realType
    up1(7)  = 0.11605614e+00_realType
    up1(8)  = 0.85352379e-01_realType
    up1(9)  = 0.61235043e-01_realType
    up1(10) = 0.42691639e-01_realType
    up1(11) = 0.29366174e-01_realType
    up1(12) = 0.20326381e-01_realType
    up1(13) = 0.14310141e-01_realType
    up1(14) = 0.10275905e-01_realType
    up1(15) = 0.75205965e-02_realType
    up1(16) = 0.55993913e-02_realType
    up1(17) = 0.42330072e-02_realType
    up1(18) = 0.32428406e-02_realType
    up1(19) = 0.25137042e-02_realType
    up1(20) = 0.19691199e-02_realType
    up1(21) = 0.15570310e-02_realType
    up1(22) = 0.12416035e-02_realType
    up1(23) = 0.99762939e-03_realType
    up1(24) = 0.80730082e-03_realType
    up1(25) = 0.65769508e-03_realType
    up1(26) = 0.53910966e-03_realType
    up1(27) = 0.44453711e-03_realType
    up1(28) = 0.36862857e-03_realType
    up1(29) = 0.30733926e-03_realType
    up1(30) = 0.25762621e-03_realType
    up1(31) = 0.21711632e-03_realType
    up1(32) = 0.18393679e-03_realType
    up1(33) = 0.15665505e-03_realType
    up1(34) = 0.13415441e-03_realType

    up2(1)  = -0.10506864e+01_realType
    up2(2)  = -0.17378349e+00_realType
    up2(3)  = -0.43906517e-01_realType
    up2(4)  = -0.13876317e-01_realType
    up2(5)  = -0.50197713e-02_realType
    up2(6)  = -0.20046033e-02_realType
    up2(7)  = -0.91800803e-03_realType
    up2(8)  = -0.53858650e-03_realType
    up2(9)  = -0.37015912e-03_realType
    up2(10) = -0.23581357e-03_realType
    up2(11) = -0.13214946e-03_realType
    up2(12) = -0.69676197e-04_realType
    up2(13) = -0.36527030e-04_realType
    up2(14) = -0.19485941e-04_realType
    up2(15) = -0.10680793e-04_realType
    up2(16) = -0.60059830e-05_realType
    up2(17) = -0.34636741e-05_realType
    up2(18) = -0.20504308e-05_realType
    up2(19) = -0.12351270e-05_realType
    up2(20) = -0.76062105e-06_realType
    up2(21) = -0.47546804e-06_realType
    up2(22) = -0.30260755e-06_realType
    up2(23) = -0.19542870e-06_realType
    up2(24) = -0.12770106e-06_realType
    up2(25) = -0.84214126e-07_realType
    up2(26) = -0.56643139e-07_realType
    up2(27) = -0.38016078e-07_realType
    up2(28) = -0.26134020e-07_realType
    up2(29) = -0.17887512e-07_realType
    up2(30) = -0.12500449e-07_realType
    up2(31) = -0.86706361e-08_realType
    up2(32) = -0.61786381e-08_realType
    up2(33) = -0.43440887e-08_realType
    up2(34) = -0.31212378e-08_realType

    up3(1)  =  0.30603769e+00_realType
    up3(2)  = -0.74623173e-02_realType
    up3(3)  = -0.35137590e-02_realType
    up3(4)  = -0.90241883e-03_realType
    up3(5)  = -0.23666410e-03_realType
    up3(6)  = -0.69336451e-04_realType
    up3(7)  = -0.21717508e-04_realType
    up3(8)  = -0.53427330e-05_realType
    up3(9)  = -0.12162441e-06_realType
    up3(10) =  0.66537991e-06_realType
    up3(11) =  0.40127135e-06_realType
    up3(12) =  0.17307016e-06_realType
    up3(13) =  0.68595664e-07_realType
    up3(14) =  0.26859566e-07_realType
    up3(15) =  0.10802572e-07_realType
    up3(16) =  0.44423253e-08_realType
    up3(17) =  0.18757232e-08_realType
    up3(18) =  0.83595370e-09_realType
    up3(19) =  0.37334228e-09_realType
    up3(20) =  0.17567414e-09_realType
    up3(21) =  0.82933451e-10_realType
    up3(22) =  0.40989510e-10_realType
    up3(23) =  0.20935863e-10_realType
    up3(24) =  0.10846455e-10_realType
    up3(25) =  0.54661649e-11_realType
    up3(26) =  0.31700296e-11_realType
    up3(27) =  0.15722041e-11_realType
    up3(28) =  0.97074333e-12_realType
    up3(29) =  0.50475514e-12_realType
    up3(30) =  0.32254746e-12_realType
    up3(31) =  0.16247920e-12_realType
    up3(32) =  0.11432002e-12_realType
    up3(33) =  0.59121027e-13_realType
    up3(34) =  0.38726995e-13_realType

    ! Set the values of tuLogFit. Both for k and omega the
    ! logarithm has been fitted.

    tuLogFit(itu1) = .true.
    tuLogFit(itu2) = .true.

    ! Set the values of constants for the cubic fits of the
    ! non-dimensional k and omega values.

    ! Constants for k.

    tup0(1,itu1)  = -0.10178274e+02_realType
    tup0(2,itu1)  = -0.79134047e+01_realType
    tup0(3,itu1)  = -0.62154735e+01_realType
    tup0(4,itu1)  = -0.48268972e+01_realType
    tup0(5,itu1)  = -0.36279650e+01_realType
    tup0(6,itu1)  = -0.25597781e+01_realType
    tup0(7,itu1)  = -0.16005079e+01_realType
    tup0(8,itu1)  = -0.76521262e+00_realType
    tup0(9,itu1)  = -0.10076775e+00_realType
    tup0(10,itu1) =  0.36262719e+00_realType
    tup0(11,itu1) =  0.65553877e+00_realType
    tup0(12,itu1) =  0.83590897e+00_realType
    tup0(13,itu1) =  0.94909088e+00_realType
    tup0(14,itu1) =  0.10224941e+01_realType
    tup0(15,itu1) =  0.10717000e+01_realType
    tup0(16,itu1) =  0.11056409e+01_realType
    tup0(17,itu1) =  0.11295908e+01_realType
    tup0(18,itu1) =  0.11467673e+01_realType
    tup0(19,itu1) =  0.11591867e+01_realType
    tup0(20,itu1) =  0.11681570e+01_realType
    tup0(21,itu1) =  0.11745296e+01_realType
    tup0(22,itu1) =  0.11788734e+01_realType
    tup0(23,itu1) =  0.11815615e+01_realType
    tup0(24,itu1) =  0.11828278e+01_realType
    tup0(25,itu1) =  0.11828094e+01_realType
    tup0(26,itu1) =  0.11815707e+01_realType
    tup0(27,itu1) =  0.11791103e+01_realType
    tup0(28,itu1) =  0.11753665e+01_realType
    tup0(29,itu1) =  0.11702319e+01_realType
    tup0(30,itu1) =  0.11635476e+01_realType
    tup0(31,itu1) =  0.11550903e+01_realType
    tup0(32,itu1) =  0.11445826e+01_realType
    tup0(33,itu1) =  0.11316601e+01_realType
    tup0(34,itu1) =  0.11158659e+01_realType

    ! Constants for omega.

    tup0(1,itu2)  =  0.68385895e+01_realType
    tup0(2,itu2)  =  0.55423492e+01_realType
    tup0(3,itu2)  =  0.45364394e+01_realType
    tup0(4,itu2)  =  0.37003435e+01_realType
    tup0(5,itu2)  =  0.29762436e+01_realType
    tup0(6,itu2)  =  0.23400254e+01_realType
    tup0(7,itu2)  =  0.17897909e+01_realType
    tup0(8,itu2)  =  0.13296526e+01_realType
    tup0(9,itu2)  =  0.94313517e+00_realType
    tup0(10,itu2) =  0.59512633e+00_realType
    tup0(11,itu2) =  0.26383242e+00_realType
    tup0(12,itu2) = -0.54289357e-01_realType
    tup0(13,itu2) = -0.35764684e+00_realType
    tup0(14,itu2) = -0.64548336e+00_realType
    tup0(15,itu2) = -0.91832029e+00_realType
    tup0(16,itu2) = -0.11773601e+01_realType
    tup0(17,itu2) = -0.14240004e+01_realType
    tup0(18,itu2) = -0.16596108e+01_realType
    tup0(19,itu2) = -0.18854088e+01_realType
    tup0(20,itu2) = -0.21024564e+01_realType
    tup0(21,itu2) = -0.23116299e+01_realType
    tup0(22,itu2) = -0.25136741e+01_realType
    tup0(23,itu2) = -0.27091934e+01_realType
    tup0(24,itu2) = -0.28986818e+01_realType
    tup0(25,itu2) = -0.30825349e+01_realType
    tup0(26,itu2) = -0.32610659e+01_realType
    tup0(27,itu2) = -0.34345194e+01_realType
    tup0(28,itu2) = -0.36030725e+01_realType
    tup0(29,itu2) = -0.37668496e+01_realType
    tup0(30,itu2) = -0.39259191e+01_realType
    tup0(31,itu2) = -0.40803056e+01_realType
    tup0(32,itu2) = -0.42299856e+01_realType
    tup0(33,itu2) = -0.43749001e+01_realType
    tup0(34,itu2) = -0.45149548e+01_realType

    ! Constants for k.

    tup1(1,itu1)  =  0.10151083e+02_realType
    tup1(2,itu1)  =  0.54871316e+01_realType
    tup1(3,itu1)  =  0.33494093e+01_realType
    tup1(4,itu1)  =  0.22113000e+01_realType
    tup1(5,itu1)  =  0.15331218e+01_realType
    tup1(6,itu1)  =  0.10899838e+01_realType
    tup1(7,itu1)  =  0.77051060e+00_realType
    tup1(8,itu1)  =  0.51657998e+00_realType
    tup1(9,itu1)  =  0.31302624e+00_realType
    tup1(10,itu1) =  0.16986834e+00_realType
    tup1(11,itu1) =  0.86387987e-01_realType
    tup1(12,itu1) =  0.43725644e-01_realType
    tup1(13,itu1) =  0.22772335e-01_realType
    tup1(14,itu1) =  0.12310034e-01_realType
    tup1(15,itu1) =  0.68940825e-02_realType
    tup1(16,itu1) =  0.39792104e-02_realType
    tup1(17,itu1) =  0.23523017e-02_realType
    tup1(18,itu1) =  0.14137727e-02_realType
    tup1(19,itu1) =  0.85642296e-03_realType
    tup1(20,itu1) =  0.51672343e-03_realType
    tup1(21,itu1) =  0.30463346e-03_realType
    tup1(22,itu1) =  0.16929149e-03_realType
    tup1(23,itu1) =  0.80893185e-04_realType
    tup1(24,itu1) =  0.21762685e-04_realType
    tup1(25,itu1) = -0.18748602e-04_realType
    tup1(26,itu1) = -0.47330395e-04_realType
    tup1(27,itu1) = -0.68310393e-04_realType
    tup1(28,itu1) = -0.84371684e-04_realType
    tup1(29,itu1) = -0.97226185e-04_realType
    tup1(30,itu1) = -0.10813606e-03_realType
    tup1(31,itu1) = -0.11791513e-03_realType
    tup1(32,itu1) = -0.12717807e-03_realType
    tup1(33,itu1) = -0.13644855e-03_realType
    tup1(34,itu1) = -0.14611996e-03_realType

    ! Constants for omega.

    tup1(1,itu2)  = -0.55838269e+01_realType
    tup1(2,itu2)  = -0.31876950e+01_realType
    tup1(3,itu2)  = -0.19989037e+01_realType
    tup1(4,itu2)  = -0.13333526e+01_realType
    tup1(5,itu2)  = -0.91990459e+00_realType
    tup1(6,itu2)  = -0.63785038e+00_realType
    tup1(7,itu2)  = -0.43381141e+00_realType
    tup1(8,itu2)  = -0.29162744e+00_realType
    tup1(9,itu2)  = -0.20386405e+00_realType
    tup1(10,itu2) = -0.15257310e+00_realType
    tup1(11,itu2) = -0.11853766e+00_realType
    tup1(12,itu2) = -0.92571574e-01_realType
    tup1(13,itu2) = -0.72154025e-01_realType
    tup1(14,itu2) = -0.56291949e-01_realType
    tup1(15,itu2) = -0.44100353e-01_realType
    tup1(16,itu2) = -0.34758707e-01_realType
    tup1(17,itu2) = -0.27583190e-01_realType
    tup1(18,itu2) = -0.22041103e-01_realType
    tup1(19,itu2) = -0.17731129e-01_realType
    tup1(20,itu2) = -0.14354453e-01_realType
    tup1(21,itu2) = -0.11689595e-01_realType
    tup1(22,itu2) = -0.95711712e-02_realType
    tup1(23,itu2) = -0.78759450e-02_realType
    tup1(24,itu2) = -0.65109013e-02_realType
    tup1(25,itu2) = -0.54047659e-02_realType
    tup1(26,itu2) = -0.45036146e-02_realType
    tup1(27,itu2) = -0.37655964e-02_realType
    tup1(28,itu2) = -0.31581617e-02_realType
    tup1(29,itu2) = -0.26558406e-02_realType
    tup1(30,itu2) = -0.22385872e-02_realType
    tup1(31,itu2) = -0.18905375e-02_realType
    tup1(32,itu2) = -0.15990497e-02_realType
    tup1(33,itu2) = -0.13540430e-02_realType
    tup1(34,itu2) = -0.11473506e-02_realType

    ! Constants for k.

    tup2(1,itu1)  = -0.13708334e+02_realType
    tup2(2,itu1)  = -0.43370192e+01_realType
    tup2(3,itu1)  = -0.16256260e+01_realType
    tup2(4,itu1)  = -0.69729756e+00_realType
    tup2(5,itu1)  = -0.32831484e+00_realType
    tup2(6,itu1)  = -0.16501608e+00_realType
    tup2(7,itu1)  = -0.93271929e-01_realType
    tup2(8,itu1)  = -0.66905532e-01_realType
    tup2(9,itu1)  = -0.49449221e-01_realType
    tup2(10,itu1) = -0.27955265e-01_realType
    tup2(11,itu1) = -0.12356466e-01_realType
    tup2(12,itu1) = -0.49538360e-02_realType
    tup2(13,itu1) = -0.19816552e-02_realType
    tup2(14,itu1) = -0.82035726e-03_realType
    tup2(15,itu1) = -0.35459496e-03_realType
    tup2(16,itu1) = -0.15988154e-03_realType
    tup2(17,itu1) = -0.74920111e-04_realType
    tup2(18,itu1) = -0.36432844e-04_realType
    tup2(19,itu1) = -0.18228517e-04_realType
    tup2(20,itu1) = -0.94187249e-05_realType
    tup2(21,itu1) = -0.49803455e-05_realType
    tup2(22,itu1) = -0.26991628e-05_realType
    tup2(23,itu1) = -0.15033797e-05_realType
    tup2(24,itu1) = -0.85865308e-06_realType
    tup2(25,itu1) = -0.50010225e-06_realType
    tup2(26,itu1) = -0.30003597e-06_realType
    tup2(27,itu1) = -0.18914449e-06_realType
    tup2(28,itu1) = -0.12284699e-06_realType
    tup2(29,itu1) = -0.82382282e-07_realType
    tup2(30,itu1) = -0.60415725e-07_realType
    tup2(31,itu1) = -0.44704797e-07_realType
    tup2(32,itu1) = -0.36272580e-07_realType
    tup2(33,itu1) = -0.30797090e-07_realType
    tup2(34,itu1) = -0.27066472e-07_realType

    ! Constants for omega.

    tup2(1,itu2)  =  0.65688815e+01_realType
    tup2(2,itu2)  =  0.22942977e+01_realType
    tup2(3,itu2)  =  0.91326107e+00_realType
    tup2(4,itu2)  =  0.40527609e+00_realType
    tup2(5,itu2)  =  0.19819770e+00_realType
    tup2(6,itu2)  =  0.11119507e+00_realType
    tup2(7,itu2)  =  0.71308510e-01_realType
    tup2(8,itu2)  =  0.41419218e-01_realType
    tup2(9,itu2)  =  0.18358706e-01_realType
    tup2(10,itu2) =  0.79158389e-02_realType
    tup2(11,itu2) =  0.45072385e-02_realType
    tup2(12,itu2) =  0.29542525e-02_realType
    tup2(13,itu2) =  0.19335211e-02_realType
    tup2(14,itu2) =  0.12420728e-02_realType
    tup2(15,itu2) =  0.79078278e-03_realType
    tup2(16,itu2) =  0.50394724e-03_realType
    tup2(17,itu2) =  0.32312896e-03_realType
    tup2(18,itu2) =  0.20923598e-03_realType
    tup2(19,itu2) =  0.13683514e-03_realType
    tup2(20,itu2) =  0.90568628e-04_realType
    tup2(21,itu2) =  0.60506139e-04_realType
    tup2(22,itu2) =  0.40936800e-04_realType
    tup2(23,itu2) =  0.27920509e-04_realType
    tup2(24,itu2) =  0.19242656e-04_realType
    tup2(25,itu2) =  0.13394215e-04_realType
    tup2(26,itu2) =  0.93908861e-05_realType
    tup2(27,itu2) =  0.66475834e-05_realType
    tup2(28,itu2) =  0.47374937e-05_realType
    tup2(29,itu2) =  0.34060533e-05_realType
    tup2(30,itu2) =  0.24642067e-05_realType
    tup2(31,itu2) =  0.17969941e-05_realType
    tup2(32,itu2) =  0.13185213e-05_realType
    tup2(33,itu2) =  0.97355039e-06_realType
    tup2(34,itu2) =  0.72438303e-06_realType

    ! Constants for k.

    tup3(1,itu1)  =  0.13357255e+02_realType
    tup3(2,itu1)  =  0.27962653e+01_realType
    tup3(3,itu1)  =  0.67564983e+00_realType
    tup3(4,itu1)  =  0.18227593e+00_realType
    tup3(5,itu1)  =  0.48230769e-01_realType
    tup3(6,itu1)  =  0.69085885e-02_realType
    tup3(7,itu1)  = -0.25075980e-02_realType
    tup3(8,itu1)  =  0.15198665e-02_realType
    tup3(9,itu1)  =  0.45292571e-02_realType
    tup3(10,itu1) =  0.29768816e-02_realType
    tup3(11,itu1) =  0.11684431e-02_realType
    tup3(12,itu1) =  0.38217835e-03_realType
    tup3(13,itu1) =  0.12135735e-03_realType
    tup3(14,itu1) =  0.39609332e-04_realType
    tup3(15,itu1) =  0.13512654e-04_realType
    tup3(16,itu1) =  0.48259094e-05_realType
    tup3(17,itu1) =  0.17973362e-05_realType
    tup3(18,itu1) =  0.70093792e-06_realType
    tup3(19,itu1) =  0.28079141e-06_realType
    tup3(20,itu1) =  0.11741966e-06_realType
    tup3(21,itu1) =  0.50025034e-07_realType
    tup3(22,itu1) =  0.21729950e-07_realType
    tup3(23,itu1) =  0.96902710e-08_realType
    tup3(24,itu1) =  0.43926199e-08_realType
    tup3(25,itu1) =  0.19274189e-08_realType
    tup3(26,itu1) =  0.80093554e-09_realType
    tup3(27,itu1) =  0.33523262e-09_realType
    tup3(28,itu1) =  0.10603408e-09_realType
    tup3(29,itu1) = -0.14220954e-10_realType
    tup3(30,itu1) = -0.43245103e-10_realType
    tup3(31,itu1) = -0.71330147e-10_realType
    tup3(32,itu1) = -0.73789519e-10_realType
    tup3(33,itu1) = -0.73223977e-10_realType
    tup3(34,itu1) = -0.73679919e-10_realType

    ! Constants for omega.

    tup3(1,itu2)  = -0.58652639e+01_realType
    tup3(2,itu2)  = -0.13617293e+01_realType
    tup3(3,itu2)  = -0.34682427e+00_realType
    tup3(4,itu2)  = -0.90911692e-01_realType
    tup3(5,itu2)  = -0.21991054e-01_realType
    tup3(6,itu2)  = -0.81494825e-02_realType
    tup3(7,itu2)  = -0.84291839e-02_realType
    tup3(8,itu2)  = -0.58630201e-02_realType
    tup3(9,itu2)  = -0.18374198e-02_realType
    tup3(10,itu2) = -0.26966903e-03_realType
    tup3(11,itu2) = -0.45904311e-04_realType
    tup3(12,itu2) = -0.34368125e-04_realType
    tup3(13,itu2) = -0.25332953e-04_realType
    tup3(14,itu2) = -0.15345625e-04_realType
    tup3(15,itu2) = -0.83950168e-05_realType
    tup3(16,itu2) = -0.44072558e-05_realType
    tup3(17,itu2) = -0.22740368e-05_realType
    tup3(18,itu2) = -0.11801059e-05_realType
    tup3(19,itu2) = -0.61295741e-06_realType
    tup3(20,itu2) = -0.32633745e-06_realType
    tup3(21,itu2) = -0.17280693e-06_realType
    tup3(22,itu2) = -0.95608569e-07_realType
    tup3(23,itu2) = -0.52411878e-07_realType
    tup3(24,itu2) = -0.29366409e-07_realType
    tup3(25,itu2) = -0.16959125e-07_realType
    tup3(26,itu2) = -0.97228139e-08_realType
    tup3(27,itu2) = -0.57661683e-08_realType
    tup3(28,itu2) = -0.33988144e-08_realType
    tup3(29,itu2) = -0.20698972e-08_realType
    tup3(30,itu2) = -0.12548468e-08_realType
    tup3(31,itu2) = -0.78279227e-09_realType
    tup3(32,itu2) = -0.49051180e-09_realType
    tup3(33,itu2) = -0.30968731e-09_realType
    tup3(34,itu2) = -0.20495795e-09_realType

  end subroutine initCurveFitDataKw

  subroutine initCurveFitDataKwMod
    !
    !       initCurveFitDataKwMod contains the curve fit constants         
    !       for the wall function data for the modified k-omega turbulence 
    !       model.                                                         
    !
    use flowVarRefState
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    ! integer :: ierr

    call terminate("initCurveFitDataKwMod", &
         "Not implemented yet")

  end subroutine initCurveFitDataKwMod

  subroutine initCurveFitDataSST
    !
    !       initCurveFitDataSST contains the curve fit constants for       
    !       the wall function data for Menter's SST turbulence model.      
    !       Warning: Wall function data developed for k-omega model        
    !
    use constants
    use flowVarRefState
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    ! Set the number of data points and allocate the memory for the
    ! arrays of the curve fits.

    nFit = 34

    allocate(ypT(0:nFit), reT(0:nFit),                     &
         up0(nFit),   up1(nFit), up2(nFit), up3(nFit), &
         tup0(nFit,nt1:nt2), tup1(nFit,nt1:nt2),       &
         tup2(nFit,nt1:nt2), tup3(nFit,nt1:nt2),       &
         tuLogFit(nt1:nt2), stat=ierr)
    if(ierr /= 0)                           &
         call terminate("initCurveFitDataSST", &
         "Memory allocation failure for curve fit &
         &coefficients")

    ! Set the values of the Reynolds numbers at interval boundaries.

    reT(0)  = 0.12529547e+00_realType
    reT(1)  = 0.44996057e+00_realType
    reT(2)  = 0.11581311e+01_realType
    reT(3)  = 0.25353238e+01_realType
    reT(4)  = 0.50446282e+01_realType
    reT(5)  = 0.94194631e+01_realType
    reT(6)  = 0.16766555e+02_realType
    reT(7)  = 0.28556753e+02_realType
    reT(8)  = 0.46274930e+02_realType
    reT(9)  = 0.71021000e+02_realType
    reT(10) = 0.10383163e+03_realType
    reT(11) = 0.14621738e+03_realType
    reT(12) = 0.20028019e+03_realType
    reT(13) = 0.26868298e+03_realType
    reT(14) = 0.35467049e+03_realType
    reT(15) = 0.46212508e+03_realType
    reT(16) = 0.59566097e+03_realType
    reT(17) = 0.76073076e+03_realType
    reT(18) = 0.96373333e+03_realType
    reT(19) = 0.12121761e+04_realType
    reT(20) = 0.15147917e+04_realType
    reT(21) = 0.18817196e+04_realType
    reT(22) = 0.23247121e+04_realType
    reT(23) = 0.28572322e+04_realType
    reT(24) = 0.34947840e+04_realType
    reT(25) = 0.42551444e+04_realType
    reT(26) = 0.51584529e+04_realType
    reT(27) = 0.62277581e+04_realType
    reT(28) = 0.74889831e+04_realType
    reT(29) = 0.89716314e+04_realType
    reT(30) = 0.10708764e+05_realType
    reT(31) = 0.12737815e+05_realType
    reT(32) = 0.15100490e+05_realType
    reT(33) = 0.17843939e+05_realType
    reT(34) = 0.21020534e+05_realType

    ! Set the values of the y+ values at interval boundaries.

    ypT(0)  = 0.35397100e+00_realType
    ypT(1)  = 0.67079200e+00_realType
    ypT(2)  = 0.10761700e+01_realType
    ypT(3)  = 0.15923000e+01_realType
    ypT(4)  = 0.22463000e+01_realType
    ypT(5)  = 0.30710600e+01_realType
    ypT(6)  = 0.41063800e+01_realType
    ypT(7)  = 0.54001200e+01_realType
    ypT(8)  = 0.70095900e+01_realType
    ypT(9)  = 0.90031400e+01_realType
    ypT(10) = 0.11461900e+02_realType
    ypT(11) = 0.14481700e+02_realType
    ypT(12) = 0.18175400e+02_realType
    ypT(13) = 0.22675200e+02_realType
    ypT(14) = 0.28135500e+02_realType
    ypT(15) = 0.34735800e+02_realType
    ypT(16) = 0.42683800e+02_realType
    ypT(17) = 0.52219300e+02_realType
    ypT(18) = 0.63617800e+02_realType
    ypT(19) = 0.77194900e+02_realType
    ypT(20) = 0.93310400e+02_realType
    ypT(21) = 0.11237300e+03_realType
    ypT(22) = 0.13484800e+03_realType
    ypT(23) = 0.16125700e+03_realType
    ypT(24) = 0.19218900e+03_realType
    ypT(25) = 0.22830600e+03_realType
    ypT(26) = 0.27034500e+03_realType
    ypT(27) = 0.31913000e+03_realType
    ypT(28) = 0.37557400e+03_realType
    ypT(29) = 0.44069100e+03_realType
    ypT(30) = 0.51559800e+03_realType
    ypT(31) = 0.60152700e+03_realType
    ypT(32) = 0.69982900e+03_realType
    ypT(33) = 0.81198500e+03_realType
    ypT(34) = 0.93960800e+03_realType

    ! Set the values of constants for the cubic fits of the
    ! non-dimensional tangential velocity.

    up0(1)  = 0.35397100e+00_realType
    up0(2)  = 0.67079000e+00_realType
    up0(3)  = 0.10761600e+01_realType
    up0(4)  = 0.15922400e+01_realType
    up0(5)  = 0.22457500e+01_realType
    up0(6)  = 0.30671700e+01_realType
    up0(7)  = 0.40830500e+01_realType
    up0(8)  = 0.52881700e+01_realType
    up0(9)  = 0.66016600e+01_realType
    up0(10) = 0.78884700e+01_realType
    up0(11) = 0.90588500e+01_realType
    up0(12) = 0.10096700e+02_realType
    up0(13) = 0.11019300e+02_realType
    up0(14) = 0.11849200e+02_realType
    up0(15) = 0.12605800e+02_realType
    up0(16) = 0.13304000e+02_realType
    up0(17) = 0.13955200e+02_realType
    up0(18) = 0.14568000e+02_realType
    up0(19) = 0.15148800e+02_realType
    up0(20) = 0.15702800e+02_realType
    up0(21) = 0.16233900e+02_realType
    up0(22) = 0.16745300e+02_realType
    up0(23) = 0.17239500e+02_realType
    up0(24) = 0.17718500e+02_realType
    up0(25) = 0.18184100e+02_realType
    up0(26) = 0.18637900e+02_realType
    up0(27) = 0.19081000e+02_realType
    up0(28) = 0.19514800e+02_realType
    up0(29) = 0.19940100e+02_realType
    up0(30) = 0.20358100e+02_realType
    up0(31) = 0.20769600e+02_realType
    up0(32) = 0.21175800e+02_realType
    up0(33) = 0.21577400e+02_realType
    up0(34) = 0.21975700e+02_realType

    up1(1)  = 0.12846958e+01_realType
    up1(2)  = 0.69922936e+00_realType
    up1(3)  = 0.44186548e+00_realType
    up1(4)  = 0.30093680e+00_realType
    up1(5)  = 0.21425046e+00_realType
    up1(6)  = 0.15674045e+00_realType
    up1(7)  = 0.11605614e+00_realType
    up1(8)  = 0.85352379e-01_realType
    up1(9)  = 0.61235043e-01_realType
    up1(10) = 0.42691639e-01_realType
    up1(11) = 0.29366174e-01_realType
    up1(12) = 0.20326381e-01_realType
    up1(13) = 0.14310141e-01_realType
    up1(14) = 0.10275905e-01_realType
    up1(15) = 0.75205965e-02_realType
    up1(16) = 0.55993913e-02_realType
    up1(17) = 0.42330072e-02_realType
    up1(18) = 0.32428406e-02_realType
    up1(19) = 0.25137042e-02_realType
    up1(20) = 0.19691199e-02_realType
    up1(21) = 0.15570310e-02_realType
    up1(22) = 0.12416035e-02_realType
    up1(23) = 0.99762939e-03_realType
    up1(24) = 0.80730082e-03_realType
    up1(25) = 0.65769508e-03_realType
    up1(26) = 0.53910966e-03_realType
    up1(27) = 0.44453711e-03_realType
    up1(28) = 0.36862857e-03_realType
    up1(29) = 0.30733926e-03_realType
    up1(30) = 0.25762621e-03_realType
    up1(31) = 0.21711632e-03_realType
    up1(32) = 0.18393679e-03_realType
    up1(33) = 0.15665505e-03_realType
    up1(34) = 0.13415441e-03_realType

    up2(1)  = -0.10506864e+01_realType
    up2(2)  = -0.17378349e+00_realType
    up2(3)  = -0.43906517e-01_realType
    up2(4)  = -0.13876317e-01_realType
    up2(5)  = -0.50197713e-02_realType
    up2(6)  = -0.20046033e-02_realType
    up2(7)  = -0.91800803e-03_realType
    up2(8)  = -0.53858650e-03_realType
    up2(9)  = -0.37015912e-03_realType
    up2(10) = -0.23581357e-03_realType
    up2(11) = -0.13214946e-03_realType
    up2(12) = -0.69676197e-04_realType
    up2(13) = -0.36527030e-04_realType
    up2(14) = -0.19485941e-04_realType
    up2(15) = -0.10680793e-04_realType
    up2(16) = -0.60059830e-05_realType
    up2(17) = -0.34636741e-05_realType
    up2(18) = -0.20504308e-05_realType
    up2(19) = -0.12351270e-05_realType
    up2(20) = -0.76062105e-06_realType
    up2(21) = -0.47546804e-06_realType
    up2(22) = -0.30260755e-06_realType
    up2(23) = -0.19542870e-06_realType
    up2(24) = -0.12770106e-06_realType
    up2(25) = -0.84214126e-07_realType
    up2(26) = -0.56643139e-07_realType
    up2(27) = -0.38016078e-07_realType
    up2(28) = -0.26134020e-07_realType
    up2(29) = -0.17887512e-07_realType
    up2(30) = -0.12500449e-07_realType
    up2(31) = -0.86706361e-08_realType
    up2(32) = -0.61786381e-08_realType
    up2(33) = -0.43440887e-08_realType
    up2(34) = -0.31212378e-08_realType

    up3(1)  =  0.30603769e+00_realType
    up3(2)  = -0.74623173e-02_realType
    up3(3)  = -0.35137590e-02_realType
    up3(4)  = -0.90241883e-03_realType
    up3(5)  = -0.23666410e-03_realType
    up3(6)  = -0.69336451e-04_realType
    up3(7)  = -0.21717508e-04_realType
    up3(8)  = -0.53427330e-05_realType
    up3(9)  = -0.12162441e-06_realType
    up3(10) =  0.66537991e-06_realType
    up3(11) =  0.40127135e-06_realType
    up3(12) =  0.17307016e-06_realType
    up3(13) =  0.68595664e-07_realType
    up3(14) =  0.26859566e-07_realType
    up3(15) =  0.10802572e-07_realType
    up3(16) =  0.44423253e-08_realType
    up3(17) =  0.18757232e-08_realType
    up3(18) =  0.83595370e-09_realType
    up3(19) =  0.37334228e-09_realType
    up3(20) =  0.17567414e-09_realType
    up3(21) =  0.82933451e-10_realType
    up3(22) =  0.40989510e-10_realType
    up3(23) =  0.20935863e-10_realType
    up3(24) =  0.10846455e-10_realType
    up3(25) =  0.54661649e-11_realType
    up3(26) =  0.31700296e-11_realType
    up3(27) =  0.15722041e-11_realType
    up3(28) =  0.97074333e-12_realType
    up3(29) =  0.50475514e-12_realType
    up3(30) =  0.32254746e-12_realType
    up3(31) =  0.16247920e-12_realType
    up3(32) =  0.11432002e-12_realType
    up3(33) =  0.59121027e-13_realType
    up3(34) =  0.38726995e-13_realType

    ! Set the values of tuLogFit. Both for k and omega the
    ! logarithm has been fitted.

    tuLogFit(itu1) = .true.
    tuLogFit(itu2) = .true.

    ! Set the values of constants for the cubic fits of the
    ! non-dimensional k and omega values.

    ! Constants for k.

    tup0(1,itu1)  = -0.10178274e+02_realType
    tup0(2,itu1)  = -0.79134047e+01_realType
    tup0(3,itu1)  = -0.62154735e+01_realType
    tup0(4,itu1)  = -0.48268972e+01_realType
    tup0(5,itu1)  = -0.36279650e+01_realType
    tup0(6,itu1)  = -0.25597781e+01_realType
    tup0(7,itu1)  = -0.16005079e+01_realType
    tup0(8,itu1)  = -0.76521262e+00_realType
    tup0(9,itu1)  = -0.10076775e+00_realType
    tup0(10,itu1) =  0.36262719e+00_realType
    tup0(11,itu1) =  0.65553877e+00_realType
    tup0(12,itu1) =  0.83590897e+00_realType
    tup0(13,itu1) =  0.94909088e+00_realType
    tup0(14,itu1) =  0.10224941e+01_realType
    tup0(15,itu1) =  0.10717000e+01_realType
    tup0(16,itu1) =  0.11056409e+01_realType
    tup0(17,itu1) =  0.11295908e+01_realType
    tup0(18,itu1) =  0.11467673e+01_realType
    tup0(19,itu1) =  0.11591867e+01_realType
    tup0(20,itu1) =  0.11681570e+01_realType
    tup0(21,itu1) =  0.11745296e+01_realType
    tup0(22,itu1) =  0.11788734e+01_realType
    tup0(23,itu1) =  0.11815615e+01_realType
    tup0(24,itu1) =  0.11828278e+01_realType
    tup0(25,itu1) =  0.11828094e+01_realType
    tup0(26,itu1) =  0.11815707e+01_realType
    tup0(27,itu1) =  0.11791103e+01_realType
    tup0(28,itu1) =  0.11753665e+01_realType
    tup0(29,itu1) =  0.11702319e+01_realType
    tup0(30,itu1) =  0.11635476e+01_realType
    tup0(31,itu1) =  0.11550903e+01_realType
    tup0(32,itu1) =  0.11445826e+01_realType
    tup0(33,itu1) =  0.11316601e+01_realType
    tup0(34,itu1) =  0.11158659e+01_realType

    ! Constants for omega.

    tup0(1,itu2)  =  0.68385895e+01_realType
    tup0(2,itu2)  =  0.55423492e+01_realType
    tup0(3,itu2)  =  0.45364394e+01_realType
    tup0(4,itu2)  =  0.37003435e+01_realType
    tup0(5,itu2)  =  0.29762436e+01_realType
    tup0(6,itu2)  =  0.23400254e+01_realType
    tup0(7,itu2)  =  0.17897909e+01_realType
    tup0(8,itu2)  =  0.13296526e+01_realType
    tup0(9,itu2)  =  0.94313517e+00_realType
    tup0(10,itu2) =  0.59512633e+00_realType
    tup0(11,itu2) =  0.26383242e+00_realType
    tup0(12,itu2) = -0.54289357e-01_realType
    tup0(13,itu2) = -0.35764684e+00_realType
    tup0(14,itu2) = -0.64548336e+00_realType
    tup0(15,itu2) = -0.91832029e+00_realType
    tup0(16,itu2) = -0.11773601e+01_realType
    tup0(17,itu2) = -0.14240004e+01_realType
    tup0(18,itu2) = -0.16596108e+01_realType
    tup0(19,itu2) = -0.18854088e+01_realType
    tup0(20,itu2) = -0.21024564e+01_realType
    tup0(21,itu2) = -0.23116299e+01_realType
    tup0(22,itu2) = -0.25136741e+01_realType
    tup0(23,itu2) = -0.27091934e+01_realType
    tup0(24,itu2) = -0.28986818e+01_realType
    tup0(25,itu2) = -0.30825349e+01_realType
    tup0(26,itu2) = -0.32610659e+01_realType
    tup0(27,itu2) = -0.34345194e+01_realType
    tup0(28,itu2) = -0.36030725e+01_realType
    tup0(29,itu2) = -0.37668496e+01_realType
    tup0(30,itu2) = -0.39259191e+01_realType
    tup0(31,itu2) = -0.40803056e+01_realType
    tup0(32,itu2) = -0.42299856e+01_realType
    tup0(33,itu2) = -0.43749001e+01_realType
    tup0(34,itu2) = -0.45149548e+01_realType

    ! Constants for k.

    tup1(1,itu1)  =  0.10151083e+02_realType
    tup1(2,itu1)  =  0.54871316e+01_realType
    tup1(3,itu1)  =  0.33494093e+01_realType
    tup1(4,itu1)  =  0.22113000e+01_realType
    tup1(5,itu1)  =  0.15331218e+01_realType
    tup1(6,itu1)  =  0.10899838e+01_realType
    tup1(7,itu1)  =  0.77051060e+00_realType
    tup1(8,itu1)  =  0.51657998e+00_realType
    tup1(9,itu1)  =  0.31302624e+00_realType
    tup1(10,itu1) =  0.16986834e+00_realType
    tup1(11,itu1) =  0.86387987e-01_realType
    tup1(12,itu1) =  0.43725644e-01_realType
    tup1(13,itu1) =  0.22772335e-01_realType
    tup1(14,itu1) =  0.12310034e-01_realType
    tup1(15,itu1) =  0.68940825e-02_realType
    tup1(16,itu1) =  0.39792104e-02_realType
    tup1(17,itu1) =  0.23523017e-02_realType
    tup1(18,itu1) =  0.14137727e-02_realType
    tup1(19,itu1) =  0.85642296e-03_realType
    tup1(20,itu1) =  0.51672343e-03_realType
    tup1(21,itu1) =  0.30463346e-03_realType
    tup1(22,itu1) =  0.16929149e-03_realType
    tup1(23,itu1) =  0.80893185e-04_realType
    tup1(24,itu1) =  0.21762685e-04_realType
    tup1(25,itu1) = -0.18748602e-04_realType
    tup1(26,itu1) = -0.47330395e-04_realType
    tup1(27,itu1) = -0.68310393e-04_realType
    tup1(28,itu1) = -0.84371684e-04_realType
    tup1(29,itu1) = -0.97226185e-04_realType
    tup1(30,itu1) = -0.10813606e-03_realType
    tup1(31,itu1) = -0.11791513e-03_realType
    tup1(32,itu1) = -0.12717807e-03_realType
    tup1(33,itu1) = -0.13644855e-03_realType
    tup1(34,itu1) = -0.14611996e-03_realType

    ! Constants for omega.

    tup1(1,itu2)  = -0.55838269e+01_realType
    tup1(2,itu2)  = -0.31876950e+01_realType
    tup1(3,itu2)  = -0.19989037e+01_realType
    tup1(4,itu2)  = -0.13333526e+01_realType
    tup1(5,itu2)  = -0.91990459e+00_realType
    tup1(6,itu2)  = -0.63785038e+00_realType
    tup1(7,itu2)  = -0.43381141e+00_realType
    tup1(8,itu2)  = -0.29162744e+00_realType
    tup1(9,itu2)  = -0.20386405e+00_realType
    tup1(10,itu2) = -0.15257310e+00_realType
    tup1(11,itu2) = -0.11853766e+00_realType
    tup1(12,itu2) = -0.92571574e-01_realType
    tup1(13,itu2) = -0.72154025e-01_realType
    tup1(14,itu2) = -0.56291949e-01_realType
    tup1(15,itu2) = -0.44100353e-01_realType
    tup1(16,itu2) = -0.34758707e-01_realType
    tup1(17,itu2) = -0.27583190e-01_realType
    tup1(18,itu2) = -0.22041103e-01_realType
    tup1(19,itu2) = -0.17731129e-01_realType
    tup1(20,itu2) = -0.14354453e-01_realType
    tup1(21,itu2) = -0.11689595e-01_realType
    tup1(22,itu2) = -0.95711712e-02_realType
    tup1(23,itu2) = -0.78759450e-02_realType
    tup1(24,itu2) = -0.65109013e-02_realType
    tup1(25,itu2) = -0.54047659e-02_realType
    tup1(26,itu2) = -0.45036146e-02_realType
    tup1(27,itu2) = -0.37655964e-02_realType
    tup1(28,itu2) = -0.31581617e-02_realType
    tup1(29,itu2) = -0.26558406e-02_realType
    tup1(30,itu2) = -0.22385872e-02_realType
    tup1(31,itu2) = -0.18905375e-02_realType
    tup1(32,itu2) = -0.15990497e-02_realType
    tup1(33,itu2) = -0.13540430e-02_realType
    tup1(34,itu2) = -0.11473506e-02_realType

    ! Constants for k.

    tup2(1,itu1)  = -0.13708334e+02_realType
    tup2(2,itu1)  = -0.43370192e+01_realType
    tup2(3,itu1)  = -0.16256260e+01_realType
    tup2(4,itu1)  = -0.69729756e+00_realType
    tup2(5,itu1)  = -0.32831484e+00_realType
    tup2(6,itu1)  = -0.16501608e+00_realType
    tup2(7,itu1)  = -0.93271929e-01_realType
    tup2(8,itu1)  = -0.66905532e-01_realType
    tup2(9,itu1)  = -0.49449221e-01_realType
    tup2(10,itu1) = -0.27955265e-01_realType
    tup2(11,itu1) = -0.12356466e-01_realType
    tup2(12,itu1) = -0.49538360e-02_realType
    tup2(13,itu1) = -0.19816552e-02_realType
    tup2(14,itu1) = -0.82035726e-03_realType
    tup2(15,itu1) = -0.35459496e-03_realType
    tup2(16,itu1) = -0.15988154e-03_realType
    tup2(17,itu1) = -0.74920111e-04_realType
    tup2(18,itu1) = -0.36432844e-04_realType
    tup2(19,itu1) = -0.18228517e-04_realType
    tup2(20,itu1) = -0.94187249e-05_realType
    tup2(21,itu1) = -0.49803455e-05_realType
    tup2(22,itu1) = -0.26991628e-05_realType
    tup2(23,itu1) = -0.15033797e-05_realType
    tup2(24,itu1) = -0.85865308e-06_realType
    tup2(25,itu1) = -0.50010225e-06_realType
    tup2(26,itu1) = -0.30003597e-06_realType
    tup2(27,itu1) = -0.18914449e-06_realType
    tup2(28,itu1) = -0.12284699e-06_realType
    tup2(29,itu1) = -0.82382282e-07_realType
    tup2(30,itu1) = -0.60415725e-07_realType
    tup2(31,itu1) = -0.44704797e-07_realType
    tup2(32,itu1) = -0.36272580e-07_realType
    tup2(33,itu1) = -0.30797090e-07_realType
    tup2(34,itu1) = -0.27066472e-07_realType

    ! Constants for omega.

    tup2(1,itu2)  =  0.65688815e+01_realType
    tup2(2,itu2)  =  0.22942977e+01_realType
    tup2(3,itu2)  =  0.91326107e+00_realType
    tup2(4,itu2)  =  0.40527609e+00_realType
    tup2(5,itu2)  =  0.19819770e+00_realType
    tup2(6,itu2)  =  0.11119507e+00_realType
    tup2(7,itu2)  =  0.71308510e-01_realType
    tup2(8,itu2)  =  0.41419218e-01_realType
    tup2(9,itu2)  =  0.18358706e-01_realType
    tup2(10,itu2) =  0.79158389e-02_realType
    tup2(11,itu2) =  0.45072385e-02_realType
    tup2(12,itu2) =  0.29542525e-02_realType
    tup2(13,itu2) =  0.19335211e-02_realType
    tup2(14,itu2) =  0.12420728e-02_realType
    tup2(15,itu2) =  0.79078278e-03_realType
    tup2(16,itu2) =  0.50394724e-03_realType
    tup2(17,itu2) =  0.32312896e-03_realType
    tup2(18,itu2) =  0.20923598e-03_realType
    tup2(19,itu2) =  0.13683514e-03_realType
    tup2(20,itu2) =  0.90568628e-04_realType
    tup2(21,itu2) =  0.60506139e-04_realType
    tup2(22,itu2) =  0.40936800e-04_realType
    tup2(23,itu2) =  0.27920509e-04_realType
    tup2(24,itu2) =  0.19242656e-04_realType
    tup2(25,itu2) =  0.13394215e-04_realType
    tup2(26,itu2) =  0.93908861e-05_realType
    tup2(27,itu2) =  0.66475834e-05_realType
    tup2(28,itu2) =  0.47374937e-05_realType
    tup2(29,itu2) =  0.34060533e-05_realType
    tup2(30,itu2) =  0.24642067e-05_realType
    tup2(31,itu2) =  0.17969941e-05_realType
    tup2(32,itu2) =  0.13185213e-05_realType
    tup2(33,itu2) =  0.97355039e-06_realType
    tup2(34,itu2) =  0.72438303e-06_realType

    ! Constants for k.

    tup3(1,itu1)  =  0.13357255e+02_realType
    tup3(2,itu1)  =  0.27962653e+01_realType
    tup3(3,itu1)  =  0.67564983e+00_realType
    tup3(4,itu1)  =  0.18227593e+00_realType
    tup3(5,itu1)  =  0.48230769e-01_realType
    tup3(6,itu1)  =  0.69085885e-02_realType
    tup3(7,itu1)  = -0.25075980e-02_realType
    tup3(8,itu1)  =  0.15198665e-02_realType
    tup3(9,itu1)  =  0.45292571e-02_realType
    tup3(10,itu1) =  0.29768816e-02_realType
    tup3(11,itu1) =  0.11684431e-02_realType
    tup3(12,itu1) =  0.38217835e-03_realType
    tup3(13,itu1) =  0.12135735e-03_realType
    tup3(14,itu1) =  0.39609332e-04_realType
    tup3(15,itu1) =  0.13512654e-04_realType
    tup3(16,itu1) =  0.48259094e-05_realType
    tup3(17,itu1) =  0.17973362e-05_realType
    tup3(18,itu1) =  0.70093792e-06_realType
    tup3(19,itu1) =  0.28079141e-06_realType
    tup3(20,itu1) =  0.11741966e-06_realType
    tup3(21,itu1) =  0.50025034e-07_realType
    tup3(22,itu1) =  0.21729950e-07_realType
    tup3(23,itu1) =  0.96902710e-08_realType
    tup3(24,itu1) =  0.43926199e-08_realType
    tup3(25,itu1) =  0.19274189e-08_realType
    tup3(26,itu1) =  0.80093554e-09_realType
    tup3(27,itu1) =  0.33523262e-09_realType
    tup3(28,itu1) =  0.10603408e-09_realType
    tup3(29,itu1) = -0.14220954e-10_realType
    tup3(30,itu1) = -0.43245103e-10_realType
    tup3(31,itu1) = -0.71330147e-10_realType
    tup3(32,itu1) = -0.73789519e-10_realType
    tup3(33,itu1) = -0.73223977e-10_realType
    tup3(34,itu1) = -0.73679919e-10_realType

    ! Constants for omega.

    tup3(1,itu2)  = -0.58652639e+01_realType
    tup3(2,itu2)  = -0.13617293e+01_realType
    tup3(3,itu2)  = -0.34682427e+00_realType
    tup3(4,itu2)  = -0.90911692e-01_realType
    tup3(5,itu2)  = -0.21991054e-01_realType
    tup3(6,itu2)  = -0.81494825e-02_realType
    tup3(7,itu2)  = -0.84291839e-02_realType
    tup3(8,itu2)  = -0.58630201e-02_realType
    tup3(9,itu2)  = -0.18374198e-02_realType
    tup3(10,itu2) = -0.26966903e-03_realType
    tup3(11,itu2) = -0.45904311e-04_realType
    tup3(12,itu2) = -0.34368125e-04_realType
    tup3(13,itu2) = -0.25332953e-04_realType
    tup3(14,itu2) = -0.15345625e-04_realType
    tup3(15,itu2) = -0.83950168e-05_realType
    tup3(16,itu2) = -0.44072558e-05_realType
    tup3(17,itu2) = -0.22740368e-05_realType
    tup3(18,itu2) = -0.11801059e-05_realType
    tup3(19,itu2) = -0.61295741e-06_realType
    tup3(20,itu2) = -0.32633745e-06_realType
    tup3(21,itu2) = -0.17280693e-06_realType
    tup3(22,itu2) = -0.95608569e-07_realType
    tup3(23,itu2) = -0.52411878e-07_realType
    tup3(24,itu2) = -0.29366409e-07_realType
    tup3(25,itu2) = -0.16959125e-07_realType
    tup3(26,itu2) = -0.97228139e-08_realType
    tup3(27,itu2) = -0.57661683e-08_realType
    tup3(28,itu2) = -0.33988144e-08_realType
    tup3(29,itu2) = -0.20698972e-08_realType
    tup3(30,itu2) = -0.12548468e-08_realType
    tup3(31,itu2) = -0.78279227e-09_realType
    tup3(32,itu2) = -0.49051180e-09_realType
    tup3(33,itu2) = -0.30968731e-09_realType
    tup3(34,itu2) = -0.20495795e-09_realType

  end subroutine initCurveFitDataSST

  subroutine initCurveFitDataSa
    !
    !       initCurveFitDataSa contains the curve fit constants for        
    !       the wall function data for the Spalart-Allmaras turbulence     
    !       model.                                                         
    !
    use constants
    use flowVarRefState
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    ! Set the number of data points and allocate the memory for the
    ! arrays of the curve fits.

    nFit = 34

    allocate(ypT(0:nFit), reT(0:nFit),                     &
         up0(nFit),   up1(nFit), up2(nFit), up3(nFit), &
         tup0(nFit,nt1:nt2), tup1(nFit,nt1:nt2),       &
         tup2(nFit,nt1:nt2), tup3(nFit,nt1:nt2),       &
         tuLogFit(nt1:nt2), stat=ierr)
    if(ierr /= 0)                          &
         call terminate("initCurveFitDataSa", &
         "Memory allocation failure for curve fit &
         &coefficients")

    ! Set the values of the Reynolds numbers at interval boundaries.

    reT(0)  = 0.12361553e+00_realType
    reT(1)  = 0.44392837e+00_realType
    reT(2)  = 0.11425793e+01_realType
    reT(3)  = 0.25011739e+01_realType
    reT(4)  = 0.49762007e+01_realType
    reT(5)  = 0.92920979e+01_realType
    reT(6)  = 0.16564578e+02_realType
    reT(7)  = 0.28414621e+02_realType
    reT(8)  = 0.46909987e+02_realType
    reT(9)  = 0.73988906e+02_realType
    reT(10) = 0.11046933e+03_realType
    reT(11) = 0.15636562e+03_realType
    reT(12) = 0.21263059e+03_realType
    reT(13) = 0.28162960e+03_realType
    reT(14) = 0.36666795e+03_realType
    reT(15) = 0.47173270e+03_realType
    reT(16) = 0.60148482e+03_realType
    reT(17) = 0.76135031e+03_realType
    reT(18) = 0.95763636e+03_realType
    reT(19) = 0.11976883e+04_realType
    reT(20) = 0.14900416e+04_realType
    reT(21) = 0.18445991e+04_realType
    reT(22) = 0.22728484e+04_realType
    reT(23) = 0.27879873e+04_realType
    reT(24) = 0.34052767e+04_realType
    reT(25) = 0.41422400e+04_realType
    reT(26) = 0.50189226e+04_realType
    reT(27) = 0.60583758e+04_realType
    reT(28) = 0.72868976e+04_realType
    reT(29) = 0.87344007e+04_realType
    reT(30) = 0.10435120e+05_realType
    reT(31) = 0.12427867e+05_realType
    reT(32) = 0.14756830e+05_realType
    reT(33) = 0.17471977e+05_realType
    reT(34) = 0.20629717e+05_realType

    ! Set the values of the y+ values at interval boundaries.

    ypT(0)  = 0.35159200e+00_realType
    ypT(1)  = 0.66628400e+00_realType
    ypT(2)  = 0.10689300e+01_realType
    ypT(3)  = 0.15816000e+01_realType
    ypT(4)  = 0.22312000e+01_realType
    ypT(5)  = 0.30504200e+01_realType
    ypT(6)  = 0.40787900e+01_realType
    ypT(7)  = 0.53638300e+01_realType
    ypT(8)  = 0.69624800e+01_realType
    ypT(9)  = 0.89426300e+01_realType
    ypT(10) = 0.11384800e+02_realType
    ypT(11) = 0.14384400e+02_realType
    ypT(12) = 0.18053200e+02_realType
    ypT(13) = 0.22522800e+02_realType
    ypT(14) = 0.27946400e+02_realType
    ypT(15) = 0.34502300e+02_realType
    ypT(16) = 0.42396900e+02_realType
    ypT(17) = 0.51868400e+02_realType
    ypT(18) = 0.63190300e+02_realType
    ypT(19) = 0.76676100e+02_realType
    ypT(20) = 0.92683300e+02_realType
    ypT(21) = 0.11161800e+03_realType
    ypT(22) = 0.13394200e+03_realType
    ypT(23) = 0.16017300e+03_realType
    ypT(24) = 0.19089800e+03_realType
    ypT(25) = 0.22677200e+03_realType
    ypT(26) = 0.26852800e+03_realType
    ypT(27) = 0.31698500e+03_realType
    ypT(28) = 0.37305000e+03_realType
    ypT(29) = 0.43772900e+03_realType
    ypT(30) = 0.51213300e+03_realType
    ypT(31) = 0.59748500e+03_realType
    ypT(32) = 0.69512600e+03_realType
    ypT(33) = 0.80652800e+03_realType
    ypT(34) = 0.93329400e+03_realType

    ! Set the values of constants for the cubic fits of the
    ! non-dimensional tangential velocity.

    up0(1)  = 0.35158800e+00_realType
    up0(2)  = 0.66627500e+00_realType
    up0(3)  = 0.10689000e+01_realType
    up0(4)  = 0.15814200e+01_realType
    up0(5)  = 0.22302800e+01_realType
    up0(6)  = 0.30461700e+01_realType
    up0(7)  = 0.40611500e+01_realType
    up0(8)  = 0.52974500e+01_realType
    up0(9)  = 0.67375400e+01_realType
    up0(10) = 0.82737300e+01_realType
    up0(11) = 0.97032300e+01_realType
    up0(12) = 0.10870500e+02_realType
    up0(13) = 0.11778000e+02_realType
    up0(14) = 0.12504200e+02_realType
    up0(15) = 0.13120400e+02_realType
    up0(16) = 0.13672500e+02_realType
    up0(17) = 0.14187000e+02_realType
    up0(18) = 0.14678500e+02_realType
    up0(19) = 0.15154800e+02_realType
    up0(20) = 0.15620100e+02_realType
    up0(21) = 0.16076700e+02_realType
    up0(22) = 0.16526000e+02_realType
    up0(23) = 0.16968900e+02_realType
    up0(24) = 0.17406100e+02_realType
    up0(25) = 0.17838200e+02_realType
    up0(26) = 0.18266100e+02_realType
    up0(27) = 0.18690500e+02_realType
    up0(28) = 0.19112500e+02_realType
    up0(29) = 0.19533300e+02_realType
    up0(30) = 0.19953900e+02_realType
    up0(31) = 0.20375800e+02_realType
    up0(32) = 0.20800300e+02_realType
    up0(33) = 0.21229000e+02_realType
    up0(34) = 0.21663200e+02_realType

    up1(1)  = 0.12933934e+01_realType
    up1(2)  = 0.70396224e+00_realType
    up1(3)  = 0.44483996e+00_realType
    up1(4)  = 0.30294593e+00_realType
    up1(5)  = 0.21569230e+00_realType
    up1(6)  = 0.15799192e+00_realType
    up1(7)  = 0.11772923e+00_realType
    up1(8)  = 0.88197525e-01_realType
    up1(9)  = 0.65306126e-01_realType
    up1(10) = 0.46660172e-01_realType
    up1(11) = 0.31523107e-01_realType
    up1(12) = 0.20308775e-01_realType
    up1(13) = 0.13042058e-01_realType
    up1(14) = 0.87147691e-02_realType
    up1(15) = 0.61456125e-02_realType
    up1(16) = 0.45422630e-02_realType
    up1(17) = 0.34735457e-02_realType
    up1(18) = 0.27173826e-02_realType
    up1(19) = 0.21579599e-02_realType
    up1(20) = 0.17315757e-02_realType
    up1(21) = 0.14003478e-02_realType
    up1(22) = 0.11397448e-02_realType
    up1(23) = 0.93291395e-03_realType
    up1(24) = 0.76764242e-03_realType
    up1(25) = 0.63503654e-03_realType
    up1(26) = 0.52818280e-03_realType
    up1(27) = 0.44172235e-03_realType
    up1(28) = 0.37160904e-03_realType
    up1(29) = 0.31442159e-03_realType
    up1(30) = 0.26761137e-03_realType
    up1(31) = 0.22916141e-03_realType
    up1(32) = 0.19742184e-03_realType
    up1(33) = 0.17107081e-03_realType
    up1(34) = 0.14902380e-03_realType

    up2(1)  = -0.10722013e+01_realType
    up2(2)  = -0.17733707e+00_realType
    up2(3)  = -0.44823897e-01_realType
    up2(4)  = -0.14179933e-01_realType
    up2(5)  = -0.51548071e-02_realType
    up2(6)  = -0.20652649e-02_realType
    up2(7)  = -0.90040092e-03_realType
    up2(8)  = -0.43873478e-03_realType
    up2(9)  = -0.26153548e-03_realType
    up2(10) = -0.19975811e-03_realType
    up2(11) = -0.15375234e-03_realType
    up2(12) = -0.93708131e-04_realType
    up2(13) = -0.46732800e-04_realType
    up2(14) = -0.21598767e-04_realType
    up2(15) = -0.10173953e-04_realType
    up2(16) = -0.51044453e-05_realType
    up2(17) = -0.27591627e-05_realType
    up2(18) = -0.15948319e-05_realType
    up2(19) = -0.96856201e-06_realType
    up2(20) = -0.60909779e-06_realType
    up2(21) = -0.39147369e-06_realType
    up2(22) = -0.25632692e-06_realType
    up2(23) = -0.16958665e-06_realType
    up2(24) = -0.11394020e-06_realType
    up2(25) = -0.76500636e-07_realType
    up2(26) = -0.52236558e-07_realType
    up2(27) = -0.35697343e-07_realType
    up2(28) = -0.24471063e-07_realType
    up2(29) = -0.17096052e-07_realType
    up2(30) = -0.11859363e-07_realType
    up2(31) = -0.83689945e-08_realType
    up2(32) = -0.58800329e-08_realType
    up2(33) = -0.42032634e-08_realType
    up2(34) = -0.30343729e-08_realType

    up3(1)  =  0.31659592e+00_realType
    up3(2)  = -0.77365039e-02_realType
    up3(3)  = -0.36297270e-02_realType
    up3(4)  = -0.92844018e-03_realType
    up3(5)  = -0.23630863e-03_realType
    up3(6)  = -0.64433685e-04_realType
    up3(7)  = -0.19446239e-04_realType
    up3(8)  = -0.64919565e-05_realType
    up3(9)  = -0.20373446e-05_realType
    up3(10) = -0.14090103e-06_realType
    up3(11) =  0.45874410e-06_realType
    up3(12) =  0.34517950e-06_realType
    up3(13) =  0.14855464e-06_realType
    up3(14) =  0.50901715e-07_realType
    up3(15) =  0.16140276e-07_realType
    up3(16) =  0.50667962e-08_realType
    up3(17) =  0.16437361e-08_realType
    up3(18) =  0.57675302e-09_realType
    up3(19) =  0.22343487e-09_realType
    up3(20) =  0.97170203e-10_realType
    up3(21) =  0.45068635e-10_realType
    up3(22) =  0.23106070e-10_realType
    up3(23) =  0.11870070e-10_realType
    up3(24) =  0.70527690e-11_realType
    up3(25) =  0.36226728e-11_realType
    up3(26) =  0.22246040e-11_realType
    up3(27) =  0.12643111e-11_realType
    up3(28) =  0.64910588e-12_realType
    up3(29) =  0.42682810e-12_realType
    up3(30) =  0.21768527e-12_realType
    up3(31) =  0.13556644e-12_realType
    up3(32) =  0.63772628e-13_realType
    up3(33) =  0.35175874e-13_realType
    up3(34) =  0.21542623e-13_realType

    ! Set the values of tuLogFit to .false., because a linear
    ! fit has been used.

    tuLogFit(itu1) = .false.

    ! Set the values of constants for the cubic fits of the
    ! non-dimensional spalart-allmaras viscosity.

    tup0(1,itu1)  = 0.14399200e+00_realType
    tup0(2,itu1)  = 0.27285000e+00_realType
    tup0(3,itu1)  = 0.43767100e+00_realType
    tup0(4,itu1)  = 0.64739300e+00_realType
    tup0(5,itu1)  = 0.91283700e+00_realType
    tup0(6,itu1)  = 0.12469600e+01_realType
    tup0(7,itu1)  = 0.16651800e+01_realType
    tup0(8,itu1)  = 0.21861800e+01_realType
    tup0(9,itu1)  = 0.28347900e+01_realType
    tup0(10,itu1) = 0.36492600e+01_realType
    tup0(11,itu1) = 0.46812500e+01_realType
    tup0(12,itu1) = 0.59588800e+01_realType
    tup0(13,itu1) = 0.74961200e+01_realType
    tup0(14,itu1) = 0.93387200e+01_realType
    tup0(15,itu1) = 0.11555500e+02_realType
    tup0(16,itu1) = 0.14225700e+02_realType
    tup0(17,itu1) = 0.17436000e+02_realType
    tup0(18,itu1) = 0.21280700e+02_realType
    tup0(19,itu1) = 0.25863100e+02_realType
    tup0(20,itu1) = 0.31296700e+02_realType
    tup0(21,itu1) = 0.37704800e+02_realType
    tup0(22,itu1) = 0.45218300e+02_realType
    tup0(23,itu1) = 0.53972900e+02_realType
    tup0(24,itu1) = 0.64103200e+02_realType
    tup0(25,itu1) = 0.75735500e+02_realType
    tup0(26,itu1) = 0.88977700e+02_realType
    tup0(27,itu1) = 0.10390800e+03_realType
    tup0(28,itu1) = 0.12056400e+03_realType
    tup0(29,itu1) = 0.13892700e+03_realType
    tup0(30,itu1) = 0.15892000e+03_realType
    tup0(31,itu1) = 0.18039200e+03_realType
    tup0(32,itu1) = 0.20312000e+03_realType
    tup0(33,itu1) = 0.22680300e+03_realType
    tup0(34,itu1) = 0.25105400e+03_realType

    tup1(1,itu1)  = 0.40950260e+00_realType
    tup1(2,itu1)  = 0.40940115e+00_realType
    tup1(3,itu1)  = 0.40919529e+00_realType
    tup1(4,itu1)  = 0.40882583e+00_realType
    tup1(5,itu1)  = 0.40819638e+00_realType
    tup1(6,itu1)  = 0.40720236e+00_realType
    tup1(7,itu1)  = 0.40598943e+00_realType
    tup1(8,itu1)  = 0.40559491e+00_realType
    tup1(9,itu1)  = 0.40881860e+00_realType
    tup1(10,itu1) = 0.41753197e+00_realType
    tup1(11,itu1) = 0.42442441e+00_realType
    tup1(12,itu1) = 0.42212075e+00_realType
    tup1(13,itu1) = 0.41529539e+00_realType
    tup1(14,itu1) = 0.41032022e+00_realType
    tup1(15,itu1) = 0.40794524e+00_realType
    tup1(16,itu1) = 0.40694094e+00_realType
    tup1(17,itu1) = 0.40625126e+00_realType
    tup1(18,itu1) = 0.40527764e+00_realType
    tup1(19,itu1) = 0.40374561e+00_realType
    tup1(20,itu1) = 0.40150883e+00_realType
    tup1(21,itu1) = 0.39842138e+00_realType
    tup1(22,itu1) = 0.39429502e+00_realType
    tup1(23,itu1) = 0.38893832e+00_realType
    tup1(24,itu1) = 0.38209495e+00_realType
    tup1(25,itu1) = 0.37349660e+00_realType
    tup1(26,itu1) = 0.36290738e+00_realType
    tup1(27,itu1) = 0.35013025e+00_realType
    tup1(28,itu1) = 0.33503951e+00_realType
    tup1(29,itu1) = 0.31766382e+00_realType
    tup1(30,itu1) = 0.29813133e+00_realType
    tup1(31,itu1) = 0.27667192e+00_realType
    tup1(32,itu1) = 0.25362172e+00_realType
    tup1(33,itu1) = 0.22930211e+00_realType
    tup1(34,itu1) = 0.20402405e+00_realType

    tup2(1,itu1)  =  0.43946228e-04_realType
    tup2(2,itu1)  =  0.90566870e-04_realType
    tup2(3,itu1)  =  0.34081055e-04_realType
    tup2(4,itu1)  =  0.50034075e-04_realType
    tup2(5,itu1)  = -0.36629521e-04_realType
    tup2(6,itu1)  = -0.33730915e-03_realType
    tup2(7,itu1)  = -0.98768734e-03_realType
    tup2(8,itu1)  = -0.17750541e-02_realType
    tup2(9,itu1)  = -0.61469972e-03_realType
    tup2(10,itu1) =  0.33676511e-02_realType
    tup2(11,itu1) =  0.22772412e-02_realType
    tup2(12,itu1) = -0.68862307e-03_realType
    tup2(13,itu1) = -0.92984441e-03_realType
    tup2(14,itu1) = -0.44253269e-03_realType
    tup2(15,itu1) = -0.14333421e-03_realType
    tup2(16,itu1) = -0.25078749e-04_realType
    tup2(17,itu1) = -0.11675748e-05_realType
    tup2(18,itu1) = -0.77479474e-05_realType
    tup2(19,itu1) = -0.19425991e-04_realType
    tup2(20,itu1) = -0.28782981e-04_realType
    tup2(21,itu1) = -0.37198549e-04_realType
    tup2(22,itu1) = -0.46839769e-04_realType
    tup2(23,itu1) = -0.52777911e-04_realType
    tup2(24,itu1) = -0.61987420e-04_realType
    tup2(25,itu1) = -0.69912435e-04_realType
    tup2(26,itu1) = -0.78150190e-04_realType
    tup2(27,itu1) = -0.84976835e-04_realType
    tup2(28,itu1) = -0.91879236e-04_realType
    tup2(29,itu1) = -0.94706478e-04_realType
    tup2(30,itu1) = -0.96428722e-04_realType
    tup2(31,itu1) = -0.95007410e-04_realType
    tup2(32,itu1) = -0.91049468e-04_realType
    tup2(33,itu1) = -0.85824231e-04_realType
    tup2(34,itu1) = -0.79305036e-04_realType

    tup3(1,itu1)  = -0.43457258e-03_realType
    tup3(2,itu1)  = -0.57319466e-03_realType
    tup3(3,itu1)  = -0.51288659e-03_realType
    tup3(4,itu1)  = -0.54857331e-03_realType
    tup3(5,itu1)  = -0.46390243e-03_realType
    tup3(6,itu1)  = -0.16364047e-03_realType
    tup3(7,itu1)  =  0.43276759e-03_realType
    tup3(8,itu1)  =  0.11606902e-02_realType
    tup3(9,itu1)  =  0.94769940e-03_realType
    tup3(10,itu1) = -0.53409400e-03_realType
    tup3(11,itu1) = -0.59146449e-03_realType
    tup3(12,itu1) = -0.43895639e-04_realType
    tup3(13,itu1) =  0.55678050e-04_realType
    tup3(14,itu1) =  0.27482853e-04_realType
    tup3(15,itu1) =  0.67866418e-05_realType
    tup3(16,itu1) = -0.15708231e-05_realType
    tup3(17,itu1) = -0.35355159e-05_realType
    tup3(18,itu1) = -0.35276554e-05_realType
    tup3(19,itu1) = -0.31393462e-05_realType
    tup3(20,itu1) = -0.28177541e-05_realType
    tup3(21,itu1) = -0.25267299e-05_realType
    tup3(22,itu1) = -0.21840943e-05_realType
    tup3(23,itu1) = -0.19739075e-05_realType
    tup3(24,itu1) = -0.16910644e-05_realType
    tup3(25,itu1) = -0.14435078e-05_realType
    tup3(26,itu1) = -0.11949962e-05_realType
    tup3(27,itu1) = -0.97317619e-06_realType
    tup3(28,itu1) = -0.75009409e-06_realType
    tup3(29,itu1) = -0.58018935e-06_realType
    tup3(30,itu1) = -0.42811291e-06_realType
    tup3(31,itu1) = -0.31260995e-06_realType
    tup3(32,itu1) = -0.22863635e-06_realType
    tup3(33,itu1) = -0.16534708e-06_realType
    tup3(34,itu1) = -0.12169915e-06_realType

  end subroutine initCurveFitDataSa
  subroutine initCurveFitDataSae
    !
    !       initCurveFitDataSae contains the curve fit constants for       
    !       the wall function data for the Spalart-Allmaras (Edwards       
    !       modification) turbulence model.                                
    !
    use flowVarRefState
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    !   integer :: ierr

    call terminate("initCurveFitDataSae", &
         "Not implemented yet")

  end subroutine initCurveFitDataSae

  subroutine initCurveFitDataVf
    !
    !       initCurveFitDataVf contains the curve fit constants for        
    !       the wall function data for the v2-f turbulence model.          
    !
    use constants
    use flowVarRefState
    use inputPhysics
    use paramTurb
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    ! Determine the version of the v2-f model.

    select case (rvfN)

    case (1_intType)

       ! Version 1 of the model.

       ! Set the number of data points and allocate the memory for
       ! the arrays of the curve fits.

       nFit = 34

       allocate(ypT(0:nFit),          reT(0:nFit),          &
            up0(nFit), up1(nFit), up2(nFit), up3(nFit), &
            tup0(nFit,nt1:nt2+1), tup1(nFit,nt1:nt2+1), &
            tup2(nFit,nt1:nt2+1), tup3(nFit,nt1:nt2+1), &
            tuLogFit(nt1:nt2+1),  stat=ierr)
       if(ierr /= 0)                          &
            call terminate("initCurveFitDataVf", &
            "Memory allocation failure for curve fit &
            &coefficients")

       ! Set the values of the Reynolds numbers at interval
       ! boundaries.

       reT( 0) = 0.12795537e+00_realType
       reT( 1) = 0.45951248e+00_realType
       reT( 2) = 0.11826997e+01_realType
       reT( 3) = 0.25888649e+01_realType
       reT( 4) = 0.51497900e+01_realType
       reT( 5) = 0.96108242e+01_realType
       reT( 6) = 0.17107619e+02_realType
       reT( 7) = 0.29176338e+02_realType
       reT( 8) = 0.47402263e+02_realType
       reT( 9) = 0.73016653e+02_realType
       reT(10) = 0.10712161e+03_realType
       reT(11) = 0.15119025e+03_realType
       reT(12) = 0.20726939e+03_realType
       reT(13) = 0.27799459e+03_realType
       reT(14) = 0.36660352e+03_realType
       reT(15) = 0.47698645e+03_realType
       reT(16) = 0.61377565e+03_realType
       reT(17) = 0.78244065e+03_realType
       reT(18) = 0.98942933e+03_realType
       reT(19) = 0.12422738e+04_realType
       reT(20) = 0.15497731e+04_realType
       reT(21) = 0.19221837e+04_realType
       reT(22) = 0.23713440e+04_realType
       reT(23) = 0.29109888e+04_realType
       reT(24) = 0.35569669e+04_realType
       reT(25) = 0.43274783e+04_realType
       reT(26) = 0.52434023e+04_realType
       reT(27) = 0.63287428e+04_realType
       reT(28) = 0.76109132e+04_realType
       reT(29) = 0.91213906e+04_realType
       reT(30) = 0.10896077e+05_realType
       reT(31) = 0.12976000e+05_realType
       reT(32) = 0.15408257e+05_realType
       reT(33) = 0.18246449e+05_realType
       reT(34) = 0.21552478e+05_realType

       ! Set the values of the y+ values at interval boundaries.

       ypT( 0) = 0.35770800e+00_realType
       ypT( 1) = 0.67787500e+00_realType
       ypT( 2) = 0.10875400e+01_realType
       ypT( 3) = 0.16091400e+01_realType
       ypT( 4) = 0.22700900e+01_realType
       ypT( 5) = 0.31037100e+01_realType
       ypT( 6) = 0.41503100e+01_realType
       ypT( 7) = 0.54585000e+01_realType
       ypT( 8) = 0.70866100e+01_realType
       ypT( 9) = 0.91042300e+01_realType
       ypT(10) = 0.11593900e+02_realType
       ypT(11) = 0.14653200e+02_realType
       ypT(12) = 0.18396800e+02_realType
       ypT(13) = 0.22959200e+02_realType
       ypT(14) = 0.28497300e+02_realType
       ypT(15) = 0.35193900e+02_realType
       ypT(16) = 0.43260500e+02_realType
       ypT(17) = 0.52941300e+02_realType
       ypT(18) = 0.64517200e+02_realType
       ypT(19) = 0.78309700e+02_realType
       ypT(20) = 0.94686000e+02_realType
       ypT(21) = 0.11406400e+03_realType
       ypT(22) = 0.13691600e+03_realType
       ypT(23) = 0.16377700e+03_realType
       ypT(24) = 0.19525000e+03_realType
       ypT(25) = 0.23200900e+03_realType
       ypT(26) = 0.27481000e+03_realType
       ypT(27) = 0.32449600e+03_realType
       ypT(28) = 0.38200300e+03_realType
       ypT(29) = 0.44837100e+03_realType
       ypT(30) = 0.52474800e+03_realType
       ypT(31) = 0.61239900e+03_realType
       ypT(32) = 0.71271500e+03_realType
       ypT(33) = 0.82722200e+03_realType
       ypT(34) = 0.95759000e+03_realType

       ! Set the values of constants for the cubic fits of the
       ! non-dimensional tangential velocity.

       up0( 1) = 0.35770900e+00_realType
       up0( 2) = 0.67787200e+00_realType
       up0( 3) = 0.10875000e+01_realType
       up0( 4) = 0.16088500e+01_realType
       up0( 5) = 0.22685400e+01_realType
       up0( 6) = 0.30965600e+01_realType
       up0( 7) = 0.41220100e+01_realType
       up0( 8) = 0.53451200e+01_realType
       up0( 9) = 0.66889900e+01_realType
       up0(10) = 0.80200800e+01_realType
       up0(11) = 0.92394800e+01_realType
       up0(12) = 0.10317900e+02_realType
       up0(13) = 0.11266600e+02_realType
       up0(14) = 0.12108200e+02_realType
       up0(15) = 0.12864500e+02_realType
       up0(16) = 0.13553100e+02_realType
       up0(17) = 0.14187900e+02_realType
       up0(18) = 0.14779400e+02_realType
       up0(19) = 0.15335900e+02_realType
       up0(20) = 0.15863600e+02_realType
       up0(21) = 0.16367500e+02_realType
       up0(22) = 0.16851800e+02_realType
       up0(23) = 0.17319700e+02_realType
       up0(24) = 0.17774100e+02_realType
       up0(25) = 0.18217500e+02_realType
       up0(26) = 0.18652200e+02_realType
       up0(27) = 0.19080100e+02_realType
       up0(28) = 0.19503300e+02_realType
       up0(29) = 0.19923700e+02_realType
       up0(30) = 0.20343400e+02_realType
       up0(31) = 0.20764400e+02_realType
       up0(32) = 0.21188800e+02_realType
       up0(33) = 0.21619100e+02_realType
       up0(34) = 0.22057500e+02_realType

       up1( 1) = 0.12712722e+01_realType
       up1( 2) = 0.69191267e+00_realType
       up1( 3) = 0.43721180e+00_realType
       up1( 4) = 0.29770939e+00_realType
       up1( 5) = 0.21186537e+00_realType
       up1( 6) = 0.15500054e+00_realType
       up1( 7) = 0.11492466e+00_realType
       up1( 8) = 0.84733790e-01_realType
       up1( 9) = 0.61015984e-01_realType
       up1(10) = 0.42707937e-01_realType
       up1(11) = 0.29393811e-01_realType
       up1(12) = 0.20241287e-01_realType
       up1(13) = 0.14118603e-01_realType
       up1(14) = 0.10028611e-01_realType
       up1(15) = 0.72611010e-02_realType
       up1(16) = 0.53541635e-02_realType
       up1(17) = 0.40146771e-02_realType
       up1(18) = 0.30560063e-02_realType
       up1(19) = 0.23578120e-02_realType
       up1(20) = 0.18410127e-02_realType
       up1(21) = 0.14534277e-02_realType
       up1(22) = 0.11589991e-02_realType
       up1(23) = 0.93274199e-03_realType
       up1(24) = 0.75723913e-03_realType
       up1(25) = 0.61991282e-03_realType
       up1(26) = 0.51149306e-03_realType
       up1(27) = 0.42528110e-03_realType
       up1(28) = 0.35632360e-03_realType
       up1(29) = 0.30082562e-03_realType
       up1(30) = 0.25590806e-03_realType
       up1(31) = 0.21932184e-03_realType
       up1(32) = 0.18942066e-03_realType
       up1(33) = 0.16482466e-03_realType
       up1(34) = 0.14450977e-03_realType

       up2( 1) = -0.10180856e+01_realType
       up2( 2) = -0.16838796e+00_realType
       up2( 3) = -0.42564371e-01_realType
       up2( 4) = -0.13467470e-01_realType
       up2( 5) = -0.49083413e-02_realType
       up2( 6) = -0.19435488e-02_realType
       up2( 7) = -0.87388594e-03_realType
       up2( 8) = -0.50925439e-03_realType
       up2( 9) = -0.34513516e-03_realType
       up2(10) = -0.22127804e-03_realType
       up2(11) = -0.12741049e-03_realType
       up2(12) = -0.68647301e-04_realType
       up2(13) = -0.36296584e-04_realType
       up2(14) = -0.19327079e-04_realType
       up2(15) = -0.10522599e-04_realType
       up2(16) = -0.58546470e-05_realType
       up2(17) = -0.33469044e-05_realType
       up2(18) = -0.19525994e-05_realType
       up2(19) = -0.11686013e-05_realType
       up2(20) = -0.71331900e-06_realType
       up2(21) = -0.44175896e-06_realType
       up2(22) = -0.27957847e-06_realType
       up2(23) = -0.17903765e-06_realType
       up2(24) = -0.11639426e-06_realType
       up2(25) = -0.76321282e-07_realType
       up2(26) = -0.51021408e-07_realType
       up2(27) = -0.34196379e-07_realType
       up2(28) = -0.23263805e-07_realType
       up2(29) = -0.15876970e-07_realType
       up2(30) = -0.10967048e-07_realType
       up2(31) = -0.76570908e-08_realType
       up2(32) = -0.53137237e-08_realType
       up2(33) = -0.37931015e-08_realType
       up2(34) = -0.26536226e-08_realType

       up3( 1) =  0.29032828e+00_realType
       up3( 2) = -0.71056794e-02_realType
       up3( 3) = -0.33374366e-02_realType
       up3( 4) = -0.85721208e-03_realType
       up3( 5) = -0.21895691e-03_realType
       up3( 6) = -0.64856031e-04_realType
       up3( 7) = -0.20819908e-04_realType
       up3( 8) = -0.51723681e-05_realType
       up3( 9) = -0.31864903e-06_realType
       up3(10) =  0.50988875e-06_realType
       up3(11) =  0.35651270e-06_realType
       up3(12) =  0.16711639e-06_realType
       up3(13) =  0.69583056e-07_realType
       up3(14) =  0.27917810e-07_realType
       up3(15) =  0.11383228e-07_realType
       up3(16) =  0.46713098e-08_realType
       up3(17) =  0.19959356e-08_realType
       up3(18) =  0.85688237e-09_realType
       up3(19) =  0.38661576e-09_realType
       up3(20) =  0.18015921e-09_realType
       up3(21) =  0.83166064e-10_realType
       up3(22) =  0.41131094e-10_realType
       up3(23) =  0.20294857e-10_realType
       up3(24) =  0.10424447e-10_realType
       up3(25) =  0.51615620e-11_realType
       up3(26) =  0.28812993e-11_realType
       up3(27) =  0.14918197e-11_realType
       up3(28) =  0.84315572e-12_realType
       up3(29) =  0.44502764e-12_realType
       up3(30) =  0.24764964e-12_realType
       up3(31) =  0.15033872e-12_realType
       up3(32) =  0.70583114e-13_realType
       up3(33) =  0.50327579e-13_realType
       up3(34) =  0.20760485e-13_realType

       ! Set the values of tuLogFit.All variables have been
       ! fitted linearly; the fifth variable is the eddy viscosity.

       tuLogFit(itu1) = .false.
       tuLogFit(itu2) = .false.
       tuLogFit(itu3) = .false.
       tuLogFit(itu4) = .false.
       tuLogFit(itu5) = .false.

       ! Set the values of constants for the cubic fits of the
       ! non-dimensional k, eps, v2 and f values.

       ! Constants for k.

       tup0( 1,itu1) = 0.23837900e-01_realType
       tup0( 2,itu1) = 0.78315200e-01_realType
       tup0( 3,itu1) = 0.18685000e+00_realType
       tup0( 4,itu1) = 0.38036300e+00_realType
       tup0( 5,itu1) = 0.70148100e+00_realType
       tup0( 6,itu1) = 0.12056400e+01_realType
       tup0( 7,itu1) = 0.19570000e+01_realType
       tup0( 8,itu1) = 0.29949800e+01_realType
       tup0( 9,itu1) = 0.42181500e+01_realType
       tup0(10,itu1) = 0.53378900e+01_realType
       tup0(11,itu1) = 0.61118500e+01_realType
       tup0(12,itu1) = 0.65046800e+01_realType
       tup0(13,itu1) = 0.66026900e+01_realType
       tup0(14,itu1) = 0.65082900e+01_realType
       tup0(15,itu1) = 0.63012300e+01_realType
       tup0(16,itu1) = 0.60354500e+01_realType
       tup0(17,itu1) = 0.57452900e+01_realType
       tup0(18,itu1) = 0.54518200e+01_realType
       tup0(19,itu1) = 0.51675400e+01_realType
       tup0(20,itu1) = 0.48994200e+01_realType
       tup0(21,itu1) = 0.46509700e+01_realType
       tup0(22,itu1) = 0.44234600e+01_realType
       tup0(23,itu1) = 0.42167200e+01_realType
       tup0(24,itu1) = 0.40297300e+01_realType
       tup0(25,itu1) = 0.38609000e+01_realType
       tup0(26,itu1) = 0.37082900e+01_realType
       tup0(27,itu1) = 0.35697500e+01_realType
       tup0(28,itu1) = 0.34429600e+01_realType
       tup0(29,itu1) = 0.33254900e+01_realType
       tup0(30,itu1) = 0.32147800e+01_realType
       tup0(31,itu1) = 0.31080800e+01_realType
       tup0(32,itu1) = 0.30023700e+01_realType
       tup0(33,itu1) = 0.28942800e+01_realType
       tup0(34,itu1) = 0.27799400e+01_realType

       tup1( 1,itu1) =  0.13413856e+00_realType
       tup1( 2,itu1) =  0.22335565e+00_realType
       tup1( 3,itu1) =  0.32434141e+00_realType
       tup1( 4,itu1) =  0.43518752e+00_realType
       tup1( 5,itu1) =  0.55218357e+00_realType
       tup1( 6,itu1) =  0.66775111e+00_realType
       tup1( 7,itu1) =  0.75987243e+00_realType
       tup1( 8,itu1) =  0.77006777e+00_realType
       tup1( 9,itu1) =  0.64264496e+00_realType
       tup1(10,itu1) =  0.42014159e+00_realType
       tup1(11,itu1) =  0.21027146e+00_realType
       tup1(12,itu1) =  0.72151582e-01_realType
       tup1(13,itu1) =  0.43462557e-03_realType
       tup1(14,itu1) = -0.29846047e-01_realType
       tup1(15,itu1) = -0.38647454e-01_realType
       tup1(16,itu1) = -0.37657148e-01_realType
       tup1(17,itu1) = -0.32885380e-01_realType
       tup1(18,itu1) = -0.27179666e-01_realType
       tup1(19,itu1) = -0.21775122e-01_realType
       tup1(20,itu1) = -0.17122657e-01_realType
       tup1(21,itu1) = -0.13311965e-01_realType
       tup1(22,itu1) = -0.10282974e-01_realType
       tup1(23,itu1) = -0.79200612e-02_realType
       tup1(24,itu1) = -0.60997017e-02_realType
       tup1(25,itu1) = -0.47109860e-02_realType
       tup1(26,itu1) = -0.36595023e-02_realType
       tup1(27,itu1) = -0.28688356e-02_realType
       tup1(28,itu1) = -0.22786936e-02_realType
       tup1(29,itu1) = -0.18420182e-02_realType
       tup1(30,itu1) = -0.15230656e-02_realType
       tup1(31,itu1) = -0.12949618e-02_realType
       tup1(32,itu1) = -0.11374337e-02_realType
       tup1(33,itu1) = -0.10354105e-02_realType
       tup1(34,itu1) = -0.97825421e-03_realType

       tup2( 1,itu1) =  0.58799048e-01_realType
       tup2( 2,itu1) =  0.57983435e-01_realType
       tup2( 3,itu1) =  0.55840141e-01_realType
       tup2( 4,itu1) =  0.52909760e-01_realType
       tup2( 5,itu1) =  0.50658583e-01_realType
       tup2( 6,itu1) =  0.55744447e-01_realType
       tup2( 7,itu1) =  0.69202235e-01_realType
       tup2( 8,itu1) =  0.43649374e-01_realType
       tup2( 9,itu1) = -0.20068023e-01_realType
       tup2(10,itu1) = -0.47375401e-01_realType
       tup2(11,itu1) = -0.35131867e-01_realType
       tup2(12,itu1) = -0.17682377e-01_realType
       tup2(13,itu1) = -0.72540324e-02_realType
       tup2(14,itu1) = -0.24963930e-02_realType
       tup2(15,itu1) = -0.61439197e-03_realType
       tup2(16,itu1) =  0.35707704e-04_realType
       tup2(17,itu1) =  0.20726554e-03_realType
       tup2(18,itu1) =  0.21257094e-03_realType
       tup2(19,itu1) =  0.17069050e-03_realType
       tup2(20,itu1) =  0.12477387e-03_realType
       tup2(21,itu1) =  0.86954464e-04_realType
       tup2(22,itu1) =  0.58869152e-04_realType
       tup2(23,itu1) =  0.39300193e-04_realType
       tup2(24,itu1) =  0.25976170e-04_realType
       tup2(25,itu1) =  0.17045223e-04_realType
       tup2(26,itu1) =  0.11152061e-04_realType
       tup2(27,itu1) =  0.72633824e-05_realType
       tup2(28,itu1) =  0.47173696e-05_realType
       tup2(29,itu1) =  0.30546703e-05_realType
       tup2(30,itu1) =  0.19644758e-05_realType
       tup2(31,itu1) =  0.12465154e-05_realType
       tup2(32,itu1) =  0.77547548e-06_realType
       tup2(33,itu1) =  0.46678879e-06_realType
       tup2(34,itu1) =  0.25787187e-06_realType

       tup3( 1,itu1) =  0.16768319e+00_realType
       tup3( 2,itu1) =  0.10621790e+00_realType
       tup3( 3,itu1) =  0.64437267e-01_realType
       tup3( 4,itu1) =  0.35904096e-01_realType
       tup3( 5,itu1) =  0.14921346e-01_realType
       tup3( 6,itu1) = -0.74747729e-02_realType
       tup3( 7,itu1) = -0.33280334e-01_realType
       tup3( 8,itu1) = -0.33896768e-01_realType
       tup3( 9,itu1) = -0.11588583e-01_realType
       tup3(10,itu1) =  0.13997083e-02_realType
       tup3(11,itu1) =  0.27365948e-02_realType
       tup3(12,itu1) =  0.14431329e-02_realType
       tup3(13,itu1) =  0.57506718e-03_realType
       tup3(14,itu1) =  0.20485602e-03_realType
       tup3(15,itu1) =  0.68525613e-04_realType
       tup3(16,itu1) =  0.21493188e-04_realType
       tup3(17,itu1) =  0.60206301e-05_realType
       tup3(18,itu1) =  0.12018375e-05_realType
       tup3(19,itu1) = -0.98181484e-07_realType
       tup3(20,itu1) = -0.34302113e-06_realType
       tup3(21,itu1) = -0.30271628e-06_realType
       tup3(22,itu1) = -0.20913696e-06_realType
       tup3(23,itu1) = -0.13440526e-06_realType
       tup3(24,itu1) = -0.82910301e-07_realType
       tup3(25,itu1) = -0.49744260e-07_realType
       tup3(26,itu1) = -0.29836016e-07_realType
       tup3(27,itu1) = -0.17773844e-07_realType
       tup3(28,itu1) = -0.10672954e-07_realType
       tup3(29,itu1) = -0.65469411e-08_realType
       tup3(30,itu1) = -0.41129357e-08_realType
       tup3(31,itu1) = -0.26461407e-08_realType
       tup3(32,itu1) = -0.17741710e-08_realType
       tup3(33,itu1) = -0.12646276e-08_realType
       tup3(34,itu1) = -0.92958789e-09_realType

       ! Constants for epsilon.

       tup0( 1,itu2) = 0.29354900e+00_realType
       tup0( 2,itu2) = 0.26347600e+00_realType
       tup0( 3,itu2) = 0.23131100e+00_realType
       tup0( 4,itu2) = 0.19880500e+00_realType
       tup0( 5,itu2) = 0.16843900e+00_realType
       tup0( 6,itu2) = 0.14348600e+00_realType
       tup0( 7,itu2) = 0.12800300e+00_realType
       tup0( 8,itu2) = 0.12650200e+00_realType
       tup0( 9,itu2) = 0.13424600e+00_realType
       tup0(10,itu2) = 0.14121500e+00_realType
       tup0(11,itu2) = 0.14022500e+00_realType
       tup0(12,itu2) = 0.13108600e+00_realType
       tup0(13,itu2) = 0.11710400e+00_realType
       tup0(14,itu2) = 0.10147100e+00_realType
       tup0(15,itu2) = 0.86210400e-01_realType
       tup0(16,itu2) = 0.72337500e-01_realType
       tup0(17,itu2) = 0.60231700e-01_realType
       tup0(18,itu2) = 0.49927500e-01_realType
       tup0(19,itu2) = 0.41292900e-01_realType
       tup0(20,itu2) = 0.34129200e-01_realType
       tup0(21,itu2) = 0.28223300e-01_realType
       tup0(22,itu2) = 0.23372700e-01_realType
       tup0(23,itu2) = 0.19397100e-01_realType
       tup0(24,itu2) = 0.16140700e-01_realType
       tup0(25,itu2) = 0.13472400e-01_realType
       tup0(26,itu2) = 0.11283300e-01_realType
       tup0(27,itu2) = 0.94837500e-02_realType
       tup0(28,itu2) = 0.80004600e-02_realType
       tup0(29,itu2) = 0.67736800e-02_realType
       tup0(30,itu2) = 0.57548200e-02_realType
       tup0(31,itu2) = 0.49043300e-02_realType
       tup0(32,itu2) = 0.41899500e-02_realType
       tup0(33,itu2) = 0.35852000e-02_realType
       tup0(34,itu2) = 0.30682200e-02_realType

       tup1( 1,itu2) = -0.10010771e+00_realType
       tup1( 2,itu2) = -0.85277160e-01_realType
       tup1( 3,itu2) = -0.69444251e-01_realType
       tup1( 4,itu2) = -0.53166462e-01_realType
       tup1( 5,itu2) = -0.37013322e-01_realType
       tup1( 6,itu2) = -0.21505994e-01_realType
       tup1( 7,itu2) = -0.72125328e-02_realType
       tup1( 8,itu2) =  0.21261451e-02_realType
       tup1( 9,itu2) =  0.40356801e-02_realType
       tup1(10,itu2) =  0.13265177e-02_realType
       tup1(11,itu2) = -0.18253838e-02_realType
       tup1(12,itu2) = -0.33986976e-02_realType
       tup1(13,itu2) = -0.35654948e-02_realType
       tup1(14,itu2) = -0.30586209e-02_realType
       tup1(15,itu2) = -0.23812190e-02_realType
       tup1(16,itu2) = -0.17596930e-02_realType
       tup1(17,itu2) = -0.12627202e-02_realType
       tup1(18,itu2) = -0.89095673e-03_realType
       tup1(19,itu2) = -0.62275508e-03_realType
       tup1(20,itu2) = -0.43321577e-03_realType
       tup1(21,itu2) = -0.30084493e-03_realType
       tup1(22,itu2) = -0.20900308e-03_realType
       tup1(23,itu2) = -0.14547503e-03_realType
       tup1(24,itu2) = -0.10156512e-03_realType
       tup1(25,itu2) = -0.71189471e-04_realType
       tup1(26,itu2) = -0.50133861e-04_realType
       tup1(27,itu2) = -0.35495151e-04_realType
       tup1(28,itu2) = -0.25282155e-04_realType
       tup1(29,itu2) = -0.18128274e-04_realType
       tup1(30,itu2) = -0.13095730e-04_realType
       tup1(31,itu2) = -0.95402614e-05_realType
       tup1(32,itu2) = -0.70178808e-05_realType
       tup1(33,itu2) = -0.52216476e-05_realType
       tup1(34,itu2) = -0.39388668e-05_realType

       tup2( 1,itu2) =  0.11572886e-01_realType
       tup2( 2,itu2) =  0.10868522e-01_realType
       tup2( 3,itu2) =  0.97691864e-02_realType
       tup2( 4,itu2) =  0.83476261e-02_realType
       tup2( 5,itu2) =  0.68769083e-02_realType
       tup2( 6,itu2) =  0.55834783e-02_realType
       tup2( 7,itu2) =  0.67702397e-02_realType
       tup2( 8,itu2) =  0.36737909e-02_realType
       tup2( 9,itu2) =  0.47795474e-03_realType
       tup2(10,itu2) = -0.81158579e-03_realType
       tup2(11,itu2) = -0.62510665e-03_realType
       tup2(12,itu2) = -0.22487250e-03_realType
       tup2(13,itu2) = -0.19693072e-04_realType
       tup2(14,itu2) =  0.41848783e-04_realType
       tup2(15,itu2) =  0.45878641e-04_realType
       tup2(16,itu2) =  0.34700119e-04_realType
       tup2(17,itu2) =  0.23057037e-04_realType
       tup2(18,itu2) =  0.14420758e-04_realType
       tup2(19,itu2) =  0.87404760e-05_realType
       tup2(20,itu2) =  0.52127418e-05_realType
       tup2(21,itu2) =  0.30833200e-05_realType
       tup2(22,itu2) =  0.18189328e-05_realType
       tup2(23,itu2) =  0.10729534e-05_realType
       tup2(24,itu2) =  0.63476327e-06_realType
       tup2(25,itu2) =  0.37690133e-06_realType
       tup2(26,itu2) =  0.22497460e-06_realType
       tup2(27,itu2) =  0.13510084e-06_realType
       tup2(28,itu2) =  0.81632875e-07_realType
       tup2(29,itu2) =  0.49681361e-07_realType
       tup2(30,itu2) =  0.30447133e-07_realType
       tup2(31,itu2) =  0.18796936e-07_realType
       tup2(32,itu2) =  0.11683669e-07_realType
       tup2(33,itu2) =  0.73153823e-08_realType
       tup2(34,itu2) =  0.46047384e-08_realType

       tup3( 1,itu2) =  0.24128568e-01_realType
       tup3( 2,itu2) =  0.13760338e-01_realType
       tup3( 3,itu2) =  0.74572108e-02_realType
       tup3( 4,itu2) =  0.39055025e-02_realType
       tup3( 5,itu2) =  0.19387640e-02_realType
       tup3( 6,itu2) =  0.79307164e-03_realType
       tup3( 7,itu2) = -0.16312244e-02_realType
       tup3( 8,itu2) = -0.12641914e-02_realType
       tup3( 9,itu2) = -0.37976445e-03_realType
       tup3(10,itu2) =  0.47821584e-04_realType
       tup3(11,itu2) =  0.80186189e-04_realType
       tup3(12,itu2) =  0.36078436e-04_realType
       tup3(13,itu2) =  0.10994524e-04_realType
       tup3(14,itu2) =  0.23244487e-05_realType
       tup3(15,itu2) =  0.52508614e-07_realType
       tup3(16,itu2) = -0.32196650e-06_realType
       tup3(17,itu2) = -0.26554053e-06_realType
       tup3(18,itu2) = -0.16334363e-06_realType
       tup3(19,itu2) = -0.90357098e-07_realType
       tup3(20,itu2) = -0.47678861e-07_realType
       tup3(21,itu2) = -0.24549302e-07_realType
       tup3(22,itu2) = -0.12513680e-07_realType
       tup3(23,itu2) = -0.63437513e-08_realType
       tup3(24,itu2) = -0.32238513e-08_realType
       tup3(25,itu2) = -0.16413320e-08_realType
       tup3(26,itu2) = -0.84056640e-09_realType
       tup3(27,itu2) = -0.43372978e-09_realType
       tup3(28,itu2) = -0.22528033e-09_realType
       tup3(29,itu2) = -0.11820379e-09_realType
       tup3(30,itu2) = -0.62596081e-10_realType
       tup3(31,itu2) = -0.33528145e-10_realType
       tup3(32,itu2) = -0.18147943e-10_realType
       tup3(33,itu2) = -0.99793712e-11_realType
       tup3(34,itu2) = -0.55553743e-11_realType

       ! Constants for v2.

       tup0( 1,itu3) = 0.47502000e-06_realType
       tup0( 2,itu3) = 0.89181800e-05_realType
       tup0( 3,itu3) = 0.62272000e-04_realType
       tup0( 4,itu3) = 0.29003500e-03_realType
       tup0( 5,itu3) = 0.10734400e-02_realType
       tup0( 6,itu3) = 0.34009200e-02_realType
       tup0( 7,itu3) = 0.95514600e-02_realType
       tup0( 8,itu3) = 0.23917300e-01_realType
       tup0( 9,itu3) = 0.52481200e-01_realType
       tup0(10,itu3) = 0.99458700e-01_realType
       tup0(11,itu3) = 0.16446900e+00_realType
       tup0(12,itu3) = 0.24385400e+00_realType
       tup0(13,itu3) = 0.33343500e+00_realType
       tup0(14,itu3) = 0.42968800e+00_realType
       tup0(15,itu3) = 0.52978700e+00_realType
       tup0(16,itu3) = 0.63143400e+00_realType
       tup0(17,itu3) = 0.73270600e+00_realType
       tup0(18,itu3) = 0.83194900e+00_realType
       tup0(19,itu3) = 0.92772000e+00_realType
       tup0(20,itu3) = 0.10187500e+01_realType
       tup0(21,itu3) = 0.11039500e+01_realType
       tup0(22,itu3) = 0.11823500e+01_realType
       tup0(23,itu3) = 0.12531300e+01_realType
       tup0(24,itu3) = 0.13155900e+01_realType
       tup0(25,itu3) = 0.13691300e+01_realType
       tup0(26,itu3) = 0.14132300e+01_realType
       tup0(27,itu3) = 0.14474100e+01_realType
       tup0(28,itu3) = 0.14712200e+01_realType
       tup0(29,itu3) = 0.14842200e+01_realType
       tup0(30,itu3) = 0.14859400e+01_realType
       tup0(31,itu3) = 0.14758600e+01_realType
       tup0(32,itu3) = 0.14533600e+01_realType
       tup0(33,itu3) = 0.14177100e+01_realType
       tup0(34,itu3) = 0.13680600e+01_realType

       tup1( 1,itu3) =  0.15655101e-04_realType
       tup1( 2,itu3) =  0.84672884e-04_realType
       tup1( 3,itu3) =  0.30186555e-03_realType
       tup1( 4,itu3) =  0.85507420e-03_realType
       tup1( 5,itu3) =  0.20814582e-02_realType
       tup1( 6,itu3) =  0.45090575e-02_realType
       tup1( 7,itu3) =  0.87126156e-02_realType
       tup1( 8,itu3) =  0.14620352e-01_realType
       tup1( 9,itu3) =  0.20720514e-01_realType
       tup1(10,itu3) =  0.24845927e-01_realType
       tup1(11,itu3) =  0.26022000e-01_realType
       tup1(12,itu3) =  0.24837349e-01_realType
       tup1(13,itu3) =  0.22373465e-01_realType
       tup1(14,itu3) =  0.19439830e-01_realType
       tup1(15,itu3) =  0.16489656e-01_realType
       tup1(16,itu3) =  0.13744920e-01_realType
       tup1(17,itu3) =  0.11298275e-01_realType
       tup1(18,itu3) =  0.91742368e-02_realType
       tup1(19,itu3) =  0.73635310e-02_realType
       tup1(20,itu3) =  0.58414654e-02_realType
       tup1(21,itu3) =  0.45756734e-02_realType
       tup1(22,itu3) =  0.35325598e-02_realType
       tup1(23,itu3) =  0.26801843e-02_realType
       tup1(24,itu3) =  0.19885487e-02_realType
       tup1(25,itu3) =  0.14310001e-02_realType
       tup1(26,itu3) =  0.98391151e-03_realType
       tup1(27,itu3) =  0.62700704e-03_realType
       tup1(28,itu3) =  0.34339929e-03_realType
       tup1(29,itu3) =  0.11882947e-03_realType
       tup1(30,itu3) = -0.58565974e-04_realType
       tup1(31,itu3) = -0.19862463e-03_realType
       tup1(32,itu3) = -0.30936281e-03_realType
       tup1(33,itu3) = -0.39707108e-03_realType
       tup1(34,itu3) = -0.46709546e-03_realType

       tup2( 1,itu3) = -0.11515786e-03_realType
       tup2( 2,itu3) = -0.19649780e-03_realType
       tup2( 3,itu3) = -0.28531197e-03_realType
       tup2( 4,itu3) = -0.35674865e-03_realType
       tup2( 5,itu3) = -0.35499057e-03_realType
       tup2( 6,itu3) = -0.96188825e-04_realType
       tup2( 7,itu3) =  0.68706040e-03_realType
       tup2( 8,itu3) =  0.16408011e-02_realType
       tup2( 9,itu3) =  0.17663922e-02_realType
       tup2(10,itu3) =  0.10532277e-02_realType
       tup2(11,itu3) =  0.31539555e-03_realType
       tup2(12,itu3) = -0.69677717e-04_realType
       tup2(13,itu3) = -0.19632838e-03_realType
       tup2(14,itu3) = -0.20683711e-03_realType
       tup2(15,itu3) = -0.17733514e-03_realType
       tup2(16,itu3) = -0.13942213e-03_realType
       tup2(17,itu3) = -0.10497064e-03_realType
       tup2(18,itu3) = -0.77063300e-04_realType
       tup2(19,itu3) = -0.55728557e-04_realType
       tup2(20,itu3) = -0.39733207e-04_realType
       tup2(21,itu3) = -0.28198531e-04_realType
       tup2(22,itu3) = -0.19838008e-04_realType
       tup2(23,itu3) = -0.13886459e-04_realType
       tup2(24,itu3) = -0.96805370e-05_realType
       tup2(25,itu3) = -0.67138161e-05_realType
       tup2(26,itu3) = -0.46515653e-05_realType
       tup2(27,itu3) = -0.32158971e-05_realType
       tup2(28,itu3) = -0.22162473e-05_realType
       tup2(29,itu3) = -0.15270107e-05_realType
       tup2(30,itu3) = -0.10497157e-05_realType
       tup2(31,itu3) = -0.72432186e-06_realType
       tup2(32,itu3) = -0.50175766e-06_realType
       tup2(33,itu3) = -0.34544811e-06_realType
       tup2(34,itu3) = -0.24279564e-06_realType

       tup3( 1,itu3) =  0.46422002e-03_realType
       tup3( 2,itu3) =  0.75115602e-03_realType
       tup3( 3,itu3) =  0.10424485e-02_realType
       tup3( 4,itu3) =  0.12956018e-02_realType
       tup3( 5,itu3) =  0.14483411e-02_realType
       tup3( 6,itu3) =  0.13404586e-02_realType
       tup3( 7,itu3) =  0.80055685e-03_realType
       tup3( 8,itu3) =  0.95237809e-04_realType
       tup3( 9,itu3) = -0.24584932e-03_realType
       tup3(10,itu3) = -0.21878056e-03_realType
       tup3(11,itu3) = -0.11092088e-03_realType
       tup3(12,itu3) = -0.46194706e-04_realType
       tup3(13,itu3) = -0.18290511e-04_realType
       tup3(14,itu3) = -0.71643534e-05_realType
       tup3(15,itu3) = -0.27476585e-05_realType
       tup3(16,itu3) = -0.10108088e-05_realType
       tup3(17,itu3) = -0.32593663e-06_realType
       tup3(18,itu3) = -0.66048306e-07_realType
       tup3(19,itu3) =  0.26645554e-07_realType
       tup3(20,itu3) =  0.44215925e-07_realType
       tup3(21,itu3) =  0.44161326e-07_realType
       tup3(22,itu3) =  0.34660614e-07_realType
       tup3(23,itu3) =  0.25119801e-07_realType
       tup3(24,itu3) =  0.17432238e-07_realType
       tup3(25,itu3) =  0.11470519e-07_realType
       tup3(26,itu3) =  0.75110077e-08_realType
       tup3(27,itu3) =  0.48557801e-08_realType
       tup3(28,itu3) =  0.30570649e-08_realType
       tup3(29,itu3) =  0.19141510e-08_realType
       tup3(30,itu3) =  0.11593848e-08_realType
       tup3(31,itu3) =  0.70447835e-09_realType
       tup3(32,itu3) =  0.42929522e-09_realType
       tup3(33,itu3) =  0.23103965e-09_realType
       tup3(34,itu3) =  0.13124926e-09_realType

       ! Constants for f.

       tup0( 1,itu4) = -0.33990e-02_realType
       tup0( 2,itu4) = -0.32796e-02_realType
       tup0( 3,itu4) = -0.31317e-02_realType
       tup0( 4,itu4) = -0.29501e-02_realType
       tup0( 5,itu4) = -0.27297e-02_realType
       tup0( 6,itu4) = -0.24652e-02_realType
       tup0( 7,itu4) = -0.21525e-02_realType
       tup0( 8,itu4) = -0.17911e-02_realType
       tup0( 9,itu4) = -0.13852e-02_realType
       tup0(10,itu4) = -0.95384e-03_realType
       tup0(11,itu4) = -0.52641e-03_realType
       tup0(12,itu4) = -0.13576e-03_realType
       tup0(13,itu4) =  0.21413e-03_realType
       tup0(14,itu4) =  0.52393e-03_realType
       tup0(15,itu4) =  0.79302e-03_realType
       tup0(16,itu4) =  0.10205e-02_realType
       tup0(17,itu4) =  0.12056e-02_realType
       tup0(18,itu4) =  0.13485e-02_realType
       tup0(19,itu4) =  0.14505e-02_realType
       tup0(20,itu4) =  0.15137e-02_realType
       tup0(21,itu4) =  0.15414e-02_realType
       tup0(22,itu4) =  0.15379e-02_realType
       tup0(23,itu4) =  0.15079e-02_realType
       tup0(24,itu4) =  0.14563e-02_realType
       tup0(25,itu4) =  0.13881e-02_realType
       tup0(26,itu4) =  0.13079e-02_realType
       tup0(27,itu4) =  0.12200e-02_realType
       tup0(28,itu4) =  0.11280e-02_realType
       tup0(29,itu4) =  0.10348e-02_realType
       tup0(30,itu4) =  0.94297e-03_realType
       tup0(31,itu4) =  0.85414e-03_realType
       tup0(32,itu4) =  0.76956e-03_realType
       tup0(33,itu4) =  0.68990e-03_realType
       tup0(34,itu4) =  0.61535e-03_realType

       tup1( 1,itu4) =  0.36489e-03_realType
       tup1( 2,itu4) =  0.35429e-03_realType
       tup1( 3,itu4) =  0.34223e-03_realType
       tup1( 4,itu4) =  0.32877e-03_realType
       tup1( 5,itu4) =  0.31384e-03_realType
       tup1( 6,itu4) =  0.29698e-03_realType
       tup1( 7,itu4) =  0.27691e-03_realType
       tup1( 8,itu4) =  0.25273e-03_realType
       tup1( 9,itu4) =  0.22211e-03_realType
       tup1(10,itu4) =  0.18427e-03_realType
       tup1(11,itu4) =  0.14258e-03_realType
       tup1(12,itu4) =  0.10527e-03_realType
       tup1(13,itu4) =  0.76806e-04_realType
       tup1(14,itu4) =  0.55424e-04_realType
       tup1(15,itu4) =  0.39245e-04_realType
       tup1(16,itu4) =  0.27022e-04_realType
       tup1(17,itu4) =  0.17875e-04_realType
       tup1(18,itu4) =  0.11139e-04_realType
       tup1(19,itu4) =  0.62935e-05_realType
       tup1(20,itu4) =  0.29155e-05_realType
       tup1(21,itu4) =  0.65684e-06_realType
       tup1(22,itu4) = -0.76742e-06_realType
       tup1(23,itu4) = -0.15875e-05_realType
       tup1(24,itu4) = -0.19855e-05_realType
       tup1(25,itu4) = -0.21022e-05_realType
       tup1(26,itu4) = -0.20422e-05_realType
       tup1(27,itu4) = -0.18806e-05_realType
       tup1(28,itu4) = -0.16696e-05_realType
       tup1(29,itu4) = -0.14434e-05_realType
       tup1(30,itu4) = -0.12233e-05_realType
       tup1(31,itu4) = -0.10215e-05_realType
       tup1(32,itu4) = -0.84412e-06_realType
       tup1(33,itu4) = -0.69337e-06_realType
       tup1(34,itu4) = -0.56947e-06_realType

       tup2( 1,itu4) = -0.66773e-05_realType
       tup2( 2,itu4) = -0.65637e-05_realType
       tup2( 3,itu4) = -0.61090e-05_realType
       tup2( 4,itu4) = -0.54483e-05_realType
       tup2( 5,itu4) = -0.44808e-05_realType
       tup2( 6,itu4) = -0.34229e-05_realType
       tup2( 7,itu4) = -0.36764e-05_realType
       tup2( 8,itu4) = -0.25226e-05_realType
       tup2( 9,itu4) = -0.39077e-05_realType
       tup2(10,itu4) = -0.50583e-05_realType
       tup2(11,itu4) = -0.63087e-05_realType
       tup2(12,itu4) = -0.41839e-05_realType
       tup2(13,itu4) = -0.25527e-05_realType
       tup2(14,itu4) = -0.15942e-05_realType
       tup2(15,itu4) = -0.10087e-05_realType
       tup2(16,itu4) = -0.64089e-06_realType
       tup2(17,itu4) = -0.40493e-06_realType
       tup2(18,itu4) = -0.25305e-06_realType
       tup2(19,itu4) = -0.15515e-06_realType
       tup2(20,itu4) = -0.92540e-07_realType
       tup2(21,itu4) = -0.53318e-07_realType
       tup2(22,itu4) = -0.29160e-07_realType
       tup2(23,itu4) = -0.14786e-07_realType
       tup2(24,itu4) = -0.65163e-08_realType
       tup2(25,itu4) = -0.20578e-08_realType
       tup2(26,itu4) =  0.19849e-09_realType
       tup2(27,itu4) =  0.11643e-08_realType
       tup2(28,itu4) =  0.14562e-08_realType
       tup2(29,itu4) =  0.14088e-08_realType
       tup2(30,itu4) =  0.12218e-08_realType
       tup2(31,itu4) =  0.99318e-09_realType
       tup2(32,itu4) =  0.77547e-09_realType
       tup2(33,itu4) =  0.58737e-09_realType
       tup2(34,itu4) =  0.43430e-09_realType

       tup3( 1,itu4) = -0.18797e-04_realType
       tup3( 2,itu4) = -0.12078e-04_realType
       tup3( 3,itu4) = -0.78817e-05_realType
       tup3( 4,itu4) = -0.53414e-05_realType
       tup3( 5,itu4) = -0.40995e-05_realType
       tup3( 6,itu4) = -0.36062e-05_realType
       tup3( 7,itu4) = -0.25946e-05_realType
       tup3( 8,itu4) = -0.26029e-05_realType
       tup3( 9,itu4) = -0.16491e-05_realType
       tup3(10,itu4) = -0.78717e-06_realType
       tup3(11,itu4) =  0.86868e-07_realType
       tup3(12,itu4) =  0.87323e-07_realType
       tup3(13,itu4) =  0.40496e-07_realType
       tup3(14,itu4) =  0.21155e-07_realType
       tup3(15,itu4) =  0.12152e-07_realType
       tup3(16,itu4) =  0.74005e-08_realType
       tup3(17,itu4) =  0.45651e-08_realType
       tup3(18,itu4) =  0.28215e-08_realType
       tup3(19,itu4) =  0.17171e-08_realType
       tup3(20,itu4) =  0.10180e-08_realType
       tup3(21,itu4) =  0.59166e-09_realType
       tup3(22,itu4) =  0.33318e-09_realType
       tup3(23,itu4) =  0.18290e-09_realType
       tup3(24,itu4) =  0.96740e-10_realType
       tup3(25,itu4) =  0.49927e-10_realType
       tup3(26,itu4) =  0.24485e-10_realType
       tup3(27,itu4) =  0.11513e-10_realType
       tup3(28,itu4) =  0.49778e-11_realType
       tup3(29,itu4) =  0.18786e-11_realType
       tup3(30,itu4) =  0.46227e-12_realType
       tup3(31,itu4) = -0.11281e-12_realType
       tup3(32,itu4) = -0.31872e-12_realType
       tup3(33,itu4) = -0.36442e-12_realType
       tup3(34,itu4) = -0.35534e-12_realType

       ! Constants for nut.

       tup0( 1,itu5) =  0.11403e-05_realType
       tup0( 2,itu5) =  0.22598e-04_realType
       tup0( 3,itu5) =  0.16840e-03_realType
       tup0( 4,itu5) =  0.84569e-03_realType
       tup0( 5,itu5) =  0.33954e-02_realType
       tup0( 6,itu5) =  0.11611e-01_realType
       tup0( 7,itu5) =  0.34239e-01_realType
       tup0( 8,itu5) =  0.12517e+00_realType
       tup0( 9,itu5) =  0.35370e+00_realType
       tup0(10,itu5) =  0.78688e+00_realType
       tup0(11,itu5) =  0.14743e+01_realType
       tup0(12,itu5) =  0.24581e+01_realType
       tup0(13,itu5) =  0.37861e+01_realType
       tup0(14,itu5) =  0.55146e+01_realType
       tup0(15,itu5) =  0.77089e+01_realType
       tup0(16,itu5) =  0.10444e+02_realType
       tup0(17,itu5) =  0.13801e+02_realType
       tup0(18,itu5) =  0.17873e+02_realType
       tup0(19,itu5) =  0.22754e+02_realType
       tup0(20,itu5) =  0.28545e+02_realType
       tup0(21,itu5) =  0.35342e+02_realType
       tup0(22,itu5) =  0.43237e+02_realType
       tup0(23,itu5) =  0.52301e+02_realType
       tup0(24,itu5) =  0.62582e+02_realType
       tup0(25,itu5) =  0.74081e+02_realType
       tup0(26,itu5) =  0.86742e+02_realType
       tup0(27,itu5) =  0.10042e+03_realType
       tup0(28,itu5) =  0.11489e+03_realType
       tup0(29,itu5) =  0.12977e+03_realType
       tup0(30,itu5) =  0.14455e+03_realType
       tup0(31,itu5) =  0.15854e+03_realType
       tup0(32,itu5) =  0.17088e+03_realType
       tup0(33,itu5) =  0.18051e+03_realType
       tup0(34,itu5) =  0.18619e+03_realType

       tup1( 1,itu5) =  0.38376e-04_realType
       tup1( 2,itu5) =  0.22169e-03_realType
       tup1( 3,itu5) =  0.85495e-03_realType
       tup1( 4,itu5) =  0.26396e-02_realType
       tup1( 5,itu5) =  0.69671e-02_realType
       tup1( 6,itu5) =  0.15867e-01_realType
       tup1( 7,itu5) =  0.46643e-01_realType
       tup1( 8,itu5) =  0.10523e+00_realType
       tup1( 9,itu5) =  0.17554e+00_realType
       tup1(10,itu5) =  0.24044e+00_realType
       tup1(11,itu5) =  0.29128e+00_realType
       tup1(12,itu5) =  0.32865e+00_realType
       tup1(13,itu5) =  0.35586e+00_realType
       tup1(14,itu5) =  0.37557e+00_realType
       tup1(15,itu5) =  0.38958e+00_realType
       tup1(16,itu5) =  0.39905e+00_realType
       tup1(17,itu5) =  0.40478e+00_realType
       tup1(18,itu5) =  0.40725e+00_realType
       tup1(19,itu5) =  0.40676e+00_realType
       tup1(20,itu5) =  0.40343e+00_realType
       tup1(21,itu5) =  0.39729e+00_realType
       tup1(22,itu5) =  0.38827e+00_realType
       tup1(23,itu5) =  0.37621e+00_realType
       tup1(24,itu5) =  0.36094e+00_realType
       tup1(25,itu5) =  0.34229e+00_realType
       tup1(26,itu5) =  0.32008e+00_realType
       tup1(27,itu5) =  0.29419e+00_realType
       tup1(28,itu5) =  0.26461e+00_realType
       tup1(29,itu5) =  0.23139e+00_realType
       tup1(30,itu5) =  0.19476e+00_realType
       tup1(31,itu5) =  0.15509e+00_realType
       tup1(32,itu5) =  0.11289e+00_realType
       tup1(33,itu5) =  0.68853e-01_realType
       tup1(34,itu5) =  0.23785e-01_realType

       tup2( 1,itu5) = -0.31406e-03_realType
       tup2( 2,itu5) = -0.62688e-03_realType
       tup2( 3,itu5) = -0.10783e-02_realType
       tup2( 4,itu5) = -0.15394e-02_realType
       tup2( 5,itu5) = -0.13969e-02_realType
       tup2( 6,itu5) = -0.14454e-01_realType
       tup2( 7,itu5) =  0.23424e-02_realType
       tup2( 8,itu5) =  0.12646e-01_realType
       tup2( 9,itu5) =  0.15053e-01_realType
       tup2(10,itu5) =  0.11222e-01_realType
       tup2(11,itu5) =  0.69106e-02_realType
       tup2(12,itu5) =  0.41344e-02_realType
       tup2(13,itu5) =  0.25022e-02_realType
       tup2(14,itu5) =  0.15280e-02_realType
       tup2(15,itu5) =  0.93303e-03_realType
       tup2(16,itu5) =  0.55503e-03_realType
       tup2(17,itu5) =  0.32030e-03_realType
       tup2(18,itu5) =  0.16369e-03_realType
       tup2(19,itu5) =  0.59340e-04_realType
       tup2(20,itu5) = -0.10107e-04_realType
       tup2(21,itu5) = -0.60320e-04_realType
       tup2(22,itu5) = -0.93922e-04_realType
       tup2(23,itu5) = -0.12076e-03_realType
       tup2(24,itu5) = -0.13968e-03_realType
       tup2(25,itu5) = -0.15343e-03_realType
       tup2(26,itu5) = -0.16392e-03_realType
       tup2(27,itu5) = -0.17035e-03_realType
       tup2(28,itu5) = -0.17361e-03_realType
       tup2(29,itu5) = -0.17458e-03_realType
       tup2(30,itu5) = -0.17152e-03_realType
       tup2(31,itu5) = -0.16658e-03_realType
       tup2(32,itu5) = -0.15883e-03_realType
       tup2(33,itu5) = -0.14875e-03_realType
       tup2(34,itu5) = -0.13710e-03_realType

       tup3( 1,itu5) =  0.11903e-02_realType
       tup3( 2,itu5) =  0.21637e-02_realType
       tup3( 3,itu5) =  0.33791e-02_realType
       tup3( 4,itu5) =  0.45914e-02_realType
       tup3( 5,itu5) =  0.50750e-02_realType
       tup3( 6,itu5) =  0.17668e-01_realType
       tup3( 7,itu5) =  0.95211e-02_realType
       tup3( 8,itu5) =  0.32633e-02_realType
       tup3( 9,itu5) =  0.15990e-03_realType
       tup3(10,itu5) = -0.34910e-03_realType
       tup3(11,itu5) = -0.21152e-03_realType
       tup3(12,itu5) = -0.10676e-03_realType
       tup3(13,itu5) = -0.58327e-04_realType
       tup3(14,itu5) = -0.35522e-04_realType
       tup3(15,itu5) = -0.24007e-04_realType
       tup3(16,itu5) = -0.16925e-04_realType
       tup3(17,itu5) = -0.13086e-04_realType
       tup3(18,itu5) = -0.10268e-04_realType
       tup3(19,itu5) = -0.82175e-05_realType
       tup3(20,itu5) = -0.67361e-05_realType
       tup3(21,itu5) = -0.54809e-05_realType
       tup3(22,itu5) = -0.45452e-05_realType
       tup3(23,itu5) = -0.36958e-05_realType
       tup3(24,itu5) = -0.30062e-05_realType
       tup3(25,itu5) = -0.24308e-05_realType
       tup3(26,itu5) = -0.19327e-05_realType
       tup3(27,itu5) = -0.15225e-05_realType
       tup3(28,itu5) = -0.11825e-05_realType
       tup3(29,itu5) = -0.89374e-06_realType
       tup3(30,itu5) = -0.66991e-06_realType
       tup3(31,itu5) = -0.48509e-06_realType
       tup3(32,itu5) = -0.34165e-06_realType
       tup3(33,itu5) = -0.23235e-06_realType
       tup3(34,itu5) = -0.14885e-06_realType

       !===============================================================

    case (6_intType)

       ! Version 6 of the model.

       ! Set the number of data points and allocate the memory for
       ! the arrays of the curve fits.

       nFit = 34

       allocate(ypT(0:nFit),          reT(0:nFit),          &
            up0(nFit), up1(nFit), up2(nFit), up3(nFit), &
            tup0(nFit,nt1:nt2+1), tup1(nFit,nt1:nt2+1), &
            tup2(nFit,nt1:nt2+1), tup3(nFit,nt1:nt2+1), &
            tuLogFit(nt1:nt2+1), stat=ierr)
       if(ierr /= 0)                              &
            call terminate("initCurveFitDataVf", &
            "Memory allocation failure for curve fit &
            &coefficients")

       ! Set the values of the Reynolds numbers at interval
       ! boundaries.

       reT( 0) =     0.13341_realType
       reT( 1) =     0.47908_realType
       reT( 2) =     1.2328_realType
       reT( 3) =     2.6968_realType
       reT( 4) =     5.3545_realType
       reT( 5) =     9.9485_realType
       reT( 6) =    17.557_realType
       reT( 7) =    29.618_realType
       reT( 8) =    47.649_realType
       reT( 9) =    73.015_realType
       reT(10) =   107.22_realType
       reT(11) =   152.01_realType
       reT(12) =   209.51_realType
       reT(13) =   282.32_realType
       reT(14) =   373.64_realType
       reT(15) =   487.34_realType
       reT(16) =   628.02_realType
       reT(17) =   801.15_realType
       reT(18) =  1013.1_realType
       reT(19) =  1271.3_realType
       reT(20) =  1584.4_realType
       reT(21) =  1962.4_realType
       reT(22) =  2416.6_realType
       reT(23) =  2960.4_realType
       reT(24) =  3608.7_realType
       reT(25) =  4378.7_realType
       reT(26) =  5289.8_realType
       reT(27) =  6364.3_realType
       reT(28) =  7627.3_realType
       reT(29) =  9107.0_realType
       reT(30) = 10835.0_realType
       reT(31) = 12849.0_realType
       reT(32) = 15187.0_realType
       reT(33) = 17896.0_realType
       reT(34) = 21028.0_realType

       ! Set the values of the y+ values at interval boundaries.

       ypT( 0) =   0.36526_realType
       ypT( 1) =   0.69218_realType
       ypT( 2) =   1.1105_realType
       ypT( 3) =   1.6431_realType
       ypT( 4) =   2.3179_realType
       ypT( 5) =   3.169_realType
       ypT( 6) =   4.2373_realType
       ypT( 7) =   5.5723_realType
       ypT( 8) =   7.2331_realType
       ypT( 9) =   9.2902_realType
       ypT(10) =  11.827_realType
       ypT(11) =  14.943_realType
       ypT(12) =  18.755_realType
       ypT(13) =  23.398_realType
       ypT(14) =  29.033_realType
       ypT(15) =  35.843_realType
       ypT(16) =  44.045_realType
       ypT(17) =  53.884_realType
       ypT(18) =  65.646_realType
       ypT(19) =  79.656_realType
       ypT(20) =  96.285_realType
       ypT(21) = 115.96_realType
       ypT(22) = 139.15_realType
       ypT(23) = 166.4_realType
       ypT(24) = 198.32_realType
       ypT(25) = 235.59_realType
       ypT(26) = 278.96_realType
       ypT(27) = 329.31_realType
       ypT(28) = 387.55_realType
       ypT(29) = 454.74_realType
       ypT(30) = 532.04_realType
       ypT(31) = 620.71_realType
       ypT(32) = 722.14_realType
       ypT(33) = 837.87_realType
       ypT(34) = 969.57_realType

       ! Set the values of constants for the cubic fits of the
       ! non-dimensional tangential velocity.

       up0( 1) =  0.36525e+00_realType
       up0( 2) =  0.69214e+00_realType
       up0( 3) =  0.11102e+01_realType
       up0( 4) =  0.16413e+01_realType
       up0( 5) =  0.23101e+01_realType
       up0( 6) =  0.31393e+01_realType
       up0( 7) =  0.41434e+01_realType
       up0( 8) =  0.53152e+01_realType
       up0( 9) =  0.65876e+01_realType
       up0(10) =  0.78593e+01_realType
       up0(11) =  0.90657e+01_realType
       up0(12) =  0.10173e+02_realType
       up0(13) =  0.11171e+02_realType
       up0(14) =  0.12066e+02_realType
       up0(15) =  0.12870e+02_realType
       up0(16) =  0.13596e+02_realType
       up0(17) =  0.14259e+02_realType
       up0(18) =  0.14868e+02_realType
       up0(19) =  0.15433e+02_realType
       up0(20) =  0.15960e+02_realType
       up0(21) =  0.16455e+02_realType
       up0(22) =  0.16923e+02_realType
       up0(23) =  0.17367e+02_realType
       up0(24) =  0.17791e+02_realType
       up0(25) =  0.18197e+02_realType
       up0(26) =  0.18586e+02_realType
       up0(27) =  0.18962e+02_realType
       up0(28) =  0.19327e+02_realType
       up0(29) =  0.19681e+02_realType
       up0(30) =  0.20027e+02_realType
       up0(31) =  0.20366e+02_realType
       up0(32) =  0.20700e+02_realType
       up0(33) =  0.21030e+02_realType
       up0(34) =  0.21359e+02_realType

       up1( 1) =  0.12450e+01_realType
       up1( 2) =  0.67756e+00_realType
       up1( 3) =  0.42800e+00_realType
       up1( 4) =  0.29111e+00_realType
       up1( 5) =  0.20657e+00_realType
       up1( 6) =  0.15024e+00_realType
       up1( 7) =  0.11062e+00_realType
       up1( 8) =  0.81225e-01_realType
       up1( 9) =  0.58625e-01_realType
       up1(10) =  0.41597e-01_realType
       up1(11) =  0.29283e-01_realType
       up1(12) =  0.20582e-01_realType
       up1(13) =  0.14529e-01_realType
       up1(14) =  0.10350e-01_realType
       up1(15) =  0.74648e-02_realType
       up1(16) =  0.54606e-02_realType
       up1(17) =  0.40524e-02_realType
       up1(18) =  0.30487e-02_realType
       up1(19) =  0.23228e-02_realType
       up1(20) =  0.17899e-02_realType
       up1(21) =  0.13940e-02_realType
       up1(22) =  0.10960e-02_realType
       up1(23) =  0.86932e-03_realType
       up1(24) =  0.69545e-03_realType
       up1(25) =  0.56075e-03_realType
       up1(26) =  0.45552e-03_realType
       up1(27) =  0.37282e-03_realType
       up1(28) =  0.30742e-03_realType
       up1(29) =  0.25530e-03_realType
       up1(30) =  0.21352e-03_realType
       up1(31) =  0.17992e-03_realType
       up1(32) =  0.15273e-03_realType
       up1(33) =  0.13061e-03_realType
       up1(34) =  0.11257e-03_realType

       up2( 1) = -0.95638e+00_realType
       up2( 2) = -0.15826e+00_realType
       up2( 3) = -0.40071e-01_realType
       up2( 4) = -0.12775e-01_realType
       up2( 5) = -0.47553e-02_realType
       up2( 6) = -0.19994e-02_realType
       up2( 7) = -0.91167e-03_realType
       up2( 8) = -0.51977e-03_realType
       up2( 9) = -0.33276e-03_realType
       up2(10) = -0.19526e-03_realType
       up2(11) = -0.11176e-03_realType
       up2(12) = -0.62687e-04_realType
       up2(13) = -0.34811e-04_realType
       up2(14) = -0.19255e-04_realType
       up2(15) = -0.10712e-04_realType
       up2(16) = -0.60184e-05_realType
       up2(17) = -0.34485e-05_realType
       up2(18) = -0.20101e-05_realType
       up2(19) = -0.11974e-05_realType
       up2(20) = -0.72812e-06_realType
       up2(21) = -0.44763e-06_realType
       up2(22) = -0.28271e-06_realType
       up2(23) = -0.17940e-06_realType
       up2(24) = -0.11595e-06_realType
       up2(25) = -0.75702e-07_realType
       up2(26) = -0.50362e-07_realType
       up2(27) = -0.33465e-07_realType
       up2(28) = -0.22609e-07_realType
       up2(29) = -0.15431e-07_realType
       up2(30) = -0.10636e-07_realType
       up2(31) = -0.73064e-08_realType
       up2(32) = -0.51550e-08_realType
       up2(33) = -0.35885e-08_realType
       up2(34) = -0.25596e-08_realType

       up3( 1) =  0.26155e+00_realType
       up3( 2) = -0.64433e-02_realType
       up3( 3) = -0.30415e-02_realType
       up3( 4) = -0.78504e-03_realType
       up3( 5) = -0.19960e-03_realType
       up3( 6) = -0.52972e-04_realType
       up3( 7) = -0.16968e-04_realType
       up3( 8) = -0.39545e-05_realType
       up3( 9) = -0.75520e-07_realType
       up3(10) =  0.29768e-06_realType
       up3(11) =  0.21782e-06_realType
       up3(12) =  0.11651e-06_realType
       up3(13) =  0.55958e-07_realType
       up3(14) =  0.25252e-07_realType
       up3(15) =  0.11133e-07_realType
       up3(16) =  0.48030e-08_realType
       up3(17) =  0.21171e-08_realType
       up3(18) =  0.93600e-09_realType
       up3(19) =  0.42750e-09_realType
       up3(20) =  0.20397e-09_realType
       up3(21) =  0.94244e-10_realType
       up3(22) =  0.48769e-10_realType
       up3(23) =  0.23922e-10_realType
       up3(24) =  0.12409e-10_realType
       up3(25) =  0.63795e-11_realType
       up3(26) =  0.36413e-11_realType
       up3(27) =  0.18837e-11_realType
       up3(28) =  0.10413e-11_realType
       up3(29) =  0.59232e-12_realType
       up3(30) =  0.35361e-12_realType
       up3(31) =  0.18262e-12_realType
       up3(32) =  0.12163e-12_realType
       up3(33) =  0.63587e-13_realType
       up3(34) =  0.41966e-13_realType

       ! Set the values of tuLogFit. All variables have been
       ! fitted linearly; the fifth variable is the eddy viscosity.

       tuLogFit(itu1) = .false.
       tuLogFit(itu2) = .false.
       tuLogFit(itu3) = .false.
       tuLogFit(itu4) = .false.
       tuLogFit(itu5) = .false.

       ! Set the values of constants for the cubic fits of the
       ! non-dimensional k, eps, v2 and f values.

       ! Constants for k.

       tup0( 1,itu1) =  0.20232e-01_realType
       tup0( 2,itu1) =  0.66720e-01_realType
       tup0( 3,itu1) =  0.15994e+00_realType
       tup0( 4,itu1) =  0.32733e+00_realType
       tup0( 5,itu1) =  0.60629e+00_realType
       tup0( 6,itu1) =  0.10415e+01_realType
       tup0( 7,itu1) =  0.16715e+01_realType
       tup0( 8,itu1) =  0.24989e+01_realType
       tup0( 9,itu1) =  0.34343e+01_realType
       tup0(10,itu1) =  0.42683e+01_realType
       tup0(11,itu1) =  0.48591e+01_realType
       tup0(12,itu1) =  0.51750e+01_realType
       tup0(13,itu1) =  0.52611e+01_realType
       tup0(14,itu1) =  0.51878e+01_realType
       tup0(15,itu1) =  0.50187e+01_realType
       tup0(16,itu1) =  0.48000e+01_realType
       tup0(17,itu1) =  0.45620e+01_realType
       tup0(18,itu1) =  0.43226e+01_realType
       tup0(19,itu1) =  0.40922e+01_realType
       tup0(20,itu1) =  0.38760e+01_realType
       tup0(21,itu1) =  0.36764e+01_realType
       tup0(22,itu1) =  0.34942e+01_realType
       tup0(23,itu1) =  0.33290e+01_realType
       tup0(24,itu1) =  0.31798e+01_realType
       tup0(25,itu1) =  0.30456e+01_realType
       tup0(26,itu1) =  0.29248e+01_realType
       tup0(27,itu1) =  0.28163e+01_realType
       tup0(28,itu1) =  0.27186e+01_realType
       tup0(29,itu1) =  0.26304e+01_realType
       tup0(30,itu1) =  0.25505e+01_realType
       tup0(31,itu1) =  0.24775e+01_realType
       tup0(32,itu1) =  0.24103e+01_realType
       tup0(33,itu1) =  0.23475e+01_realType
       tup0(34,itu1) =  0.22880e+01_realType

       tup1( 1,itu1) =  0.11194e+00_realType
       tup1( 2,itu1) =  0.18748e+00_realType
       tup1( 3,itu1) =  0.27407e+00_realType
       tup1( 4,itu1) =  0.36966e+00_realType
       tup1( 5,itu1) =  0.46803e+00_realType
       tup1( 6,itu1) =  0.55499e+00_realType
       tup1( 7,itu1) =  0.60642e+00_realType
       tup1( 8,itu1) =  0.58843e+00_realType
       tup1( 9,itu1) =  0.47590e+00_realType
       tup1(10,itu1) =  0.31013e+00_realType
       tup1(11,itu1) =  0.16039e+00_realType
       tup1(12,itu1) =  0.58026e-01_realType
       tup1(13,itu1) =  0.15163e-02_realType
       tup1(14,itu1) = -0.23590e-01_realType
       tup1(15,itu1) = -0.31162e-01_realType
       tup1(16,itu1) = -0.30422e-01_realType
       tup1(17,itu1) = -0.26463e-01_realType
       tup1(18,itu1) = -0.21749e-01_realType
       tup1(19,itu1) = -0.17331e-01_realType
       tup1(20,itu1) = -0.13569e-01_realType
       tup1(21,itu1) = -0.10516e-01_realType
       tup1(22,itu1) = -0.81061e-02_realType
       tup1(23,itu1) = -0.62327e-02_realType
       tup1(24,itu1) = -0.47901e-02_realType
       tup1(25,itu1) = -0.36855e-02_realType
       tup1(26,itu1) = -0.28427e-02_realType
       tup1(27,itu1) = -0.22005e-02_realType
       tup1(28,itu1) = -0.17118e-02_realType
       tup1(29,itu1) = -0.13404e-02_realType
       tup1(30,itu1) = -0.10584e-02_realType
       tup1(31,itu1) = -0.84476e-03_realType
       tup1(32,itu1) = -0.68362e-03_realType
       tup1(33,itu1) = -0.56284e-03_realType
       tup1(34,itu1) = -0.47348e-03_realType

       tup2( 1,itu1) =  0.46642e-01_realType
       tup2( 2,itu1) =  0.46777e-01_realType
       tup2( 3,itu1) =  0.47060e-01_realType
       tup2( 4,itu1) =  0.48495e-01_realType
       tup2( 5,itu1) =  0.50649e-01_realType
       tup2( 6,itu1) =  0.49396e-01_realType
       tup2( 7,itu1) =  0.43509e-01_realType
       tup2( 8,itu1) =  0.22224e-01_realType
       tup2( 9,itu1) = -0.22230e-01_realType
       tup2(10,itu1) = -0.32317e-01_realType
       tup2(11,itu1) = -0.23975e-01_realType
       tup2(12,itu1) = -0.13064e-01_realType
       tup2(13,itu1) = -0.57706e-02_realType
       tup2(14,itu1) = -0.20813e-02_realType
       tup2(15,itu1) = -0.52362e-03_realType
       tup2(16,itu1) =  0.28081e-04_realType
       tup2(17,itu1) =  0.17226e-03_realType
       tup2(18,itu1) =  0.17433e-03_realType
       tup2(19,itu1) =  0.13808e-03_realType
       tup2(20,itu1) =  0.99662e-04_realType
       tup2(21,itu1) =  0.68603e-04_realType
       tup2(22,itu1) =  0.46166e-04_realType
       tup2(23,itu1) =  0.30634e-04_realType
       tup2(24,itu1) =  0.20237e-04_realType
       tup2(25,itu1) =  0.13310e-04_realType
       tup2(26,itu1) =  0.87532e-05_realType
       tup2(27,itu1) =  0.57793e-05_realType
       tup2(28,itu1) =  0.38107e-05_realType
       tup2(29,itu1) =  0.25226e-05_realType
       tup2(30,itu1) =  0.16733e-05_realType
       tup2(31,itu1) =  0.11109e-05_realType
       tup2(32,itu1) =  0.73765e-06_realType
       tup2(33,itu1) =  0.49068e-06_realType
       tup2(34,itu1) =  0.32170e-06_realType

       tup3( 1,itu1) =  0.14048e+00_realType
       tup3( 2,itu1) =  0.90412e-01_realType
       tup3( 3,itu1) =  0.53426e-01_realType
       tup3( 4,itu1) =  0.24095e-01_realType
       tup3( 5,itu1) =  0.34292e-03_realType
       tup3( 6,itu1) = -0.15803e-01_realType
       tup3( 7,itu1) = -0.25093e-01_realType
       tup3( 8,itu1) = -0.22520e-01_realType
       tup3( 9,itu1) = -0.58532e-02_realType
       tup3(10,itu1) =  0.73777e-03_realType
       tup3(11,itu1) =  0.16152e-02_realType
       tup3(12,itu1) =  0.98832e-03_realType
       tup3(13,itu1) =  0.44037e-03_realType
       tup3(14,itu1) =  0.16675e-03_realType
       tup3(15,itu1) =  0.56569e-04_realType
       tup3(16,itu1) =  0.17339e-04_realType
       tup3(17,itu1) =  0.45579e-05_realType
       tup3(18,itu1) =  0.76381e-06_realType
       tup3(19,itu1) = -0.18121e-06_realType
       tup3(20,itu1) = -0.31622e-06_realType
       tup3(21,itu1) = -0.24874e-06_realType
       tup3(22,itu1) = -0.16604e-06_realType
       tup3(23,itu1) = -0.10191e-06_realType
       tup3(24,itu1) = -0.61291e-07_realType
       tup3(25,itu1) = -0.35829e-07_realType
       tup3(26,itu1) = -0.20759e-07_realType
       tup3(27,itu1) = -0.12251e-07_realType
       tup3(28,itu1) = -0.71234e-08_realType
       tup3(29,itu1) = -0.42082e-08_realType
       tup3(30,itu1) = -0.25154e-08_realType
       tup3(31,itu1) = -0.15209e-08_realType
       tup3(32,itu1) = -0.93515e-09_realType
       tup3(33,itu1) = -0.60259e-09_realType
       tup3(34,itu1) = -0.38327e-09_realType

       ! Constants for epsilon.

       tup0( 1,itu2) =  0.24164e+00_realType
       tup0( 2,itu2) =  0.22008e+00_realType
       tup0( 3,itu2) =  0.19751e+00_realType
       tup0( 4,itu2) =  0.17571e+00_realType
       tup0( 5,itu2) =  0.15723e+00_realType
       tup0( 6,itu2) =  0.14536e+00_realType
       tup0( 7,itu2) =  0.14328e+00_realType
       tup0( 8,itu2) =  0.15220e+00_realType
       tup0( 9,itu2) =  0.16784e+00_realType
       tup0(10,itu2) =  0.17712e+00_realType
       tup0(11,itu2) =  0.17387e+00_realType
       tup0(12,itu2) =  0.15954e+00_realType
       tup0(13,itu2) =  0.13909e+00_realType
       tup0(14,itu2) =  0.11713e+00_realType
       tup0(15,itu2) =  0.96528e-01_realType
       tup0(16,itu2) =  0.78571e-01_realType
       tup0(17,itu2) =  0.63547e-01_realType
       tup0(18,itu2) =  0.51265e-01_realType
       tup0(19,itu2) =  0.41354e-01_realType
       tup0(20,itu2) =  0.33411e-01_realType
       tup0(21,itu2) =  0.27064e-01_realType
       tup0(22,itu2) =  0.21997e-01_realType
       tup0(23,itu2) =  0.17947e-01_realType
       tup0(24,itu2) =  0.14703e-01_realType
       tup0(25,itu2) =  0.12100e-01_realType
       tup0(26,itu2) =  0.10003e-01_realType
       tup0(27,itu2) =  0.83083e-02_realType
       tup0(28,itu2) =  0.69337e-02_realType
       tup0(29,itu2) =  0.58144e-02_realType
       tup0(30,itu2) =  0.48991e-02_realType
       tup0(31,itu2) =  0.41476e-02_realType
       tup0(32,itu2) =  0.35278e-02_realType
       tup0(33,itu2) =  0.30142e-02_realType
       tup0(34,itu2) =  0.25866e-02_realType

       tup1( 1,itu2) = -0.70667e-01_realType
       tup1( 2,itu2) = -0.59210e-01_realType
       tup1( 3,itu2) = -0.46663e-01_realType
       tup1( 4,itu2) = -0.33366e-01_realType
       tup1( 5,itu2) = -0.19886e-01_realType
       tup1( 6,itu2) = -0.72638e-02_realType
       tup1( 7,itu2) =  0.28452e-02_realType
       tup1( 8,itu2) =  0.81949e-02_realType
       tup1( 9,itu2) =  0.67022e-02_realType
       tup1(10,itu2) =  0.13143e-02_realType
       tup1(11,itu2) = -0.31089e-02_realType
       tup1(12,itu2) = -0.50215e-02_realType
       tup1(13,itu2) = -0.50170e-02_realType
       tup1(14,itu2) = -0.41408e-02_realType
       tup1(15,itu2) = -0.30979e-02_realType
       tup1(16,itu2) = -0.21970e-02_realType
       tup1(17,itu2) = -0.15136e-02_realType
       tup1(18,itu2) = -0.10274e-02_realType
       tup1(19,itu2) = -0.69277e-03_realType
       tup1(20,itu2) = -0.46639e-03_realType
       tup1(21,itu2) = -0.31444e-03_realType
       tup1(22,itu2) = -0.21272e-03_realType
       tup1(23,itu2) = -0.14458e-03_realType
       tup1(24,itu2) = -0.98812e-04_realType
       tup1(25,itu2) = -0.67941e-04_realType
       tup1(26,itu2) = -0.47013e-04_realType
       tup1(27,itu2) = -0.32748e-04_realType
       tup1(28,itu2) = -0.22968e-04_realType
       tup1(29,itu2) = -0.16220e-04_realType
       tup1(30,itu2) = -0.11536e-04_realType
       tup1(31,itu2) = -0.82628e-05_realType
       tup1(32,itu2) = -0.59619e-05_realType
       tup1(33,itu2) = -0.43340e-05_realType
       tup1(34,itu2) = -0.31753e-05_realType

       tup2( 1,itu2) =  0.82009e-02_realType
       tup2( 2,itu2) =  0.78044e-02_realType
       tup2( 3,itu2) =  0.72279e-02_realType
       tup2( 4,itu2) =  0.66236e-02_realType
       tup2( 5,itu2) =  0.61231e-02_realType
       tup2( 6,itu2) =  0.54757e-02_realType
       tup2( 7,itu2) =  0.46058e-02_realType
       tup2( 8,itu2) =  0.31013e-02_realType
       tup2( 9,itu2) = -0.57395e-03_realType
       tup2(10,itu2) = -0.13230e-02_realType
       tup2(11,itu2) = -0.82054e-03_realType
       tup2(12,itu2) = -0.27325e-03_realType
       tup2(13,itu2) = -0.29977e-05_realType
       tup2(14,itu2) =  0.73241e-04_realType
       tup2(15,itu2) =  0.70951e-04_realType
       tup2(16,itu2) =  0.50202e-04_realType
       tup2(17,itu2) =  0.31490e-04_realType
       tup2(18,itu2) =  0.18675e-04_realType
       tup2(19,itu2) =  0.10781e-04_realType
       tup2(20,itu2) =  0.61493e-05_realType
       tup2(21,itu2) =  0.34945e-05_realType
       tup2(22,itu2) =  0.19882e-05_realType
       tup2(23,itu2) =  0.11363e-05_realType
       tup2(24,itu2) =  0.65315e-06_realType
       tup2(25,itu2) =  0.37828e-06_realType
       tup2(26,itu2) =  0.22089e-06_realType
       tup2(27,itu2) =  0.13005e-06_realType
       tup2(28,itu2) =  0.77269e-07_realType
       tup2(29,itu2) =  0.46321e-07_realType
       tup2(30,itu2) =  0.28022e-07_realType
       tup2(31,itu2) =  0.17109e-07_realType
       tup2(32,itu2) =  0.10536e-07_realType
       tup2(33,itu2) =  0.65500e-08_realType
       tup2(34,itu2) =  0.41065e-08_realType

       tup3( 1,itu2) =  0.19008e-01_realType
       tup3( 2,itu2) =  0.11466e-01_realType
       tup3( 3,itu2) =  0.65783e-02_realType
       tup3( 4,itu2) =  0.33228e-02_realType
       tup3( 5,itu2) =  0.10124e-02_realType
       tup3( 6,itu2) = -0.46456e-03_realType
       tup3( 7,itu2) = -0.12995e-02_realType
       tup3( 8,itu2) = -0.14253e-02_realType
       tup3( 9,itu2) = -0.23840e-03_realType
       tup3(10,itu2) =  0.11860e-03_realType
       tup3(11,itu2) =  0.10989e-03_realType
       tup3(12,itu2) =  0.47896e-04_realType
       tup3(13,itu2) =  0.13977e-04_realType
       tup3(14,itu2) =  0.22846e-05_realType
       tup3(15,itu2) = -0.47051e-06_realType
       tup3(16,itu2) = -0.69412e-06_realType
       tup3(17,itu2) = -0.45961e-06_realType
       tup3(18,itu2) = -0.25229e-06_realType
       tup3(19,itu2) = -0.12854e-06_realType
       tup3(20,itu2) = -0.63372e-07_realType
       tup3(21,itu2) = -0.30806e-07_realType
       tup3(22,itu2) = -0.14923e-07_realType
       tup3(23,itu2) = -0.72529e-08_realType
       tup3(24,itu2) = -0.35415e-08_realType
       tup3(25,itu2) = -0.17443e-08_realType
       tup3(26,itu2) = -0.86787e-09_realType
       tup3(27,itu2) = -0.43581e-09_realType
       tup3(28,itu2) = -0.22141e-09_realType
       tup3(29,itu2) = -0.11372e-09_realType
       tup3(30,itu2) = -0.59087e-10_realType
       tup3(31,itu2) = -0.31082e-10_realType
       tup3(32,itu2) = -0.16513e-10_realType
       tup3(33,itu2) = -0.88934e-11_realType
       tup3(34,itu2) = -0.48436e-11_realType

       ! Constants for v2.

       tup0( 1,itu3) =  0.86118e-05_realType
       tup0( 2,itu3) =  0.76561e-04_realType
       tup0( 3,itu3) =  0.38909e-03_realType
       tup0( 4,itu3) =  0.14605e-02_realType
       tup0( 5,itu3) =  0.44732e-02_realType
       tup0( 6,itu3) =  0.11704e-01_realType
       tup0( 7,itu3) =  0.26786e-01_realType
       tup0( 8,itu3) =  0.54333e-01_realType
       tup0( 9,itu3) =  0.97821e-01_realType
       tup0(10,itu3) =  0.15774e+00_realType
       tup0(11,itu3) =  0.23273e+00_realType
       tup0(12,itu3) =  0.32009e+00_realType
       tup0(13,itu3) =  0.41667e+00_realType
       tup0(14,itu3) =  0.51933e+00_realType
       tup0(15,itu3) =  0.62515e+00_realType
       tup0(16,itu3) =  0.73157e+00_realType
       tup0(17,itu3) =  0.83660e+00_realType
       tup0(18,itu3) =  0.93915e+00_realType
       tup0(19,itu3) =  0.10385e+01_realType
       tup0(20,itu3) =  0.11342e+01_realType
       tup0(21,itu3) =  0.12256e+01_realType
       tup0(22,itu3) =  0.13125e+01_realType
       tup0(23,itu3) =  0.13944e+01_realType
       tup0(24,itu3) =  0.14708e+01_realType
       tup0(25,itu3) =  0.15412e+01_realType
       tup0(26,itu3) =  0.16053e+01_realType
       tup0(27,itu3) =  0.16624e+01_realType
       tup0(28,itu3) =  0.17120e+01_realType
       tup0(29,itu3) =  0.17534e+01_realType
       tup0(30,itu3) =  0.17862e+01_realType
       tup0(31,itu3) =  0.18095e+01_realType
       tup0(32,itu3) =  0.18227e+01_realType
       tup0(33,itu3) =  0.18250e+01_realType
       tup0(34,itu3) =  0.18158e+01_realType

       tup1( 1,itu3) =  0.13148e-03_realType
       tup1( 2,itu3) =  0.51056e-03_realType
       tup1( 3,itu3) =  0.14555e-02_realType
       tup1( 4,itu3) =  0.33824e-02_realType
       tup1( 5,itu3) =  0.67133e-02_realType
       tup1( 6,itu3) =  0.11625e-01_realType
       tup1( 7,itu3) =  0.17737e-01_realType
       tup1( 8,itu3) =  0.23712e-01_realType
       tup1( 9,itu3) =  0.27813e-01_realType
       tup1(10,itu3) =  0.29364e-01_realType
       tup1(11,itu3) =  0.28719e-01_realType
       tup1(12,itu3) =  0.26553e-01_realType
       tup1(13,itu3) =  0.23565e-01_realType
       tup1(14,itu3) =  0.20284e-01_realType
       tup1(15,itu3) =  0.17054e-01_realType
       tup1(16,itu3) =  0.14085e-01_realType
       tup1(17,itu3) =  0.11506e-01_realType
       tup1(18,itu3) =  0.93475e-02_realType
       tup1(19,itu3) =  0.75674e-02_realType
       tup1(20,itu3) =  0.61075e-02_realType
       tup1(21,itu3) =  0.49132e-02_realType
       tup1(22,itu3) =  0.39366e-02_realType
       tup1(23,itu3) =  0.31369e-02_realType
       tup1(24,itu3) =  0.24817e-02_realType
       tup1(25,itu3) =  0.19440e-02_realType
       tup1(26,itu3) =  0.15022e-02_realType
       tup1(27,itu3) =  0.11385e-02_realType
       tup1(28,itu3) =  0.83880e-03_realType
       tup1(29,itu3) =  0.59169e-03_realType
       tup1(30,itu3) =  0.38785e-03_realType
       tup1(31,itu3) =  0.21987e-03_realType
       tup1(32,itu3) =  0.81692e-04_realType
       tup1(33,itu3) = -0.31589e-04_realType
       tup1(34,itu3) = -0.12412e-03_realType

       tup2( 1,itu3) = -0.45876e-03_realType
       tup2( 2,itu3) = -0.56213e-03_realType
       tup2( 3,itu3) = -0.48447e-03_realType
       tup2( 4,itu3) = -0.12712e-03_realType
       tup2( 5,itu3) =  0.51606e-03_realType
       tup2( 6,itu3) =  0.12754e-02_realType
       tup2( 7,itu3) =  0.20365e-02_realType
       tup2( 8,itu3) =  0.19980e-02_realType
       tup2( 9,itu3) =  0.11629e-02_realType
       tup2(10,itu3) =  0.48125e-03_realType
       tup2(11,itu3) =  0.38876e-04_realType
       tup2(12,itu3) = -0.17127e-03_realType
       tup2(13,itu3) = -0.23448e-03_realType
       tup2(14,itu3) = -0.22707e-03_realType
       tup2(15,itu3) = -0.19347e-03_realType
       tup2(16,itu3) = -0.15335e-03_realType
       tup2(17,itu3) = -0.11100e-03_realType
       tup2(18,itu3) = -0.78059e-04_realType
       tup2(19,itu3) = -0.54121e-04_realType
       tup2(20,itu3) = -0.37687e-04_realType
       tup2(21,itu3) = -0.26073e-04_realType
       tup2(22,itu3) = -0.18195e-04_realType
       tup2(23,itu3) = -0.12730e-04_realType
       tup2(24,itu3) = -0.89271e-05_realType
       tup2(25,itu3) = -0.63091e-05_realType
       tup2(26,itu3) = -0.44582e-05_realType
       tup2(27,itu3) = -0.31881e-05_realType
       tup2(28,itu3) = -0.22705e-05_realType
       tup2(29,itu3) = -0.16359e-05_realType
       tup2(30,itu3) = -0.11754e-05_realType
       tup2(31,itu3) = -0.85145e-06_realType
       tup2(32,itu3) = -0.61411e-06_realType
       tup2(33,itu3) = -0.44454e-06_realType
       tup2(34,itu3) = -0.32118e-06_realType

       tup3( 1,itu3) =  0.21178e-02_realType
       tup3( 2,itu3) =  0.26960e-02_realType
       tup3( 3,itu3) =  0.28709e-02_realType
       tup3( 4,itu3) =  0.25635e-02_realType
       tup3( 5,itu3) =  0.18560e-02_realType
       tup3( 6,itu3) =  0.98939e-03_realType
       tup3( 7,itu3) =  0.10046e-03_realType
       tup3( 8,itu3) = -0.30638e-03_realType
       tup3( 9,itu3) = -0.25468e-03_realType
       tup3(10,itu3) = -0.15987e-03_realType
       tup3(11,itu3) = -0.82673e-04_realType
       tup3(12,itu3) = -0.38599e-04_realType
       tup3(13,itu3) = -0.17054e-04_realType
       tup3(14,itu3) = -0.70519e-05_realType
       tup3(15,itu3) = -0.23971e-05_realType
       tup3(16,itu3) = -0.31545e-06_realType
       tup3(17,itu3) =  0.88126e-07_realType
       tup3(18,itu3) =  0.13543e-06_realType
       tup3(19,itu3) =  0.96125e-07_realType
       tup3(20,itu3) =  0.71267e-07_realType
       tup3(21,itu3) =  0.42315e-07_realType
       tup3(22,itu3) =  0.27372e-07_realType
       tup3(23,itu3) =  0.17327e-07_realType
       tup3(24,itu3) =  0.10546e-07_realType
       tup3(25,itu3) =  0.68294e-08_realType
       tup3(26,itu3) =  0.40887e-08_realType
       tup3(27,itu3) =  0.27989e-08_realType
       tup3(28,itu3) =  0.17076e-08_realType
       tup3(29,itu3) =  0.11817e-08_realType
       tup3(30,itu3) =  0.76583e-09_realType
       tup3(31,itu3) =  0.54352e-09_realType
       tup3(32,itu3) =  0.36625e-09_realType
       tup3(33,itu3) =  0.25793e-09_realType
       tup3(34,itu3) =  0.18203e-09_realType

       ! Constants for f.

       tup0( 1,itu4) =  0.36932e-03_realType
       tup0( 2,itu4) =  0.69498e-03_realType
       tup0( 3,itu4) =  0.11057e-02_realType
       tup0( 4,itu4) =  0.16193e-02_realType
       tup0( 5,itu4) =  0.22555e-02_realType
       tup0( 6,itu4) =  0.30335e-02_realType
       tup0( 7,itu4) =  0.39685e-02_realType
       tup0( 8,itu4) =  0.50643e-02_realType
       tup0( 9,itu4) =  0.63046e-02_realType
       tup0(10,itu4) =  0.76520e-02_realType
       tup0(11,itu4) =  0.90512e-02_realType
       tup0(12,itu4) =  0.10438e-01_realType
       tup0(13,itu4) =  0.11748e-01_realType
       tup0(14,itu4) =  0.12920e-01_realType
       tup0(15,itu4) =  0.13898e-01_realType
       tup0(16,itu4) =  0.14631e-01_realType
       tup0(17,itu4) =  0.15089e-01_realType
       tup0(18,itu4) =  0.15297e-01_realType
       tup0(19,itu4) =  0.15285e-01_realType
       tup0(20,itu4) =  0.15080e-01_realType
       tup0(21,itu4) =  0.14711e-01_realType
       tup0(22,itu4) =  0.14207e-01_realType
       tup0(23,itu4) =  0.13596e-01_realType
       tup0(24,itu4) =  0.12904e-01_realType
       tup0(25,itu4) =  0.12154e-01_realType
       tup0(26,itu4) =  0.11369e-01_realType
       tup0(27,itu4) =  0.10567e-01_realType
       tup0(28,itu4) =  0.97653e-02_realType
       tup0(29,itu4) =  0.89762e-02_realType
       tup0(30,itu4) =  0.82104e-02_realType
       tup0(31,itu4) =  0.74760e-02_realType
       tup0(32,itu4) =  0.67787e-02_realType
       tup0(33,itu4) =  0.61223e-02_realType
       tup0(34,itu4) =  0.55092e-02_realType

       tup1( 1,itu4) =  0.10014e-02_realType
       tup1( 2,itu4) =  0.98811e-03_realType
       tup1( 3,itu4) =  0.97209e-03_realType
       tup1( 4,itu4) =  0.95225e-03_realType
       tup1( 5,itu4) =  0.92674e-03_realType
       tup1( 6,itu4) =  0.89247e-03_realType
       tup1( 7,itu4) =  0.84501e-03_realType
       tup1( 8,itu4) =  0.77980e-03_realType
       tup1( 9,itu4) =  0.69604e-03_realType
       tup1(10,itu4) =  0.59785e-03_realType
       tup1(11,itu4) =  0.49284e-03_realType
       tup1(12,itu4) =  0.38934e-03_realType
       tup1(13,itu4) =  0.29358e-03_realType
       tup1(14,itu4) =  0.20915e-03_realType
       tup1(15,itu4) =  0.13748e-03_realType
       tup1(16,itu4) =  0.79309e-04_realType
       tup1(17,itu4) =  0.36899e-04_realType
       tup1(18,itu4) =  0.90827e-05_realType
       tup1(19,itu4) = -0.84239e-05_realType
       tup1(20,itu4) = -0.18715e-04_realType
       tup1(21,itu4) = -0.24036e-04_realType
       tup1(22,itu4) = -0.26019e-04_realType
       tup1(23,itu4) = -0.25851e-04_realType
       tup1(24,itu4) = -0.24377e-04_realType
       tup1(25,itu4) = -0.22183e-04_realType
       tup1(26,itu4) = -0.19673e-04_realType
       tup1(27,itu4) = -0.17109e-04_realType
       tup1(28,itu4) = -0.14652e-04_realType
       tup1(29,itu4) = -0.12395e-04_realType
       tup1(30,itu4) = -0.10383e-04_realType
       tup1(31,itu4) = -0.86269e-05_realType
       tup1(32,itu4) = -0.71206e-05_realType
       tup1(33,itu4) = -0.58456e-05_realType
       tup1(34,itu4) = -0.47776e-05_realType

       tup2( 1,itu4) = -0.73381e-05_realType
       tup2( 2,itu4) = -0.67793e-05_realType
       tup2( 3,itu4) = -0.58449e-05_realType
       tup2( 4,itu4) = -0.49347e-05_realType
       tup2( 5,itu4) = -0.41793e-05_realType
       tup2( 6,itu4) = -0.40232e-05_realType
       tup2( 7,itu4) = -0.54838e-05_realType
       tup2( 8,itu4) = -0.91450e-05_realType
       tup2( 9,itu4) = -0.12056e-04_realType
       tup2(10,itu4) = -0.13441e-04_realType
       tup2(11,itu4) = -0.12747e-04_realType
       tup2(12,itu4) = -0.10761e-04_realType
       tup2(13,itu4) = -0.84210e-05_realType
       tup2(14,itu4) = -0.62599e-05_realType
       tup2(15,itu4) = -0.45849e-05_realType
       tup2(16,itu4) = -0.34480e-05_realType
       tup2(17,itu4) = -0.19626e-05_realType
       tup2(18,itu4) = -0.10949e-05_realType
       tup2(19,itu4) = -0.59189e-06_realType
       tup2(20,itu4) = -0.30260e-06_realType
       tup2(21,itu4) = -0.14033e-06_realType
       tup2(22,itu4) = -0.51286e-07_realType
       tup2(23,itu4) = -0.65167e-08_realType
       tup2(24,itu4) =  0.14883e-07_realType
       tup2(25,itu4) =  0.22563e-07_realType
       tup2(26,itu4) =  0.23657e-07_realType
       tup2(27,itu4) =  0.21434e-07_realType
       tup2(28,itu4) =  0.18140e-07_realType
       tup2(29,itu4) =  0.14649e-07_realType
       tup2(30,itu4) =  0.11482e-07_realType
       tup2(31,itu4) =  0.88134e-08_realType
       tup2(32,itu4) =  0.66562e-08_realType
       tup2(33,itu4) =  0.49705e-08_realType
       tup2(34,itu4) =  0.36802e-08_realType

       tup3( 1,itu4) = -0.26409e-04_realType
       tup3( 2,itu4) = -0.19723e-04_realType
       tup3( 3,itu4) = -0.15989e-04_realType
       tup3( 4,itu4) = -0.13800e-04_realType
       tup3( 5,itu4) = -0.12497e-04_realType
       tup3( 6,itu4) = -0.11352e-04_realType
       tup3( 7,itu4) = -0.94581e-05_realType
       tup3( 8,itu4) = -0.64516e-05_realType
       tup3( 9,itu4) = -0.38268e-05_realType
       tup3(10,itu4) = -0.19060e-05_realType
       tup3(11,itu4) = -0.82598e-06_realType
       tup3(12,itu4) = -0.31494e-06_realType
       tup3(13,itu4) = -0.96327e-07_realType
       tup3(14,itu4) = -0.11829e-07_realType
       tup3(15,itu4) =  0.30746e-07_realType
       tup3(16,itu4) =  0.70108e-07_realType
       tup3(17,itu4) =  0.37205e-07_realType
       tup3(18,itu4) =  0.19880e-07_realType
       tup3(19,itu4) =  0.10689e-07_realType
       tup3(20,itu4) =  0.57169e-08_realType
       tup3(21,itu4) =  0.30479e-08_realType
       tup3(22,itu4) =  0.15779e-08_realType
       tup3(23,itu4) =  0.82116e-09_realType
       tup3(24,itu4) =  0.40692e-09_realType
       tup3(25,itu4) =  0.19883e-09_realType
       tup3(26,itu4) =  0.90676e-10_realType
       tup3(27,itu4) =  0.39255e-10_realType
       tup3(28,itu4) =  0.14131e-10_realType
       tup3(29,itu4) =  0.32481e-11_realType
       tup3(30,itu4) = -0.10665e-11_realType
       tup3(31,itu4) = -0.24047e-11_realType
       tup3(32,itu4) = -0.24412e-11_realType
       tup3(33,itu4) = -0.20525e-11_realType
       tup3(34,itu4) = -0.15829e-11_realType

       ! Constants for nut.

       tup0( 1,itu5) =  0.23125e-04_realType
       tup0( 2,itu5) =  0.21542e-03_realType
       tup0( 3,itu5) =  0.11557e-02_realType
       tup0( 4,itu5) =  0.45993e-02_realType
       tup0( 5,itu5) =  0.14891e-01_realType
       tup0( 6,itu5) =  0.40523e-01_realType
       tup0( 7,itu5) =  0.93406e-01_realType
       tup0( 8,itu5) =  0.19626e+00_realType
       tup0( 9,itu5) =  0.44036e+00_realType
       tup0(10,itu5) =  0.83627e+00_realType
       tup0(11,itu5) =  0.14308e+01_realType
       tup0(12,itu5) =  0.22842e+01_realType
       tup0(13,itu5) =  0.34675e+01_realType
       tup0(14,itu5) =  0.50606e+01_realType
       tup0(15,itu5) =  0.71506e+01_realType
       tup0(16,itu5) =  0.98324e+01_realType
       tup0(17,itu5) =  0.13213e+02_realType
       tup0(18,itu5) =  0.17422e+02_realType
       tup0(19,itu5) =  0.22609e+02_realType
       tup0(20,itu5) =  0.28946e+02_realType
       tup0(21,itu5) =  0.36629e+02_realType
       tup0(22,itu5) =  0.45870e+02_realType
       tup0(23,itu5) =  0.56903e+02_realType
       tup0(24,itu5) =  0.69976e+02_realType
       tup0(25,itu5) =  0.85344e+02_realType
       tup0(26,itu5) =  0.10326e+03_realType
       tup0(27,itu5) =  0.12397e+03_realType
       tup0(28,itu5) =  0.14767e+03_realType
       tup0(29,itu5) =  0.17452e+03_realType
       tup0(30,itu5) =  0.20457e+03_realType
       tup0(31,itu5) =  0.23779e+03_realType
       tup0(32,itu5) =  0.27396e+03_realType
       tup0(33,itu5) =  0.31270e+03_realType
       tup0(34,itu5) =  0.35337e+03_realType

       tup1( 1,itu5) =  0.37001e-03_realType
       tup1( 2,itu5) =  0.15197e-02_realType
       tup1( 3,itu5) =  0.46103e-02_realType
       tup1( 4,itu5) =  0.11376e-01_realType
       tup1( 5,itu5) =  0.23543e-01_realType
       tup1( 6,itu5) =  0.40906e-01_realType
       tup1( 7,itu5) =  0.64800e-01_realType
       tup1( 8,itu5) =  0.11582e+00_realType
       tup1( 9,itu5) =  0.17214e+00_realType
       tup1(10,itu5) =  0.21559e+00_realType
       tup1(11,itu5) =  0.25612e+00_realType
       tup1(12,itu5) =  0.29399e+00_realType
       tup1(13,itu5) =  0.32838e+00_realType
       tup1(14,itu5) =  0.35836e+00_realType
       tup1(15,itu5) =  0.38343e+00_realType
       tup1(16,itu5) =  0.40383e+00_realType
       tup1(17,itu5) =  0.42066e+00_realType
       tup1(18,itu5) =  0.43496e+00_realType
       tup1(19,itu5) =  0.44719e+00_realType
       tup1(20,itu5) =  0.45758e+00_realType
       tup1(21,itu5) =  0.46621e+00_realType
       tup1(22,itu5) =  0.47302e+00_realType
       tup1(23,itu5) =  0.47789e+00_realType
       tup1(24,itu5) =  0.48067e+00_realType
       tup1(25,itu5) =  0.48112e+00_realType
       tup1(26,itu5) =  0.47894e+00_realType
       tup1(27,itu5) =  0.47384e+00_realType
       tup1(28,itu5) =  0.46551e+00_realType
       tup1(29,itu5) =  0.45364e+00_realType
       tup1(30,itu5) =  0.43792e+00_realType
       tup1(31,itu5) =  0.41811e+00_realType
       tup1(32,itu5) =  0.39406e+00_realType
       tup1(33,itu5) =  0.36565e+00_realType
       tup1(34,itu5) =  0.33291e+00_realType

       tup2( 1,itu5) = -0.15144e-02_realType
       tup2( 2,itu5) = -0.21673e-02_realType
       tup2( 3,itu5) = -0.22502e-02_realType
       tup2( 4,itu5) = -0.80447e-03_realType
       tup2( 5,itu5) =  0.27762e-02_realType
       tup2( 6,itu5) =  0.17677e-02_realType
       tup2( 7,itu5) = -0.10701e-01_realType
       tup2( 8,itu5) =  0.22380e-01_realType
       tup2( 9,itu5) =  0.85050e-02_realType
       tup2(10,itu5) =  0.62055e-02_realType
       tup2(11,itu5) =  0.49141e-02_realType
       tup2(12,itu5) =  0.39356e-02_realType
       tup2(13,itu5) =  0.30516e-02_realType
       tup2(14,itu5) =  0.22492e-02_realType
       tup2(15,itu5) =  0.15552e-02_realType
       tup2(16,itu5) =  0.10053e-02_realType
       tup2(17,itu5) =  0.70099e-03_realType
       tup2(18,itu5) =  0.50284e-03_realType
       tup2(19,itu5) =  0.36969e-03_realType
       tup2(20,itu5) =  0.27228e-03_realType
       tup2(21,itu5) =  0.20061e-03_realType
       tup2(22,itu5) =  0.14330e-03_realType
       tup2(23,itu5) =  0.98897e-04_realType
       tup2(24,itu5) =  0.62505e-04_realType
       tup2(25,itu5) =  0.33283e-04_realType
       tup2(26,itu5) =  0.68299e-05_realType
       tup2(27,itu5) = -0.14485e-04_realType
       tup2(28,itu5) = -0.32453e-04_realType
       tup2(29,itu5) = -0.47992e-04_realType
       tup2(30,itu5) = -0.61021e-04_realType
       tup2(31,itu5) = -0.71877e-04_realType
       tup2(32,itu5) = -0.79636e-04_realType
       tup2(33,itu5) = -0.86125e-04_realType
       tup2(34,itu5) = -0.90063e-04_realType

       tup3( 1,itu5) =  0.66740e-02_realType
       tup3( 2,itu5) =  0.93418e-02_realType
       tup3( 3,itu5) =  0.10767e-01_realType
       tup3( 4,itu5) =  0.97000e-02_realType
       tup3( 5,itu5) =  0.58162e-02_realType
       tup3( 6,itu5) =  0.58753e-02_realType
       tup3( 7,itu5) =  0.14886e-01_realType
       tup3( 8,itu5) = -0.21765e-02_realType
       tup3( 9,itu5) =  0.66623e-03_realType
       tup3(10,itu5) =  0.46817e-03_realType
       tup3(11,itu5) =  0.24852e-03_realType
       tup3(12,itu5) =  0.10079e-03_realType
       tup3(13,itu5) =  0.25392e-04_realType
       tup3(14,itu5) = -0.29448e-05_realType
       tup3(15,itu5) = -0.56548e-05_realType
       tup3(16,itu5) =  0.17092e-05_realType
       tup3(17,itu5) =  0.17225e-05_realType
       tup3(18,itu5) =  0.97328e-06_realType
       tup3(19,itu5) =  0.61893e-07_realType
       tup3(20,itu5) = -0.52008e-06_realType
       tup3(21,itu5) = -0.93452e-06_realType
       tup3(22,itu5) = -0.10961e-05_realType
       tup3(23,itu5) = -0.11719e-05_realType
       tup3(24,itu5) = -0.11584e-05_realType
       tup3(25,itu5) = -0.11184e-05_realType
       tup3(26,itu5) = -0.10094e-05_realType
       tup3(27,itu5) = -0.90341e-06_realType
       tup3(28,itu5) = -0.79513e-06_realType
       tup3(29,itu5) = -0.68439e-06_realType
       tup3(30,itu5) = -0.57895e-06_realType
       tup3(31,itu5) = -0.47935e-06_realType
       tup3(32,itu5) = -0.39672e-06_realType
       tup3(33,itu5) = -0.31868e-06_realType
       tup3(34,itu5) = -0.25374e-06_realType

    end select

  end subroutine initCurveFitDataVf
end module turbCurveFits
