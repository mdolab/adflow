!
!      ******************************************************************
!      *                                                                *
!      * File:          initCurveFitDataSST.f90                         *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 08-21-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initCurveFitDataSST
!
!      ******************************************************************
!      *                                                                *
!      * initCurveFitDataSST contains the curve fit constants for       *
!      * the wall function data for Menter's SST turbulence model.      *
!      *                                                                *
!      * Warning: Wall function data developed for k-omega model        *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use paramTurb
       implicit none
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
