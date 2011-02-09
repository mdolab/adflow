!
!      ******************************************************************
!      *                                                                *
!      * File:          initCurveFitDataSa.f90                          *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 07-27-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initCurveFitDataSa
!
!      ******************************************************************
!      *                                                                *
!      * initCurveFitDataSa contains the curve fit constants for        *
!      * the wall function data for the Spalart-Allmaras turbulence     *
!      * model.                                                         *
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
