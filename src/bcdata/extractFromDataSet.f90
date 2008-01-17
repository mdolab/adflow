!
!      ******************************************************************
!      *                                                                *
!      * File:          extractFromDataSet.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-24-2004                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine extractFromDataSet(blockFaceID)
!
!      ******************************************************************
!      *                                                                *
!      * extractFromDataSet tries to extract and interpolate the        *
!      * variables in bcVarNames from the cgns data set.                *
!      * If successful the corresponding entry of bcVarPresent is       *
!      * set to .true., otherwise it is set to .false.                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: blockFaceID
!
!      Local variables.
!
       integer(kind=intType) :: k, l, m, n
       integer(kind=intType) :: nInter, nDim, nVarPresent, nCoor

       integer(kind=intType), dimension(3) :: dataDim, coor
       integer(kind=intType), dimension(2,3) :: indCoor
       integer(kind=intType), dimension(2,nbcVar) :: ind

       character(len=maxStringLen) :: errorMessage

       logical :: xPresent, yPresent, zPresent, rPresent
       logical :: sameInterpol, firstVar
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether the variables are specified and if so,
       ! where they are located in the data set. As the number of
       ! variables specified is usually not so big, a linear search
       ! algorithm is perfectly okay. At the moment only the Dirichlet
       ! arrays are checked.

       nVarPresent = 0

       do m=1,nbcVar
         bcVarPresent(m) = .false.

         dataSetLoop: do k=1,nDataSet
           do l=1,dataSet(k)%nDirichletArrays
             if(dataSet(k)%dirichletArrays(l)%arrayName == &
                bcVarNames(m)) then

               ! Variable is present. Store the indices, update
               ! nVarPresent and set bcVarPresent(m) to .True.

               ind(1,m) = k; ind(2,m) = l

               nVarPresent      = nVarPresent + 1
               bcVarPresent(m) = .true.

               ! Set the units for this variable.

               mass(m)   = dataSet(k)%dirichletArrays(l)%mass
               length(m) = dataSet(k)%dirichletArrays(l)%len
               time(m)   = dataSet(k)%dirichletArrays(l)%time
               temp(m)   = dataSet(k)%dirichletArrays(l)%temp
               angle(m)  = dataSet(k)%dirichletArrays(l)%angle

               ! Exit the search loop, as the variable was found.

               exit dataSetLoop

             endif
           enddo
         enddo dataSetLoop
       enddo

       ! Return if none of the variables are present.

       if(nVarPresent == 0) return

       ! Find out whether the given data points are equal for every
       ! variable or that every variable must be interpolated
       ! differently.

       sameInterpol = .true.
       firstVar     = .true.

       do m=1,nbcVar
         if( bcVarPresent(m) ) then
           k = ind(1,m)
           l = ind(2,m)

           if( firstVar ) then
             nDim = dataSet(k)%dirichletArrays(l)%nDimensions
             firstVar = .false.

             do n=1,nDim
               dataDim(n) = dataSet(k)%dirichletArrays(l)%dataDim(n)
             enddo
           else
             if(nDim == dataSet(k)%dirichletArrays(l)%nDimensions) then
               do n=1,nDim
                 if(dataSet(k)%dirichletArrays(l)%dataDim(n) /= &
                    dataDim(n)) sameInterpol = .false.
               enddo
             else
               sameInterpol = .false.
             endif
           endif

         endif
       enddo

       ! Determine the situation we are dealing with here.

       testSameInterpol: if( sameInterpol ) then

         ! The interpolation is the same for all variables.
         ! First determine the number of interpolation points.

         nInter = dataDim(1)
         do m=2,nDim
           nInter = nInter*dataDim(m)
         enddo

         ! If nInter == 1 then the prescribed data is constant
         ! everywhere and the variables can be determined easily.

         testConstant1: if(nInter == 1) then

           ! Data is constant for this subface. Set the data.

           do m=1,nbcVar
             if( bcVarPresent(m) ) then
               k = ind(1,m)
               l = ind(2,m)

               bcVarArray(:,:,m) = &
                  dataSet(k)%dirichletArrays(l)%dataArr(1)
             endif
           enddo

         else testConstant1

           ! Data varies over the interface and must be interpolated.
           ! Determine the indices of the coordinates in the dataset.

           rPresent = .false.
           xPresent = .false.
           yPresent = .false.
           zPresent = .false.

           do k=1,nDataSet
             do l=1,dataSet(k)%nDirichletArrays

               if(dataSet(k)%dirichletArrays(l)%arrayName == &
                  cgnsCoorr) then
                 indCoor(1,1) = k; indCoor(2,1) = l
                 rPresent = .true.
                 exit
               endif

               if(dataSet(k)%dirichletArrays(l)%arrayName == &
                  cgnsCoorx) then
                 indCoor(1,1) = k; indCoor(2,1) = l
                 xPresent = .true.
               endif

               if(dataSet(k)%dirichletArrays(l)%arrayName == &
                  cgnsCoory) then
                 indCoor(1,2) = k; indCoor(2,2) = l
                 yPresent = .true.
               endif

               if(dataSet(k)%dirichletArrays(l)%arrayName == &
                  cgnsCoorz) then
                 indCoor(1,3) = k; indCoor(2,3) = l
                 zPresent = .true.
               endif

             enddo
           enddo

           ! Check if a radial coordinate is present.

           if( rPresent ) then

             ! Radial coordinate is present. Use radial interpolation
             ! for the given variable.

             call radialInterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                                        nbcVar, cgnsBoco,       &
                                        blockFaceID, indCoor,   &
                                        ind, bcVarPresent,      &
                                        bcVarArray, axAssumed)

           else if(xPresent .or. yPresent .or. zPresent) then

             ! Cartesian interpolation will be performed. Determine
             ! which coordinates are present.

             nCoor = 0
             if( xPresent ) then
               nCoor = nCoor + 1; coor(nCoor) = 1
             endif
             if( yPresent ) then
               nCoor = nCoor + 1; coor(nCoor) = 2
             endif
             if( zPresent ) then
               nCoor = nCoor + 1; coor(nCoor) = 3
             endif

             ! The number of dimensions cannot be larger than the
             ! number of coordinates. Check this.

             if(nDim > nCoor) then
               write(errorMessage,100) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
               call terminate("extractFromDataSet", errorMessage)
             endif

             ! Check what kind of interpolation must be used.

             select case (nCoor)

               case (1_intType)

                 ! 1D line interpolation.

                 call cart1D_InterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                                             nbcVar, cgnsBoco,       &
                                             blockFaceID, coor(1),   &
                                             indCoor, ind,           &
                                             bcVarPresent,           &
                                             bcVarArray)

               !=======================================================

               case (2_intType, 3_intType)

                  call terminate("extractFromDataSet", &
                                 "Multi-D Cartesian interpolation &
                                 &not implemented yet")

             end select

           else

             ! Neither the radial nor the cartesian coordinates are
             ! present. So there is not enough information available
             ! for the interpolation. Print an error message and exit.

             write(errorMessage,101) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName)
               call terminate("extractFromDataSet", errorMessage)

           endif

         endif testConstant1

       else testSameInterpol

         ! Different interpolation must be used for the different
         ! variables. Loop over the number of variables and test
         ! whether they are present.

         bcVarLoop: do m=1,nbcVar
           testBcVarPresent: if( bcVarPresent(m) ) then

             ! Store the indices of the corresponding data set and
             ! dirichlet array a bit easier. Abbreviate the number
             ! of dimensions a bit easier.

             k    = ind(1,m)
             l    = ind(2,m)
             nDim = dataSet(k)%dirichletArrays(l)%nDimensions

             ! Determine the number of interpolation points.

             nInter = dataSet(k)%dirichletArrays(l)%dataDim(1)
             do n=2,nDim
               nInter = nInter &
                      * dataSet(k)%dirichletArrays(l)%dataDim(n)
             enddo

             ! If nInter == 1 then the prescribed data is constant
             ! everywhere and the variable can be determined easily.

             testConstant2: if(nInter == 1) then

               ! Data is constant for this subface. Set it

               bcVarArray(:,:,m) = &
                   dataSet(k)%dirichletArrays(l)%dataArr(1)

             else testConstant2

               ! Data varies over the interface and must be
               ! interpolated. Determine the indices of the
               ! coordinates in the dirichlet arrays of the given
               ! data set. Note that the coordinates now have to be
               ! specified in the same dataset as the variable.

               rPresent = .false.
               xPresent = .false.
               yPresent = .false.
               zPresent = .false.

               do l=1,dataSet(k)%nDirichletArrays

                 if(dataSet(k)%dirichletArrays(l)%arrayName == &
                    cgnsCoorr) then
                   indCoor(1,1) = k; indCoor(2,1) = l
                   rPresent = .true.
                   exit
                 endif

                 if(dataSet(k)%dirichletArrays(l)%arrayName == &
                    cgnsCoorx) then
                   indCoor(1,1) = k; indCoor(2,1) = l
                   xPresent = .true.
                 endif

                 if(dataSet(k)%dirichletArrays(l)%arrayName == &
                    cgnsCoory) then
                   indCoor(1,2) = k; indCoor(2,2) = l
                   yPresent = .true.
                 endif

                 if(dataSet(k)%dirichletArrays(l)%arrayName == &
                    cgnsCoorz) then
                   indCoor(1,3) = k; indCoor(2,3) = l
                   zPresent = .true.
                 endif

               enddo

               ! Check if a radial coordinate is present.

               if( rPresent ) then

                 ! Radial coordinate is present. Use radial
                 ! interpolation for the given variable.

                 call radialInterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                                            1_intType, cgnsBoco,    &
                                            blockFaceID, indCoor,   &
                                            ind(1,m),               &
                                            bcVarPresent(m),        &
                                            bcVarArray(1,1,m),      &
                                            axAssumed)

               else if(xPresent .or. yPresent .or. zPresent) then

                 ! Cartesian interpolation will be performed. Determine
                 ! which coordinates are present.

                 nCoor = 0
                 if( xPresent ) then
                   nCoor = nCoor + 1; coor(nCoor) = 1
                 endif
                 if( yPresent ) then
                   nCoor = nCoor + 1; coor(nCoor) = 2
                 endif
                 if( zPresent ) then
                   nCoor = ncoor + 1; coor(nCoor) = 3
                 endif

                 ! The number of dimensions cannot be larger than the
                 ! number of coordinates. Check this.

                 if(nDim > nCoor) then
                   write(errorMessage,200) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName), &
            trim(bcVarNames(m))
                   call terminate("extractFromDataSet", errorMessage)
                 endif

                 ! Check what kind of interpolation must be used.

                 select case (nCoor)

                   case (1_intType)

                     ! 1D line interpolation.

                     call cart1D_InterpolSubface(iBeg, iEnd, jBeg, jEnd, &
                                                 1_intType, cgnsBoco,    &
                                                 blockFaceID, coor(1),   &
                                                 indCoor, ind,           &
                                                 bcVarPresent(m),        &
                                                 bcVarArray(1,1,m))

                   !===================================================

                   case (2_intType, 3_intType)

                      call terminate("extractFromDataSet", &
                                     "Multi-D Cartesian interpolation &
                                     &not implemented yet")

                 end select

               else

                 ! Neither the radial nor the cartesian coordinates
                 ! are present. So there is not enough information
                 ! available for the interpolation. Print an error
                 ! message and exit.

                 write(errorMessage,201) &
            trim(cgnsDoms(nbkGlobal)%zonename), &
            trim(cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%bocoName), &
            trim(bcVarNames(m))
                   call terminate("extractFromDataSet", &
                                  errorMessage)
               endif

             endif testConstant2

           endif testBcVarPresent
         enddo bcVarLoop

       endif testSameInterpol

       ! Format statements.

 100   format("Zone",1X,A,", subface",1X,A,": Number of dimensions &
              &is larger than number of coordinates.")
 101   format("Zone",1X,A,", subface",1X,A,": No coordinates &
              &are present for the interpolation.")
 200   format("Zone",1X,A,", subface",1X,A,", variable",1X,A, &
              ": Number of dimensions is larger than number of &
              &coordinates.")
 201   format("Zone",1X,A,", subface",1X,A,": No coordinates &
              &are present for the interpolation of",1X,A,".")

       end subroutine extractFromDataSet
