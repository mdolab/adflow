!
!      ******************************************************************
!      *                                                                *
!      * File:          determineBleedFlowRegions.f90                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-11-2005                                      *
!      * Last modified: 08-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineBleedFlowRegions
!
!      ******************************************************************
!      *                                                                *
!      * determineBleedFlowAreas determines the number of inflow and    *
!      * outflow bleed regions as well as the mapping to and from the   *
!      * corresponding family ID.                                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use bleedFlows
       use cgnsGrid
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, ii, jj
!
!      Function definitions
!
       real(kind=realType) :: getMassFlux
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of inflow and outflow bleed regions.

       nInflowBleeds  = 0
       nOutflowBleeds = 0

       do nn=1,cgnsNFamilies
         select case (cgnsFamilies(nn)%BCType)
           case (MassBleedInflow)
             nInflowBleeds = nInflowBleeds + 1

           case (MassBleedOutflow)
             nOutflowBleeds = nOutflowBleeds + 1

         end select
       enddo

       ! Allocate the memory for the arrays to store the bleed flow
       ! information.

       allocate( inflowBleeds( nInflowBleeds), &
                outflowBleeds(nOutflowBleeds), stat=ierr)
       if(ierr /= 0)                                 &
         call terminate("determineBleedFlowRegions", &
                        "Memory allocation failure for inflowBleeds &
                        &and outflowBleeds")

       ! Loop again over the families but now store some info;
       ! ii is the counter for the inflow and jj for the outflow.

       ii = 0
       jj = 0

       do nn=1,cgnsNFamilies
         select case (cgnsFamilies(nn)%BCType)
           case (MassBleedInflow)

             ! Inflow bleed. Store the mapping to and from the
             ! corresponding family and determine the mass
             ! flux for this region.

             ii = ii + 1
             cgnsFamilies(nn)%bleedRegionID = ii
             inflowBleeds(ii)%familyID      = nn

             inflowBleeds(ii)%massFlux = getMassFlux(nn)

           !=============================================================

           case (MassBleedOutflow)

             ! Outflow bleed. Store the mapping to and from the
             ! corresponding family and determine the relative mass
             ! flux for this region.

             jj = jj + 1
             cgnsFamilies(nn)%bleedRegionID = jj
             outflowBleeds(jj)%familyID     = nn

             outflowBleeds(jj)%massFlux = getMassFlux(nn)

         end select
       enddo

       end subroutine determineBleedFlowRegions

!      ==================================================================

       function getMassFlux(famID)
!
!      ******************************************************************
!      *                                                                *
!      * getMassFlux extracts the relative mass flux from the           *
!      * prescribed data for the given family. If not present           *
!      * processor 0 will print an error message and the computation    *
!      * terminates.                                                    *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use cgnsNames
       use communication
       implicit none
!
!      Function type.
!
       real(kind=realType) :: getMassFlux
!
!      Function arguments.
!
       integer(kind=intType), intent(in) :: famID
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm
       integer(kind=intType) :: nDataSet

       character(len=maxStringLen) :: errorMessage

       type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the number of boundary condition data sets and the 
       ! data sets a bit easier.

       nDataSet = cgnsFamilies(famID)%nDataSet
       dataSet => cgnsFamilies(famID)%dataSet

       ! Loop over the number of data sets and Dirichlet arrays and
       ! find out if the relative mass flux has been prescribed.

       dataSetLoop: do nn=1,nDataSet
         do mm=1,dataSet(nn)%nDirichletArrays

           if(dataSet(nn)%dirichletArrays(mm)%arrayName == &
               cgnsMassFlow) exit dataSetLoop
         enddo
       enddo dataSetLoop

       ! Check if the mass flux was prescribed.
       ! If not, processor 0 prints an error message and the
       ! program terminates.

       if(nn > nDataSet) then
         write(errorMessage,100) trim(cgnsFamilies(famID)%familyName)
 100     format("Family ", a, ": No mass flow prescribed for &
                &a bleed flow region")
         if(myID == 0) call terminate("getMassFlux", errorMessage)
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Mass flux is present in the dataset for the family.
       ! This is an integral quantity and thus constant for the entire
       ! bleed region, even if e.g. the velocity varies. Therefore
       ! set it to the first value of the Dirichlet array.

       getMassFlux = dataSet(nn)%dirichletArrays(mm)%dataArr(1)

       end function getMassFlux
