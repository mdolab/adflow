!
!      ******************************************************************
!      *                                                                *
!      * File:          coeffConvegenceCheck.f90                        *
!      * Author:        Eran Arad, Yuval Dagan                          *
!      * Starting date: 04-08-2010                                      *
!      * Last modified: 13-05-2010                                                *
!      *                                                                *
!      ******************************************************************
!
 subroutine coeffConvergenceCheck(iConv,iterTot,sps)
!
!      ******************************************************************
!      *   Evaluate the level of monitored coefficients convergence     *
!      *   Activated by input parameter   eran-coeffConv                *
!      *                                                                *
!      *   "coefficients convergence criterion" > 0   (epsCoefConv)     *
!      *                                                                *
!      *  Basically, over interation window coeffConvIterWin, the mean  *
!      *  and std dev. (sigma) of each monitored coefficient is computed*
!      *  if  sigma/mean < epsCoefConv then set converged = true, even  *
!      *  if residual convergence criterin was not fullfilled yet.      *
!      *                                                                *
!      *  Triple checks are performed, where to window size and         *
!      *  convergence criterion are multiplied by 1,10 and 100.         *
!      *  The triple check was implemented to handle cases where        *
!      *  non-decaying oscillations of coefficients develop, over       *
!      *  a large period of iterations.                                 *
!      *  The program returns also the convergenceQuality parameter:    *
!      *                                                                *
!      *     convergenceQuality = 0   No convergence                    *
!      *     convergenceQuality = 2   Coefficient-convergence X100      *
!      *     convergenceQuality = 4   Coefficient-convergence X10       *
!      *     convergenceQuality = 6   Coefficient-convergence X1        *
!      *     convergenceQuality = 10  residue convergence               *
!      *                                                                *
!      *  This scheme is required for the automatic polar sweep.        *
!      *  Use it very carefully, since it might lead to unconverged     *
!      *  solutions. If reaching residual convegence criterion is       *
!      *  possible, it is much safer to do that.                        *
!      *                                                                *
!      ******************************************************************
!
   use monitor
   use cgnsNames
   use constants
   use inputIteration
   implicit none
! routine variables
   integer(kind=intType) :: iConv,sps,iterTot
!
! parameters
!
   integer, parameter :: numLevelsOfConvCheck = 3 ! (that is times 1,10,100)
!
! Global variables
!
   integer(kind=intType) ,save :: levelsConvCheckWindowSize(numLevelsOfConvCheck)
   integer(kind=intType) ,save :: nMonCoeff
   integer(kind=intType) ,save,allocatable :: mcoeff(:)
   real(kind=realType), save,allocatable :: varMean(:,:), varStd(:,:)
   real(kind=realType), save :: roWinSize(numLevelsOfConvCheck),&
        levelsEpsCoefConv(numLevelsOfConvCheck)
!
! Local variables
!
   integer(kind=intType) :: mm, i, level, testConv, actualNumLevelsOfConvCheck
   real(kind=realType)  :: aVarMean, varStdMax
!
!----------------------------------------------------------------------
!    
   if(iterTot == minIterNum)then
!
! First time that you are here:
!-----------------------------------
! Allocate arrays
! Go through the monitored variables and locate the relevant coefficients
! Set inital values to global parameters
!
      nMonCoeff = 0
      allocate(mcoeff(nMon),varMean(nMon,numLevelsOfConvCheck),&
           varStd(nMon,numLevelsOfConvCheck))
      do level=1,numLevelsOfConvCheck
         levelsConvCheckWindowSize(level) = 10**(level-1)*ConvCheckWindowSize
         roWinSize(level)    = one/float(levelsConvCheckWindowSize(level))
         levelsEpsCoefConv(level) = 10.0**(level-1)*epsCoefConv
      end do ! level

      do mm=1,nMon
      
         select case (monNames(mm))

         case (cgnsCl,  cgnsClp,cgnsCd,cgnsCdp,&
              cgnsCfx, cgnsCfy,cgnsCfz,&
              cgnsCmx, cgnsCmy,cgnsCmz)
            nMonCoeff         = nMonCoeff + 1
            mcoeff(nMonCoeff) = mm
            varMean(nMonCoeff,:)=zero
            varStd(nMonCoeff,:) =zero
         end select
      end do ! mm
     

   end if !  First time
!--------------------------------------------------------------------------------------
!
!-------------------compute mean value
!
   do level=1,numLevelsOfConvCheck
      do mm=1,nMonCoeff
         varMean(mm,level) =varMean(mm,level) + roWinSize(level)*convArray(iConv,sps,mcoeff(mm))
      end do
   end do

   if (iterTot == minIterNum)return ! no need to update convergenceQuality since
                                    ! convergence test was not performed yet

   varStdMax = 1000.0  ! default large value
   
! find out if for present iteration, STD should be computed, at any level
 
   testConv = 1
   do level=1,numLevelsOfConvCheck
      testConv = testConv*mod(iterTot-minIterNum,levelsConvCheckWindowSize(level))
   end do

   if (testConv == 0 )then
!
! Set the acual level of convergence check, that might be smaller than 
! numLevelsOfConvCheck since not enough iterations have been performed
! and the lebvel window size is too large
!
      actualNumLevelsOfConvCheck = numLevelsOfConvCheck
      do level=numLevelsOfConvCheck,1,-1
         if((iterTot-minIterNum) < levelsConvCheckWindowSize(level))&
              actualNumLevelsOfConvCheck = actualNumLevelsOfConvCheck -1
      end do
      if (actualNumLevelsOfConvCheck == 0)return ! not enough iterations even
!                                                  for the smallest window size
! ---------compute STD
!
      do level=1,actualNumLevelsOfConvCheck
         if(mod(iterTot-minIterNum,levelsConvCheckWindowSize(level)) == 0)then
            do mm=1,nMonCoeff
               do i=iConv-levelsConvCheckWindowSize(level),iConv
                  varStd(mm,level) = varStd(mm,level) +  (convArray(i,sps,mcoeff(mm)) - &
                       varMean(mm,level))**2 
               end do
               varStd(mm,level) = sqrt( varStd(mm,level)/&
                    float(levelsConvCheckWindowSize(level)-1))
               aVarMean = abs(varMean(mm,level))
               if(aVarMean  > eps )varStd(mm,level) = varStd(mm,level)/aVarMean
            end do
!
! Locate the largest STD
!
            varStdMax = maxval(varStd(:,level))
!
! Check for convergence
!
            if(varStdMax <= epsCoefConv)convergenceQuality = 8 - 2*level
! 
! now prepare for next cycle
!
            varMean(:,level) = zero
            varStd(:,level)  = zero
         end if ! mod
      end do ! level 

   end if ! testConv

   return
 end subroutine coeffConvergenceCheck
