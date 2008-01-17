!
!      ******************************************************************
!      *                                                                *
!      * File:          nullifyFlowDomPointers.f90                      *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 08-16-2004                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine nullifyFlowDomPointers(nn,level,sps)
!
!      ******************************************************************
!      *                                                                *
!      * nullifyFlowDomPointers nullifies all the pointers of the       *
!      * given block.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use block
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, level, sps
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nullify(flowDoms(nn,level,sps)%BCType)
       nullify(flowDoms(nn,level,sps)%BCFaceID)
       nullify(flowDoms(nn,level,sps)%cgnsSubface)

       nullify(flowDoms(nn,level,sps)%inBeg)
       nullify(flowDoms(nn,level,sps)%jnBeg)
       nullify(flowDoms(nn,level,sps)%knBeg)
       nullify(flowDoms(nn,level,sps)%inEnd)
       nullify(flowDoms(nn,level,sps)%jnEnd)
       nullify(flowDoms(nn,level,sps)%knEnd)

       nullify(flowDoms(nn,level,sps)%dinBeg)
       nullify(flowDoms(nn,level,sps)%djnBeg)
       nullify(flowDoms(nn,level,sps)%dknBeg)
       nullify(flowDoms(nn,level,sps)%dinEnd)
       nullify(flowDoms(nn,level,sps)%djnEnd)
       nullify(flowDoms(nn,level,sps)%dknEnd)

       nullify(flowDoms(nn,level,sps)%icBeg)
       nullify(flowDoms(nn,level,sps)%jcBeg)
       nullify(flowDoms(nn,level,sps)%kcBeg)
       nullify(flowDoms(nn,level,sps)%icEnd)
       nullify(flowDoms(nn,level,sps)%jcEnd)
       nullify(flowDoms(nn,level,sps)%kcEnd)

       nullify(flowDoms(nn,level,sps)%neighBlock)
       nullify(flowDoms(nn,level,sps)%neighProc)
       nullify(flowDoms(nn,level,sps)%l1)
       nullify(flowDoms(nn,level,sps)%l2)
       nullify(flowDoms(nn,level,sps)%l3)
       nullify(flowDoms(nn,level,sps)%groupNum)

       nullify(flowDoms(nn,level,sps)%ibndry)
       nullify(flowDoms(nn,level,sps)%idonor)
       nullify(flowDoms(nn,level,sps)%overint)
       nullify(flowDoms(nn,level,sps)%neighProcOver)
       nullify(flowDoms(nn,level,sps)%neighBlockOver)
       nullify(flowDoms(nn,level,sps)%iblank)

       nullify(flowDoms(nn,level,sps)%BCData)
       nullify(flowDoms(nn,level,sps)%viscSubface)

       nullify(flowDoms(nn,level,sps)%viscIminPointer)
       nullify(flowDoms(nn,level,sps)%viscImaxPointer)
       nullify(flowDoms(nn,level,sps)%viscJminPointer)
       nullify(flowDoms(nn,level,sps)%viscJmaxPointer)
       nullify(flowDoms(nn,level,sps)%viscKminPointer)
       nullify(flowDoms(nn,level,sps)%viscKmaxPointer)

       nullify(flowDoms(nn,level,sps)%x)
       nullify(flowDoms(nn,level,sps)%xOld)
       nullify(flowDoms(nn,level,sps)%si)
       nullify(flowDoms(nn,level,sps)%sj)
       nullify(flowDoms(nn,level,sps)%sk)
       nullify(flowDoms(nn,level,sps)%vol)
       nullify(flowDoms(nn,level,sps)%volOld)

       nullify(flowDoms(nn,level,sps)%pori)
       nullify(flowDoms(nn,level,sps)%porj)
       nullify(flowDoms(nn,level,sps)%pork)

       nullify(flowDoms(nn,level,sps)%indFamilyI)
       nullify(flowDoms(nn,level,sps)%indFamilyJ)
       nullify(flowDoms(nn,level,sps)%indFamilyK)

       nullify(flowDoms(nn,level,sps)%factFamilyI)
       nullify(flowDoms(nn,level,sps)%factFamilyJ)
       nullify(flowDoms(nn,level,sps)%factFamilyK)

       nullify(flowDoms(nn,level,sps)%rotMatrixI)
       nullify(flowDoms(nn,level,sps)%rotMatrixJ)
       nullify(flowDoms(nn,level,sps)%rotMatrixK)

       nullify(flowDoms(nn,level,sps)%sFaceI)
       nullify(flowDoms(nn,level,sps)%sFaceJ)
       nullify(flowDoms(nn,level,sps)%sFaceK)

       nullify(flowDoms(nn,level,sps)%w)
       nullify(flowDoms(nn,level,sps)%wOld)
       nullify(flowDoms(nn,level,sps)%p)
       nullify(flowDoms(nn,level,sps)%gamma)
       nullify(flowDoms(nn,level,sps)%rlv)
       nullify(flowDoms(nn,level,sps)%rev)
       nullify(flowDoms(nn,level,sps)%s)

       nullify(flowDoms(nn,level,sps)%dw)
       nullify(flowDoms(nn,level,sps)%fw)

       nullify(flowDoms(nn,level,sps)%dwOldRK)

       nullify(flowDoms(nn,level,sps)%p1)
       nullify(flowDoms(nn,level,sps)%w1)
       nullify(flowDoms(nn,level,sps)%wr)

       nullify(flowDoms(nn,level,sps)%mgIFine)
       nullify(flowDoms(nn,level,sps)%mgJFine)
       nullify(flowDoms(nn,level,sps)%mgKFine)

       nullify(flowDoms(nn,level,sps)%mgIWeight)
       nullify(flowDoms(nn,level,sps)%mgJWeight)
       nullify(flowDoms(nn,level,sps)%mgKWeight)

       nullify(flowDoms(nn,level,sps)%mgICoarse)
       nullify(flowDoms(nn,level,sps)%mgJCoarse)
       nullify(flowDoms(nn,level,sps)%mgKCoarse)

       nullify(flowDoms(nn,level,sps)%ico)
       nullify(flowDoms(nn,level,sps)%jco)
       nullify(flowDoms(nn,level,sps)%kco)

       nullify(flowDoms(nn,level,sps)%wn)
       nullify(flowDoms(nn,level,sps)%pn)
       nullify(flowDoms(nn,level,sps)%dtl)
       nullify(flowDoms(nn,level,sps)%radI)
       nullify(flowDoms(nn,level,sps)%radJ)
       nullify(flowDoms(nn,level,sps)%radK)

       nullify(flowDoms(nn,level,sps)%d2Wall)

       nullify(flowDoms(nn,level,sps)%bmti1)
       nullify(flowDoms(nn,level,sps)%bmti2)
       nullify(flowDoms(nn,level,sps)%bmtj1)
       nullify(flowDoms(nn,level,sps)%bmtj2)
       nullify(flowDoms(nn,level,sps)%bmtk1)
       nullify(flowDoms(nn,level,sps)%bmtk2)

       nullify(flowDoms(nn,level,sps)%bvti1)
       nullify(flowDoms(nn,level,sps)%bvti2)
       nullify(flowDoms(nn,level,sps)%bvtj1)
       nullify(flowDoms(nn,level,sps)%bvtj2)
       nullify(flowDoms(nn,level,sps)%bvtk1)
       nullify(flowDoms(nn,level,sps)%bvtk2)

       end subroutine nullifyFlowDomPointers
