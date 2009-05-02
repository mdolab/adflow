!
!     ******************************************************************
!     *                                                                *
!     * File:          updateFlow.f90                                  *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 24-02-2009                                      *
!     * Last modified: 24-02-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine updateFlow
!
!     ******************************************************************
!     *                                                                *
!     *  reruns the initialization routines to update AOA and other    *
!     *  flow variables after a design change                          *
!     *                                                                *
!     ******************************************************************
!
        implicit none

        call referenceState
       
        call setFlowInfinityState

      end subroutine updateFlow
