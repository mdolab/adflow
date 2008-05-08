MODULE TAPENADE_DD_TEST
   IMPLICIT NONE
   PUBLIC
   INTEGER,PUBLIC,PARAMETER :: wwp=selected_real_kind(12, 307)
   integer, PUBLIC :: ddphase = 0
   real*8, PUBLIC :: ddeps = 1.d-4
   real*8, PUBLIC, PARAMETER :: epszero = 1.d-5
   CHARACTER(len=16), PUBLIC, PARAMETER :: rwTmpFilename='rwTmpFile'
   integer, PUBLIC, PARAMETER :: rwTmpFile = 36
   CHARACTER(len=16), PUBLIC, PARAMETER :: epsValuesFilename='epsvalues'
   integer, PUBLIC, PARAMETER :: ddepsValuesFile = 37
   CHARACTER(len=16), PUBLIC, PARAMETER :: ddValuesFilename='ddvalues'
   integer, PUBLIC, PARAMETER :: ddddValuesFile = 38
   CHARACTER(len=16), PUBLIC, PARAMETER :: testResultFilename='ddtestresult'
   integer, PUBLIC, PARAMETER :: ddtestResultFile = 39
   CHARACTER(len=16), PUBLIC, PARAMETER :: eps2ValuesFilename='eps2values'
   integer, PUBLIC, PARAMETER :: ddeps2ValuesFile = 40
CONTAINS

      subroutine tracesetphase(phase)
      integer phase
      if (phase.eq.0) then
! phase that only prints values.         
         ddphase = 0
      else if (phase.eq.1) then
! phase that stores the intermediate values on "epsValuesFile"
         ddphase = 1
         OPEN(ddepsValuesFile, FILE=epsValuesFilename)
      else if (phase.eq.2) then
! phase that reads the values from "epsValuesFile", uses them to run
! divided diffs with current value, and compares with AD diffs
! Puts results into "testResultFile".
         CLOSE(ddepsValuesFile)
         ddphase = 2
         OPEN(ddepsValuesFile, FILE=epsValuesFilename)
         OPEN(ddeps2ValuesFile, FILE=eps2ValuesFilename)
         OPEN(ddtestResultFile, FILE=testResultFilename)
      else if (phase.eq.3) then
! phase that reads the values from "epsValuesFile", uses them to run
! centered divided diffs with current value,
! and stores the divided diff on "ddValuesFile".
         CLOSE(ddepsValuesFile)
         ddphase = 3
         OPEN(ddepsValuesFile, FILE=epsValuesFilename)
         OPEN(ddeps2ValuesFile, FILE=eps2ValuesFilename)
         OPEN(ddddValuesFile, FILE=ddValuesFilename)
      else if (phase.eq.4) then
! phase that reads the divided diffs from "ddValuesFile" and
! compares with AD diffs. Puts results into "testResultFile".
         CLOSE(ddepsValuesFile)
         CLOSE(ddddValuesFile)
         ddphase = 4
         OPEN(ddddValuesFile, FILE=ddValuesFilename)
         OPEN(ddtestResultFile, FILE=testResultFilename)
      else
! phase that does nothing         
         ddphase = phase
      endif
      end subroutine

      subroutine traceactivescalarreal(varname, var, vard)
      real var, vard, ddvar, dd, diff
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         WRITE(ddepsValuesFile, *) var
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',&
&                         ' read ', ddvarname, ' expecting ', varname
            endif
            READ(ddepsValuesFile, *) ddvar
            dd = (ddvar-var)/ddeps
            diff = (abs(vard-dd)*100.0)/  &
&                         max(abs(vard),abs(dd))
            if (diff.gt.10.0) then
               diffstr = 'DIFFERENCE!!'
            else
               diffstr = '            '
            endif
            if ((abs(vard).gt.epszero).or.(abs(dd).gt.epszero)) then
               write (*, 300) varname, vard, dd, diffstr
            endif
         else
            write (*, 350) varname, vard
         endif
      endif
 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
 350  format("   ", a,":",e8.2)
      end subroutine

      subroutine traceactivearray1real(varname, var, vard, lower1, upper1)
      integer lower1,upper1
      real var(lower1:upper1)
      real vard(lower1:upper1)
      real ddvar, dd, diff
      real valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         do i1=lower1,upper1
            WRITE(ddepsValuesFile, *) var(i1)
         enddo
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.0) dd = 0.0
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
         endif
         notprintedheader=.true.
         ibuf = 1
         do i1=lower1,upper1
            if (ddphase.eq.2) then
               READ(ddepsValuesFile, *) ddvar
               dd = (ddvar-var(i1))/ddeps
            endif
            if ((abs(vard(i1)).gt.epszero).or.(abs(dd).gt.epszero)) then
               valbuf(ibuf) = vard(i1)
               ddbuf(ibuf) = dd
               indexbuf1(ibuf) = i1
               if (ddphase.eq.2) then
                  diff = (abs(vard(i1)-dd)*100.0)/  &
&                      max(abs(vard(i1)),abs(dd))
                  if (diff.gt.10.0) then
                     diffbuf(ibuf) = '  DIFFERENCE!!'
                  else
                     diffbuf(ibuf) = '              '
                  endif
               endif
               ibuf = ibuf+1
            endif
            if (ibuf.gt.10) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)  &
&                   ((indexbuf1(j),valbuf(j)), j=1,10)
               if (ddphase.eq.2) then
                  write (*, 250) ((ddbuf(j)), j=1,10)
                  write (*, 270) ((diffbuf(j)), j=1,10)
               endif
               ibuf = 1
            endif
         end do
         if (ibuf.gt.1) then
            if (notprintedheader) then
               write(*,300) varname
               notprintedheader=.false.
            endif
            write (*, 200)  &
&                ((indexbuf1(j),valbuf(j)), j=1,ibuf-1)
            if (ddphase.eq.2) then
               write (*, 250) ((ddbuf(j)), j=1,ibuf-1)
               write (*, 270) ((diffbuf(j)), j=1,ibuf-1)
            endif
         endif
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end subroutine

      subroutine traceactivearray2real(varname, var, vard, lower1, upper1, lower2, upper2)
      integer lower1,upper1,lower2,upper2
      real var(lower1:upper1, lower2:upper2)
      real vard(lower1:upper1, lower2:upper2)
      real ddvar, dd, diff
      real valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10),indexbuf2(10)
      character*19 diffbuf(10)
      integer i1,i2,j
      integer ibuf
      logical notprintedheader

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         do i2=lower2,upper2
            do i1=lower1,upper1
               WRITE(ddepsValuesFile, *) var(i1,i2)
            enddo
         enddo
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.0) dd = 0.0
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
         endif
         notprintedheader=.true.
         ibuf = 1
         do i2=lower2,upper2
            do i1=lower1,upper1
               if (ddphase.eq.2) then
                  READ(ddepsValuesFile, *) ddvar
                  dd = (ddvar-var(i1,i2))/ddeps
               endif
            if ((abs(vard(i1,i2)).gt.epszero).or.  &
&               (abs(dd).gt.epszero)) then
                  valbuf(ibuf) = vard(i1,i2)
                  ddbuf(ibuf) = dd
                  indexbuf1(ibuf) = i1
                  indexbuf2(ibuf) = i2
                  if (ddphase.eq.2) then
                     diff = (abs(vard(i1,i2)-dd)*100.0)/  &
&                         max(abs(vard(i1,i2)),abs(dd))
                     if (diff.gt.10.0) then
                        diffbuf(ibuf) = '       DIFFERENCE!!'
                     else
                        diffbuf(ibuf) = '                   '
                     endif
                  endif
                  ibuf = ibuf+1
               endif
               if (ibuf.gt.10) then
                  if (notprintedheader) then
                     write(*,300) varname
                     notprintedheader=.false.
                  endif
                  write (*, 200)  &
&                      ((indexbuf1(j),indexbuf2(j),valbuf(j)), j=1,10)
                  if (ddphase.eq.2) then
                     write (*, 250) ((ddbuf(j)), j=1,10)
                     write (*, 270) ((diffbuf(j)), j=1,10)
                  endif
                  ibuf = 1
               endif
            end do
         end do
         if (ibuf.gt.1) then
            if (notprintedheader) then
               write(*,300) varname
               notprintedheader=.false.
            endif
            write (*, 200)  &
&                ((indexbuf1(j),indexbuf2(j),valbuf(j)), j=1,ibuf-1)
            if (ddphase.eq.2) then
               write (*, 250) ((ddbuf(j)), j=1,ibuf-1)
               write (*, 270) ((diffbuf(j)), j=1,ibuf-1)
            endif
         endif
      endif
 200  format('          ', 10(i4,',',i4,'->',e8.2))
 250  format('      (dd:)', 10('         (',e8.2,')'))
 270  format('          ', 10(a19))
 300  format('   ', a, ':')
      end subroutine

      subroutine traceactivescalarreal8(varname, var, vard)
      real*8 var, vard, ddvar, dd, diff, varrw, vardrw
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname

      if (ddphase.eq.0) then
         write (*, 350) varname, vard
      else if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         WRITE(ddepsValuesFile, *) var
      else if (ddphase.eq.2) then
         READ(ddepsValuesFile, '(a)') ddvarname
         if (ddvarname.ne.varname) then
            write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
         endif
         READ(ddepsValuesFile, *) ddvar
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         WRITE(rwTmpFile, *) var
         CLOSE(rwTmpFile)
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         READ(rwTmpFile, *) varrw
         CLOSE(rwTmpFile)
         WRITE(ddeps2ValuesFile, '(a)') varname
         WRITE(ddeps2ValuesFile, *) varrw
         dd = (ddvar-varrw)/ddeps
         diff = (abs(vard-dd)*100.0)/  &
&                         max(abs(vard),abs(dd))
!            if (diff.gt.10.0) then
         if (diff.gt.0.1) then
            diffstr = 'DIFFERENCE!!'
         else
            diffstr = '            '
         endif
!            if ((abs(vard).gt.epszero).or.(abs(dd).gt.epszero)) then
!               write (*,*) ddvar, var, diff, varrw, (ddvar-varrw)/ddeps
         write (ddtestResultFile, 300) varname, vard, dd, diffstr
!            endif
      else if (ddphase.eq.3) then
         READ(ddepsValuesFile, '(a)') ddvarname
         if (ddvarname.ne.varname) then
            write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
         endif
         READ(ddepsValuesFile, *) ddvar
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         WRITE(rwTmpFile, *) var
         CLOSE(rwTmpFile)
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         READ(rwTmpFile, *) varrw
         CLOSE(rwTmpFile)
         WRITE(ddeps2ValuesFile, '(a)') varname
         WRITE(ddeps2ValuesFile, *) varrw
         dd = (ddvar-varrw)/(2.d0*ddeps)
         WRITE(ddddValuesFile, '(a)') varname
         WRITE(ddddValuesFile, *) dd
      else if (ddphase.eq.4) then
         READ(ddddValuesFile, '(a)') ddvarname
         if (ddvarname.ne.varname) then
            write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
         endif
         READ(ddddValuesFile, *) dd
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         WRITE(rwTmpFile, *) vard
         CLOSE(rwTmpFile)
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         READ(rwTmpFile, *) vardrw
         CLOSE(rwTmpFile)
         diff = (abs(vardrw-dd)*100.0)/  &
&                         max(abs(vardrw),abs(dd))
!            if (diff.gt.10.0) then
         if (diff.gt.0.1) then
            diffstr = 'DIFFERENCE!!'
         else
            diffstr = '            '
         endif
         write (ddtestResultFile, 300) varname, vardrw, dd, diffstr
      endif
! 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
! 350  format("   ", a,":",e8.2)
 300  format("   ", a,":",e17.10," (dd:",e17.10,")   ",a12)
 350  format("   ", a,":",e17.10)
      end subroutine

      subroutine traceactivescalarrealWP(varname, var, vard)
      real(wwp) :: var, vard, ddvar, dd, diff, varrw, vardrw
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname

      if (ddphase.eq.0) then
         write (*, 350) varname, vard
      else if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         WRITE(ddepsValuesFile, *) var
      else if (ddphase.eq.2) then
         READ(ddepsValuesFile, '(a)') ddvarname
         if (ddvarname.ne.varname) then
            write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
         endif
         READ(ddepsValuesFile, *) ddvar
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         WRITE(rwTmpFile, *) var
         CLOSE(rwTmpFile)
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         READ(rwTmpFile, *) varrw
         CLOSE(rwTmpFile)
         WRITE(ddeps2ValuesFile, '(a)') varname
         WRITE(ddeps2ValuesFile, *) varrw
         dd = (ddvar-varrw)/ddeps
         diff = (abs(vard-dd)*100.0_wwp)/  &
&                         max(abs(vard),abs(dd))
!            if (diff.gt.10.0) then
         if (diff.gt.0.1) then
            write (ddtestResultFile, 310) varname, vard, dd, NINT(diff)
         else
            write (ddtestResultFile, 320) varname, vard, dd
         endif
      else if (ddphase.eq.3) then
         READ(ddepsValuesFile, '(a)') ddvarname
         if (ddvarname.ne.varname) then
            write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
         endif
         READ(ddepsValuesFile, *) ddvar
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         WRITE(rwTmpFile, *) var
         CLOSE(rwTmpFile)
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         READ(rwTmpFile, *) varrw
         CLOSE(rwTmpFile)
         WRITE(ddeps2ValuesFile, '(a)') varname
         WRITE(ddeps2ValuesFile, *) varrw
         dd = (ddvar-varrw)/(2.d0*ddeps)
         WRITE(ddddValuesFile, '(a)') varname
         WRITE(ddddValuesFile, *) dd
      else if (ddphase.eq.4) then
         READ(ddddValuesFile, '(a)') ddvarname
         if (ddvarname.ne.varname) then
            write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
         endif
         READ(ddddValuesFile, *) dd
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         WRITE(rwTmpFile, *) vard
         CLOSE(rwTmpFile)
         OPEN(rwTmpFile, FILE=rwTmpFilename)
         READ(rwTmpFile, *) vardrw
         CLOSE(rwTmpFile)
         diff = (abs(vardrw-dd)*100.0_wwp)/  &
&                         max(abs(vardrw),abs(dd))
!            if (diff.gt.10.0) then
         if (diff.gt.0.1) then
            write (ddtestResultFile, 310) varname, vard, dd, NINT(diff)
         else
            write (ddtestResultFile, 320) varname, vard, dd
         endif
      endif
! 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
! 350  format("   ", a,":",e8.2)
 300  format("   ", a,":",e20.13," (dd:",e20.13,")   ",a12)
 310  format("   ", a,":",e20.13," (dd:",e20.13,")   DIFFERENCE!! ",i4,"%")
 320  format("   ", a,":",e20.13," (dd:",e20.13,")   ")
 350  format("   ", a,":",e20.13)
      end subroutine

      subroutine traceactivearray1real8(varname, var, vard, lower1, upper1)
      integer lower1,upper1
      real*8 var(lower1:upper1)
      real*8 vard(lower1:upper1)
      real*8 ddvar, dd, diff
      real*8 valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         do i1=lower1,upper1
            WRITE(ddepsValuesFile, *) var(i1)
         enddo
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.0) dd = 0.0
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
         endif
         notprintedheader=.true.
         ibuf = 1
         do i1=lower1,upper1
            if (ddphase.eq.2) then
               READ(ddepsValuesFile, *) ddvar
               dd = (ddvar-var(i1))/ddeps
            endif
            if ((abs(vard(i1)).gt.epszero).or.(abs(dd).gt.epszero)) then
               valbuf(ibuf) = vard(i1)
               ddbuf(ibuf) = dd
               indexbuf1(ibuf) = i1
               if (ddphase.eq.2) then
                  diff = (abs(vard(i1)-dd)*100.0)/  &
&                      max(abs(vard(i1)),abs(dd))
                  if (diff.gt.10.0) then
                     diffbuf(ibuf) = '  DIFFERENCE!!'
                  else
                     diffbuf(ibuf) = '              '
                  endif
               endif
               ibuf = ibuf+1
            endif
            if (ibuf.gt.10) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)  &
&                   ((indexbuf1(j),valbuf(j)), j=1,10)
               if (ddphase.eq.2) then
                  write (*, 250) ((ddbuf(j)), j=1,10)
                  write (*, 270) ((diffbuf(j)), j=1,10)
               endif
               ibuf = 1
            endif
         end do
         if (ibuf.gt.1) then
            if (notprintedheader) then
               write(*,300) varname
               notprintedheader=.false.
            endif
            write (*, 200)  &
&                ((indexbuf1(j),valbuf(j)), j=1,ibuf-1)
            if (ddphase.eq.2) then
               write (*, 250) ((ddbuf(j)), j=1,ibuf-1)
               write (*, 270) ((diffbuf(j)), j=1,ibuf-1)
            endif
         endif
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end subroutine

      subroutine traceactivearray2real8  &
&          (varname, var, vard, lower1, upper1, lower2, upper2)
      integer lower1,upper1,lower2,upper2
      real*8 var(lower1:upper1, lower2:upper2)
      real*8 vard(lower1:upper1, lower2:upper2)
      real*8 ddvar, dd, diff
      real*8 valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10),indexbuf2(10)
      character*19 diffbuf(10)
      integer i1,i2,j
      integer ibuf
      logical notprintedheader

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         do i2=lower2,upper2
            do i1=lower1,upper1
               WRITE(ddepsValuesFile, *) var(i1,i2)
            enddo
         enddo
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.0) dd = 0.0
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
         endif
         notprintedheader=.true.
         ibuf = 1
         do i2=lower2,upper2
            do i1=lower1,upper1
               if (ddphase.eq.2) then
                  READ(ddepsValuesFile, *) ddvar
                  dd = (ddvar-var(i1,i2))/ddeps
               endif
            if ((abs(vard(i1,i2)).gt.epszero).or.  &
&               (abs(dd).gt.epszero)) then
                  valbuf(ibuf) = vard(i1,i2)
                  ddbuf(ibuf) = dd
                  indexbuf1(ibuf) = i1
                  indexbuf2(ibuf) = i2
                  if (ddphase.eq.2) then
!                   diff = (abs(vard(i1,i2)-dd)*100.0)/max(abs(vard(i1,i2)),abs(dd))
                    diff = dabs(dabs(vard(i1,i2))-dabs(dd))
                     if (diff.gt.10.e-5) then
                        diffbuf(ibuf) = '       DIFFERENCE!!'
                     else
                        diffbuf(ibuf) = '                   '
                     endif
                  endif
                  ibuf = ibuf+1
               endif
               if (ibuf.gt.10) then
                  if (notprintedheader) then
                     write(*,300) varname
                     notprintedheader=.false.
                  endif
                  write (*, 200)  &
&                      ((indexbuf1(j),indexbuf2(j),valbuf(j)), j=1,10)
                  if (ddphase.eq.2) then
                     write (*, 250) ((ddbuf(j)), j=1,10)
                     write (*, 270) ((diffbuf(j)), j=1,10)
                  endif
                  ibuf = 1
               endif
            end do
         end do
         if (ibuf.gt.1) then
            if (notprintedheader) then
               write(*,300) varname
               notprintedheader=.false.
            endif
            write (*, 200)  &
&                ((indexbuf1(j),indexbuf2(j),valbuf(j)), j=1,ibuf-1)
            if (ddphase.eq.2) then
               write (*, 250) ((ddbuf(j)), j=1,ibuf-1)
               write (*, 270) ((diffbuf(j)), j=1,ibuf-1)
            endif
         endif
      endif
 200  format('          ', 10(i4,',',i4,'->',e8.2))
 250  format('      (dd:)', 10('         (',e8.2,')'))
 270  format('          ', 10(a19))
 300  format('   ', a, ':')
      end subroutine

      subroutine traceactivescalardouble(varname, var, vard)
      double precision var, vard, ddvar, dd, diff
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         WRITE(ddepsValuesFile, *) var
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
            READ(ddepsValuesFile, *) ddvar
            dd = (ddvar-var)/ddeps
            diff = (abs(vard-dd)*100.0)/  &
&                         max(abs(vard),abs(dd))
            if (diff.gt.10.0) then
               diffstr = 'DIFFERENCE!!'
            else
               diffstr = '            '
            endif
            if ((abs(vard).gt.epszero).or.(abs(dd).gt.epszero)) then
               write (*, 300) varname, vard, dd, diffstr
            endif
         else
            write (*, 350) varname, vard
         endif
      endif
 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
 350  format("   ", a,":",e8.2)
      end subroutine

      subroutine traceactivearray1double(varname, var, vard, lower1, upper1)
      integer lower1,upper1
      double precision var(lower1:upper1)
      double precision vard(lower1:upper1)
      double precision ddvar, dd, diff
      double precision valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         do i1=lower1,upper1
            WRITE(ddepsValuesFile, *) var(i1)
         enddo
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.0) dd = 0.0
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
         endif
         notprintedheader=.true.
         ibuf = 1
         do i1=lower1,upper1
            if (ddphase.eq.2) then
               READ(ddepsValuesFile, *) ddvar
               dd = (ddvar-var(i1))/ddeps
            endif
            if ((abs(vard(i1)).gt.epszero).or.(abs(dd).gt.epszero)) then
               valbuf(ibuf) = vard(i1)
               ddbuf(ibuf) = dd
               indexbuf1(ibuf) = i1
               if (ddphase.eq.2) then
                  diff = (abs(vard(i1)-dd)*100.0)/  &
&                      max(abs(vard(i1)),abs(dd))
                  if (diff.gt.10.0) then
                     diffbuf(ibuf) = '  DIFFERENCE!!'
                  else
                     diffbuf(ibuf) = '              '
                  endif
               endif
               ibuf = ibuf+1
            endif
            if (ibuf.gt.10) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)  &
&                   ((indexbuf1(j),valbuf(j)), j=1,10)
               if (ddphase.eq.2) then
                  write (*, 250) ((ddbuf(j)), j=1,10)
                  write (*, 270) ((diffbuf(j)), j=1,10)
               endif
               ibuf = 1
            endif
         end do
         if (ibuf.gt.1) then
            if (notprintedheader) then
               write(*,300) varname
               notprintedheader=.false.
            endif
            write (*, 200)  &
&                ((indexbuf1(j),valbuf(j)), j=1,ibuf-1)
            if (ddphase.eq.2) then
               write (*, 250) ((ddbuf(j)), j=1,ibuf-1)
               write (*, 270) ((diffbuf(j)), j=1,ibuf-1)
            endif
         endif
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end subroutine

      subroutine traceactivearray2double  &
&          (varname, var, vard, lower1, upper1, lower2, upper2)
      integer lower1,upper1,lower2,upper2
      double precision var(lower1:upper1, lower2:upper2)
      double precision vard(lower1:upper1, lower2:upper2)
      double precision ddvar, dd, diff
      double precision valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10),indexbuf2(10)
      character*19 diffbuf(10)
      integer i1,i2,j
      integer ibuf
      logical notprintedheader

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         do i2=lower2,upper2
            do i1=lower1,upper1
               WRITE(ddepsValuesFile, *) var(i1,i2)
            enddo
         enddo
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.0) dd = 0.0
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
         endif
         notprintedheader=.true.
         ibuf = 1
         do i2=lower2,upper2
            do i1=lower1,upper1
               if (ddphase.eq.2) then
                  READ(ddepsValuesFile, *) ddvar
                  dd = (ddvar-var(i1,i2))/ddeps
               endif
            if ((abs(vard(i1,i2)).gt.epszero).or.  &
&               (abs(dd).gt.epszero)) then
                  valbuf(ibuf) = vard(i1,i2)
                  ddbuf(ibuf) = dd
                  indexbuf1(ibuf) = i1
                  indexbuf2(ibuf) = i2
                  if (ddphase.eq.2) then
                     diff = (abs(vard(i1,i2)-dd)*100.0)/  &
&                         max(abs(vard(i1,i2)),abs(dd))
                     if (diff.gt.10.0) then
                        diffbuf(ibuf) = '       DIFFERENCE!!'
                     else
                        diffbuf(ibuf) = '                   '
                     endif
                  endif
                  ibuf = ibuf+1
               endif
               if (ibuf.gt.10) then
                  if (notprintedheader) then
                     write(*,300) varname
                     notprintedheader=.false.
                  endif
                  write (*, 200)  &
&                      ((indexbuf1(j),indexbuf2(j),valbuf(j)), j=1,10)
                  if (ddphase.eq.2) then
                     write (*, 250) ((ddbuf(j)), j=1,10)
                     write (*, 270) ((diffbuf(j)), j=1,10)
                  endif
                  ibuf = 1
               endif
            end do
         end do
         if (ibuf.gt.1) then
            if (notprintedheader) then
               write(*,300) varname
               notprintedheader=.false.
            endif
            write (*, 200)  &
&                ((indexbuf1(j),indexbuf2(j),valbuf(j)), j=1,ibuf-1)
            if (ddphase.eq.2) then
               write (*, 250) ((ddbuf(j)), j=1,ibuf-1)
               write (*, 270) ((diffbuf(j)), j=1,ibuf-1)
            endif
         endif
      endif
 200  format('          ', 10(i4,',',i4,'->',e8.2))
 250  format('      (dd:)', 10('         (',e8.2,')'))
 270  format('          ', 10(a19))
 300  format('   ', a, ':')
      end subroutine

      subroutine traceactivearray3double  &
&          (varname, var, vard, lower1, upper1, lower2, upper2,  &
&           lower3, upper3)
      integer lower1,upper1,lower2,upper2,lower3,upper3
      double precision var(lower1:upper1, lower2:upper2, lower3:upper3)
      double precision vard(lower1:upper1, lower2:upper2, lower3:upper3)
      double precision ddvar, dd, diff
      double precision valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10),indexbuf2(10),indexbuf3(10)
      character*19 diffbuf(10)
      integer i1,i2,i3,j
      integer ibuf
      logical notprintedheader

      if (ddphase.eq.0) return
      if (ddphase.eq.1) then
         WRITE(ddepsValuesFile, '(a)') varname
         do i3=lower3,upper3
            do i2=lower2,upper2
               do i1=lower1,upper1
                  WRITE(ddepsValuesFile, *) var(i1,i2,i3)
               enddo
            enddo
         enddo
      else if ((ddphase.eq.2).or.(ddphase.eq.0)) then
         if (ddphase.eq.0) dd = 0.0
         if (ddphase.eq.2) then
            READ(ddepsValuesFile, '(a)') ddvarname
            if (ddvarname.ne.varname) then
               write(*,*) 'ERROR: mismatch in DD program control !!!',  &
&                         ' read ', ddvarname, ' expecting ', varname
            endif
         endif
         notprintedheader=.true.
         ibuf = 1
         do i3=lower3,upper3
          do i2=lower2,upper2
           do i1=lower1,upper1
               if (ddphase.eq.2) then
                  READ(ddepsValuesFile, *) ddvar
                  dd = (ddvar-var(i1,i2,i3))/ddeps
               endif
            if ((abs(vard(i1,i2,i3)).gt.epszero).or.  &
&               (abs(dd).gt.epszero)) then
                  valbuf(ibuf) = vard(i1,i2,i3)
                  ddbuf(ibuf) = dd
                  indexbuf1(ibuf) = i1
                  indexbuf2(ibuf) = i2
                  indexbuf3(ibuf) = i3
                  if (ddphase.eq.2) then
                     diff = (abs(vard(i1,i2,i3)-dd)*100.0)/  &
&                         max(abs(vard(i1,i2,i3)),abs(dd))
                     if (diff.gt.10.0) then
                        diffbuf(ibuf) = '       DIFFERENCE!!'
                     else
                        diffbuf(ibuf) = '                   '
                     endif
                  endif
                  ibuf = ibuf+1
               endif
               if (ibuf.gt.10) then
                  if (notprintedheader) then
                     write(*,300) varname
                     notprintedheader=.false.
                  endif
                  write (*, 200)  &
&          ((indexbuf1(j),indexbuf2(j),indexbuf3(j),valbuf(j)), j=1,10)
                  if (ddphase.eq.2) then
                     write (*, 250) ((ddbuf(j)), j=1,10)
                     write (*, 270) ((diffbuf(j)), j=1,10)
                  endif
                  ibuf = 1
               endif
            end do
           end do
         end do
         if (ibuf.gt.1) then
            if (notprintedheader) then
               write(*,300) varname
               notprintedheader=.false.
            endif
            write (*, 200)  &
&       ((indexbuf1(j),indexbuf2(j),indexbuf3(j),valbuf(j)), j=1,ibuf-1)
            if (ddphase.eq.2) then
               write (*, 250) ((ddbuf(j)), j=1,ibuf-1)
               write (*, 270) ((diffbuf(j)), j=1,ibuf-1)
            endif
         endif
      endif
 200  format('          ', 10(i4,',',i4,',',i4,'->',e8.2))
 250  format('      (dd:)', 10('              (',e8.2,')'))
 270  format('          ', 10(a24))
 300  format('   ', a, ':')
      end subroutine

END MODULE TAPENADE_DD_TEST
