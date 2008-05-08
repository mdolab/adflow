C
C $Id: ddTest.f 367 2006-11-29 00:47:02Z acmarta $
C

c Utility subroutines for validation of tangent AD derivatives
c compared to Divided Differences.

      subroutine DDREAL4(varname, var, vard)
      IMPLICIT NONE
      real*4 var, vard, ddvar, dd, diff, varwr
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname
      integer ddphase, ddfile
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         WRITE(ddfile, *) var
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         call DDWRITEANDREAD4(var, varwr)
         READ(ddfile, *) ddvar
         dd = (ddvar-varwr)/ddeps
         if ((abs(vard).gt.epszero).or.(abs(dd).gt.epszero)) then
            diff = (abs(vard-dd)*100.0)/
     +           max(abs(vard),abs(dd))
            if (diff.gt.1.0) then
               diffstr = 'DIFFERENCE!!'
            else
               diffstr = '            '
            endif
            write (*, 300) varname, vard, dd, diffstr
         endif
      endif
 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
      end

      subroutine DDREAL4ARRAY
     +     (varname, var, vard, length)
      IMPLICIT NONE
      integer length
      real*4 var(length)
      real*4 vard(length)
      real*4 ddvar, dd, diff, varwr
      real*4 valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader
      integer ddphase, ddfile
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         do i1=1,length
            WRITE(ddfile, *) var(i1)
         enddo
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         notprintedheader=.true.
         ibuf = 1
         do i1=1,length
            call DDWRITEANDREAD4(var(i1), varwr)
            READ(ddfile, *) ddvar
            dd = (ddvar-varwr)/ddeps
            if ((abs(vard(i1)).gt.epszero).or.(abs(dd).gt.epszero)) then
               valbuf(ibuf) = vard(i1)
               ddbuf(ibuf) = dd
               indexbuf1(ibuf) = i1
               diff = (abs(vard(i1)-dd)*100.0)/
     +                 max(abs(vard(i1)),abs(dd))
               if (diff.gt.1.0) then
                  diffbuf(ibuf) = '  DIFFERENCE!!'
               else
                  diffbuf(ibuf) = '              '
               endif
               ibuf = ibuf+1
            endif
            if(ibuf.gt.10.or.i1.eq.length) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)
     +              (indexbuf1(j),valbuf(j), j=1,ibuf-1)
               write (*, 250) (ddbuf(j), j=1,ibuf-1)
               write (*, 270) (diffbuf(j), j=1,ibuf-1)
               ibuf = 1
            endif
         end do
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end

      subroutine DDWRITEANDREAD4(var, varwr)
      IMPLICIT NONE
      real*4 var,varwr

      OPEN(38, FILE='ddwrfile')
      WRITE(38, *) var
      REWIND(38)
      READ(38, *) varwr
      CLOSE(38)
      end



      subroutine DDREAL8(varname, var, vard)
      IMPLICIT NONE
      real*8 var, vard, ddvar, dd, diff, varwr
      character varname*(*)
      character*12 diffstr
      character*50 ddvarname
      integer ddphase, ddfile
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         WRITE(ddfile, *) var
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         call DDWRITEANDREAD8(var, varwr)
         READ(ddfile, *) ddvar
         dd = (ddvar-varwr)/ddeps
         if ((abs(vard).gt.epszero).or.(abs(dd).gt.epszero)) then
            diff = (abs(vard-dd)*100.0)/
     +           max(abs(vard),abs(dd))
            if (diff.gt.1.0) then
               diffstr = 'DIFFERENCE!!'
            else
               diffstr = '            '
            endif
            write (*, 300) varname, vard, dd, diffstr
         endif
      endif
 300  format("   ", a,":",e8.2," (dd:",e8.2,")   ",a12)
      end

      subroutine DDREAL8ARRAY
     +     (varname, var, vard, length)
      IMPLICIT NONE
      integer length
      real*8 var(length)
      real*8 vard(length)
      real*8 ddvar, dd, diff, varwr
      real*8 valbuf(10),ddbuf(10)
      character varname*(*)
      character*50 ddvarname
      integer indexbuf1(10)
      character*14 diffbuf(10)
      integer i1,j
      integer ibuf
      logical notprintedheader
      integer ddphase, ddfile
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      if (ddphase.eq.1) then
         WRITE(ddfile, '(a)') varname
         do i1=1,length
            WRITE(ddfile, *) var(i1)
         enddo
      else if (ddphase.eq.2) then
         call DDCHECKVARNAME(varname)
         notprintedheader=.true.
         ibuf = 1
         do i1=1,length
            call DDWRITEANDREAD8(var(i1), varwr)
            READ(ddfile, *) ddvar
            dd = (ddvar-varwr)/ddeps
            if ((abs(vard(i1)).gt.epszero).or.(abs(dd).gt.epszero)) then
               valbuf(ibuf) = vard(i1)
               ddbuf(ibuf) = dd
               indexbuf1(ibuf) = i1
               diff = (abs(vard(i1)-dd)*100.0)/
     +                 max(abs(vard(i1)),abs(dd))
               if (diff.gt.1.0) then
                  diffbuf(ibuf) = '  DIFFERENCE!!'
               else
                  diffbuf(ibuf) = '              '
               endif
               ibuf = ibuf+1
            endif
            if(ibuf.gt.10.or.i1.eq.length) then
               if (notprintedheader) then
                  write(*,300) varname
                  notprintedheader=.false.
               endif
               write (*, 200)
     +              (indexbuf1(j),valbuf(j), j=1,ibuf-1)
               write (*, 250) (ddbuf(j), j=1,ibuf-1)
               write (*, 270) (diffbuf(j), j=1,ibuf-1)
               ibuf = 1
            endif
         end do
      endif
 200  format('          ', 10(i4,'->',e8.2))
 250  format('      (dd:)', 10('    (',e8.2,')'))
 270  format('          ', 10(a14))
 300  format('   ', a, ':')
      end

      subroutine DDWRITEANDREAD8(var, varwr)
      IMPLICIT NONE
      real*8 var,varwr

      OPEN(38, FILE='ddwrfile')
      WRITE(38, *) var
      REWIND(38)
      READ(38, *) varwr
      CLOSE(38)
      end


      subroutine DDCHECKVARNAME(varname)
      IMPLICIT NONE
      character varname*(*)
      character*50 ddvarname
      integer ddphase, ddfile
      real*8 ddeps, epszero
      common /DDTESTCOMMON/ ddphase, ddfile, ddeps, epszero

      READ(ddfile, '(a)') ddvarname
      if (ddvarname.ne.varname) then
         write(*,*) 'ERROR: mismatch in DD program control !!!',
     +        ' read ', ddvarname, ' expecting ', varname
      endif
      end

