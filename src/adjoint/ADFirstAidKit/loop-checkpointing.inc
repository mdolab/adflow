c Benjamin Dauvergne November 2005
c
c This file implements functions to efficiently checkpoint loops
c You must change your big loop in a recursive function as
c prototype at the end of this file and change the choose_median
c parameter to fit your needs

      subroutine beta(s,t,res)
       integer s,t,res,i
       res = 1
       do i=1,s
        res = res * (t+i)
        res = res / i
       enddo
      end

      subroutine max_t_as_beta_under(s,l,res)
       integer s,l,res,t

       res = 0
       t = 0
       do while (res.lt.l)
        t = t + 1
        call beta(s,t,res)
       enddo
       if (t.eq.0) then
        call beta(s,t,res)
       else
        call beta(s,t-1,res)
       endif
      end
      
c Choose a good point to put a checkpoint
c ntype determine the algorithm used:
c nfrom = starting iteration (in)
c nto = ending iteration (in)
c median = where to checkpoint (out)
c ntype = see below (in)
c s = max number of checkpoints
c minl, maxl = min and max on the size of the iteration domain to revert on another
c              technic than treeverse
c ntype2 = the checkpointing technic to use when nto - nfrom "< minl" or "> maxl"
c 
c Type kind:
c 1 dichotomy    ckp = ( from + to ) / 2
c 2 dichohomy/ra ckp = ( from + to ) / 2 if to-from > maxl else ckp = to-1
c 3 RA ckp = to-1
c 4 SA ckp = from
c 5 Treeverse ckp = from+binomial(s+t-1,t-1) where t = min { t / binomial(s+t,t)  > to-from } - 1
c 6 2-level ckp = if to-from > maxl then to - maxl else from
      SUBROUTINE choose_median(nfrom, nto, median, ntype, s, minl 
     +,maxl,ntype2)
      integer nfrom,nto,ntype,s,l
      integer maxl,minl,ntype2

       if (ntype.eq.1) then
c---------Dichotomy
          median = (nfrom+nto)/2
       else if (ntype.eq.2) then
c---------Dichotomy over 20 steps, recompute all under
          if ((nto-nfrom).gt.maxl) then
            median = (nfrom+nto)/2
          else
            median = nto-1
          endif
c---------Recompute all
       else if (ntype.eq.3) then 
         median = nto-1
       else if (ntype.eq.4) then
c---------Save all
         median = nfrom
       else if (ntype.eq.5) then
c---------Treeverse over 20 steps, recompute all under
         if ((nto-nfrom.gt.minl).and.
     +(nto-nfrom.lt.maxl)) then
           call max_t_as_beta_under(s, nto-nfrom, l)
           median = nfrom+l-1
         else
c Remove the limits to protect from infinite recursivity
           call choose_median(nfrom, nto, median, ntype2, s, 0, 2000000, 
     + ntype2)
         endif
       else if (ntype.eq.6) then
c---------2-level checkpoints
         if (nto-nfrom.gt.maxl) then
           median = nto-maxl
         else
           median = nfrom
         endif
       endif
      END

c Prototype for the recursive function
c Iterate OneBigLooptIter from nfrom to nto inclusive
c      subroutine RecursiveBigLoopIters(nfrom,nto)
c
c      if (nfrom.eq.nto) then
cc $AD NOCHECKPOINT
c      call oneBigLoopIter()
c      goto 1000
c      endif
c
c      call choose_median(nfrom, nto, median,4,5,0,10000000,0)
c                                            \- ntype
c                                              \= max number of checkpoints
c                                                \= minl
c                                                  \= maxl
c                                                           \= ntype2
cc $AD NOCHECKPOINT
c      call RecursiveBigLoopIters(nfrom,median)
c      call RecursiveBigLoopIters(median+1,nto) 
c1000  END
