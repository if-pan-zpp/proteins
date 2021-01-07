      subroutine get_chains(nres,ch,nch,ichain)
      implicit none
      integer nres		! NUMBER OF RESIDUES
      integer nch
      integer ichain(nch,2)
      integer ires,ich
      character*1 ch(nres)
      character*80 line
     
      ich=1
      ichain(ich,1)=1
     
      do ires=1,nres-1
       if (ch(ires).ne.ch(ires+1)) then
        ichain(ich,2)=ires
        ich=ich+1
        ichain(ich,1)=ires+1
       endif
      enddo
      ichain(ich,2)=nres
C TEST:      if (ich.ne.nch) stop
      return
      end
