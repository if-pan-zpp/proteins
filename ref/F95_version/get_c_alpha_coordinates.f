      subroutine get_c_alpha_coordinates(pdbfile,nres,x,y,z,
     $           seq,iseq,ch,nch)
C THIS SURROUTINE READS IN C-ALPHA COORDINTES FROM THE PDB FILE
C HERE X BECOME XN AND SO ON
      implicit none
      character (len=*) pdbfile ! PDB FILE
      integer nres		! NUMBER OF RESIDUES
      real*8 x(nres)
      real*8 y(nres)
      real*8 z(nres)
      integer iseq(nres)
      character*3 seq(nres)
      character*1 ch(nres)
      integer ires,nch
      character*80 line
C OPEN PDB FILE

      open(20, file=trim(pdbfile),
     $         form='formatted', status='old')
      ires=0
 11   read(20, '(a)', end = 12, err = 12) line
      if (line(1:4) .eq. 'ATOM' ) then
        if (line (14:16) .eq. 'CA ' ) then
          ires=ires+1
          read( line(18:20), '(a3)' ) seq(ires)
          read( line(22:22), '(a1)' ) ch(ires)
          read( line(24:26), '(i3)' ) iseq(ires)
          read( line(31:54), '(3f8.3)' ) x(ires),y(ires),z(ires)
        endif
      endif
      goto 11
 12   close (20)
C CLOSE PDB FILE
      nch=1
      do ires=1,nres-1
       if (ch(ires).ne.ch(ires+1)) then
        nch=nch+1
       endif
      enddo
      return
      end
