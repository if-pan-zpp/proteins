      subroutine get_nz(pot_file,nz)
      implicit none
      character (len=*) pot_file
      integer nz
      character*3 a1,a2
      character*80 line
      integer i,j,k
      i=0
      j=1
      k=0
      nz=0
      open(23,file=trim(pot_file),form='formatted',status='old')
11    read(23,'(a)',end=12,err=12) line
      read(line(1:3), '(a3)') a1
      i=i+1
      k=k+1
      if (i.gt.1) then
       if (a1.ne.a2) then
        j=j+1
        if (k.gt.nz) nz=k
        k=0
       endif
      endif
      a2=a1
      goto 11
12    close(23)
      if (k.gt.nz) nz=k
      return
      end
