      subroutine get_ss_bonds(pdbfile,nres,iseq,seq,ch,icmap,
     $           kcont,kk_map)
     
      implicit none
      
      character (len=*) pdbfile
      integer kk_map
      integer nres
      integer icmap(nres,nres)
      integer iseq(nres)
      character*3 seq(nres)
      character*1 ch(nres)
      real*8 xn(nres),yn(nres),zn(nres)
      integer kcont
      integer ksb(kcont,2)
      character*80 line,outputfile
      character*1 chain1,chain2
      integer cys1,cys2,nssb
      integer i1,i2,k,j
      
      nssb=0
      open(5, file=trim(pdbfile),form='formatted',status='old')
15    read(5,'(a)',end=20,err=20) line
      if(line(1:6).ne.'SSBOND') goto 15
      read (line(16:16), '(a1)' ) chain1
      read (line(17:21), '(i5)' ) cys1
      read (line(30:30), '(a1)' ) chain2
      read (line(31:35), '(i5)' ) cys2
      nssb=nssb+1
      do j=1,nres
       if (ch(j).eq.chain1 .and. iseq(j).eq.cys1) then
        ksb(nssb,1)=j
       endif
       if (ch(j).eq.chain2 .and. iseq(j).eq.cys2) then
        ksb(nssb,2)=j
       endif
      enddo
      goto 15
20    close(5)

      if (nssb.gt.0) then
       write(kk_map,'(/,a)') '# SS BONDS'
       do k=1,nssb
        i1=ksb(k,1)
        i2=ksb(k,2)
        if(seq(i1).ne.'CYS') then
         write(kk_map,*)'RESIDUE ',ksb(k,1),' IS NOT A CYSTEINE.
     $    PLEASE CHECK!'
         stop
        endif
        if(seq(i2).ne.'CYS') then
         write(kk_map,*)'RESIDUE ',ksb(k,2),' IS NOT A CYSTEINE. 
     $    PLEASE CHECK!'
         stop
        endif
!40      format(a8,i4,2x,a3,2x,i4,1x,a1,6x,a3,2x,i4,1x,a1)
40        format(a1,1x,i6,2x,i4,2x,a3,1x,a1,1x,i3, 5x,
     $         i4,2x,a3,1x,a1,1x,i3)
        if(icmap(i1,i2).eq.1) then
         write(kk_map,40)'S',k,
     $    i1,seq(i1),ch(i1),iseq(i1),i2,seq(i2),ch(i2),iseq(i2)
         icmap(i1,i2)=3
        endif
       enddo
      endif
      return
      end
