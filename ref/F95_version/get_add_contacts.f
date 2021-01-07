      subroutine get_add_contacts(pdbfile,nres,iseq,seq,ch,icmap,
     $           kcont,kk_map,rr,arg,outputfile,kk_out,iresn,nssb2)
     
      implicit none
      
      character (len=*) pdbfile
      integer kk_map,rr,kk_out
      integer nres,test1,test2
      integer icmap(nres,nres)
      integer iseq(nres),iresn(nres)
      character*3 seq(nres)
      character*1 ch(nres)
      real*8 xn(nres),yn(nres),zn(nres)
      integer kcont
      integer ksb(nssb2,2)
      character*80 line,outputfile,arg
      character*1 chain1,chain2
      integer cys1,cys2,nssb,nssb2
      integer i1,i2,k,j
      
      nssb=0
      test1=0
      test2=0
      open(rr,file=trim(arg),status='old')
15    read(rr,'(a)',end=20,err=20) line
      if(line(1:9).ne.'ADCONTACT') goto 15
      read(rr,*) chain1,cys1,chain2,cys2      
      nssb=nssb+1
      do j=1,nres
       if (ch(j).eq.chain1.and.iresn(j).eq.cys1) then
        ksb(nssb,1)=j
        test1=1
       endif
       if (ch(j).eq.chain2.and.iresn(j).eq.cys2) then
        ksb(nssb,2)=j
        test2=1
       endif
      enddo
      if ((test1.eq.0).or.(test2.eq.0)) then
      open(kk_out,file=trim(outputfile),
     $            form='formatted',status='unknown')
                 
      write(kk_out,'(a37)')'ADDITIONAL BONDS DEFINED INCORRECTLY!'
      close(kk_out)
      stop
      endif
      test1=0
      test2=0
      goto 15
20    close(rr)
      if (nssb.gt.0) then
       write(kk_map,'(/,a)') '# ADDITIONAL CONTACTS'
       do k=1,nssb
        i1=ksb(k,1)
        i2=ksb(k,2)
!40      format(a8,i4,2x,a3,2x,i4,1x,a1,6x,a3,2x,i4,1x,a1)
40        format(a1,1x,i6,2x,i4,2x,a3,1x,a1,1x,i3, 5x, 
     $         i4,2x,a3,1x,a1,1x,i3)
     
        if(icmap(i1,i2).eq.0) then
         write(kk_map,40)'N',k,
     $   i1,seq(i1),ch(i1),iresn(i1),i2,seq(i2),ch(i2),iresn(i2)
         icmap(i1,i2)=1
         kcont=kcont+1
        endif
        
       enddo
      endif
      return
      end
