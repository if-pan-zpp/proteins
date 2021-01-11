      subroutine get_knot_limits(nres,iseq,seq,ch,icmap,
     & kcont,rr,arg,kb_t,ke_t,knumbers,kk_out)
     
      implicit none
      
      integer nres,kcont,kk_out
      integer iseq(nres)  
      integer icmap(nres,nres) 
      character*3 seq(nres)
      character*1 ch(nres)            
      integer kk_map,rr,knumbers
      integer ksb(kcont,2)
      integer kb_t(knumbers,2)
      integer ke_t(knumbers,2)
      character*80 line,outputfile,arg
      character*1 chain1,chain2
      integer cys1,cys2,nssb,test1,test2
      integer i1,i2,k,j
      
      nssb=0
      test1=0
      test2=0
      open(rr,file=trim(arg),status='old')
15    read(rr,'(a)',end=20,err=20) line
      if(line(1:8).ne.'KKLIMITS') goto 15
      nssb=nssb+1
      read(rr,*) chain1,cys1,chain2,cys2      
      do j=1,nres
       if (ch(j).eq.chain1.and.iseq(j).eq.cys1) then
        ksb(nssb,1)=j
        kb_t(nssb,1)=j
        test1=1
       endif   
       if (ch(j).eq.chain2.and.iseq(j).eq.cys2) then
        ksb(nssb,2)=j
        ke_t(nssb,1)=j
        test2=1
       endif   
       enddo          
      if ((test1.eq.0).or.(test2.eq.0)) then
!      open(kk_out,file=trim(outputfile),
!     $            form='formatted',status='unknown')
                                  
      write(kk_out,'(a32)')'KNOT LIMITS DEFINED INCORRECTLY!'
!      close(kk_out)
      stop
      endif
      test1=0
      test2=0     
      
      goto 15
20    close(rr)

      return
      end
