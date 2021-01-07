      subroutine get_knot_number(pdbfile,nres,iseq,seq,ch,icmap,
     $           kcont,kk_map,rr,arg,nssb)
     
      implicit none
      
      character (len=*) pdbfile
      integer kk_map,rr
      integer nres
      integer icmap(nres,nres)
      integer iseq(nres)
      character*3 seq(nres)
      character*1 ch(nres)
      real*8 xn(nres),yn(nres),zn(nres)
      integer kcont
      integer ksb(kcont,2)
      character*80 line,outputfile,arg
      character*1 chain1,chain2
      integer cys1,cys2,nssb
      integer i1,i2,k,j
      
      nssb=0
      open(rr,file=trim(arg),status='old')
15    read(rr,'(a)',end=20,err=20) line
      if(line(1:8).ne.'KKLIMITS') goto 15      
      nssb=nssb+1
      goto 15
20    close(rr)

      return
      end
