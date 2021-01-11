      subroutine get_add_bonds_numb(pdbfile,nres,iseq,seq,ch,icmap,
     $           kcont,kk_map,rr,arg,outputfile,kk_out,iresn,
     $           nssb1,nssb2,nssb3)
     
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
      integer ksb(kcont,2)
      character*80 line,outputfile,arg
      character*1 chain1,chain2
      integer nssb1,nssb2,nssb3,cys1,cys2
      integer i1,i2,k,j
      
      nssb1=0
      nssb2=0
      nssb3=0
      
      open(rr,file=trim(arg),status='old')
15    read(rr,'(a)',end=20,err=20) line
      if(line(1:6).ne.'ADBOND') goto 15 
      nssb1=nssb1+1
      goto 15
20    close(rr)

      open(rr,file=trim(arg),status='old')
30    read(rr,'(a)',end=40,err=40) line
      if(line(1:9).ne.'ADCONTACT') goto 30      
      nssb2=nssb2+1
      goto 30
40    close(rr) 

      open(rr,file=trim(arg),status='old')
50    read(rr,'(a)',end=60,err=60) line
      if(line(1:13).ne.'REMOVECONTACT') goto 50      
      nssb3=nssb3+1
      goto 50
60    close(rr)       
      
      return
      end
