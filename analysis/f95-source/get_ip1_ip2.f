      subroutine get_ip1_ip2(nres,ch,snt,sct,nch,ichain,ip1,ip2,kk,
     $ snch,scch)
      implicit none
      integer nres,nch
      integer kk,snch_i,scch_i
      character*1 ch(nres)
      character*1 snt,sct,snch1,scch1
      integer ichain(nch,2)
      integer ip1,ip2,i1,i2,i
      character*6 snch,scch
      
      snch1=snch
      scch1=scch
      
      ip1=0
      ip2=0
      do i=1,nch
      if (trim(snch).eq.'N') i1=ichain(i,1)
      if (trim(snch).eq.'C') i1=ichain(i,2)
      if ((ichar(snch1).gt.48).and.(ichar(snch1).lt.58)) then
      read(snch,*)snch_i
      i1=ichain(i,1)+snch_i-1
      endif      
      if (trim(scch).eq.'N') i2=ichain(i,1)
      if (trim(scch).eq.'C') i2=ichain(i,2)
      if ((ichar(scch1).gt.48).and.(ichar(scch1).lt.58)) then
      read(scch,*)scch_i
      i2=ichain(i,1)+scch_i-1
      endif
      
      if ((trim(snch).ne.'N').and.
     $     (trim(snch).ne.'C').and.
     $     ((ichar(snch1).lt.49).or.(ichar(snch1).gt.57)).and.
     $     (ichar(snch1).ne.67).and.
     $     (ichar(snch1).ne.78)) then
       write(kk,'(/,a)') 'CHECK THE FIRST TERMINUS DEFINITION!'
       stop
      endif
       
      if ((trim(scch).ne.'N').and.
     $     (trim(scch).ne.'C').and.
     $     ((ichar(scch1).lt.49).or.(ichar(scch1).gt.57)).and.
     $     (ichar(scch1).ne.67).and.
     $     (ichar(scch1).ne.78)) then
       write(kk,'(/,a)') 'CHECK THE SECOND TERMINUS DEFINITION!'
       stop
      endif 
              
       if (ch(i1).eq.snt) then
        ip1=i1
       endif
       if (ch(i2).eq.sct) then
        ip2=i2
       endif
      enddo
      if (ip1.eq.0 .or. ip2.eq.0) then
       write(kk,'(/,a)') 'PROBLEM WITH IDENTIFYING AA CHEINS!'
       write(kk,'(a,/)') 'CHECK THE PDB AND INPUT FILES' 
       stop
      endif
      return
      end