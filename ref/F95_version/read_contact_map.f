      subroutine read_contact_map(mapfile,nres,icmap,kcont)
     
      implicit none
      
      character (len=*) mapfile
      character*80 line
      integer rr2,nres,kcont
      integer icmap(nres,nres)
      integer i1,i2
      
      kcont=0
      
      do i1=1,nres
       do i2=1,nres
        icmap(i1,i2)=0
       enddo
      enddo 
      
      rr2=59
      
      open(rr2,file=trim(mapfile),status='old')
15    read(rr2,'(a)',end=20,err=20) line
      if (line(1:2) .eq. 'M ' ) then
         read( line(11:14), '(i4)' ) i1
         read( line(31:34), '(i4)' ) i2
         icmap(i1,i2)=1
         kcont=kcont+1               
      endif
      if (line(1:2) .eq. 'B ' ) then
         read( line(11:14), '(i4)' ) i1
         read( line(31:34), '(i4)' ) i2
         if (icmap(i1,i2).eq.0) then
         icmap(i1,i2)=3
         kcont=kcont+1
         endif
         if (icmap(i1,i2).eq.1) then
         icmap(i1,i2)=3
         endif               
      endif
      if (line(1:2) .eq. 'S ' ) then
         read( line(11:14), '(i4)' ) i1
         read( line(31:34), '(i4)' ) i2
         icmap(i1,i2)=3               
      endif
      if (line(1:2) .eq. 'N ' ) then
         read( line(11:14), '(i4)' ) i1
         read( line(31:34), '(i4)' ) i2
         if (icmap(i1,i2).eq.0) then
         icmap(i1,i2)=1
         kcont=kcont+1
         endif                      
      endif
      if (line(1:2) .eq. 'R ' ) then
         read( line(11:14), '(i4)' ) i1
         read( line(31:34), '(i4)' ) i2
         if (icmap(i1,i2).ne.0) then
         icmap(i1,i2)=0
         kcont=kcont-1
         endif               
      endif
      goto 15
20    close(rr2)

      return
      end