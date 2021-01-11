      subroutine get_eij(pot_file,nres,icmap,seq,eij)
      implicit none
      character (len=*) pot_file
      integer nres
      integer icmap(nres,nres)
      real*8 eij(nres,nres)
      character*3 seq(nres)
      character*3 aa1,aa2
      character*80 line
      real*8 epsi
      integer i,j
      
      do i=1,nres
       do j=1,nres
        eij(i,j)=0.d0
       enddo
      enddo
      open(20,file=trim(pot_file),form='formatted', status='old')
 11   read(20,'(a)',end=12,err=12) line
      read(line(1:3),'(a3)') aa1
      read(line(5:7),'(a3)') aa2
      read(line(9:13),'(f5.3)') epsi
      do i=1,nres
       do j=1,nres
        if (icmap(i,j).eq.1) then
         if (aa1.eq.seq(i) .and. aa2.eq.seq(j) .or.
     $      aa1.eq.seq(j) .and. aa2.eq.seq(i)) then
          eij(i,j)=epsi
         endif
        endif
       enddo
      enddo
      goto 11
 12   close(20)
      return
      end
      
