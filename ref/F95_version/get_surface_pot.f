      subroutine get_surface_pot(pot_file,nres,seq,nz,wz,wp,wf,
     $           unit,temp)
      implicit none
      character (len=*) pot_file
      integer nres,nz
      character*3 seq(nres)
      real*8 wz(nres,nz)
      real*8 wp(nres,nz)
      real*8 wf(nres,nz)
      real*8 unit,temp
      real*8 v(3)
      integer icount(nres)
      character*3 aa
      character*80 line
      integer i,j,k
      
      do i=1,nres
       icount(i)=0
       do j=1,nz
        wz(i,j)=0.d0
        wp(i,j)=0.d0
        wf(i,j)=0.d0
       enddo
      enddo
      open(24,file=trim(pot_file),form='formatted',status='old')
11    read(24,'(a)',end=12,err=12) line
      read(line(1:3), '(a3)') aa
      read(line(4:33), '(3f10.4)') (v(k),k=1,3)
      do i=1,nres
       if (seq(i).eq.aa) then
        icount(i)=icount(i)+1
        j=icount(i)
        wz(i,j)=v(1)/unit
        wp(i,j)=v(2)*temp
        wf(i,j)=v(3)*temp*unit
       endif
      enddo
      goto 11
12    close(24)
      do i=1,nres
       j=icount(i)
       if (j.lt.nz) then
        do k=j+1,nz
         wz(i,k)=wz(i,j)+(k-j)*0.1d0/unit
         wp(i,k)=0.d0
         wf(i,k)=0.d0
        enddo
       endif
      enddo
C TEST
c      do i=1,nres
c       do j=1,nz
c        write(6,'(i6,2x,a3,2x,3f10.3)') i,seq(i),wz(i,j),wp(i,j),wf(i,j)
c       enddo
c      enddo
c      stop
C
      return
      end
