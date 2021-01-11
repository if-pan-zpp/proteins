      subroutine line_config(nres,ntrj,xu,yu,zu,nch,ichain,unit)
      implicit none
      integer nres,ntrj,nch
      real*8 xu(nres,ntrj),yu(nres,ntrj),zu(nres,ntrj)
      integer ichain(nch,2)
      real*8 bond,unit
      parameter (bond=3.8d0)
      integer i,j,i1,i2,k
      
      do k=1,ntrj
       do i=1,nch
        i1=ichain(i,1)
        i2=ichain(i,2)
        do j=i1,i2
         xu(j,k)=bond*i/unit
         yu(j,k)=0.d0
         zu(j,k)=bond*(j-i1)/unit
        enddo
       enddo
      enddo
      return
      end
