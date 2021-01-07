      subroutine center(nres,x0,y0,z0)
      implicit none
      integer nres,i
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 xcm,ycm,zcm
      xcm=0.d0
      ycm=0.d0
      zcm=0.d0
      do i=1,nres
       xcm=xcm+x0(i)
       ycm=ycm+y0(i)
       zcm=zcm+z0(i)
      enddo
      xcm=xcm/nres
      ycm=ycm/nres
      zcm=zcm/nres
      do i=1,nres
       x0(i)=x0(i)-xcm
       y0(i)=y0(i)-ycm
       z0(i)=z0(i)-zcm
      enddo
      return
      end
