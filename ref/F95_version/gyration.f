      subroutine gyration(nres,x0,y0,z0,rg)
C THIS SUBROUTINE COMPUTES THE RADIUS OF GYRATION
      implicit none
      integer nres
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 rg
      real*8 xcm,ycm,zcm
      real*8 dx,dy,dz
      integer i
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

      rg=0.d0
      do i=1,nres
       dx=x0(i)-xcm
       dy=y0(i)-ycm
       dz=z0(i)-zcm
       rg=rg+dx*dx+dy*dy+dz*dz
      enddo
      rg=rg/nres
      rg=sqrt(rg)

      return
      end
 
