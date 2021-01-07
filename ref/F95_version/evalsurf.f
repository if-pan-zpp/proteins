      subroutine evalsurf(nres,z0,fz,nz,wz,wp,wf,enewall,boxl)
      implicit none
      integer nres,nz
      real*8 z0(nres),fz(nres)
      real*8 wz(nres,nz),wp(nres,nz),wf(nres,nz)
      real*8 enewall,boxl
      real*8 ene,df,z1,z2,zi,y
      integer i,j
      
      enewall=0.d0
C   SURFACE AT Z=-BOXL/2
      do i=1,nres
       zi=0.5d0*boxl+z0(i)
       if (zi.le.wz(i,nz)) then
        do j=1,nz-1
         z1=wz(i,j)
         z2=wz(i,j+1)
         if ((zi.gt.z1).and.(zi.le.z2)) then
          ene=wp(i,j)+(wp(i,j+1)-wp(i,j))*(zi-z1)/(z2-z1)
           df=wf(i,j)+(wf(i,j+1)-wf(i,j))*(zi-z1)/(z2-z1)
          enewall=enewall+ene
          fz(i)=fz(i)+df
         endif
        enddo
       endif
      enddo
C   SOFT WALL AT Z=BOXL/2
      do i=1,nres
       zi=0.5d0*boxl-z0(i)
       y=sqrt(3.d0)/zi**3
       if (y.ge.1.d0) then
        df=-45.d0*(y-1.d0)*(y+1.d0)*y/zi
        fz(i)=fz(i)+df
       endif
      enddo
C TEST OF THE SOFT WALL
c      do i=1,nres
c       zi=0.5d0*boxl+z0(i)
c       y=sqrt(3.d0)/zi**3
c       if (y.ge.1.d0) then
c        df=45.d0*(y-1.d0)*(y+1.d0)*y/zi
c        fz(i)=fz(i)+df
c       endif
c      enddo
C
      return
      end
      