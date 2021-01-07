      subroutine get_zmin(nres,z0,boxl,unit,zmin,zcm,zmax,icw,icwmap)
      implicit none
      integer nres
      real*8 z0(nres)
      real*8 unit,boxl,zcut
      real*8 zmin,zcm,zmax,dz
      integer icwmap(nres)
      integer i,icw,icw0
      
      zcut=5.d0/unit
      zmax=-boxl
      zmin=boxl
      zcm=0.d0
      icw0=0
      do i=1,nres
       zmin=min(zmin,z0(i))
       zmax=max(zmax,z0(i))
       zcm=zcm+z0(i)
       dz=0.5d0*boxl+z0(i)
       if (dz.lt.zcut) then
        icw0=icw0+1
        icwmap(i)=icwmap(i)+1
       endif
      enddo
      zcm=zcm/nres
      icw=icw0
      return
      end
