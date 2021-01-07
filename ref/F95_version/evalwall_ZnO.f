      subroutine evalwall_ZnO(nres,z0,fz,nch,ichain,wpd,wpm,boxl,
     $           enewall,icw)
      implicit none
      integer nres
      real*8 z0(nres)
      real*8 fz(nres)
      integer nch
      integer ichain(nch,2)
      real*8 wpd(nres),wpm(nres)
      real*8 boxl,enewall
      real*8 a,a1,a2,ene,dfz,z,y
      integer i,i1,i2,j,icw
      a=8.3d0
      a1=(a+1.d0)/(9.d0-a)
      a2=10.d0*exp(a)/(9.d0-a)
      enewall=0.d0
      do i=1,nch
       i1=ichain(i,1)
       i2=ichain(i,2)
       icw=0
       do j=i1,i2
C HARD WALL AT Z=BOXL/2
        z=0.5d0*boxl-z0(j)
        y=sqrt(3.d0)/z**3
        if (y.ge.1.d0) then
         ene=0.5d0*(y**3-3.d0*y)+1.d0
         dfz=-4.5d0*(y-1.d0)*(y+1.d0)*y/z
         fz(j)=fz(j)+dfz
        endif
C ZnO SURFACE AT Z=-BOXL/2
        z=(0.5d0*boxl+z0(j))/wpm(j)
        ene=(a1/z**10-a2*exp(-a*z)/z)*wpd(j)
        if (ene.lt.0) icw=icw+1
        dfz=(10.d0*a1/z**11-a2*(1.d0+a*z)*exp(-a*z)/z**2)*wpd(j)/wpm(j)
        fz(j)=fz(j)+dfz
        enewall=enewall+ene
       enddo
      enddo
      return
      end
