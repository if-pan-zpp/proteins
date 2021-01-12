      subroutine evalwall_Au(nres,z0,fz,nch,ichain,wpd,wpm,boxl,
     $                       enewall)
      implicit none
      integer nres
      real*8 z0(nres)
      real*8 fz(nres)
      integer nch
      integer ichain(nch,2)
      real*8 wpd(nres),wpm(nres)
      real*8 boxl,enewall
      real*8 ene,dfz,z,y
      integer i,i1,i2,j
      enewall=0.d0
      do i=1,nch
       i1=ichain(i,1)
       i2=ichain(i,2)
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
        z=0.5d0*boxl+z0(j)
        y=(wpm(j)/z)**6
        ene=wpd(j)*y*(y-2.d0)
        enewall=enewall+ene
        dfz=12.d0*wpd(j)*y*(y-1.d0)/z
        fz(j)=fz(j)+dfz
       enddo
      enddo
      return
      end
