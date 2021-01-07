      subroutine intvel3d(nres,x1,y1,z1,temp,delta,unit)
C  ASSIGN INITIAL VELOCITIES TO ATOMS
      implicit none
      integer nres,i
      real*8 x1(nres),y1(nres),z1(nres)
      real*8 vx(nres),vy(nres),vz(nres)
      real*8 temp,delta,unit
      real*8 sumx,sumy,sumz,xx,yy,zz,xyz,av2,ran2

      sumx=0.d0
      sumy=0.d0
      sumz=0.d0
      do i=1,nres
       xx=2.d0*(ran2(0)-0.5d0)
       yy=2.d0*(ran2(0)-0.5d0)
       zz=2.d0*(ran2(0)-0.5d0)
       xyz=1.d0/sqrt(xx*xx+yy*yy+zz*zz)
       vx(i)=xx*xyz
       vy(i)=yy*xyz
       vz(i)=zz*xyz
       sumx=sumx+vx(i)
       sumy=sumy+vy(i)
       sumz=sumz+vz(i)
      enddo
C SCALE VELOCITIES SO THAT TOTAL MOMENTUM = ZERO
      av2=0.d0
      do i=1,nres
       vx(i)=vx(i)-sumx/nres
       vy(i)=vy(i)-sumy/nres
       vz(i)=vz(i)-sumz/nres
       av2=av2+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
      enddo
      av2=av2/nres
C SCALE VELOCITIES TO DISIRED TEMPERATURE
      xyz=sqrt(3.d0*temp/av2)
      do i=1,nres
       vx(i)=vx(i)*xyz
       vy(i)=vy(i)*xyz
       vz(i)=vz(i)*xyz
      enddo
      do i=1,nres
       x1(i)=vx(i)*delta/unit
       y1(i)=vy(i)*delta/unit
       z1(i)=vz(i)*delta/unit
      enddo
      return
      end