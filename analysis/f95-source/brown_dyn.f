      subroutine brown_dyn(nres,x0,y0,z0,fx,fy,fz,gamma,delta,temp)
      implicit none
      integer nres,i
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 gamma
      real*8 delta
      real*8 temp
      real*8 twopi
      real*8 r1,r2,ran2,gdrn,dl,dt
      real*8 xs,ys,zs
      
      twopi=2.d0*dacos(-1.d0)
      dl=sqrt(2.d0*temp*delta/gamma)
      dt=delta/gamma
      
C X-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gdrn=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       xs=x0(i)+fx(i)*dt+gdrn*dl
       x0(i)=xs
      enddo
C Y-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gdrn=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       ys=y0(i)+fy(i)*dt+gdrn*dl
       y0(i)=ys
      enddo
C Z-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gdrn=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       zs=z0(i)+fz(i)*dt+gdrn*dl
       z0(i)=zs
      enddo
      return
      end
