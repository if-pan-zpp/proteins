      subroutine leapfrog_dyn(nres,x0,y0,z0,x1,y1,z1,fx,fy,fz,
     $           gamma,delta,temp)
      implicit none
      integer nres,i
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 x1(nres),y1(nres),z1(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 gamma
      real*8 delta
      real*8 temp
      real*8 twopi
      real*8 xs,ys,zs
      real*8 r1,r2,ran2,gdrn,rn,ek
      real*8 dp,a1,a2,c1,c2
C USEFUL CONSTANTS      
      twopi=2.d0*dacos(-1.d0)
      dp=sqrt(2.d0*gamma*temp/delta)
      a1=1.d0-0.5d0*gamma*delta
      a2=1.d0+0.5d0*gamma*delta
      c1=a1/a2
      c2=delta*delta/a2
C X-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gdrn=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       xs=x0(i)+x1(i)*c1+(fx(i)+gdrn*dp)*c2
       x1(i)=xs-x0(i)
       x0(i)=xs
       enddo
C Y-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gdrn=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)  
       ys=y0(i)+y1(i)*c1+(fy(i)+gdrn*dp)*c2
       y1(i)=ys-y0(i)
       y0(i)=ys
      enddo
C Z-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gdrn=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       zs=z0(i)+z1(i)*c1+(fz(i)+gdrn*dp)*c2
       z1(i)=zs-z0(i)
       z0(i)=zs
      enddo
      return
      end
