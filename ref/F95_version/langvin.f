      subroutine langvin(nres,x1,y1,z1,fx,fy,fz,gamma,delta,temp)
      implicit none
      integer nres
      double precision x1(nres), y1(nres), z1(nres)	! VELOCITIES
      double precision fx(nres),fy(nres),fz(nres)	! FORCES
      double precision gamma,delta,temp
      integer i
      double precision r1,r2,ran2
      double precision gam,const2,gamma2
      double precision twopi
      
      twopi=2.0d0*dacos(-1.d0)
      const2=sqrt(2.d0*temp*gamma*delta)*delta
      gamma2=gamma/delta
C X-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       x1(i)=x1(i)+const2*gam
       fx(i)=fx(i)-gamma2*x1(i)
      enddo
C Y-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       y1(i)=y1(i)+const2*gam
       fy(i)=fy(i)-gamma2*y1(i)
      enddo
C Z-COMPONENT
      do i=1,nres
       r1=ran2(0)
       r2=ran2(0)
       gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
       z1(i)=z1(i)+const2*gam
       fz(i)=fz(i)-gamma2*z1(i)
      enddo
      return
      end