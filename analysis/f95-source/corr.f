      subroutine corr(nres,x0,y0,z0,x1,y1,z1,x2,y2,z2,
     $ x3,y3,z3,x4,y4,z4,x5,y5,z5,fx,fy,fz,delta)
C  CORRECT PREDICTED POSITIONS AND THEIR DERIVATIVES
      implicit none
      integer nres
      double precision x0(nres),y0(nres),z0(nres)
      double precision x1(nres),y1(nres),z1(nres)
      double precision x2(nres),y2(nres),z2(nres)
      double precision x3(nres),y3(nres),z3(nres)
      double precision x4(nres),y4(nres),z4(nres)
      double precision x5(nres),y5(nres),z5(nres)
      double precision fx(nres),fy(nres),fz(nres)
      double precision delta
      double precision f02,f12,f32,f42,f52
      integer i
      double precision xerr,yerr,zerr
      double precision deltsq
      deltsq=0.5d0*delta*delta
      f02=dble(3./16.)
      f12=dble(251./360.)
      f32=dble(11./18.) 
      f42=dble(1./6.)
      f52=dble(1./60.)
      do i=1,nres
        xerr=x2(i)-deltsq*fx(i)
        yerr=y2(i)-deltsq*fy(i)
        zerr=z2(i)-deltsq*fz(i)
        x0(i)=x0(i)-xerr*f02
        x1(i)=x1(i)-xerr*f12
        x2(i)=x2(i)-xerr
        x3(i)=x3(i)-xerr*f32
        x4(i)=x4(i)-xerr*f42
        x5(i)=x5(i)-xerr*f52
        y0(i)=y0(i)-yerr*f02
        y1(i)=y1(i)-yerr*f12
        y2(i)=y2(i)-yerr
        y3(i)=y3(i)-yerr*f32
        y4(i)=y4(i)-yerr*f42
        y5(i)=y5(i)-yerr*f52
        z0(i)=z0(i)-zerr*f02
        z1(i)=z1(i)-zerr*f12
        z2(i)=z2(i)-zerr
        z3(i)=z3(i)-zerr*f32
        z4(i)=z4(i)-zerr*f42
        z5(i)=z5(i)-zerr*f52      
      enddo
      return
      end
