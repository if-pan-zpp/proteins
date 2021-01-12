      subroutine predct(nres,x0,y0,z0,x1,y1,z1,x2,y2,z2,
     $ x3,y3,z3,x4,y4,z4,x5,y5,z5)
C USE FIFTH-ORDER TAYLOR SERIES TO PREDICT POSITIONS & THEIR
C DERIVATIVES AT NEXT TIME-STEP
      implicit none
      integer nres
      double precision x0(nres),y0(nres),z0(nres)
      double precision x1(nres),y1(nres),z1(nres)
      double precision x2(nres),y2(nres),z2(nres)
      double precision x3(nres),y3(nres),z3(nres)
      double precision x4(nres),y4(nres),z4(nres)
      double precision x5(nres),y5(nres),z5(nres)
      integer i
      do i=1,nres
       x0(i)=x0(i)+x1(i)+x2(i)+x3(i)+x4(i)+x5(i)
       y0(i)=y0(i)+y1(i)+y2(i)+y3(i)+y4(i)+y5(i)
       z0(i)=z0(i)+z1(i)+z2(i)+z3(i)+z4(i)+z5(i)
       x1(i)=x1(i)+2.d0*x2(i)+3.d0*x3(i)+4.d0*x4(i)+5.d0*x5(i)
       y1(i)=y1(i)+2.d0*y2(i)+3.d0*y3(i)+4.d0*y4(i)+5.d0*y5(i)
       z1(i)=z1(i)+2.d0*z2(i)+3.d0*z3(i)+4.d0*z4(i)+5.d0*z5(i)
       x2(i)=x2(i)+3.d0*x3(i)+6.d0*x4(i)+10.d0*x5(i)
       y2(i)=y2(i)+3.d0*y3(i)+6.d0*y4(i)+10.d0*y5(i)
       z2(i)=z2(i)+3.d0*z3(i)+6.d0*z4(i)+10.d0*z5(i)
       x3(i)=x3(i)+4.d0*x4(i)+10.d0*x5(i)
       y3(i)=y3(i)+4.d0*y4(i)+10.d0*y5(i)
       z3(i)=z3(i)+4.d0*z4(i)+10.d0*z5(i)
       x4(i)=x4(i)+5.d0*x5(i)
       y4(i)=y4(i)+5.d0*y5(i)
       z4(i)=z4(i)+5.d0*z5(i)
      enddo
      return
      end
 
