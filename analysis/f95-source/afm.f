      subroutine afm(nres,x0,y0,z0,xn,yn,zn,ip1,ip2,fx,fy,fz,fresist,
     $     HH1,HH2,afx,afy,afz,xpull,ypull,zpull,aforce)
      implicit none
      integer nres,ip1,ip2
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 xn(nres),yn(nres),zn(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 dx,dy,dz,rsq,r,ene1,fce,repx,repy,repz
      real*8 xpull,ypull,zpull,aforce,fax,fay,faz
      real*8 afx,afy,afz,fresist,HH1,HH2
C STICK THE N-TERMINUS WITH A HARMONIC POTENTIAL
      dx=x0(ip1)-xn(ip1)
      dy=y0(ip1)-yn(ip1)
      dz=z0(ip1)-zn(ip1)
      rsq=dx*dx+dy*dy+dz*dz
      r=sqrt(rsq)
      ene1=HH1*rsq+HH2*rsq*rsq
      fce=(2*HH1+4*HH2*rsq)
      repx=-fce*dx
      repy=-fce*dy
      repz=-fce*dz
      fx(ip1)=fx(ip1)+repx
      fy(ip1)=fy(ip1)+repy
      fz(ip1)=fz(ip1)+repz
C CONSTANT VELOCITY IN THE DIRECTION PARRALLEL TO THE LINE
C CONNECTING FIRST AND LAST MONOMERS IN THE NATIVE STATE
      
      fax=aforce*afx
      fay=aforce*afy
      faz=aforce*afz 
      
      fx(ip2)=fx(ip2)+fax
      fy(ip2)=fy(ip2)+fay
      fz(ip2)=fz(ip2)+faz
      
      return
      end

