      subroutine evalgo1012(nres,x0,y0,z0,fx,fy,fz,b,icmap,eij,
     $           kcont,klist,dij,krist,kront,epot,
     $           H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,
     $           jstep,cutparam,nr_l)
     
      implicit none
      
      integer nres
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 b(nres)
      integer icmap(nres,nres),imap(nres,nres)
      integer kcont,kront,nch
      integer conttime(kcont)
      integer klist(kcont,2)
      real*8 dij(kcont)
      real*8 eij(nres,nres)
      integer krist(nres*nres,2)
      integer ichain(nch,2)
      real*8 epot,H1,H2,rcut,rcut2,sigma0,cutparam
      integer i,j,k,i1,i2,inc,jstep,nr_l
      real*8 dx,dy,dz,dxn,dyn,dzn,r,r2,rb,rn,rb2
      real*8 ene,fce,repx,repy,repz,dcut
      real*8 dis,d6,dm12,dm10,dm6,dm5,dm3,dm2,eps
      
      epot=0.d0
      inc=0
      do i=1,nres
       do j=1,nres
        imap(i,j)=0
       enddo
      enddo 
      
      do k=1,nch
       i1=ichain(k,1)
       i2=ichain(k,2)
       do i=i1,i2-1
        j=i+1
        dx=x0(i)-x0(j)
        dy=y0(i)-y0(j)
        dz=z0(i)-z0(j)
        r2=dx*dx+dy*dy+dz*dz
        r=sqrt(r2)
        rb=r-b(i)
        rb2=rb*rb
        ene=H1*rb2+H2*rb2*rb2
        fce=(2*H1+4*H2*rb2)*rb
        epot=epot+ene
        repx=-fce*dx/r
        repy=-fce*dy/r
        repz=-fce*dz/r
        fx(i)=fx(i)+repx
        fx(j)=fx(j)-repx
        fy(i)=fy(i)+repy
        fy(j)=fy(j)-repy
        fz(i)=fz(i)+repz
        fz(j)=fz(j)-repz
       enddo
      enddo
      do k=1,kcont
       i=klist(k,1)
       j=klist(k,2)
       dx=x0(i)-x0(j)
       dy=y0(i)-y0(j)
       dz=z0(i)-z0(j)
       r2=dx*dx+dy*dy+dz*dz
       r=sqrt(r2)
       
       dcut=dij(k)*cutparam
       if (r.le.dcut) then
         imap(i,j)=1
         inc=inc+1
        if (nr_l.eq.3) then  
         conttime(k)=jstep
        endif
        if (nr_l.eq.2) then 
         if (conttime(k).eq.0) conttime(k)=jstep
        endif
       endif        
              
       if (icmap(i,j).eq.1) then
        eps=eij(i,j)
        dis=dij(k)/r
        dm2=dis*dis
        dm3=dm2*dis
        dm5=dm2*dm3
        dm6=dm3*dm3
        dm10=dm5*dm5
        dm12=dm6*dm6
        ene=eps*(5.d0*dm12-6.d0*dm10)
        epot=epot+ene
        fce=-60.d0*eps*(dm12-dm10)/r2
        repx=-fce*dx
        repy=-fce*dy
        repz=-fce*dz
        fx(i)=fx(i)+repx
        fx(j)=fx(j)-repx
        fy(i)=fy(i)+repy
        fy(j)=fy(j)-repy
        fz(i)=fz(i)+repz
        fz(j)=fz(j)-repz
       endif
       if (icmap(i,j).eq.3) then
        rn=dij(k)
        rb=r-rn
        rb2=rb*rb
        ene=H1*rb2+H2*rb2*rb2
        fce=(2*H1+4*H2*rb2)*rb
        epot=epot+ene
        repx=-fce*dx/r
        repy=-fce*dy/r
        repz=-fce*dz/r
        fx(i)=fx(i)+repx
        fx(j)=fx(j)-repx
        fy(i)=fy(i)+repy
        fy(j)=fy(j)-repy
        fz(i)=fz(i)+repz
        fz(j)=fz(j)-repz
       endif
      enddo
      rcut2=rcut**2
      do k=1,kront
        i=krist(k,1)
        j=krist(k,2)
        dx=x0(i)-x0(j)
        dy=y0(i)-y0(j)
        dz=z0(i)-z0(j)
        r2=dx*dx+dy*dy+dz*dz
        if (r2.le.rcut2) then
           r=sqrt(r2)        
           dis=sigma0/r
           d6=dis**6
           ene=4.d0*d6*(d6-1.d0)+1.d0
           epot=epot+ene
           fce=24.d0*d6*(1.d0-2.d0*d6)/r2
           repx=-fce*dx
           repy=-fce*dy
           repz=-fce*dz
           fx(i)=fx(i)+repx
           fx(j)=fx(j)-repx
           fy(i)=fy(i)+repy
           fy(j)=fy(j)-repy
           fz(i)=fz(i)+repz
           fz(j)=fz(j)-repz
        endif
      enddo
      return
      end