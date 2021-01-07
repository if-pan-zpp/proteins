      subroutine eval_angles(nres,x0,y0,z0,fx,fy,fz,nch,ichain,
     $           kphi1,kphi2,phin,enephi,ktheta,thetan,eneth)
      implicit none
      integer nres,nch
      integer ichain(nch,2)
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 phin(nres),thetan(nres)
      real*8 kphi1,kphi2,enephi,ktheta,eneth
      real*8 phi,theta,dvdp
      real*8 fi(3),fj(3),fk(3),fl(3)
      real*8 rij(3),rkj(3),rkl(3)
      integer i1,i2,i3,i4,j,j1,j2
      
      enephi=0.d0
      eneth=0.d0
      do j=1,nch
       j1=ichain(j,1)
       j2=ichain(j,2)
       do i1=j1,j2-2
        i2=i1+1
        i3=i1+2
        i4=i1+3
        rij(1)=x0(i1)-x0(i2)
        rij(2)=y0(i1)-y0(i2)
        rij(3)=z0(i1)-z0(i2)
        rkj(1)=x0(i3)-x0(i2)
        rkj(2)=y0(i3)-y0(i2)
        rkj(3)=z0(i3)-z0(i2)
        call force_bending(rij,rkj,fi,fj,fk,theta)
        eneth=eneth+ktheta*(theta-thetan(i1))*(theta-thetan(i1))
        dvdp=2.d0*ktheta*(theta-thetan(i1))
        fx(i1)=fx(i1)-dvdp*fi(1)
        fy(i1)=fy(i1)-dvdp*fi(2)
        fz(i1)=fz(i1)-dvdp*fi(3)
        fx(i2)=fx(i2)-dvdp*fj(1)
        fy(i2)=fy(i2)-dvdp*fj(2)
        fz(i2)=fz(i2)-dvdp*fj(3)
        fx(i3)=fx(i3)-dvdp*fk(1)
        fy(i3)=fy(i3)-dvdp*fk(2)
        fz(i3)=fz(i3)-dvdp*fk(3)
        if (i4.le.j2) then
         rkl(1)=x0(i3)-x0(i4)
         rkl(2)=y0(i3)-y0(i4)
         rkl(3)=z0(i3)-z0(i4)
         call force_torsion(rij,rkj,rkl,fi,fj,fk,fl,phi)
         enephi=enephi+kphi1*(1.d0-cos(phi-phin(i1)))
     $         +kphi2*(1.d0-cos(3.d0*(phi-phin(i1))))
         dvdp=kphi1*sin(phi-phin(i1))
     $       +3.d0*kphi2*sin(3.d0*(phi-phin(i1)))
         fx(i1)=fx(i1)-dvdp*fi(1)
         fy(i1)=fy(i1)-dvdp*fi(2)
         fz(i1)=fz(i1)-dvdp*fi(3)
         fx(i2)=fx(i2)-dvdp*fj(1)
         fy(i2)=fy(i2)-dvdp*fj(2)
         fz(i2)=fz(i2)-dvdp*fj(3)
         fx(i3)=fx(i3)-dvdp*fk(1)
         fy(i3)=fy(i3)-dvdp*fk(2)
         fz(i3)=fz(i3)-dvdp*fk(3)
         fx(i4)=fx(i4)-dvdp*fl(1)
         fy(i4)=fy(i4)-dvdp*fl(2)
         fz(i4)=fz(i4)-dvdp*fl(3)
        endif
       enddo
      enddo
      return
      end
