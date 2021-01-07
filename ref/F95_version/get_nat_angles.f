      subroutine get_nat_angles(nres,x,y,z,nch,ichain,phi,theta)
      implicit none
      integer nres,nch
      integer ichain(nch,2)
      real*8 x(nres),y(nres),z(nres)
      real*8 phi(nres),theta(nres)
      real*8 rij(3),rkj(3),rkl(3),m(3),n(3)
      real*8 dn,dm,cosphi,rijrkj,d2ij,d2kj,costh
      integer i1,i2,i3,i4,k,j,j1,j2
      
      do i1=1,nres
       phi(i1)=0.d0
       theta(i1)=0.d0
      enddo
C DIHEDRAL ANGLES      
      do j=1,nch
      j1=ichain(j,1)
      j2=ichain(j,2)
       do i1=j1,j2-3
        i2=i1+1
        i3=i1+2
        i4=i1+3
        rij(1)=x(i1)-x(i2)
        rij(2)=y(i1)-y(i2)
        rij(3)=z(i1)-z(i2)
        rkj(1)=x(i3)-x(i2)
        rkj(2)=y(i3)-y(i2)
        rkj(3)=z(i3)-z(i2)
        rkl(1)=x(i3)-x(i4)
        rkl(2)=y(i3)-y(i4)
        rkl(3)=z(i3)-z(i4)
        m(1)=rij(2)*rkj(3)-rij(3)*rkj(2)
        m(2)=rij(3)*rkj(1)-rij(1)*rkj(3)
        m(3)=rij(1)*rkj(2)-rij(2)*rkj(1)
        n(1)=rkj(2)*rkl(3)-rkj(3)*rkl(2)
        n(2)=rkj(3)*rkl(1)-rkj(1)*rkl(3)
        n(3)=rkj(1)*rkl(2)-rkj(2)*rkl(1)
        dn=0.d0
        dm=0.d0
        do k=1,3
         dn=dn+n(k)*n(k)
         dm=dm+m(k)*m(k)
        enddo
        if (dn.gt.0.d0 .and. dm.gt.0.d0) then
         do k=1,3
          n(k)=n(k)/sqrt(dn)
          m(k)=m(k)/sqrt(dm)
         enddo
         cosphi=0.d0
         do k=1,3
          cosphi=cosphi+n(k)*m(k)
         enddo
         cosphi=min(cosphi,1.d0)
         cosphi=max(cosphi,-1.d0)
         dn=0.d0
         do k=1,3
          dn=dn+rij(k)*n(k)
         enddo
         if (dn.lt.0) then
           phi(i1)=-1.d0*acos(cosphi)
         else
           phi(i1)=acos(cosphi)
         endif
        endif
       enddo
      enddo
C BOND ANGLES
      do j=1,nch
       j1=ichain(j,1)
       j2=ichain(j,2)
       do i1=j1,j2-2
        i2=i1+1
        i3=i1+2
        rij(1)=x(i1)-x(i2)
        rij(2)=y(i1)-y(i2)
        rij(3)=z(i1)-z(i2)
        rkj(1)=x(i3)-x(i2)
        rkj(2)=y(i3)-y(i2)
        rkj(3)=z(i3)-z(i2)       
        rijrkj=0.d0
        d2ij=0.d0
        d2kj=0.d0
        do k=1,3
         rijrkj=rijrkj+rij(k)*rkj(k)
         d2ij=d2ij+rij(k)*rij(k)
         d2kj=d2kj+rkj(k)*rkj(k)
        enddo
        costh=rijrkj/sqrt(d2ij*d2kj)
        theta(i1)=acos(costh)
       enddo
      enddo
      return
      end
