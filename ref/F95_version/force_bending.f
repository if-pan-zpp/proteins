      subroutine force_bending(rij,rkj,fi,fj,fk,theta)
C SWOPE, W. C. AND FERGUSON, D. M. (1992)
C ALTERNATIVE EXPRESSIONS FOR ENERGIES AND FORCES DUE TO ANGLE BENDING AND TORSIONAL ENERGY
C J. COMPUT. CHEM. 13: 585-594. DOI: 10.1002/JCC.540130508
      implicit none
      real*8 ri(3),rj(3),rk(3)
      real*8 fi(3),fj(3),fk(3)
      real*8 rij(3),rkj(3)
      real*8 rp(3),qij(3),qkj(3),qp(3),a(3),b(3)
      real*8 theta,rijrkj,d2ij,d2kj,dij,dkj,drp,costh
      integer k
      
      rijrkj=0.d0
      d2ij=0.d0
      d2kj=0.d0
      do k=1,3
       rijrkj=rijrkj+rij(k)*rkj(k)
       d2ij=d2ij+rij(k)*rij(k)
       d2kj=d2kj+rkj(k)*rkj(k)
      enddo
      dij=sqrt(d2ij)
      dkj=sqrt(d2kj)
      costh=rijrkj/(dij*dkj)
      costh=min(costh,1.d0)
      costh=max(costh,-1.d0)
      theta=acos(costh)
      if (costh.eq.0.d0) then
       do k=1,3
        fi(k)=0.d0
        fj(k)=0.d0
        fk(k)=0.d0
       enddo
      else
       rp(1)=rij(2)*rkj(3)-rij(3)*rkj(2)
       rp(2)=rij(3)*rkj(1)-rij(1)*rkj(3)
       rp(3)=rij(1)*rkj(2)-rij(2)*rkj(1)
       drp=sqrt(rp(1)*rp(1)+rp(2)*rp(2)+rp(3)*rp(3))
       do k=1,3
        qij(k)=rij(k)/dij
        qkj(k)=rkj(k)/dkj
        qp(k)=rp(k)/drp
       enddo
       a(1)=qij(2)*qp(3)-qij(3)*qp(2)
       a(2)=qij(3)*qp(1)-qij(1)*qp(3)
       a(3)=qij(1)*qp(2)-qij(2)*qp(1)
       b(1)=qkj(2)*qp(3)-qkj(3)*qp(2)
       b(2)=qkj(3)*qp(1)-qkj(1)*qp(3)
       b(3)=qkj(1)*qp(2)-qkj(2)*qp(1)
       do k=1,3
        fi(k)=a(k)/dij
        fk(k)=-b(k)/dkj
        fj(k)=-fi(k)-fk(k)
       enddo
      endif
      return
      end
