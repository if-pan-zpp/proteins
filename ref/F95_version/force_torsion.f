      subroutine force_torsion(rij,rkj,rkl,fi,fj,fk,fl,phi)
C B EKKER, H., BERENDSEN, H. J. C. AND VAN GUNSTEREN, W. F. (1995)
C FORCE AND VIRIAL OF TORSIONAL-ANGLE-DEPENDENT POTENTIALS
C J. COMPUT. CHEM., 16: 527-533. DOI: 10.1002/JCC.540160502
      implicit none
      real*8 fi(3),fj(3),fk(3),fl(3)
      real*8 rij(3),rkj(3),rkl(3),m(3),n(3)
      real*8 d2m,d2n,dmn,d2rkj,rijn,drkj,rijrkj,rklrkj,cosphi,df,phi
      integer i,k
      
      m(1)=rij(2)*rkj(3)-rij(3)*rkj(2)
      m(2)=rij(3)*rkj(1)-rij(1)*rkj(3)
      m(3)=rij(1)*rkj(2)-rij(2)*rkj(1)
      n(1)=rkj(2)*rkl(3)-rkj(3)*rkl(2)
      n(2)=rkj(3)*rkl(1)-rkj(1)*rkl(3)
      n(3)=rkj(1)*rkl(2)-rkj(2)*rkl(1)
      d2n=0.d0
      d2m=0.d0
      d2rkj=0.d0
      do k=1,3
       d2n=d2n+n(k)*n(k)
       d2m=d2m+m(k)*m(k)
       d2rkj=d2rkj+rkj(k)*rkj(k)
      enddo
      if (d2n.eq.0.d0 .or. d2m.eq.0.d0) then
       phi=0.d0
       do k=1,3
        fi(k)=0.d0
        fj(k)=0.d0
        fk(k)=0.d0
        fl(k)=0.d0
       enddo
      else
       dmn=0.d0
       do k=1,3
         dmn=dmn+n(k)*m(k)
       enddo
       cosphi=dmn/sqrt(d2n*d2m)
       cosphi=min(cosphi,1.d0)
       cosphi=max(cosphi,-1.d0)
       rijn=0.d0
       do k=1,3
        rijn=rijn+rij(k)*n(k)
       enddo
       if (rijn.lt.0.d0) then
        phi=-1.d0*acos(cosphi)
       else
        phi=acos(cosphi)
       endif
       drkj=sqrt(d2rkj)
       do k=1,3
        fi(k)=m(k)*drkj/d2m
        fl(k)=-n(k)*drkj/d2n
       enddo
       rijrkj=0.d0
       rklrkj=0.d0
       do k=1,3
        rijrkj=rijrkj+rij(k)*rkj(k)
        rklrkj=rklrkj+rkl(k)*rkj(k)
       enddo
       do k=1,3
        df=(fi(k)*rijrkj-fl(k)*rklrkj)/d2rkj
        fj(k)=-fi(k)+df
        fk(k)=-fl(k)-df
       enddo
      endif
      return
      end
