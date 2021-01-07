      subroutine kabsch(ns,q,q_trg,d,qkk)
      implicit none
      integer ns
      real*8 q(ns,3),q_trg(ns,3),qk(ns,3),qkk(ns,3)
      real*8 r(3,3),s(3,3),u(3,3),a(3,3),b(3,3),c(3,3)
      real*8 q_trg_cm(3),q_cm(3),aumu(3),v1(3),v2(3),v3(3)
      real*8 e_0,sigma3,bdotc,e,d,vmod
      integer k,i,l,m,nrot
C*********************************************************************
C               NO TRANSLATIONS
C*********************************************************************
      e_0=0.E0
      do k=1,3
       q_trg_cm(k) =  0.E0
       q_cm(k) = 0.E0
       do i=1,ns
        q_trg_cm(k) = q_trg_cm(k) + q_trg(i,k)
        q_cm(k)     = q_cm(k) + q(i,k)
       enddo
       q_trg_cm(k) = q_trg_cm(k)/ns
       q_cm(k)     = q_cm(k)/ns
       do i=1,ns
        qk(i,k)    = q(i,k) - q_cm(k)
        q_trg(i,k) = q_trg(i,k) - q_trg_cm(k)
C NOW THEY HAVE THE SAME CENTER OF MASS: (0,0,0)
        e_0=e_0+qk(i,k)*qk(i,k)+q_trg(i,k)*q_trg(i,k)
       enddo
      enddo
      e_0=e_0/(1.0*ns)
C NO MORE NEED FOR Q, JUST QK!!!  CALCOLO MATRICE R
      do k=1,3
       do l=1,3
          r(k,l)=0.E0
          do i=1,ns
             r(k,l) = r(k,l) + q_trg(i,k)*qk(i,l)
          enddo
          r(k,l)=2.E0*r(k,l)/(1.0*ns)
       enddo
      enddo
C*****************************************************
C    CALCULATE MATRIX S = R_T*R
C*****************************************************
      do k=1,3
        do l=1,3
          s(k,l)=0.E0
          do m=1,3
             s(k,l)=s(k,l)+r(m,k)*r(m,l)
          enddo
        enddo
      enddo
      call jacobi(s,3,3,aumu,a,nrot)
      call eigsrt(aumu,a,3,3)
C CALCOLO A3 COME A1XA2
      do i=1,3
       v1(i)=a(i,1)
       v2(i)=a(i,2)
      enddo
      call pvector (v1,v2,v3,vmod)
      do k=1,3
       a(k,3)=v3(k)
      enddo
C CALCOLO I VETTORI B
      do k=1,3
       do l=1,3
          b(k,l)=0.E0
          do m=1,3
             b(k,l)=b(k,l)+r(k,m)*a(m,l)
          enddo
          c(k,l)=b(k,l)
       enddo
      enddo
C********************************************************************
C      NORMALIZATION
C********************************************************************
      call norma(b)
      do i=1,3
        v1(i)=b(i,1)
        v2(i)=b(i,2)
      enddo
      call pvector (v1,v2,v3,vmod)
      do i=1,3
        b(i,3)=v3(i)
      end do
C CONTROLLO I PRODOTTI SCALARI DEGLI A_K CON I B_K
      bdotc=0.E0
      do l=1,3
        bdotc=bdotc+b(l,3)*c(l,3)
      end do
C SE IL PRODOTTO SCALARE DEL TERZO E'<0, ALLORA CAMBIA IL SEGNO
      sigma3 = 1.0
      if (bdotc.lt.0.E0) sigma3=-1.0
C************************************************************************
C           MATRICE U
C************************************************************************
      do l=1,3
       do m=1,3
          u(l,m)=0.E0
          do k=1,3
             u(l,m)=u(l,m)+b(l,k)*a(m,k)
          enddo
       enddo
      enddo
C******************** ROTAZIONE DI KABSCH   *****************************
      do i=1,ns
        do k=1,3
          qkk(i,k)=0.E0
          do l=1,3
            qkk(i,k)=qkk(i,k)+u(k,l)*qk(i,l)
          enddo
        enddo
      enddo
C**************************************************************************
C        DISTANZA DI KABSCH
C**************************************************************************
      d=0.E0
      do i=1,ns
C PRINT *,(QKK(I,K),K=1,3)  ! SERVE PER SCRIVERE LE COORDINATE
        do k =1,3
           d=d+(q_trg(i,k)-qkk(i,k))**2
         enddo
      enddo
      d=sqrt(d/ns)
      e = sqrt(e_0-sqrt(aumu(1))-sqrt(aumu(2))-sigma3*sqrt(aumu(3)))
      return
      end
      