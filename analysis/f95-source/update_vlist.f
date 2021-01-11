      subroutine update_vlist(nres,x0,y0,z0,icmap,krist,kront,isi,
     $           temp,delta,gamma,kupdate,rcut,unit)
     
      implicit none

      integer nres,kront
      real*8 x0(nres),y0(nres),z0(nres)
      integer icmap(nres,nres)
      integer krist(nres*nres,2)
      real*8 temp,delta,gamma,rcut,unit
      real*8 dl,r2v,rv,r2
      integer kupdate,i,j,isi,in
      
      do i=1,nres*nres
        krist(i,1)=0
        krist(i,2)=0
      enddo
      in=2
c
c      dl=sqrt(6*temp*delta/gamma)
c      rv=2*kupdate*dl+rcut
c      rv=rv/unit
c
      rv=12.d0/unit
      r2v=rv*rv
      kront=0
      do i=1,nres-in
        do j=i+in,nres
          if (icmap(i,j).eq.0) then
            r2=(x0(i)-x0(j))**2
     $        +(y0(i)-y0(j))**2
     $        +(z0(i)-z0(j))**2
            if (r2.le.r2v) then
              kront=kront+1
              krist(kront,1)=i
              krist(kront,2)=j
            endif
          endif
        enddo
      enddo
      end
