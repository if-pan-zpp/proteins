      subroutine contact_fluct(nres,kcont,klist,dij,x0,y0,z0,drij,dr2ij)
      implicit none
      integer nres
      integer kcont
      integer klist(kcont,2)
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 dij(kcont)
      real*8 drij(kcont)
      real*8 dr2ij(kcont)
      real*8 d2,d1
      integer i,j,k
      
      do k=1,kcont
       i=klist(k,1)
       j=klist(k,2)
       d2=(x0(i)-x0(j))**2
     $   +(y0(i)-y0(j))**2
     $   +(z0(i)-z0(j))**2
       d1=sqrt(d2)
       drij(k)=drij(k) + d1-dij(k)
       dr2ij(k)=dr2ij(k) + (d1-dij(k))**2
      enddo
      return
      end
