      subroutine gopotential612(nres,xn,yn,zn,icmap,b,bond,unit,
     $           kcont,klist,sij,dij,rcut,sigma0,nch,ichain,rcut3)
      implicit none
      integer nres
      real*8 xn(nres),yn(nres),zn(nres)
      integer icmap(nres,nres)
      real*8 b(nres)
      integer kcont
      integer klist(kcont,2)
      integer nch
      integer ichain(nch,2)
      real*8 sij(kcont)
      real*8 dij(kcont)
      integer i,j,k,i1,i2
      real*8 dx,dy,dz,d2,d1,sigmaij,rcut,sigma0,unit
      real*8 bond,rcut3
      
      bond=3.8d0/unit
      do i=1,nres
         b(i)=bond
      enddo
      bond=0.0d0
      j=0
      do k=1,nch
       i1=ichain(k,1)
       i2=ichain(k,2)
       do i=i1,i2-1
         j=j+1
         dx=xn(i)-xn(i+1)
         dy=yn(i)-yn(i+1)
         dz=zn(i)-zn(i+1)
         d2=dx*dx+dy*dy+dz*dz
         d1=dsqrt(d2)
         b(i)=d1
         bond=bond+d1
       enddo
      enddo
      bond=bond/j
      
      k=0
      do i=1,nres
        do j=1,nres
          if (icmap(i,j).ne.0) then
            k=k+1
            klist(k,1)=i
            klist(k,2)=j
          endif
        enddo
      enddo
      !if (k.ne.kcont) write(*,*)'ERROR'
      do k=1,kcont
         i=klist(k,1)
         j=klist(k,2)
         dx=xn(i)-xn(j)
         dy=yn(i)-yn(j)
         dz=zn(i)-zn(j)
         d2=dx*dx+dy*dy+dz*dz
         d1=sqrt(d2)
         dij(k)=d1
         sigmaij=d1*0.5d0**(1.d0/6.d0)
         sij(k)=sigmaij
      enddo
      rcut=rcut3/unit
      sigma0=rcut*0.5d0**(1.d0/6.d0)
      return
      end
