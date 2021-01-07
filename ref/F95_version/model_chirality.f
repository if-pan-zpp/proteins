      subroutine model_chirality(nres,xn,yn,zn,bond,chirn,nch,ichain)
      implicit none
      integer nres,nch
      real*8 xn(nres),yn(nres),zn(nres)
      real*8 chirn(nres)
      integer ichain(nch,2)
      real*8 xs(nres),ys(nres),zs(nres)
      real*8 xa,ya,za,as,aa
      real*8 bond
      integer ib,ib1,ib2,k,i1,i2
      
      do k=1,nch
       i1=ichain(k,1)
       i2=ichain(k,2)
       do ib=i1,i2-1
        ib1=ib+1
        xs(ib)=xn(ib1)-xn(ib)
        ys(ib)=yn(ib1)-yn(ib)
        zs(ib)=zn(ib1)-zn(ib)
       enddo
      enddo
      do k=1,nch
       i1=ichain(k,1)
       i2=ichain(k,2)
       do ib=i1,i2-3
        ib1=ib+1
        ib2=ib+2
        xa=ys(ib)*zs(ib1)-zs(ib)*ys(ib1)
        ya=zs(ib)*xs(ib1)-xs(ib)*zs(ib1)
        za=xs(ib)*ys(ib1)-ys(ib)*xs(ib1)
        as=xa*xs(ib2)+ya*ys(ib2)+za*zs(ib2)
        aa=bond**3
        chirn(ib)=as/aa
       enddo
      enddo
      return
      end

