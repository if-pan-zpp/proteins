      subroutine calc_w_rg(nres,x0,y0,z0,rg,w)
      implicit none
      integer nres
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 w,rg
      real*8 mit(3,3),vmit(3,3),ra(3)
      real*8 xcm,ycm,zcm,dx,dy,dz,rb,dr
      integer i,j,k,nrot
      
      xcm=0.d0
      ycm=0.d0
      zcm=0.d0
      do i=1,nres
       xcm=xcm+x0(i)
       ycm=ycm+y0(i)
       zcm=zcm+z0(i)
      enddo
      xcm=xcm/nres
      ycm=ycm/nres
      zcm=zcm/nres
      do i=1,3
       do j=1,3
        mit(i,j)=0.d0
       enddo
      enddo
      rg=0.d0
      do i=1,nres
       dx=x0(i)-xcm
       dy=y0(i)-ycm
       dz=z0(i)-zcm
       rg=rg+dx*dx+dy*dy+dz*dz
       mit(1,1)=mit(1,1)+dy*dy+dz*dz
       mit(2,2)=mit(2,2)+dz*dz+dx*dx
       mit(3,3)=mit(3,3)+dx*dx+dy*dy
       mit(1,2)=mit(1,2)-dx*dy
       mit(2,3)=mit(2,3)-dy*dz
       mit(3,1)=mit(3,1)-dz*dx
      enddo
      mit(2,1)=mit(1,2)
      mit(3,2)=mit(2,3)
      mit(1,3)=mit(3,1)
      rg=sqrt(rg/nres)
      call jacobi(mit,3,3,ra,vmit,nrot)
      do i=1,3
       ra(i)=sqrt(ra(i)/nres)
      enddo
      call sort_set(3,ra)
      rb=(ra(1)+ra(3))*0.5d0
      dr=ra(2)-rb
      w=dr/rb
      return
      end
