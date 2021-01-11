      subroutine compute_rmsd(nres,x0,y0,z0,xn,yn,zn,rms,kw,kw2)
      implicit none
C THIS SUBROUTINES COMPUTE THE RMSD OF THE C ALPHA BACKBONE TO THE
C NATIVE BACKBONE TAKEN FROM PDB
      integer nres
      real*8 x0(nres),y0(nres),z0(nres)		! COORDINATES
      real*8 xn(nres),yn(nres),zn(nres)		! NATIVE COORDINATES
      real*8 r(nres,3),rn(nres,3),qll(nres,3)
      real*8 kw(nres,3),kw2(nres,3)
      real*8 rms
      integer ib,j
      do ib=1,nres
       r(ib,1)=x0(ib)
       r(ib,2)=y0(ib)
       r(ib,3)=z0(ib)
       rn(ib,1)=xn(ib)
       rn(ib,2)=yn(ib)
       rn(ib,3)=zn(ib)
      enddo
      call kabsch(nres,r,rn,rms,qll)
      do ib=1,nres
       do j=1,3
       kw(ib,j)=kw(ib,j)+qll(ib,j)
       kw2(ib,j)=kw2(ib,j)+qll(ib,j)*qll(ib,j)
       enddo
      enddo
      return
      end
