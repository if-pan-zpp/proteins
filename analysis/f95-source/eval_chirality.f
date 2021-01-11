      subroutine eval_chirality(nres,x0,y0,z0,fx,fy,fz,
     &           bond,chirn,echi,enechi,nch,ichain)

C THIS SUBROUTINE COMPUTES THE ENERGY AND FORCE GIVEN BY THE CHIRALITY POTENTIALS

      implicit none
      
      integer nres,nch
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 chir(nres),chirn(nres)
      integer ichain(nch,2)
      real*8 xs(nres),ys(nres),zs(nres)
      real*8 bond,enechi,echi,xa,ya,za,aa,cc
      real*8 repx,repy,repz,fchi
      integer ib,ib1,ib2,ib3,k,i1,i2
      
      do k=1,nch
       i1=ichain(k,1)
       i2=ichain(k,2)
       do ib=i1,i2-1
        ib1=ib+1
        xs(ib)=x0(ib1)-x0(ib)
        ys(ib)=y0(ib1)-y0(ib)
        zs(ib)=z0(ib1)-z0(ib)
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
        aa=xa*xs(ib2)+ya*ys(ib2)+za*zs(ib2)
        chir(ib)=aa/bond**3
       enddo
      enddo
      
      enechi=0.d0
      do k=1,nch
       i1=ichain(k,1)
       i2=ichain(k,2)
       do ib=i1,i2-3
        ib1=ib+1
        ib2=ib+2
        ib3=ib+3
        cc=chir(ib)*chirn(ib)
        if (cc.le.0.d0) then
C COMPUTE ENERGY
         aa=(chir(ib)-chirn(ib))*(chir(ib)-chirn(ib))
         enechi=enechi+echi*0.5d0*aa
C COMPUTE FORCES
         fchi=-echi*(chir(ib)-chirn(ib))/bond**3
C FIRST RESIDUE
         repx=ys(ib2)*zs(ib1)-zs(ib2)*ys(ib1)
         repy=zs(ib2)*xs(ib1)-xs(ib2)*zs(ib1)
         repz=xs(ib2)*ys(ib1)-ys(ib2)*xs(ib1)
         fx(ib)=fx(ib)+fchi*repx
         fy(ib)=fy(ib)+fchi*repy
         fz(ib)=fz(ib)+fchi*repz
C SECOND RESIDUE
         xa=xs(ib)+xs(ib1)
         ya=ys(ib)+ys(ib1)
         za=zs(ib)+zs(ib1)
         repx=ya*zs(ib2)-za*ys(ib2)
         repy=za*xs(ib2)-xa*zs(ib2)
         repz=xa*ys(ib2)-ya*xs(ib2)
         fx(ib1)=fx(ib1)+fchi*repx
         fy(ib1)=fy(ib1)+fchi*repy
         fz(ib1)=fz(ib1)+fchi*repz
C THIRD RESIDUE
         xa=xs(ib1)+xs(ib2)
         ya=ys(ib1)+ys(ib2)
         za=zs(ib1)+zs(ib2)
         repx=ya*zs(ib)-za*ys(ib)
         repy=za*xs(ib)-xa*zs(ib)
         repz=xa*ys(ib)-ya*xs(ib)
         fx(ib2)=fx(ib2)+fchi*repx
         fy(ib2)=fy(ib2)+fchi*repy
         fz(ib2)=fz(ib2)+fchi*repz
C FOURTH RESIDUE
         repx=ys(ib)*zs(ib1)-zs(ib)*ys(ib1)
         repy=zs(ib)*xs(ib1)-xs(ib)*zs(ib1)
         repz=xs(ib)*ys(ib1)-ys(ib)*xs(ib1)
         fx(ib3)=fx(ib3)+fchi*repx
         fy(ib3)=fy(ib3)+fchi*repy
         fz(ib3)=fz(ib3)+fchi*repz
        endif
       enddo
      enddo
      return
      end