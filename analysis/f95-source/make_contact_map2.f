      subroutine make_contact_map2(nres,xn,yn,zn,rcm,
     $           icmap,kcont,isi,kk_map,seq,ch,iresn)
      implicit none
      integer nres
      real*8 xn(nres),yn(nres),zn(nres)
      integer iresn(nres)
      character*3 seq(nres)
      character*1 ch(nres)
      real*8 dca(nres,nres)
      real*8 rcm,rcm2,r2
      integer icmap(nres,nres)
      integer kcont,isi,kk_map,i1,i2
      
      rcm2=rcm*rcm
      do i1=1,nres
       do i2=1,nres
        icmap(i1,i2)=0
        dca(i1,i2)=0.d0
       enddo
      enddo
      do i1=1,nres-isi
       do i2=i1+isi,nres
        r2=(xn(i1)-xn(i2))**2
     $    +(yn(i1)-yn(i2))**2
     $    +(zn(i1)-zn(i2))**2
        if (r2.le.rcm2) then
         icmap(i1,i2)=1
         dca(i1,i2)=sqrt(r2)
        endif
       enddo
      enddo
      
      write(kk_map,'(/,a13,i1)') 'CONTACT MAP M',isi
      write(kk_map,'(a24,f8.3,/)') 'C-ALPHA DISTANCE CUT-OFF',rcm
!      write(kk_map,'(10x,a,13x,a,10x,a)') 'I1  AA1','I2  AA2',
!     $  'C-ALPHA DISTANCE'
      write(kk_map,'(12x,a)')
     $ 'I1  AA  C I(PDB)    I2  AA  C I(PDB) C-ALPHA DISTANCE'
      write(kk_map,'(a,a)')
     $ '==================================',
     $ '==============================='
                           

      kcont=0
      do i1=1,nres
       do i2=1,nres
        if (icmap(i1,i2).eq.1) then
         kcont=kcont+1
         write (kk_map, '(a1,1x,i6,2x,i4,2x,a3,1x,a1,1x,i3, 5x, 
     $                    i4,2x,a3,1x,a1,1x,i3, 5x, f8.4)' ) 
     $         'M',kcont,i1,seq(i1),ch(i1),iresn(i1),
     $         i2,seq(i2),ch(i2),iresn(i2),dca(i1,i2)
        endif
       enddo
      enddo

      return
      end
