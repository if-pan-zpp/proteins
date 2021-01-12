      subroutine make_contact_map(pdb_name,nres,nat,icmap,alpha,
     $           kcont,isi,kk_map,iresn)
C THIS SUBROUTINE CALCULATES THE CONTACT MAP
      implicit none
      character (len=*) pdb_name
      integer kk_map
      integer nres			! NUMBER OF RESIDUES
      integer nat			! MAX NUMBER OF ATOMS IN A RESIDUE
      integer kcont			! NUMBER OF CONTACTS
      integer isi 			! CONTACT MAP TYPE 
      integer icmap(nres, nres)		! CONTACT MAP
      character*3 seq(nres)		! AMINO ACID SEQUENCE
      integer iresn(nres)		! RESIDUR NUMBERS IN PDB FILE
      integer na(nres)			! NUMBER OF ATOMS IN RESIDUES
      character*1 at(nres, nat)		! HEAVY ATOMS N, S, O, C
      character*1 ch(nres)		! CHAIN SYMBOL
      real*8 x(nres, nat)		! X-COORDINATES OF ATOMS
      real*8 y(nres, nat)		! Y-COORDINATES OF ATOMS
      real*8 z(nres, nat)		! Z-COORDINATES OF ATOMS
      real*8 vrad(nres, nat)		! VAN DER WAALS RADII
      real*8 dca(nres, nres)		! C-ALPHA DISTANCES
      integer i1,j1,i2,j2,in
      real*8 r2, rmin, dist, rc
      character*80 outputfile
      real*8 alpha
C VAN DER WAALS RADII ARE SCALED BY FACTOR ALPHA

      call read_all_atom(pdb_name,nres,nat,seq,na,x,y,z,
     $     at,ch,iresn,kk_map)
      
      call assign_vdw_radii(pdb_name,nres,nat,seq,na,at,vrad,kk_map)    
      
      do i1=1,nres
       do i2=1,nres
        icmap(i1,i2) = 0
       enddo
      enddo
      do i1=1,nres-isi
       do i2=i1+isi,nres
        do j1=1,na(i1)
         do j2=1,na(i2)
          r2=(x(i1,j1)-x(i2,j2))**2
     $      +(y(i1,j1)-y(i2,j2))**2
     $      +(z(i1,j1)-z(i2,j2))**2
          dist=sqrt(r2)
          rc=vrad(i1,j1)+vrad(i2,j2)
          if (dist.lt.rc*alpha) then
           icmap(i1,i2) = 1
          endif
         enddo
        enddo
       enddo
      enddo

      do i1=1,nres
       do i2=1,nres
C C-ALPHA IS GIVEN AS THE SECOND ATOM OF A RESIDUES IN PDB
        j1=2
        j2=2
        r2=(x(i1,j1)-x(i2,j2))**2
     $    +(y(i1,j1)-y(i2,j2))**2
     $    +(z(i1,j1)-z(i2,j2))**2
        dist=sqrt(r2)
        dca(i1,i2)=dist
       enddo
      enddo
C
C WRITE OUT THE CONTACT MAP
C
      write(kk_map,'(/,a13,i1,/)') 'CONTACT MAP M',isi
      write(kk_map,'(12x,a)') 
     $ 'I1  AA  C I(PDB)    I2  AA  C I(PDB) C-ALPHA DISTANCE'
      write(kk_map,'(a,a)') 
     $ '==================================',
     $ '==============================='

      kcont=0
      do i1 = 1, nres
       do i2 = 1, nres
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
