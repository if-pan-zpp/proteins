      subroutine write_pdb(iun,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
C THIS SUBROUTINE PRINTS THE CHAIN COORDINATES IN PDB FORMAT
      implicit none
      integer iun		! SPECIFICATION OF A LOGICAL UNIT FOR THE PDB FILE
      integer nres		! NUMBER OF RESIDUES
      real*8 x0(nres)
      real*8 y0(nres)
      real*8 z0(nres)		! COORDINATES
      character*3 seq(nres)	! SEQUENCE
      character*1 ch(nres)	! CHAIN SYMBOL
      integer iseq(nres)	! RESIDUE NUMBER
      real*8 time,energy,rms
      real*8 xmin,ymin,zmin,unit
      integer ib,iatom,imodel
      xmin=0.d0
      ymin=0.d0
      zmin=0.d0
      do ib=1,nres
        xmin=xmin+x0(ib)
        ymin=ymin+y0(ib)
        zmin=zmin+z0(ib)
      enddo
      xmin=xmin/nres
      ymin=ymin/nres
      zmin=zmin/nres
      
C      write(iun,'(/,a,f10.1,4x,a,f10.4,4x,a,f8.4)')
C     +  'REMARK    TIME =',time,'ENERGY =',energy,'RMSD =',rms

      write(iun,'(a5,i5)') 'MODEL',imodel
      iatom=0
1200  format(a4,i7,2x,a4,a3,1x,a1,i4,4x,3f8.3)
      do ib=1,nres
       iatom=iatom+1
       write(iun,1200) 'ATOM',iatom,'CA  ',seq(ib),ch(ib),iseq(ib),
     $  (x0(ib)-xmin)*unit,(y0(ib)-ymin)*unit,(z0(ib)-zmin)*unit
       if (ib.lt.nres) then
        if (ch(ib).ne.ch(ib+1)) then
         write(iun,'(a3)') 'TER'
        endif
       endif
      enddo
      write(iun,'(a3)') 'TER'
      write(iun,'(a5)') 'ENDMDL'
      return
      end
