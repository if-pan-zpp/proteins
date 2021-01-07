      subroutine read_all_atom (pdbfile,nres,nat,seq,na,x,y,z,
     $           at,ch,iresn,kk_map)
C THIS SURROUTINE READS IN THE STRUCTURE FROM THE PDB FILE
C IT IS USED IN THE CONTACT MAP CALCULATION
      implicit none
      character (len=*) pdbfile ! PDB FILE
      integer kk_map
      integer nres		! NUMBER OF RESIDUES
      integer nat		! MAX NUMBER OF ATOMS IN A RESIDUE
      character*3 seq(nres)	! AMINO ACID SEQUENCE
      integer na(nres)		! NUMBER OF ATOMS IN RESIDUES
      character*1 at(nres, nat)	! HEAVY ATOMS N, S, O, C
      character*1 ch(nres)	! CHAIN SYMBOL
      real*8 x(nres, nat)	! X-COORDINATES OF ATOMS
      real*8 y(nres, nat)	! Y-COORDINATES OF ATOMS
      real*8 z(nres, nat)	! Z-COORDINATES OF ATOMS
      integer iresn(nres)	! RESIDUR NUMBERS IN PDB FILE
      logical struc_error	! PROBLEMS WITH PDB FILE
      integer ires,iat
      integer i,j
      character*80 line, outputfile
      character*1 atm
C OPEN PDB FILE

      open(20, file=trim(pdbfile),
     $         form='formatted', status='old')
      ires = 0
 1    read(20, '(a)', end = 100, err = 100) line
      if (line(1:4) .eq. 'ATOM' ) then
       if (line (14:16) .eq. 'N  ' ) then
        ires = ires + 1
        read( line(18:20), '(a3)' ) seq(ires)
        read( line(22:22), '(a1)' ) ch(ires)
        read( line(24:26), '(i3)' ) iresn(ires)
        iat = 0
       endif
       read( line(14:14), '(a1)' ) atm
       if (atm.eq.'N' .or. atm.eq.'S' .or. 
     $     atm.eq.'O' .or. atm.eq.'C' ) then
        iat = iat + 1
        na(ires) = iat
        at(ires, iat) = atm
        read( line(31:54), '(3f8.3)' ) x(ires, iat), 
     $   y(ires, iat), z(ires, iat)
       endif
      endif
      goto 1
 100  close (20)
C CLOSE PDB FILE
C CHECK THE BACKBONE
    
      struc_error=.false.
      do i = 1, nres
       if (at(i,1).ne.'N') then
        write(kk_map,*)'BACKBONE PROBLEM: NO NITROGEN IN RESIDUE ',
     $   iresn(i)
        struc_error=.true.
       endif
       if (at(i,2).ne.'C') then
        write(kk_map,*) 'BACKBONE PROBLEM: NO C-ALPHA IN RESIDUE ',
     $   iresn(i)
        struc_error=.true.
       endif
       if (at(i,3).ne.'C') then
        write(kk_map,*) 'BACKBONE PROBLEM: NO CARBONE IN RESIDUE ',
     $   iresn(i)
        struc_error=.true.
       endif
       if (at(i,4).ne.'O') then
        write(kk_map,*) 'BACKBONE PROBLEM: NO OXYGEN IN RESIDUE ',
     $   iresn(i)
        struc_error=.true.
       endif
      enddo
C CHECK THE AMINO ACID SEQUENCE
      do i = 1, nres-1
       j=iresn(i+1)-iresn(i)
       if (j.ne.1 .and. ch(i+1).eq.ch(i)) then
        write(kk_map, *) 'MISSING RESIDUE ', iresn(i+1)
        struc_error=.true.
       endif
      enddo
      if (struc_error) then
       close(kk_map)
       stop
      else
       write(kk_map, '(/,a,a)') 'NO PROBLEMS DETECTED IN STRUCTURE  ',
     $   pdbfile
      endif
      return
      end 
