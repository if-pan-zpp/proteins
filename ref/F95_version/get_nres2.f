      subroutine get_nres2 (pdbfile, nres, nat)
C THIS SUBROUTINE READS THE PDB FILE AND FINDS 
C THE NUMBER OF RESIDUES (nres) AS WELL AS THE MAXIMAL 
C NUMBER OF HEAVY ATOMS IN A RESIDUE (nat) WHICH IS USUALLY 14
      implicit none
      character (len=*) pdbfile 	! PDB FILE
      integer nres		! NUMBER OR RESIDUES
      integer nat
      integer iat
      character*1 a
      character*80 line
      open (20, file=trim(pdbfile),form='formatted',status='old')
      
      nres = 0
      nat = 0
      iat = 0
 1    read (20, '(a)', end = 100, err = 100) line
      if (line(1:4) .eq. 'ATOM') then
       if (line (14:16) .eq. 'CA  ') then
         nres = nres + 1
    !     if (iat.gt.nat) then
    !      nat = iat
    !     endif
    !     iat = 0
       endif
    !  a=line(14:14)
    !  if (a.eq.'N' .or. a.eq.'S' .or. a.eq.'O' .or. a.eq. 'C') then
    !   iat = iat + 1
    !  endif
      endif
      goto 1
 100  close (20)
      if (nat.lt.14 ) then
        nat=14
      endif
      return
      end
