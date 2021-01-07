      subroutine get_initial_coordinates(pdbfile,nres,ntrj,x,y,z,
     $           unit,seq,kk_out)
C THIS SURROUTINE READS IN C-ALPHA COORDINTES FROM THE PDB FILE
      implicit none
      character (len=*) pdbfile ! PDB FILE
      integer nres		! NUMBER OF RESIDUES
      integer ntrj		! NUMBER OF STRUCTURES
      integer kk_out
      real*8 x(nres,ntrj)
      real*8 y(nres,ntrj)
      real*8 z(nres,ntrj)
      character*3 seq(nres),seq1(nres)
      integer ires,imod
      real*8 unit
      character*80 line
C CHECK THE NUMBER OF INITIAL UNFOLDED STRUCTURES
      open(20, file=trim(pdbfile),
     $         form='formatted', status='old')
      imod=0
 8    read(20, '(a)', end = 9, err = 9) line
       if (line(1:5) .eq. 'MODEL' ) then
        imod=imod+1
       endif
      goto 8
 9    close(20)
      if (imod.ne.ntrj) then
       write(kk_out,'(a)') 'PROBLEM WITH INITIAL UNFOLDED STRUCTURES!'
       write(kk_out,'(a,i6)') 'NUMBER OF INITIAL STRUCTURES: ', imod
       write(kk_out,'(a,i6)') 'NUMBER OF TRAJECTORIES: ',ntrj
       stop
      endif
C READ THE STRUCTURES IN
      imod=0
      open(20, file=trim(pdbfile),
     $         form='formatted', status='old')
 11   read(20, '(a)', end = 12, err = 12) line
      if (line(1:5) .eq. 'MODEL' ) then
       ires=0
       imod=imod+1
      endif
      if (line(1:4) .eq. 'ATOM' ) then
        if (line (14:16) .eq. 'CA ' ) then
          ires=ires+1
          read( line(18:20), '(a3)' ) seq1(ires)
          read( line(31:54), '(3f8.3)' ) x(ires,imod),
     $                                   y(ires,imod),
     $                                   z(ires,imod)
        endif
      endif
      if (line(1:5) .eq. 'ENDMD') then
C CHECK THE NUMBER OF RESIDUES
       if (ires.ne.nres) then
        write(kk_out,'(a)') 'PROBLEM WITH INITIAL UNFOLDED STRUCTURES!'
        write(kk_out,'(a,i6)') 'NUMBER OF RESIDUES IN STRUCTURE ',imod
        stop
       endif
C CHECK THE SEQUENCE
       do ires=1,nres
        if (seq1(ires).ne.seq(ires)) then
         write(kk_out,'(a)') 'PROBLEM WITH UNFOLDED STRUCTURES!'
         write(kk_out,'(a)') 'WRONG SEQUENCE IN STRUCTURE ',imod
         stop
        endif
       enddo
      endif
      goto 11
C CLOSE PDB FILE      
 12   close (20)
      do imod=1,ntrj
       do ires=1,nres
        x(ires,imod)=x(ires,imod)/unit
        y(ires,imod)=y(ires,imod)/unit
        z(ires,imod)=z(ires,imod)/unit
       enddo
      enddo
      return
      end
