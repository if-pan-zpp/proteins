      subroutine assign_vdw_radii (pdb_name,nres,nat,seq,na,at,
     $           vrad,kk_map)
C THIS SUBROUTINE ASSIGNS VAN DER WAALS RADII
      character (len=*) pdb_name
      integer kk_map
      integer nres		! MAX NUMBER OF RESIDUES
      integer nat		! MAX NUMBER OF ATOMS IN A RESIDUE
      character*3 seq(nres)	! AMINO ACID SEQUENCE
      integer na(nres)		! NUMBER OF ATOMS IN RESIDUES
      character*1 at(nres, nat)	! HEAVY ATOMS N, S, O, C
      real*8 vrad(nres, nat)	! VAN DER WAALS RADII
      integer ib, iat, j
      integer nb(nres)
      character*3 ares
      character*1 ana
      character*80 outputfile
      real*8 rad
      do ib = 1, nres
       vrad(ib,1)=1.64d0           ! BACKBONE NITROGEN
       vrad(ib,2)=1.88d0           ! BACKBONE C-ALPHA
       vrad(ib,3)=1.61d0           ! BACKBONE C'
       vrad(ib,4)=1.42d0           ! BACKBONE OXYGEN
       ares=seq(ib)
       if(ares.eq.'GLY') then
        nb(ib)=4
       else if(ares.eq.'ALA') then
        nb(ib)=5
        vrad(ib,5)=1.88d0
       else if(ares.eq.'VAL') then
        nb(ib)=7
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.88d0
       else if(ares.eq.'PHE') then
        nb(ib)=11
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.61d0
        vrad(ib,8)=1.76d0
        vrad(ib,9)=1.76d0
        vrad(ib,10)=1.76d0
        vrad(ib,11)=1.76d0
       else if(ares.eq.'PRO') then
        nb(ib)=7
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.88d0
       else if(ares.eq.'MET') then
        nb(ib)=8
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.77d0
        vrad(ib,8)=1.88d0
       else if(ares.eq.'ILE'.or.ares.eq.'LEU') then
        nb(ib)=8
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.88d0
        vrad(ib,8)=1.88d0
       else if(ares.eq.'ASP') then
        nb(ib)=8
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.61d0
        vrad(ib,7)=1.46d0
        vrad(ib,8)=1.42d0
       else if(ares.eq.'GLU') then
        nb(ib)=9
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.61d0
        vrad(ib,8)=1.46d0
        vrad(ib,9)=1.42d0
       else if(ares.eq.'LYS') then
        nb(ib)=9
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.88d0
        vrad(ib,8)=1.88d0
        vrad(ib,9)=1.64d0
       else if(ares.eq.'ARG') then
        nb(ib)=11
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.88d0
        vrad(ib,8)=1.64d0
        vrad(ib,9)=1.61d0
        vrad(ib,10)=1.64d0
        vrad(ib,11)=1.64d0
       else if(ares.eq.'SER') then
        nb(ib)=6
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.46d0
       else if(ares.eq.'THR') then
        nb(ib)=7
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.46d0
        vrad(ib,7)=1.88d0
       else if(ares.eq.'TYR') then
        nb(ib)=12
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.61d0
        vrad(ib,7)=1.76d0
        vrad(ib,8)=1.76d0
        vrad(ib,9)=1.76d0
        vrad(ib,10)=1.76d0
        vrad(ib,11)=1.61d0
        vrad(ib,12)=1.46d0
       else if(ares.eq.'HIS') then
        nb(ib)=10
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.61d0
        vrad(ib,7)=1.64d0
        vrad(ib,8)=1.76d0
        vrad(ib,9)=1.76d0
        vrad(ib,10)=1.64d0
       else if(ares.eq.'CYS') then
        nb(ib)=6
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.77d0
       else if(ares.eq.'ASN') then
        nb(ib)=8
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.61d0
        vrad(ib,7)=1.42d0
        vrad(ib,8)=1.64d0
       else if(ares.eq.'GLN') then
        nb(ib)=9
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.88d0
        vrad(ib,7)=1.61d0
        vrad(ib,8)=1.42d0
        vrad(ib,9)=1.64d0
       else if(ares.eq.'TRP') then
        nb(ib)=14
        vrad(ib,5)=1.88d0
        vrad(ib,6)=1.61d0
        vrad(ib,7)=1.76d0
        vrad(ib,8)=1.61d0
        vrad(ib,9)=1.64d0
        vrad(ib,10)=1.61d0
        vrad(ib,11)=1.76d0
        vrad(ib,12)=1.76d0
        vrad(ib,13)=1.76d0
        vrad(ib,14)=1.76d0
       endif
      enddo
C CHECK AND CORRECTION
      write(kk_map,'(a)') 'VAN DER WAALS RADII ASSIGNED TO  N, C, O, S'
      do ib = 1, nres
       iat=na(ib)
       if(iat.lt.nb(ib)) then
        write(kk_map,'(a)')'AMINO ACID HAS FEWER ATOMS THAN SHOULD BE:'
        write(kk_map,'(i5,2x,a3,2x,2i6)') ib,seq(ib),iat,nb(ib)
       else if(iat.gt.nb(ib)) then
        write(kk_map,'(a)')'AMINO ACID HAS MORE ATOMS THAN SHOULD BE:'
        write(kk_map,'(i5,2x,a3,2x,2i6)') ib,seq(ib),iat,nb(ib)
c         do j = 1, iat
c          write(1,'(a3,2x,i3,2x,a3)') seq(ib),iat,at(ib,j)
c         enddo
        iat=nb(ib)
        na(ib)=nb(ib)
       endif
       do j=1,iat
        ares=seq(ib)
        ana=at(ib,j)
        rad=vrad(ib,j)
        if(ana.eq.'N'.and.rad.ne.1.64d0) then
         vrad(ib,j)=1.64d0
        else if(ana.eq.'S'.and.rad.ne.1.77d0) then
         vrad(ib,j)=1.77d0
        else if(ana.eq.'O'.and.rad.ne.1.42d0.and.rad.ne.1.46d0) then
         vrad(ib,j)=1.46d0
        else if(ana.eq.'C'.and.rad.ne.1.88d0.and.rad.ne.1.76d0.
     $          and.rad.ne.1.61d0) then
         vrad(ib,j)=1.88d0
        endif
        if(rad.eq.0) then
         write(kk_map,'(/,a,/)') 'ATOM ERROR:'
         write(kk_map,'(a3,2x,i3,2x,a3)') ib,seq(ib),at(ib,j)
         stop
        endif
       enddo
      enddo
      return
      end

