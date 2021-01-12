      program Go_sim
      implicit none

C DECLARATION OF VARIABLES
      
      integer i_arg,rr,kk_out,kk_pdb,kk_map,kk_knot
      integer nres,nat,kcont,nch,iseed,ntrj,kupdate
      integer i,j,in,isi,nz,ip1,ip2,itbs,knumbers,ktight
      integer nssb1,nssb2,nssb3
      
      integer, dimension(:,:), allocatable :: icmap
      integer, dimension(:,:), allocatable :: ichain
      integer, dimension(:), allocatable :: iseq,iresn
      real*8, dimension(:,:), allocatable :: eij
      real*8, dimension(:), allocatable :: xn
      real*8, dimension(:), allocatable :: yn
      real*8, dimension(:), allocatable :: zn
      real*8, dimension(:,:), allocatable :: xu
      real*8, dimension(:,:), allocatable :: yu
      real*8, dimension(:,:), allocatable :: zu
      real*8, dimension(:,:), allocatable :: wz
      real*8, dimension(:,:), allocatable :: wp
      real*8, dimension(:,:), allocatable :: wf
      character*3, dimension(:), allocatable :: seq
      character*1, dimension(:), allocatable :: ch
      
      character*80 arg,desc,pdb_name,output_name,outputfile
      character*80 pot_file1,pot_file2,unf_name,cont_map_name
      character*1 snt,sct,knot_ask
      character*6 snch,scch
      real*8 alpha,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta
      real*8 tsim,tsav,tpdb,teql,tfave,boxl,vpull,velo,HH1,HH2,tbs
      real*8 forcon,c_param,a_param,tstop,rcm,s612,rcut3
      real*8 r2,ran2,unit,max_pull_f,cutparam,tstart,tend,tstep,tksav
      logical lpc,llf,lbd,lchir,langl,l612,l612s,l1012,l61012
      logical lthrm,lfold,lpull,lsurf,lmj,stop_pulling,c_velo,c_for
      logical lcm1,lcm2,lcm3,stop_pulling2
      
      rr=1 ! SPECIFICATION OF A LOGICAL UNIT FOR A FILE
      
C GET THE INPUT FILE NAME
      
      do i_arg=1,iargc()
         call getarg(i_arg,arg)
      enddo
      
C DEFAULT PARAMETERS

      alpha=1.244455099105835d0
      unit=5.d0
      gamma=2.d0
      delta=0.005d0
      temp=0.3d0
      tstart=0.3d0
      tend=0.3d0
      tstep=0.d0
      H1=50.d0
      H2=0.d0
      HH1=0.06d0
      HH2=0.0d0
      ktight=40
      echi=1.d0
      kphi1=1.d0
      kphi2=0.5d0
      ktheta=20.d0
      velo=0.005d0
      forcon=1.0d0
      teql=100
      tsav=100
      tpdb=1000
      tsim=30000
      tstop=0d0
      tbs=100
      ntrj=5
      kupdate=200
      iseed=8140
      isi=3
      boxl=50.d0
      max_pull_f=25.0d0
      cutparam=1.5d0
      snt='A'
      sct='A'
      snch='N'
      scch='C'
      knot_ask='N'
      tksav=500
      c_param=0.85d0
      a_param=3.8d0
      rcm=7.5d0
      rcut3=4.0d0
      
      c_velo=.false.
      c_for=.false.
      lchir=.false.
      langl=.false.
      l612=.false.
      l612s=.false.
      l1012=.false.
      l61012=.false.
      lpc=.false.
      llf=.false.
      lbd=.false.
      lthrm=.false.
      lfold=.false.
      lpull=.false.
      lsurf=.false.
      lmj=.false.
      lcm1=.false.
      lcm2=.false.
      lcm3=.false.
      stop_pulling=.false.
      stop_pulling2=.false.
      unf_name='none'
      output_name='none'
      cont_map_name='none'
      
      
      
C READING PARAMETERS FROM THE INPUT FILE

      open(rr,file=trim(arg),status='old') 
 11     read(rr,'(a)',end=12, err=12) desc
        if (desc(1:28) .eq. 'all_atom_overlap_contact_map' ) then
         lcm1=.true.
        endif
        if (desc(1:35) .eq. 'contact_map_C_alpha_distane_cutoff:' ) then
         lcm2=.true.
         read(desc(36:), *) rcm
        endif
        if (desc(1:14) .eq. 'thermodynamics' ) then
         lthrm=.true.
        endif
        if (desc(1:9) .eq. 'pulling  ' ) then
         lpull=.true.
        endif
        if (desc(1:17) .eq. 'constant_velocity' ) then
         c_velo=.true.
        endif
        if (desc(1:14) .eq. 'constant_force' ) then
         c_for=.true.
        endif        
        if (desc(1:6) .eq. 'iseed:' ) then
         read(desc(7: ), *) iseed
        endif
        if (desc(1:23) .eq. 'contact_length_cut-off:' ) then
         read(desc(24: ), *) cutparam
        endif
        if (desc(1:18) .eq. 'max_pulling_force:' ) then
          read(desc(19: ), *) max_pull_f
        endif
        if (desc(1:23) .eq. 'stop_pulling_when_ICN=0' ) then
         stop_pulling=.true.
        endif
        if (desc(1:34) .eq. 'stop_pulling_when_D(1,N)>c*A*(N-1)' ) then
         stop_pulling2=.true.
        endif
        if (desc(1:7) .eq. 'folding' ) then
         lfold=.true.
        endif
        if (desc(1:20) .eq. 'unfolded_structures:' ) then
         read(desc(21: ), *) unf_name
        endif
        if (desc(1:15) .eq. 'chirality_model' ) then
         lchir=.true.
        endif
        if (desc(1:13) .eq. 'angular_model' ) then
         langl=.true.
        endif
        if (desc(1:14) .eq. '6_12_potential' ) then
         l612=.true.
        endif
        if (desc(1:28) .eq. '6_12_type_potential_shifted:' ) then
         l612s=.true.
         read(desc(29: ), *) s612
        endif
        if (desc(1:15) .eq. '10_12_potential' ) then
         l1012=.true.
        endif
        if (desc(1:17) .eq. '6_10_12_potential' ) then
         l61012=.true.
        endif
        if (desc(1:28) .eq. 'langevin_predictor_corrector' ) then
         lpc=.true.
        endif
        if (desc(1:18) .eq. 'langevin_leap_frog' ) then
         llf=.true.
        endif
        if (desc(1:27) .eq. 'brownian_dynamics_leap_frog' ) then
         lbd=.true.
        endif
        if (desc(1:14) .eq. 'PDB_file_name:' ) then
         read(desc(15: ), *) pdb_name
        endif
        if (desc(1:17) .eq. 'output_file_name:' ) then
         read(desc(18: ), *) output_name
        endif
        if (desc(1:22) .eq. 'contact_map_file_name:' ) then
         read(desc(23: ), *) cont_map_name
        endif
        if (desc(1:21) .eq. 'contact_map_|i-j|.ge.' ) then
         read(desc(22: ), *) isi
        endif
        if (desc(1:26) .eq. 'non-native_length_cut-off:' ) then
         read(desc(27: ), *) rcut3
        endif
        if (desc(1:16) .eq. 'alpha_parameter:' ) then
         read(desc(17: ), *) alpha
        endif
        if (desc(1:25) .eq. 'friction_parameter_gamma:' ) then
         read(desc(26: ), *) gamma
        endif
        if (desc(1:23) .eq. 'stop_criterion_c_param:' ) then
         read(desc(24: ), *) c_param
        endif
        if (desc(1:23) .eq. 'stop_criterion_A_param:' ) then
         read(desc(24: ), *) a_param
        endif
        if (desc(1:12) .eq. 'length_unit:' ) then
         read(desc(13: ), *) unit
        endif
        if (desc(1:16) .eq. 'time_step_delta:' ) then
         read(desc(17: ), *) delta
        endif
        if (desc(1:12) .eq. 'temperature:' ) then
         read(desc(13: ), *) temp
        endif
        if (desc(1:18) .eq. 'start_temperature:' ) then
         read(desc(19: ), *) tstart
        endif
        if (desc(1:16) .eq. 'end_temperature:' ) then
         read(desc(17: ), *) tend
        endif
        if (desc(1:38) .eq. 'knot_is_tightened_when_consists_of_AA:' ) 
     $   then
         read(desc(39: ), *) ktight
        endif
        if (desc(1:17) .eq. 'temperature_step:' ) then
         read(desc(18: ), *) tstep
        endif
        if (desc(1:18) .eq. 'harmonic_coeff_H1:' ) then
         read(desc(19: ), *) H1
        endif
        if (desc(1:20) .eq. 'anharmonic_coeff_H2:' ) then
         read(desc(21: ), *) H2
        endif
        if (desc(1:20) .eq. 'chirality_parameter:' ) then
         read(desc(21: ), *) echi
        endif
        if (desc(1:21) .eq. 'torsion_potential_K1:' ) then
         read(desc(22: ), *) kphi1
        endif
        if (desc(1:21) .eq. 'torsion_potential_K2:' ) then
         read(desc(22: ), *) kphi2
        endif
        if (desc(1:23) .eq. 'bond_angle_potential_K:' ) then
         read(desc(24: ), *) ktheta
        endif                
        if (desc(1:20) .eq. 'data_save_intervals:' ) then
         read(desc(21: ), *) tsav
        endif
        if (desc(1:25) .eq. 'structure_save_intervals:' ) then
         read(desc(26: ), *) tpdb
        endif
        if (desc(1:29) .eq. 'knot_position_save_intervals:' ) then
         read(desc(30: ), *) tksav
        endif
        if (desc(1:19) .eq. 'equilibration_time:' ) then
         read(desc(20: ), *) teql
        endif
        if (desc(1:16) .eq. 'simulation_time:' ) then
         read(desc(17: ), *) tsim
        endif
        if (desc(1:24) .eq. 'stop_pulling_after_time:' ) then
         read(desc(25: ), *) tstop
        endif                         
        if (desc(1:23) .eq. 'number_of_trajectories:' ) then
         read(desc(24: ), *) ntrj
        endif
        if (desc(1:14) .eq. 'pulling_speed:' ) then
         read(desc(15: ), *) velo
        endif
        if (desc(1:14) .eq. 'pulling_force:' ) then
         read(desc(15: ), *) forcon
        endif
        if (desc(1:28) .eq. 'pulling_spring_constant_HH1:' ) then
         read(desc(29: ), *) HH1
        endif
        if (desc(1:28) .eq. 'pulling_spring_constant_HH2:' ) then
         read(desc(29: ), *) HH2
        endif
        if (desc(1:37).eq.'pull_on_[N/AA_No.]_terminus_of_chain:' ) then
         read(desc(38:39), *) snt
         read(desc(40: ),*) snch
        endif
        if (desc(1:24).eq.'knots_examination_[Y/N]:' ) then
         read(desc(25:26), *) knot_ask
        endif
        if (desc(1:37).eq.'pull_on_[C/AA_No.]_terminus_of_chain:' ) then
         read(desc(38:39), *) sct
         read(desc(40: ),*) scch
        endif
        if (desc(1:24) .eq. 'force_average_intervals:' ) then
         read(desc(25: ), *) tfave
        endif
        if (desc(1:9) .eq. 'box_size:' ) then
         read(desc(10: ), *) boxl
        endif
        if (desc(1:17) .eq. 'bound_state_time:' ) then
         read(desc(18: ), *) tbs
        endif
        if (desc(1:23) .eq. 'surface_potential_file:' ) then
         read(desc(24: ), *) pot_file2
         lsurf=.true.
        endif
        if (desc(1:18) .eq. 'MJ_potential_file:' ) then
         lmj=.true.
         read(desc(19: ), *) pot_file1
        endif
        goto 11
 12   close(rr)
      
C OPEN THE CONTACT MAP FILE   
     
      kk_map=13
      if (output_name.eq.'none') then
       in=index(pdb_name,'.')
       if (in.gt.0) then
        outputfile=pdb_name(1:in-1)//'.map'
       else
        outputfile=pdb_name//'.map'
       endif
      else
       in=index(output_name,' ')
       outputfile=output_name(1:in-1)//'.map'
      endif
      
      if (cont_map_name.eq.'none') then
      open(kk_map,file=trim(outputfile),
     $            form='formatted',status='unknown')
      endif

      r2=ran2(iseed)
      
      if (cont_map_name.eq.'none') then
      call get_nres(trim(pdb_name),nres,nat)
      else
      call get_nres2(trim(pdb_name),nres,nat) 
      endif
      
      allocate(icmap(nres,nres))      
      allocate(xn(nres))
      allocate(yn(nres))
      allocate(zn(nres))
      allocate(seq(nres))
      allocate(iseq(nres))
      allocate(iresn(nres))
      allocate(ch(nres))     
      allocate(eij(nres,nres))
      allocate(xu(nres,ntrj))
      allocate(yu(nres,ntrj))
      allocate(zu(nres,ntrj))
 
      call get_c_alpha_coordinates(trim(pdb_name),nres,xn,yn,zn,
     $     seq,iseq,ch,nch)
      
      
      allocate(ichain(nch,2))

       
      if (cont_map_name.eq.'none') then
      if (lcm1) then
       call make_contact_map(trim(pdb_name),nres,nat,icmap,alpha,
     $      kcont,isi,kk_map,iresn)
      else
       if (lcm2) then
        call make_contact_map2(nres,xn,yn,zn,rcm,
     $       icmap,kcont,isi,kk_map,seq,ch,iseq)
       endif
      endif
      else
      call read_contact_map(cont_map_name,nres,icmap,kcont)
      endif

      kk_out=10
      if (output_name.eq.'none') then
       in=index(pdb_name,'.')
       if (in.gt.0) then
        outputfile=pdb_name(1:in-1)//'.out'
       else
        outputfile=pdb_name//'.out'
       endif
      else
       in=index(output_name,' ')
       outputfile=output_name(1:in-1)//'.out'
      endif      

      
      call get_chains(nres,ch,nch,ichain)

      if (cont_map_name.eq.'none') then

      call get_add_bonds_numb(trim(pdb_name),nres,iseq,seq,ch,icmap,
     $    kcont,kk_map,rr,arg,outputfile,kk_out,iresn,nssb1,nssb2,nssb3)
      
      call get_ss_bonds(trim(pdb_name),nres,iseq,seq,ch,icmap,
     $     kcont,kk_map)
   
      call get_add_bonds(trim(pdb_name),nres,iseq,seq,ch,icmap,
     $     kcont,kk_map,rr,arg,outputfile,kk_out,iseq,nssb1)
     
      call get_add_contacts(trim(pdb_name),nres,iseq,seq,ch,icmap,
     $     kcont,kk_map,rr,arg,outputfile,kk_out,iseq,nssb2)
     
      call get_rem_contacts(trim(pdb_name),nres,iseq,seq,ch,icmap,
     $     kcont,kk_map,rr,arg,outputfile,kk_out,iseq,nssb3)
     
      endif
         
      close(kk_map)
     
      call get_knot_number(trim(pdb_name),nres,iseq,seq,ch,icmap,
     $     kcont,kk_map,rr,arg,knumbers)
            
      if (lsurf) then
       call get_nz(trim(pot_file2),nz)
       allocate(wz(nres,nz))
       allocate(wp(nres,nz))
       allocate(wf(nres,nz))
      else
       allocate(wz(nres,1))
       allocate(wp(nres,1))
       allocate(wf(nres,1))
      endif
      
C SCALE

      do i=1,nres
       xn(i)=xn(i)/unit
       yn(i)=yn(i)/unit
       zn(i)=zn(i)/unit
      enddo
      H1=H1*unit*unit
      H2=H2*unit**4
      HH1=HH1*unit*unit
      HH2=HH2*unit**4
      vpull=velo/unit
      boxl=boxl/unit
      forcon=forcon*unit
      s612=s612/unit
      
C OPEN OUTPUT FILE

      open(kk_out,file=trim(outputfile),
     $            form='formatted',status='unknown')
      
      write(kk_out, '(a10,a,/)') 'PDB FILE: ', pdb_name
      if (cont_map_name.eq.'none') then
      write(kk_out, '(a14,i1,/)') 'CONTACT MAP: M',isi
      else
      write(kk_out, '(a32,a,/)') 'CONTACT MAP READ FROM THE FILE: ',
     $ trim(cont_map_name)
      endif
      if (lchir) then
       write(kk_out,'(a)') 'CHIRALITY MODEL'
       write(kk_out,'(a32,f8.4)') 'CHIRALITY POTENTIAL PARAMETER:  ',
     $   echi
      endif
      if (langl) then
       write(kk_out,'(a)') 'BOND-ANGLE AND DIHEDRAL-ANGLE POTENTIAL'
       write(kk_out,'(a8,f8.4)') 'K_THETA:  ', ktheta
       write(kk_out,'(a8,f8.4)') 'K1_PHI:   ', kphi1
       write(kk_out,'(a8,f8.4)') 'K2_PHI:   ', kphi2
      endif
      if (l612) then
       write(kk_out,'(/,a,/)') '12-6 POTENTIAL'
      endif
      if (l612s) then
       write(kk_out,'(/,a,f8.3,/)') '12-6-TYPE POTENTIAL WITH s =',
     $       s612*unit
      endif
      if (l1012) then
       write(kk_out,'(/,a,/)') '12-10 POTENTIAL'
      endif
      if (l61012) then
       write(kk_out,'(/,a,/)') '12-10-6 POTENTIAL'
      endif
      
      if (lmj) then
       call get_eij(trim(pot_file1),nres,icmap,seq,eij)
       write(kk_out,'(a,a,/)') 'POTENTIAL DEPTH E_ij: ',pot_file1
      else
       do i=1,nres
        do j=1,nres
         if (icmap(i,j).eq.1) then
          eij(i,j)=1.d0
         else
          eij(i,j)=0.d0
         endif
        enddo
       enddo
       write(kk_out,'(a,/)') 'UNIFORM POTENTIAL DEPTH: E_ij=1'
      endif
      
      write(kk_out,'(a28,f12.4)') 'BOND HARMONIC COEFFICIENT:  ', 
     $      H1/unit**2
      write(kk_out,'(a28,f12.4)') 'BOND AHARMONIC COEFFICIENT: ',
     $      H2/unit**4
      
      if (lpc) then
       write(kk_out,'(/,a)') 'LANGEVIN DYNAMICS'
       write(kk_out,'(a,/)') 'PREDICTOR-CORRECTOR ALGORITHM'
      endif
      if (llf) then
       write(kk_out,'(/,a)') 'LANGEVIN DYNAMICS'
       write(kk_out,'(a,/)') 'LEAP-FROG ALGORITHM'
      endif
      if (lbd) then
       write(kk_out,'(/,a)') 'BROWNIAN DYNAMICS'
       write(kk_out,'(a,/)') 'LEAP-FROG ALGORITHM'
      endif
      
      write(kk_out,'(a30,i12)') 'RANDOM NUMBER GENERATOR SEED  ',iseed
      write(kk_out,'(a24,f8.4)') 'INTEGRATION TIME STEP:  ', delta
      write(kk_out,'(a20,f8.4)') 'FRICTION PARAMETER: ', gamma
      write(kk_out,'(a13,f8.4)') 'LENGTH UNIT: ', unit
      
      if (lthrm) then
       if (tstep.eq.0.d0) then
        write(kk_out,'(a12,f8.4)') 'TEMPERATURE:', temp
       else
        write(kk_out,'(a13,f7.4,a3,f7.4,a13,f7.4)') 'TEMPERATURES ',
     $        tstart,' - ',tend,' WITH A STEP ',tstep        
       endif
      else
       write(kk_out,'(a12,f8.4)') 'TEMPERATURE:', temp
      endif
      
      if (lthrm) then
       write(kk_out,'(/,a)') 'THERMODYNAMICS'
      endif
      
      if (lfold) then
       write(kk_out,'(/,a)') 'FOLDING'
       if (unf_name.eq.'none') then
        call line_config(nres,ntrj,xu,yu,zu,nch,ichain,unit)
        write(kk_out,'(a)') 'STRAIGHT INITIAL CONFIGURATION'
       else
        call get_initial_coordinates(trim(unf_name),nres,ntrj,
     $       xu,yu,zu,unit,seq,kk_out)
        write(kk_out,'(a,a)') 'INITIAL UNFOLDED STRUCTURES: ',unf_name
       endif
      endif
      
      if (lpull) then
       write(kk_out,'(/,a)') 'PULLING'
       if (c_velo) then
       write(kk_out,'(a,f8.4)') 'PULLING SPEED: ',velo
       endif
       if (c_for) then
       write(kk_out,'(a,f8.4)') 'PULLING FORCE: ',forcon/unit
       endif
       write(kk_out,'(a,f8.4)') 'SPRING CONSTANT HH1: ',HH1/unit**2
       write(kk_out,'(a,f8.4)') 'SPRING CONSTANT HH2: ',HH2/unit**4
       call get_ip1_ip2(nres,ch,snt,sct,nch,ichain,ip1,ip2,kk_out,snch,
     $  scch)
       write(kk_out,'(a)') 'SPRINGS ATTACHED TO'
       write(kk_out,'(a,a,a,a1,a,i6)') '-  ',trim(snch),
     $  ' TERMINUS OF CHAIN '
     $  ,snt,', WHICH IS BEAD ',ip1
       write(kk_out,'(a,a,a,a1,a,i6)') '-  ',trim(scch),
     $  ' TERMINUS OF CHAIN '
     $  ,sct,', WHICH IS BEAD ',ip2
      endif

      if (lsurf) then
       write(kk_out,'(/,a)') 'PROTEIN-SURFACE INTERACTION'
       call get_surface_pot(trim(pot_file2),nres,seq,nz,wz,wp,wf,
     $      unit,temp)
       write(kk_out,'(a24,a)') 'SURFACE POTENTIAL FILE: ', pot_file2
       write(kk_out,'(a10,f8.2)') 'BOX SIZE: ', boxl*unit
       write(kk_out,'(a18,f8.2)') 'BOUND STATE TIME: ', tbs
       itbs=tbs/tsav
      endif
        call flush(kk_out)
C OPEN PDB OUTPUT FILE

      if ((lpull.and.(c_velo.and.c_for)).or.
     $ (lpull.and.(.not.c_velo.and.(.not.c_for)))) then 
       write(kk_out, '(a35/)') '***********************************'
       write(kk_out, '(a35/)') 'PULLING MANNER DEFINED INCORRECTLY!'
       stop
      endif

      kk_pdb=11
      if (output_name.eq.'none') then
       in=index(pdb_name,'.')
       if (in.gt.0) then
        outputfile=pdb_name(1:in-1)//'.sav'
       else
        outputfile=pdb_name//'.sav'
       endif
      else
       in=index(output_name,' ')
       outputfile=output_name(1:in-1)//'.sav'
      endif

      open(kk_pdb,file=trim(outputfile),
     $            form='formatted',status='unknown')
     

      if (knot_ask.eq.'Y') then
      kk_knot=15
      if (output_name.eq.'none') then
       in=index(pdb_name,'.')
       if (in.gt.0) then
        outputfile=pdb_name(1:in-1)//'.knot'
       else
        outputfile=pdb_name//'.knot'
       endif
      else
       in=index(output_name,' ')
       outputfile=output_name(1:in-1)//'.knot'
      endif

      open(kk_knot,file=trim(outputfile),
     $            form='formatted',status='unknown')
1223  format(a10,8x,a7,8x,a7/)
        write(kk_knot,1223) '#     TIME','N_LIMIT','C_LIMIT'
      endif
    
C MAIN LOOP
          
      if (lthrm) then
       call main_loop_th(nres,kcont,xn,yn,zn,seq,iseq,unit,
     $  ch,icmap,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta,isi,
     $  tsim,tsav,tpdb,teql,kupdate,ntrj,nch,ichain,kk_out,kk_pdb,
     $  eij,lchir,langl,l612,l612s,l1012,l61012,lpc,llf,lbd,
     $  tstart,tend,tstep,cutparam,s612,1,rcut3)
      endif
      
      if (lfold) then
       call main_loop_folding(nres,kcont,xn,yn,zn,seq,iseq,unit,
     $  ch,icmap,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta,isi,
     $  tsim,tsav,tpdb,kupdate,ntrj,nch,ichain,kk_out,kk_pdb,
     $  eij,lchir,langl,l612,l612s,l1012,l61012,lpc,llf,lbd,xu,yu,zu,
     $  cutparam,s612,2,rcut3)
      endif
      
      if (lpull) then
       call main_loop_pulling(nres,kcont,xn,yn,zn,seq,iseq,unit,
     $  ch,icmap,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta,isi,
     $  tsim,tsav,tpdb,kupdate,ntrj,nch,ichain,kk_out,kk_pdb,
     $  eij,lchir,langl,l612,l612s,l1012,l61012,lpc,llf,lbd,
     $  vpull,HH1,HH2,ip1,ip2,tfave,max_pull_f,stop_pulling,cutparam,
     $  s612,3,knot_ask,kk_knot,tksav,knumbers,rr,arg,c_velo,c_for,
     $  forcon,c_param,a_param,stop_pulling2,tstop,ktight,rcut3)
      endif
      
      if (lsurf) then
       call main_loop_su(nres,kcont,xn,yn,zn,seq,iseq,unit,
     $  ch,icmap,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta,isi,
     $  tsim,tsav,tpdb,teql,kupdate,ntrj,nch,ichain,kk_out,kk_pdb,
     $  eij,lchir,langl,l612,l1012,l61012,lpc,llf,lbd,
     $  boxl,nz,wz,wp,wf,itbs,tbs,cutparam,4,rcut3)
      endif
     
C CLOSE FLIES 

      close(kk_out)
      close(kk_pdb)
      close(kk_knot)

      end
