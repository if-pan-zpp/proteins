c       The i,i+2 contacts purely repulsive

        program cg
        implicit double precision(a-h,o-z)

        parameter(len=10000) !maximum number of all residues together
        character aseq*3,pdbfile*32,outfile*32,savfile*32
        character seqfile*32,arg*32,stafile*32,filname*32
        character buffer*128,paramfile*32,mapfile*32
        
        logical lconftm
        logical lallatom,langle
        logical l3rdcn(len)
        logical lfrompdb(len),ldi,lsim,lstartpdb,lpbc
        logical lparam,lpbcx,lpbcy,lpbcz
        logical lnatend,ldisimp
        logical lii4
        logical lcospid,lsldh
        logical lcontin,lconect(len),lpdb
        
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/parm/f02,f12,f32,f42,f52
        common/sig/sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/verl/oxv(3,len),vrcut2sq
        common/cmap/klont,klist(3,len*50)
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(2,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/hhar/H1,H2,tolerance
        common/misc/ad,lconftm,lsim
        common/bas/unit,men,lpdb
        common/kier/lpbcx,lpbcy,lpbcz
        common/mass/rmas(len),rsqmas(len)
        common/pid/Pi,c216
        common/masscenter/xcm,ycm,zcm,xmcm,ymcm,zmcm
        common/gyr/cgyr(len),cend(len),xyzcm(3,len),w(len),knts(2,len)
        common/chiral/CDH,langle,ldisimp
        common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
        common/angnat/the0(len),phi0(len),lfrompdb
        common/angtemp/thetemp(len),phitemp(len),chir(len)
        common/restart/delta,stafile,filname,klenstr
        common/ssb/knct3rd(len),l3rdcn
        common/respul/ip1,ip2
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/xyforces/xuforce,xdforce,yuforce,ydforce
        common/neigh/nei(2,len),rnei,rneisq
        dimension ttab(150) ! table of temperatures
        dimension kbt(len*50),fbt(len*50),fbt2(len*50)
        dimension kut(len*50),fut(len*50),fut2(len*50)
        dimension mtraj(len*20)
        dimension caac(len*50),kaak(len*50)

        pdbfile='1ubq.pdb'         ! PDB INPUT FILEAME
        seqfile='glut.txt'         ! SEQUENCE INPUT FILE
        outfile='saw.out'          ! OUTPUT FILENAME
        savfile='saw.pdb'          ! FILENAME FOR SAVING CONFORMATIONS
        paramfile='parameters.txt' ! PARAMETER FILE
        stafile='(a,0,a)'          ! FORMAT STRING FOR OTHER FILES
        ! SIMULATION MODE ----------------------------------------------
        lpdb=.false.      ! USING A PDB FILE FOR GO-MODEL SIMULATIONS
        ! SIMULATION OPTIONS -------------------------------------------
        lstartpdb=.false. ! START FROM PDB FILE
        lnatend=.false.   ! END WHEN ALL NATIVE CONTACTS PRESENT
        lconftm=.false.   ! CALCULATE TIME FOR CONTACT FORMATION
        lallatom=.false.  ! GO MODEL WITH ALL-ATOM CONTACT MAP
        lpbc=.false.      ! PERIODIC BOUNDARY COND-S IN EVERY DIRECTION
        lpbcx=.false.     ! PERIODIC BOUNDARY CONDITIONS IN X DIRECTION
        lpbcy=.false.     ! PERIODIC BOUNDARY CONDITIONS IN Y DIRECTION
        lpbcz=.false.     ! PERIODIC BOUNDARY CONDITIONS IN Z DIRECTION
        ! SIMULATION ITSELF --------------------------------------------
        lparam=.false.     ! IF NO PARAMFILE IS SPECIFIED, LSIM=.TRUE.
        lsim=.false.      ! FOR MODEL WITH SIMPLIFIED CONTACT POTENTIAL
        ldisimp=.false.   ! HARMONIC DIHEDRAL POTENTIAL
        lcospid=.false.! COSINE FUNCTION FOR THE PID POTENTIAL
        langle=.true.  ! FOR MODEL WITH ANGLE POTENTIALS
        ldi=.true.     ! FOR DIHEDRAL TERM IN POTENTIALS, IF LANGLE IS T
        lii4=.true.    ! i+4 CONTACTS ARE INCLUDED BY DEFAULT, EXCEPT BB
        lsldh=.false.  ! SLOWLY TURN ON DIHEDRAL POTENTIAL IN THE START
        ! SIMULATION PARAMETERS - GENERAL ------------------------------
        c216=2.0d0**(1.d0/6.d0)
        Pi=dacos(-1.d0)! PI
        twopi=2.d0*Pi  ! 2PI
        unit=5.d0      ! LENGTH UNIT
        iseed=448      ! RANDOM SEED
        ntraj=1        ! NUMBER OF TRAJECTORIES
        mskip=0        ! SKIPPING STEPS [tau]
        mstep=3000000  ! TOTAL SIMULATION TIME [tau]
        kwrite=100     ! HOW OFTEN TO PRINT ENERGY, 0: DON'T PRINT
        ksave= 1000    ! HOW OFTEN TO SAVE CONFORMATION, 0: DON'T SAVE
        tstart=0.35    ! STARTING TEMPERATURE
        tend=0.35      ! ENDING TEMPERATURE
        tstep=0.05     ! TEMPERATURE STEP
        nen1=1         ! FIRST SIMULATED INDEX (INDICES<nen1 ARE FROZEN)
        ! SIMULATION PARAMETERS - DISULFIDE BRIDGES
        nssb=0         ! NUMBER OF NATIVE DISULFIDE BRIDGE, 0: NO BRIDGE
        ! SIMULATION PARAMETERS - ANGLE POTENTIALS
        CBA=30.d0!20.0COEFF FOR BOND ANGLE POTENTIALS,   20 for Clementi
        CDA=0.66 !1.0 COEFF FOR DIHED. ANGLE POT. (K1),   1 for Clementi
        CDB=0.66 !0.5 COEFF FOR DIHED. ANGLE POT. (K3), 0.5 for Clementi
        CDH=3.33 !1.0 COEFF FOR DIHED. ANGLE POT. (HARMONIC), 1 for Clem
        ! SIMULATION PARAMETERS - PULLING AND STRETCHING
        dar=1.d-10     ! DISPLACEMENT FOR COMPUTING THE FORCES
        coef=0.01      ! COEFFICIENT FOR CONST FORCE  [epsilon/(a*tau)]
        veldist=12.d0  ! DISTANCE BEFORE ACQUIRING FULL VELOCITY
        naver= 100     ! AVERAGING INSTATANEOUS F OVER MD STEPS (tau)
        ! SIMULATION PARAMETERS - SIMULATION BOX
        tdens=0.001    ! TARGET DENSITY (IN RESIDUES/ANSTREM^3)
        sdens=0.0001   ! INITIAL DENSITY (THE BOX WILL SQUEEZE TO TDENS)
        kbperiodmax=6  ! NUMBER OF OSCILLATIONS AFTER WHICH IS PULLING
        ! SIMULATION PARAMETERS - TECHNICAL
        rnei=7.5       ! CUT-OFF DISTANCE FOR NEIGHBOUR RES. [Angstrem]
        ad=2000        ! TIME FOR ADIABATIC POTENTIAL TURNOFF [steps]
        cntfct=1.3d0   ! MULTIPLIES SIGMA TO CHECK IF CONTACT PRESENT
        confcut=4.56   ! MIN. DISTANCE FOR CHAINS GENERATED BY SAW [A]
        bond=3.8       ! BOND LENGTH [Angstrem]
        dnat=0.0       ! CUT-OFF DISTANCE FOR NATIVE CONTACTS [Angstrem]
        delta=0.005    ! INTEGRATION TIME STEP [tau]
        gamma=2.0      ! LANGEVIN PARAMETER [m/tau]
        H1 = 50.0      ! HARMONIC COEFF. [epsilon/A^2]
        H2 = 0.0       ! ANHARMONIC COEFF. [epsilon/A^4]
        verlcut=10.0   ! CUT-OFF DISTANCE FOR VERLET LIST[Angstrem]
        tolerance=0.d0 ! HOW FAR FROM THE MINIMUM CONTACT TURNS ON [%]
        cut=5.0        ! CUT-OFF DISTANCE FOR REPULSIVE TERM [Angstrem]
        rcut=18.d0     ! CUT-OFF DISTANCE FOR THE REST [Angstrem]
        ! READING VARIABLES --------------------------------------------
        if(iargc().gt.0) then ! reading parameters from external file
        call getarg(1,arg) ! for more than 1 args, do i_arg=1,iargc()
        open(7,file=arg,status='old') ! no spaces in filenames allowed
        write(*,*) 'RUNNING INPUTFILE ',arg
11      read(7,'(a)',end=12, err=12) buffer
        if(buffer(1:7).eq.'pdbfile') then
            read(buffer(8:),*) pdbfile 
        elseif(buffer(1:8).eq.'lallatom') then
            read(buffer(9:),*) lallatom
        elseif(buffer(1:6).eq.'cntfct') then
            read(buffer(7:),*) cntfct
        elseif(buffer(1:7).eq.'lconftm') then
            read(buffer(8:),*) lconftm
        elseif(buffer(1:7).eq.'lnatend') then
            read(buffer(8:),*) lnatend
        elseif(buffer(1:6).eq.'lparam') then
            read(buffer(7:),*) lparam
        elseif(buffer(1:5).eq.'ntraj') then
            read(buffer(6:),*) ntraj
        elseif(buffer(1:5).eq.'mstep') then
            read(buffer(6:),*) mstep
        elseif(buffer(1:4).eq.'lpdb') then
            read(buffer(5:),*) lpdb
        elseif(buffer(1:3).eq.'cut') then
            read(buffer(4:),*) cut
        elseif(buffer(1:7).eq.'klenstr') then
            read(buffer(8:),*) klenstr
        elseif(buffer(1:4).eq.'file') then
            read(buffer(5:),*) filname
            write(stafile,*) '(a',klenstr,',a)'
            write(outfile,stafile) filname,'.out'
            write(mapfile,stafile) filname,'.map'
            write(savfile,stafile) filname,'.pdb'
        else ! writing to console, unless file indexes 5 or 6 are in use
            write(*,*) 'UNRECOGNIZED OPTION: ',buffer
        endif
        goto 11
12      close(7)
        endif
        do i=1,len      ! LFROMPDB MUST BE ZEROED BEFORE LOADING CMAPS
            the0(i)=-1.0
            phi0(i)=0.0
            lfrompdb(i)=.false.
        enddo
        if(lpbc) then
            lpbcx=.true.
            lpbcy=.true.
            lpbcz=.true.
        endif
        klont=0
        open(1,file=outfile,status='unknown')
        if(ksave.ne.0) open(2,file=savfile,status='unknown')

        write(1,*)'#I,I+2 CONTACTS PURELY REPULSIVE'
        
        ! SCALE LENGTHS
        cut=cut/unit  ! REPULSIVE INTERACTIONS CUT-OFF
        sigma0=cut/c216 !*0.5d0**(1.d/6.d)! REPULSIVE INTERACTIONS SIGMA
        rcut = rcut/unit ! OTHER POTENTIALS CUT-OFF
        rnei=rnei/unit
        rcutsq=rcut*rcut
        cutsq=cut*cut
        rneisq=rnei*rnei
        verlcut=verlcut/unit
        confcut=confcut/unit
        vrcut2sq=verlcut*verlcut/4.d0
        dnat=dnat/unit
        H1=H1*unit*unit
        H2=H2*unit**4
        bond=bond/unit
        sdens=sdens*unit**3

        ! SETUP TABLE OF TEMPERATURE
        nt=nint(abs(tstart-tend)/tstep)+1
        if(nt.gt.150) nt=150
        ttstep=tstep
        if(tstart.gt.tend) ttstep=-tstep
        do i=1,nt
        ttab(i)=tstart+(i-1)*ttstep
        enddo
        
!       LOAD PROTEIN CONFORMATION
        write(1,'(/,a,2x,a,/)')'#PDB FILE =',pdbfile
        call load_protein(pdbfile)
        do i=1,men
           x0(i)=xn(i)
           y0(i)=yn(i)
           z0(i)=zn(i)
        enddo
        
        if(lallatom) then
            write(1,'(a)')'#CONSTRUCT CONTACT-MAP BASED ON ALL-ATOM'
            call compute_contact_map(pdbfile)
        elseif(klont.eq.0) then
        write(1,'(a,f6.2,a)')'#CONSTRUCT CONTACTMAP BASED ON CUT-OFF =',
     +      dnat*unit,' ANGSTROM'
            call compute_cmap(dnat)
        else
            write(1,'(a)')'#CONTACT-MAP BASED ON THE SEQUENCE FILE'
        endif
        
        
        ip1=1
        ip2=men
        do i=1,men
           rmas(i)=1.d0
        enddo
        jq=1
        part=men-nen1+1
        xmin=x0(1)
        ymin=y0(1)
        zmin=z0(1)
        xmax=x0(1)
        ymax=y0(1)
        zmax=z0(1)
        do ib=2,men
           if(x0(ib).lt.xmin) xmin=x0(ib)
           if(y0(ib).lt.ymin) ymin=y0(ib)
           if(z0(ib).lt.zmin) zmin=z0(ib)
           if(x0(ib).gt.xmax) xmax=x0(ib)
           if(y0(ib).gt.ymax) ymax=y0(ib)
           if(z0(ib).gt.zmax) zmax=z0(ib)
        enddo
        xdown=xmin-2*bond
        xup=xmax+2*bond
        ydown=ymin-2*bond
        yup=ymax+2*bond
        zdown=zmin-2*bond
        zup=zmax+2*bond
        xsep=xup-xdown
        ysep=yup-ydown
        zsep=zup-zdown
        xinv=1.d0/xsep
        yinv=1.d0/ysep
        zinv=1.d0/zsep
        call update_verlet_list(verlcut,nen1)

        write(1,'(a)') '#NO NON-REPULSIVE CONTACTS'

        if(klont.gt.0) then
            flush(1)
            corder=0
            do k=1,klont
            corder=corder + abs(klist(1,k)-klist(2,k))
            enddo
            corder=corder/klont
            corder=corder/men
            write(1,'(a,f8.4)') '#RELATIVE NATIVE CONTACT ORDER', corder
        else
            write(1,'(a)') '#NO NATIVE CONTACTS'
        endif
        
        do 47 ks=1,nssb
        i1=ksb(1,ks)
        i2=ksb(2,ks)
        if(aseq(i1).ne.'CYS') then
        write(1,*)'#RESIDUE ',ksb(ks,1),' IS NOT A CYSTEINE. PLS CHECK!' 
        stop
        endif
        if(aseq(i2).ne.'CYS') then
        write(1,*)'#RESIDUE ',ksb(ks,2),' IS NOT A CYSTEINE. PLS CHECK!' 
        stop
        endif
        icheck=0
        do k=1,klont
        ki1=klist(1,k)
        ki2=klist(2,k)
        if((ki1.eq.i1.and.ki2.eq.i2).or.(ki1.eq.i2.and.ki2.eq.i1)) then
           klist(3,k)=sign(631,klist(3,k)) ! SS BONDS HAVE VALUE +-631
           icheck=1
           write(1,'(a,i4,5x,2(a3,i4,3x))')'#NATIVE SS BOND',ks,
     +          aseq(i1),iseq(i1),aseq(i2),iseq(i2)
        endif
      enddo
      if(icheck.eq.0) then
      write(1,'(a,i4,5x,2(a3,i4,3x))')'#SS BOND COULD NOT BE MADE',ks,
     +        aseq(i1),iseq(i1),aseq(i2),iseq(i2)
      endif
 47   continue
        
        ! PUT CUSTOM POTENTIALS HERE
        
        
        write(1,'(a)')'#USING HARMONIC POTENTIALS FOR NATIVE SS BONDS!'
        olno=ran2(iseed) ! FOR LEGACY (COMPARE OLD VERSION, SAME SEED)
        write(1,'(a,2(a,f6.2))')'#USING ANHARMONIC POTENTIAL',
     +  '   H1 =',H1/unit/unit,'   H2 =',H2/unit**4

        if(langle) then
        if(ldi) then
          write(1,'(a)')'#USING POTENTIALS FOR BOND AND DIHEDRAL ANGLES'
        else
          write(1,'(a)')'#USING POTENTIALS ONLY FOR BOND ANGLES'
        endif
        endif
        
        if(lpdb) call compute_native_angles
        
        write(1,'(a)')'#DISABLE NATIVE CONTACTS (I,I+2)' ! AND (I,I+3)'
        km=0
        do k=1,klont
        i=klist(1,k)
        j=klist(2,k)
        ijdiff=j-i
        if(.not.lconect(i)) ijdiff=5
        if(ijdiff.ge.3) then ! 4) then
        km=km+1
        klist(1,km)=i
        klist(2,km)=j
        klist(3,km)=klist(3,k)
        endif
        enddo
        klont=km
        
        ngln=0
        do ib=1,men
            x0(ib)=xn(ib)     ! assign native coordinates to actual ones
            y0(ib)=yn(ib)
            z0(ib)=zn(ib)
            nei(1,ib)=0       ! set neighbour counter to zero
            nei(2,ib)=0
            if(aseq(ib).eq.'GLN') ngln=ngln+1 ! count glutamines
        enddo
        
        bond=0.d0 ! length of the Ca-Ca bond is averaged over all bonds
        do ic=1,nchains
            do i=menchain(ic)+1,menchain(ic+1)-1 
                dx=xn(i)-xn(i+1)
                dy=yn(i)-yn(i+1)
                dz=zn(i)-zn(i+1)
                dal=dx*dx+dy*dy+dz*dz
                dal=dsqrt(dal)
                b(i)=dal
                bond = bond + dal
            enddo
        enddo
        bond=bond/(men-1) ! this averaging previously was in gopotential
        
        write(1,'(a)')'#USING THE GO-LIKE 6-12 LJ POTENTIALS'
        if(lpdb) call gopotential(asigma)
        
        call prepare(edsg)
        call evalgo(edsg,chi)
        call evalcpot(edsg)
        if(langle) call evalangles(edsg,lsldh,1.d0)

        write(1,'(/,a,i10)') '#TOTAL PROTEIN LENGTH      ',men
        write(1,'(a,i10)')  '#NUMBER OF NATIVE CONTACTS ',klont
        if(lpdb) then
          write(1,'(a,f10.4)')'#AVERAGE LENGTH OF CONTACTS',asigma*unit
        else
          write(1,'(a)')'#NO PDB FILE USED FOR GO MODEL CONSTRUCTION'
        endif
        write(1,'(a,f10.2)')'#ENERGY OF NATIVE STATE    ',edsg
        write(1,*)
        if(naver.ne.0) 
     +  write(1,'(a,i8,a)') '#FORCE AVERAGED OVER ',naver,' TAU'
        write(1,'(a,f7.2)') '#VERLET LIST CUTOFF ',verlcut*unit
        if(lconftm) then
        write(1,'(a)')'#COMPUTING AVERAGED TIMES FOR CONTACT FORMATION'
        endif
        write(1,'(/,a,f10.3)')'#DELTA    =',delta
        write(1,'(a,f10.3)')'#GAMMA    =',gamma
        write(1,'(/,a,i10)') '#NUMBER OF TRAJECTORIES ',ntraj
        write(1,'(a,i10)') '#SIMULATION TIME        ',mstep
        write(1,'(a,i10)') '#SKIPPING STEPS         ',mskip
        write(1,'(a,i10)')'#RANDOM SEED            ',iseed
        write(1,'(/,a,7x,3f7.3)')'#TSTART TEND TSTEP',tstart,tend,tstep
        write(1,'(a,i6)')'#NUMBER OF TEMPERATURE STEPS',nt

        ! RESCALE TIME PARAMETERS
        kunit=nint(1.d0/delta) ! if delta = 0.005, this is 200
        naver=naver*kunit
        mskip=mskip*kunit
        mstep=mstep*kunit
        kwrite=kwrite*kunit
        ksave=ksave*kunit
        ! SCALE FACTORS FOR VELOCITIES DURING EQUILIBRATION
        delsq=delta*delta
        deltsq=0.5d0*delsq
        ! SET PARAMETERS IN PREDICTOR-CORRECTOR METHOD
        f02=dble(3./16.)
        f12=dble(251./360.)
        f32=dble(11./18.)
        f42=dble(1./6.)
        f52=dble(1./60.)
        dfr=ran2(iseed) ! dfr is not used anywhere, just here (legacy)
        
C ===============================================
        ! LOOP OVER TEMPERATURES
        do 2000 it=1,nt
        tr = ttab(it)
        write(1,'(/,a,f7.3)')'#TEMPERATURE ',tr

        ! LANGEVIN PARAMETERS 
        gamma2=gamma/delta
        const2=2.d0*tr*gamma*delta     ! assume xmas=1
        const2=dsqrt(const2)*delta
        aheat=delsq*part*3*tr

        if(lconftm) then
        do i=1,klont
        fbt(i)=0.d0
        fbt2(i)=0.d0
        mtraj(i)=0
        enddo
        endif

        inot=0



C =========================================
        ! LOOP OVER STARTING CONFIGURATIONS
        iterate=0
        do 1000 itraj=1,ntraj
        
        ! LANGEVIN PARAMETERS 
        
        lcontin=.true.  !true: simulation proceeds; false: it stops
        fresist = 0.d0
        cofdih=0.0 !    COEFFICIENT FOR SLOWLY TURNING ON DIHEDRAL POT.
        kqont=0    !    NUMBER OF NON-NATIVE CONTACTS IS ALS0 0
        jq=1       !    VERLET LIST HAS 2 COPIES, JQ=1 AND JQ=2

        iterate=iterate+1

        do i=1,klont
        kbt(i)=0
        kut(i)=0
        imap(i)=0
        enddo
        

        ! STARTING CONFORMATION
        if(lconftm) then   ! A STRAIGHT-LINE 
            do i=1,men
                x0(i)=0.d0
                y0(i)=0.d0
                z0(i)=(i-1)*bond
            enddo
            call confstart(confcut)
        endif
        
            startvolume=(men)/(sdens) !*unit**3)
            startboxsize=0.5*startvolume**(1.0/3.0)
            xmin=-startboxsize
            ymin=-startboxsize
            zmin=-startboxsize
            xmax=startboxsize
            ymax=startboxsize
            zmax=startboxsize
            do ib=1,men
                if(x0(ib)-2*bond.lt.xmin) xmin=x0(ib)-2*bond
                if(y0(ib)-2*bond.lt.ymin) ymin=y0(ib)-2*bond
                if(z0(ib)-2*bond.lt.zmin) zmin=z0(ib)-2*bond
                if(x0(ib)+2*bond.gt.xmax) xmax=x0(ib)+2*bond
                if(y0(ib)+2*bond.gt.ymax) ymax=y0(ib)+2*bond
                if(z0(ib)+2*bond.gt.zmax) zmax=z0(ib)+2*bond
            enddo
            xdown=xmin
            xup=xmax
            ydown=ymin
            yup=ymax
            zdown=zmin
            zup=zmax
            xsep=xup-xdown
            ysep=yup-ydown
            zsep=zup-zdown
            xinv=1.d0/xsep
            yinv=1.d0/ysep
            zinv=1.d0/zsep
        oldxup=xup
        oldxdown=xdown
        oldyup=yup
        oldydown=ydown
        oldzup=zup
        oldzdown=zdown
        


        ! LOAD INITIAL VELOCITIES OF PARTICLES
        call intvel3d(aheat,part,nen1)
        do i=1,men ! ZERO THE CONNECTION TABLES
        l3rdcn(i)=.false. ! true if residue i forms a disulfide bond
        knct3rd(i)=0        ! index of residue that is bonded to i
        nei(1,ib)=0
        nei(2,ib)=0
        enddo
        do ks=1,nssb
           i1=ksb(1,ks)
           i2=ksb(2,ks)
           l3rdcn(i1)=.true.
           l3rdcn(i2)=.true.
           knct3rd(i1)=i2
           knct3rd(i2)=i1
        enddo
        
        xsep=xup-xdown
        ysep=yup-ydown
        zsep=zup-zdown
        xinv=1.d0/xsep
        yinv=1.d0/ysep
        zinv=1.d0/zsep
        call update_verlet_list(verlcut,nen1) ! set up Verlet list
        ! ASSIGN INITIAL ACCELERATION BASED ON INITIAL POSITIONS
        call prepare(epot)
        call evalgo(epot,chi)
        call evalcpot(epot)
        if(langle) call evalangles(epot,lsldh,0.d0)

        ! SCALE ACCELERATIONS
        do 530 i=1,men
            x2(i)=fx(i)*deltsq
            y2(i)=fy(i)*deltsq
            z2(i)=fz(i)*deltsq
530     continue

        ! TABLE HEADING
        if(ksave.ne.0) write(2,'(a,i9)')'MODEL',iterate
        if(kwrite.ne.0) then
        write(1,'(//,a,i4)')'#TRAJECTORY',iterate
        write(1,'(a,a,a)')'#    TIME          EPOT          ETOT',
     +  '   ICN B1-B2 S1-S2 B1-S2 B1-B1 S1-S1 B1-S1',
     +  '      RG       L    RMSD NCORD     W CORDR KNOTS KNOTE'
        endif


c -----------------------------------------
c        ENTER MAIN LOOP OF SIMULATION
c -----------------------------------------

        kb=0
533     continue
        kb=kb+1

        call lang(twopi,gamma2,const2,nen1)
        call corr(deltsq,nen1)
        call predct(nen1)
        call prepare(epot)
        call evalgo(epot,chi)
        call evalcpot(epot)
        !endif
        if(langle) then
            if(lsldh.and.kb.le.ad) cofdih=kb*1.d0/ad
            call evalangles(epot,lsldh,cofdih)
        endif

        
        
        if(lnatend) lcontin=(icn.lt.klont)
        if(kb.lt.mskip) goto 533             ! SKIPPING STEPS


        ! CALCULATE KINETIC ENERGY AND MEAN COORDINATION NUMBER
        sumvel=0.d0
        do 540 i=1,men
            sumvel=sumvel+(x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i))*rmas(i)
540     continue

        ekin=sumvel/(2*delsq)           ! kinetic energy
        etot=epot+ekin                  ! total energy

        ! FIND CONTACTS WHICH APPEAR FOR THE FIRST TIME DURING FOLDING
        if(lconftm) then
           do k=1,klont
              i=klist(1,k)
              j=klist(2,k)
              if(kbt(k).eq.0) then
                 if(imap(k).eq.1) kbt(k)=kb ! CONTACT APPEARS
              endif
           enddo
        endif
        
        ! UPDATE VERLET LIST IF CERTAIN DISTANCE IS CROSSED
        vrcut2sqlim = vrcut2sq
        if(lpbcx) vrcut2sqlim = vrcut2sqlim-max(xdown-oldxdown,0.0)
        if(lpbcx) vrcut2sqlim = vrcut2sqlim-max(oldxup-xup,0.0)
        if(lpbcy) vrcut2sqlim = vrcut2sqlim-max(ydown-oldydown,0.0)
        if(lpbcy) vrcut2sqlim = vrcut2sqlim-max(oldyup-yup,0.0)
        if(lpbcz) vrcut2sqlim = vrcut2sqlim-max(zdown-oldzdown,0.0)
        if(lpbcz) vrcut2sqlim = vrcut2sqlim-max(oldzup-zup,0.0)
        do 1525 i=1,men
        if((x0(i)-oxv(1,i))**2+(y0(i)-oxv(2,i))**2+(z0(i)-oxv(3,i))**2
     +  .gt.vrcut2sqlim) then
            call update_verlet_list(verlcut,nen1)
            oldxup=xup
            oldxdown=xdown
            oldyup=yup
            oldydown=ydown
            oldzup=zup
            oldzdown=zdown
            goto 1526
        endif
1525    continue
1526    continue
        
        ! PRINTING TO FILE EVERY KWRITE STEPS AND AT THE END
        if(kwrite.ne.0) then
        if(kb.eq.1.or.mod(kb,kwrite).eq.0 .or. .not. lcontin) then
        ktime=nint(kb*delta)
        call gyration(rg)               ! RADIUS OF GYRATION
        call cgyration()                ! Rg and other parameters like W
        call compute_rmsd(rms)          ! RMSD

        dx=x0(ip1)-x0(ip2) ! END TO END DISTANCE
        dy=y0(ip1)-y0(ip2)
        dz=z0(ip1)-z0(ip2)
        ree=sqrt(dx*dx+dy*dy+dz*dz)

        write(1,'(i9,2f14.3,7i6,2f8.2,f8.3,3f6.2,2i6)') ktime,epot,etot,
     +  icn,0,0,0,0,0,0,
     +  rg*unit,ree*unit,rms*unit,0.0,w(1),0.0,knts(1,1),knts(2,1)

        call flush(1)

        endif
        endif
        if(ksave.ne.0.and.(((mod(kb,ksave).eq.0)).or.(kb.eq.1))) then
            time=kb*delta
            if(mod(ksave,kwrite).ne.0) then
                call compute_rmsd(rms)   ! RMSD
                call cgyration() ! RG FOR INDIVIDUAL CHAINS
            endif
            call print_conformation(2,time,epot,rms)
            call flush(2)
        endif

        
        if(abs(etot).gt.99999999.9) lcontin=.false. ! stop if it blew up
        
        if(kb.lt.mstep+mskip .AND. lcontin) goto 533

c -----------------------------------------
        ! END LOOP OF SIMULATION
c -----------------------------------------

        if(ksave.ne.0.and.mod(kb,ksave).ne.0) then
        time=kb*delta
        call compute_rmsd(rms)       ! RMSD
        call cgyration()
        call print_conformation(2,time,epot,rms)
        call flush(2)
        endif

        ! ACCUMULATE CONTACT BREAKING TIMES
        if(lconftm) then
        do k=1,klont
        if(kbt(k).ne.0) then
        aa=kbt(k)*delta
        fbt(k)=fbt(k)+aa
        fbt2(k)=fbt2(k)+aa*aa
        bb=kut(k)*delta
        fut(k)=fut(k)+bb
        fut2(k)=fut2(k)+bb*bb
        mtraj(k)=mtraj(k)+1
        endif
        enddo
        endif

c       print distances in the contacts

1000    continue
        ! END LOOP OVER CONFIGURATIONS

        if(lconftm) then
        write(1,'(/,a)')'#AVERAGE TIME NEEDED FOR FORMING EACH CONTACT'
        write(1,'(a)')'# ICN    I    J  J-I       t0    DISP.'
        do k=1,klont
        mntraj=mtraj(k)
        if(mntraj.ne.0) then
        fbt(k)=fbt(k)/mntraj
        fbt2(k)=fbt2(k)/mntraj
        fut(k)=fut(k)/mntraj
        fut2(k)=fut2(k)/mntraj
        endif
        aa=fbt(k)
        bb=sqrt(abs(fbt2(k)-aa*aa))/2
        i=klist(1,k)
        j=klist(2,k)
        caac(k)=aa ! time of contact break
        kaak(k)=k  ! index of contact
        enddo
            call sort2(klont,caac,kaak)
            do k=1,klont
            fcaa=caac(k)
            if(fcaa.gt.0) then
            ll=kaak(k)
            i=klist(1,ll)
            j=klist(2,ll)
            ise=iseq(i)
            jse=iseq(j)
            if(fcaa.gt.0.01) write(1,'(4i5,2f11.2)')ll,ise,jse,j-i,fcaa
            endif
            enddo
        endif


2000    continue
        ! END LOOP OVER TEMPERATURES
C ===============================================

        close(1)
c        write(2,'(a)') 'ENDMDL'
        close(2)
        close(22)
c ----------------------------
c         END OF EXECUTION
c ----------------------------
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine predct(nen1)
        implicit double precision(a-h,o-z)

        ! USE FIFTH-ORDER TAYLOR SERIES TO PREDICT POSITIONS & THEIR
        ! DERIVATIVES AT NEXT TIME-STEP

        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/bas/unit,men,lpdb
!$omp parallel default(shared), private(i)
!$omp do
        do 300 i=nen1,men
        x0(i)=x0(i)+x1(i)+x2(i)+x3(i)+x4(i)+x5(i)
        y0(i)=y0(i)+y1(i)+y2(i)+y3(i)+y4(i)+y5(i)
        z0(i)=z0(i)+z1(i)+z2(i)+z3(i)+z4(i)+z5(i)
        x1(i)=x1(i)+2.d0*x2(i)+3.d0*x3(i)+4.d0*x4(i)+5.d0*x5(i)
        y1(i)=y1(i)+2.d0*y2(i)+3.d0*y3(i)+4.d0*y4(i)+5.d0*y5(i)
        z1(i)=z1(i)+2.d0*z2(i)+3.d0*z3(i)+4.d0*z4(i)+5.d0*z5(i)
        x2(i)=x2(i)+3.d0*x3(i)+6.d0*x4(i)+10.d0*x5(i)
        y2(i)=y2(i)+3.d0*y3(i)+6.d0*y4(i)+10.d0*y5(i)
        z2(i)=z2(i)+3.d0*z3(i)+6.d0*z4(i)+10.d0*z5(i)
        x3(i)=x3(i)+4.d0*x4(i)+10.d0*x5(i)
        y3(i)=y3(i)+4.d0*y4(i)+10.d0*y5(i)
        z3(i)=z3(i)+4.d0*z4(i)+10.d0*z5(i)
        x4(i)=x4(i)+5.d0*x5(i)
        y4(i)=y4(i)+5.d0*y5(i)
        z4(i)=z4(i)+5.d0*z5(i)
300     continue
!$omp enddo nowait
!$omp end parallel 
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine corr(deltsq,nen1)
        implicit double precision(a-h,o-z)

        ! CORRECT PREDICTED POSITIONS AND THEIR DERIVATIVES

        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/parm/f02,f12,f32,f42,f52
        common/bas/unit,men,lpdb
!$omp parallel default(shared), private(i,xerr,yerr,zerr)
!$omp do
        do 670 i=nen1,men
          xerr=x2(i)-deltsq*fx(i)
          yerr=y2(i)-deltsq*fy(i)
          zerr=z2(i)-deltsq*fz(i)
          x0(i)=x0(i)-xerr*f02
          x1(i)=x1(i)-xerr*f12
          x2(i)=x2(i)-xerr
          x3(i)=x3(i)-xerr*f32
          x4(i)=x4(i)-xerr*f42
          x5(i)=x5(i)-xerr*f52
          y0(i)=y0(i)-yerr*f02
          y1(i)=y1(i)-yerr*f12
          y2(i)=y2(i)-yerr
          y3(i)=y3(i)-yerr*f32
          y4(i)=y4(i)-yerr*f42
          y5(i)=y5(i)-yerr*f52
          z0(i)=z0(i)-zerr*f02
          z1(i)=z1(i)-zerr*f12
          z2(i)=z2(i)-zerr
          z3(i)=z3(i)-zerr*f32
          z4(i)=z4(i)-zerr*f42
          z5(i)=z5(i)-zerr*f52
670       continue
!$omp enddo nowait
!$omp end parallel 
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE PREPARES VARIABLES FOR COMPUTING ENERGY AND FORCE

        subroutine prepare(epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/bas/unit,men,lpdb
        common/cmapi/cntfct,imap(len*50),icn
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/misc/ad,lconftm,lsim
        common/xyforces/xuforce,xdforce,yuforce,ydforce
        common/neigh/nei(2,len),rnei,rneisq
        
        epot=0.d0
        icn=0
        
        yuforce=0.d0
        ydforce=0.d0
        xuforce=0.d0
        xdforce=0.d0
!$omp parallel default(shared), private(i,j)
!$omp do
        do i=1,men
            nei(2,i)=nei(1,i)
            nei(1,i)=0
            fx(i) = 0.d0
            fy(i) = 0.d0
            fz(i) = 0.d0
            j=i+1
            v(1,i)=x0(j)-x0(i)
            v(2,i)=y0(j)-y0(i)
            v(3,i)=z0(j)-z0(i)
            j=i+2
            v(4,i)=x0(i)-x0(j)
            v(5,i)=y0(i)-y0(j)
            v(6,i)=z0(i)-z0(j)
        enddo
!$omp enddo nowait
!$omp end parallel
!$omp parallel default(shared), private(i,j,vxvnrm)
!$omp do
        do i=1,men-1
            j=i+1
            vxv(4,j)=v(3,j)*v(5,i)-v(2,j)*v(6,i)
            vxv(5,j)=v(1,j)*v(6,i)-v(3,j)*v(4,i)
            vxv(6,j)=v(2,j)*v(4,i)-v(1,j)*v(5,i)
            vxv(1,j)=v(3,i)*v(2,j)-v(2,i)*v(3,j)
            vxv(2,j)=v(1,i)*v(3,j)-v(3,i)*v(1,j)
            vxv(3,j)=v(2,i)*v(1,j)-v(1,i)*v(2,j)
            vxvnrm=vxv(1,j)*vxv(1,j)+vxv(2,j)*vxv(2,j)+vxv(3,j)*vxv(3,j)
            vxvnrm=dsqrt(vxvnrm)
            vxv(1,j)=vxv(1,j)/vxvnrm
            vxv(2,j)=vxv(2,j)/vxvnrm
            vxv(3,j)=vxv(3,j)/vxvnrm
            vnrm(j)=vxvnrm
        enddo
!$omp enddo nowait
!$omp end parallel
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTE THE ENERGY AND FORCE OF THE CUSTOM POTENTIAL

        subroutine evalcpot(epot)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        logical lconftm
        logical lfrompdb(len),lconect(len),l3rdcn(len),lsim
        logical lpbcx,lpbcy,lpbcz
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/sig/sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(2,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn
        common/angnat/the0(len),phi0(len),lfrompdb
        common/kier/lpbcx,lpbcy,lpbcz
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/hhar/H1,H2,tolerance
        common/misc/ad,lconftm,lsim
        common/bas/unit,men,lpdb
        common/ssb/knct3rd(len),l3rdcn
        common/pid/Pi,c216
        common/neigh/nei(2,len),rnei,rneisq
        
!$omp parallel default(shared), private(i,j,k,l,dx,dy,dz,r,rsq,icm,rsig,
!$omp& ene,fce,repx,repy,repz,kqistabs,iconttype,bbir,bbjr,bbij,
!$omp& vxvinrm,vxvjnrm,sdchni,sdchnj,i0,j0,sdir,sdirnrm,sdjrnrm,screenr,
!$omp& adiab,sdjr,r6,enec,expc,coulpotcoeff,totcoeff,phi,rb,rb2,rsi,
!$omp& ksdchnsi,ksdchnsj,kmaxbckbi,kmaxbckbj,iconttype2,it2,
!$omp& jt2,icheck,lssn,kss,lss,kqadabs)
!$omp do reduction(+:epot)
        do 466 k=1,kqont ! non-native attractive contacts
            rsig=1.5 
            i=kqist(1,k,jq) ! kqist 1 or 2 = residue numbers
            j=kqist(2,k,jq) ! jq = number of verlet list (1 or 2)
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.lt.rneisq) then
                nei(1,i)=nei(1,i)+1
                nei(1,j)=nei(1,j)+1
            endif
            if(rsq.gt.rcutsq) then
                goto 465
            endif
            r = sqrt(rsq)

            ! END OF DEFINING CONTACT TYPES, NOW COMPUTING FORCES
            if(r.gt.cut) then
               goto 465
            else
               rsi=sigma0/r
               r6=rsi**6
               ene=4.d0*r6*(r6-1.d0)+1.d0
               fce= 24.d0*r6*(1.d0-2.d0*r6)/r
            endif
            
            epot=epot+ene
            if(fce.gt.1.d+3) fce=1.d+3
            if(fce.lt.-1.d+3) fce=-1.d+3
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
465         continue
466     continue
!$omp enddo nowait
!$omp end parallel
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES THE ENERGY AND FORCE OF THE GO MODEL

        subroutine evalgo(epot,chi)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3
        logical lconftm
        logical lfrompdb(len),lconect(len),l3rdcn(len)
        logical lsim
        logical lpbcx,lpbcy,lpbcz
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/sig/sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmap/klont,klist(3,len*50)
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(2,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn
        common/angnat/the0(len),phi0(len),lfrompdb
        common/kier/lpbcx,lpbcy,lpbcz
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/hhar/H1,H2,tolerance
        common/misc/ad,lconftm,lsim
        common/bas/unit,men,lpdb
        common/ssb/knct3rd(len),l3rdcn
        common/neigh/nei(2,len),rnei,rneisq
        
!$omp parallel default(shared), private(i,j,xi,yi,zi,dx,dy,dz,
!$omp& rsq,r,rb,rb2,ene,fce,repx,repy,repz,rsi,r6)
!$omp do reduction(+:epot)
        do 495 i=1,men-1
        if(lconect(i)) then ! lconect is true if i and i+1 are connected
            j=i+1   ! contacts between neighbour residues
            xi=x0(i)
            yi=y0(i)
            zi=z0(i)
            dx = xi-x0(j)
            dy = yi-y0(j)
            dz = zi-z0(j)
            rsq=dx*dx+dy*dy+dz*dz
            nei(1,i)=nei(1,i)+1
            nei(1,j)=nei(1,j)+1
            r = sqrt(rsq)
            rb = r - b(i)
            rb2 = rb*rb
            ene = H1*rb2 + H2*rb2*rb2
            fce = (2*H1+4*H2*rb2)*rb ! NEGATIVE FCE MEANS REPULSION
            if(fce.gt.1.d+3) fce=1.d+3
            if(fce.lt.-1.d+3) fce=-1.d+3
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
            epot=epot+ene
            if(lconect(j)) then ! i,i+2 contacts purely repulsive
                j=i+2
                dx = xi-x0(j)
                dy = yi-y0(j)
                dz = zi-z0(j)
                rsq=dx*dx+dy*dy+dz*dz
                if(rsq.lt.rneisq) then
                    nei(1,i)=nei(1,i)+1
                    nei(1,j)=nei(1,j)+1
                endif
                if(rsq.lt.cutsq) then
                    r = sqrt(rsq)
                    rsi=sigma0/r
                    r6=rsi**6
                    ene=4.d0*r6*(r6-1.d0)+1.d0
                    fce=24.d0*r6*(1.d0-2.d0*r6)/r

                    if(fce.gt.1.d+3) fce=1.d+3
                    if(fce.lt.-1.d+3) fce=-1.d+3
                    fce=-fce/r
                    repx=fce*dx
                    repy=fce*dy
                    repz=fce*dz
                    fx(i) = fx(i) + repx
                    fx(j) = fx(j) - repx
                    fy(i) = fy(i) + repy
                    fy(j) = fy(j) - repy
                    fz(i) = fz(i) + repz
                    fz(j) = fz(j) - repz
                    epot=epot+ene
                endif
            endif
        endif
495     continue
!$omp enddo nowait
!$omp end parallel

        do 496 k=1,kcont ! all other contacts (including native)
            i=kcist(1,k)
            j=kcist(2,k)
            dx = x0(i)-x0(j)
            dy = y0(i)-y0(j)
            dz = z0(i)-z0(j)
            if(lpbcx) dx = dx-xsep*nint(dx*xinv)
            if(lpbcy) dy = dy-ysep*nint(dy*yinv)
            if(lpbcz) dz = dz-zsep*nint(dz*zinv)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq.lt.rneisq) then
                nei(1,i)=nei(1,i)+1
                nei(1,j)=nei(1,j)+1
            endif
            if(rsq.gt.rcutsq) goto 496
            r = sqrt(rsq)
            icm=kcist(3,k)
            if(abs(klist(3,icm)).lt.631) then ! L-J potential
               rsig=sont(icm)
               if(r.le.rsig*cntfct) then
                  icn=icn+1
                  if(imap(icm).eq.0) then
                     imap(icm)=1
                  endif
               else
                  if(imap(icm).eq.1) then
                     imap(icm)=0
                  endif
               endif
               rsi=rsig/r
               r6=rsi**6
               ene=4.d0*r6*(r6-1.d0)
               fce= 24.d0*r6*(1.d0-2.d0*r6)/r
            else if(abs(klist(3,icm)).eq.631) then ! for SSbonds
               icn=icn+1
               rb = r - 6.d0/unit
               rb2 = rb*rb
               ene = H1*rb2 + H2*rb2*rb2
               fce = (2*H1+4*H2*rb2)*rb
            endif
            epot=epot+ene
            if(fce.gt.1.d+3) fce=1.d+3
            if(fce.lt.-1.d+3) fce=-1.d+3
            fce=-fce/r
            repx=fce*dx
            repy=fce*dy
            repz=fce*dz
            fx(i) = fx(i) + repx
            fx(j) = fx(j) - repx
            fy(i) = fy(i) + repy
            fy(j) = fy(j) - repy
            fz(i) = fz(i) + repz
            fz(j) = fz(j) - repz
496     continue

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  ASSIGN INITIAL VELOCITIES TO ATOMS

        subroutine intvel3d(aheat,part,nen1)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/vel/x1(len),y1(len),z1(len)
        common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
        common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
        common/bas/unit,men,lpdb

        sumx=0.
        sumy=0.
        sumz=0.
        do 200 i=nen1,men
        xx=2.d0*(ran2(0)-0.5d0)
        yy=2.d0*(ran2(0)-0.5d0)
        zz=2.d0*(ran2(0)-0.5d0)
        xyz=1.d0/ dsqrt(xx*xx+yy*yy+zz*zz)
        x1(i)=xx*xyz
        y1(i)=yy*xyz
        z1(i)=zz*xyz
        x2(i)=0.d0
        y2(i)=0.d0
        z2(i)=0.d0
        x3(i)=0.d0
        y3(i)=0.d0
        z3(i)=0.d0
        x4(i)=0.d0
        y4(i)=0.d0
        z4(i)=0.d0
        x5(i)=0.d0
        y5(i)=0.d0
        z5(i)=0.d0
        sumx=sumx+x1(i)
        sumy=sumy+y1(i)
        sumz=sumz+z1(i)
200     continue

        ! SCALE VELOCITIES SO THAT TOTAL MOMENTUM = ZERO
        x=0.d0
        do 210 i=nen1,men
        x1(i)=x1(i)-sumx/part
        y1(i)=y1(i)-sumy/part
        z1(i)=z1(i)-sumz/part
        x=x+x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i)
210     continue

        ! SCALE VELOCITIES TO DESIRED TEMPERATURE
        heat= dsqrt(aheat/x)
        do 220 i=nen1,men
        x1(i)=x1(i)*heat
        y1(i)=y1(i)*heat
        z1(i)=z1(i)*heat
220     continue
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  THE LANGEVIN NOISE

        subroutine lang(twopi,gamma2,const2,nen1)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/vel/x1(len),y1(len),z1(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/bas/unit,men,lpdb

        ! X-COMPONENT
        do 10 i=nen1,men
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         x1(i)=x1(i)+const2*gam
         fx(i)=fx(i)-gamma2*x1(i)
10      continue

        ! Y-COMPONENT
        do 20 i=nen1,men
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         y1(i)=y1(i)+const2*gam
         fy(i)=fy(i)-gamma2*y1(i)
20      continue

        ! Z-COMPONENT
        do 30 i=nen1,men
         r1=ran2(0)
         r2=ran2(0)
         gam=dsqrt(-2.d0*log(r1))*dcos(twopi*r2)
         z1(i)=z1(i)+const2*gam
         fz(i)=fz(i)-gamma2*z1(i)
30      continue

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    THIS SUBROUTINE UPDATES VERLET LIST

        subroutine update_verlet_list(verlcut,nen1)
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lsamechain,lendchn
        logical lconftm,lsim,l3rdcn(len),lii4
        logical lpbcx,lpbcy,lpbcz
        character aseq*3
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/verl/oxv(3,len),vrcut2sq
        common/sig/sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/misc/ad,lconftm,lsim
        common/ssb/knct3rd(len),l3rdcn
        common/cmap/klont,klist(3,len*50)
        common/kier/lpbcx,lpbcy,lpbcz
        common/cmp2/kcont,kcist(3,len*500),kqont,kqist(2,len*1500,2),jq
        common/cmapi/cntfct,imap(len*50),icn
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lpdb
        kront = 0
        kcont = 0
        kqont2 = kqont
        kqont = 0
        kactual = 1
        kqactual = 1
        kfactual = 1
        jq2=jq
        jq=3-jq ! if jq was 2, it is 1, if it was 1, it is 2
        cutoff=verlcut+rcut
        cutoff2=cutoff**2
        
        do ib=1,men
            oxv(1,ib)=x0(ib)
            oxv(2,ib)=y0(ib)
            oxv(3,ib)=z0(ib)
        enddo
        if((lpbcx.and.xsep.lt.0.001).or.(lpbcy.and.ysep.lt.0.001)
     +  .or.(lpbcz.and.zsep.lt.0.001)) then
            write(1,*) 'BROKEN PBC. PLS CHECK!' 
            stop
        endif
        
        ic=2
        do i=1,men
            if(i.gt.menchain(ic)) ic=ic+1
                
            if(lconect(i).and.lconect(i+1)) then
                kdist=i+3
            else if(lconect(i)) then
                kdist=i+2
            else
                kdist=i+1
            endif
            do j=kdist,men
                lsamechain=j.le.menchain(ic)
                lendchn=j.eq.i+1
                if(i.lt.nen1 .and. j.lt.nen1) goto 1129
                xi=x0(i)
                yi=y0(i)
                zi=z0(i)
                dx = xi-x0(j)
                dy = yi-y0(j)
                dz = zi-z0(j)
                if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                if(lpbcy) dy = dy-ysep*nint(dy*yinv)
                if(lpbcz) dz = dz-zsep*nint(dz*zinv)
                rsq=dx*dx+dy*dy+dz*dz

                if(rsq.lt.cutoff2) then
9797                continue
                    if(kactual.le.klont) then ! find native contacts
                        k1=klist(1,kactual)
                        k2=klist(2,kactual)
                        if(k1.lt.i .or. (k1.eq.i .and. k2.lt.j)) then
                            kactual=kactual+1
                            goto 9797
                        else if(k1.eq.i .and. k2.eq.j) then
                            kcont=kcont+1
                            kcist(1,kcont)=i
                            kcist(2,kcont)=j
                            kcist(3,kcont)=kactual
                            kactual=kactual+1
                            goto 1129
                        endif
                    endif
                    ! 3579 is the label for "the rest" of contacts
                    if(lendchn) goto 1129
                    iname1=inameseq(i)
                    iname2=inameseq(j)

                    if(lsamechain) then
                        if(abs(j-i).eq.4 .and. .not. lii4) goto 1129
                    endif

9898                continue
                    if(kqactual.le.kqont2) then ! find non-native c.
                        kq1=kqist(1,kqactual,jq2)
                        kq2=kqist(2,kqactual,jq2)
                        if(kq1.lt.i .or.(kq1.eq.i.and.kq2.lt.j)) then
                            kqactual=kqactual+1
                            goto 9898
                        else if(kq1.eq.i .and. kq2.eq.j) then
                            kqont=kqont+1
                            kqist(1,kqont,jq)=i
                            kqist(2,kqont,jq)=j
                            kqactual=kqactual+1
                            goto 1129
                        endif
                    endif

                    ! make a new non-native contact
                    kqont=kqont+1
                    kqist(1,kqont,jq)=i
                    kqist(2,kqont,jq)=j
                endif
1129            continue
            enddo
        enddo
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    THIS SUBROUTINE GENERATE THE STARTING CONFIGURATION

        subroutine confstart(confcut)
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        logical lsim,lconect(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/misc/ad,lconftm,lsim
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/bas/unit,men,lpdb
        common/pid/Pi,c216
        dimension phi(len),theta(len)
        dimension T(len,3,3),R(len,3)
        
        ! this is "infinite" box, for close-packed one insert zeroes
        xmin=0.d0
        ymin=0.d0
        zmin=0.d0
        xmax=0.d0
        ymax=0.d0
        zmax=0.d0
        
        do 1000 ic=1,nchains
        kter=0
100     continue
        ! ASSIGN RANDOM VALUE FOR PHI AND THETA
        phi(menchain(ic)+1)=Pi/2
        theta(menchain(ic)+1)=0.d0
        phi(menchain(ic)+2)=0.d0
        theta(menchain(ic)+2)=ran2(0)*Pi/3
        do i=menchain(ic)+3,menchain(ic+1)-1
        phi(i)=(2.d0*ran2(0)-1.d0)*Pi
        theta(i)=ran2(0)*Pi/3        ! theta <= Pi/3
        enddo

        ! COMPUTE THE TRANSFER MATRICES
        do ib=menchain(ic)+1,menchain(ic+1)-1
        T(ib,1,1)=cos(theta(ib))
        T(ib,1,2)=sin(theta(ib))
        T(ib,1,3)=0.d0
        T(ib,2,1)=sin(theta(ib))*cos(phi(ib))
        T(ib,2,2)=-cos(theta(ib))*cos(phi(ib))
        T(ib,2,3)=sin(phi(ib))
        T(ib,3,1)=sin(theta(ib))*sin(phi(ib))
        T(ib,3,2)=-cos(theta(ib))*sin(phi(ib))
        T(ib,3,3)=-cos(phi(ib))
        enddo
        ! COMPUTE RANDOM DIRECTION OF PROTEIN CHAIN
        theta0=dacos(1.d0-2.d0*ran2(0))
        phi0=2.d0*Pi*ran2(0)
        ! COMPUTE THE BACK-BONE VECTORS
        do ib=menchain(ic)+1,menchain(ic+1)-1
        r1=bond*sin(theta0)*cos(phi0)
        r2=bond*sin(theta0)*sin(phi0)
        r3=bond*cos(theta0)
        do i=ib,menchain(ic)+1,-1
        R(ib,1)=T(i,1,1)*r1+T(i,1,2)*r2+T(i,1,3)*r3
        R(ib,2)=T(i,2,1)*r1+T(i,2,2)*r2+T(i,2,3)*r3
        R(ib,3)=T(i,3,1)*r1+T(i,3,2)*r2+T(i,3,3)*r3
        r1=R(ib,1)
        r2=R(ib,2)
        r3=R(ib,3)
        enddo
        enddo
        ranx=(xmax-xmin)*ran2(0)+xmin
        rany=(ymax-ymin)*ran2(0)+ymin
        ranz=(zmax-zmin)*ran2(0)+zmin
        ! COMPUTE THE POSITIONS OF MONOMERS
        do ib=menchain(ic)+1,menchain(ic+1)
        x0(ib)=ranx
        y0(ib)=rany
        z0(ib)=ranz
        do i=menchain(ic)+1,ib-1
        x0(ib)=x0(ib)+R(i,1)
        y0(ib)=y0(ib)+R(i,2)
        z0(ib)=z0(ib)+R(i,3)
        enddo
        enddo
        ! CHECK THE CONFORMATION
        do i=1,menchain(ic+1)-3
        do j=i+3,menchain(ic+1)
        dx = x0(j) - x0(i)
        dy = y0(j) - y0(i)
        dz = z0(j) - z0(i)

        rs2 = dx*dx+dy*dy+dz*dz
        if(rs2.lt.confcut*confcut) then
          kter = kter+1
          if(kter.lt.9000) then
            goto 100
          else
            write(2,*)'Confstart: FAILURE'
          stop
          endif
        endif
        enddo
        enddo
            do ib=1,menchain(ic+1)
            if(x0(ib).lt.xmin) xmin=x0(ib)
            if(y0(ib).lt.ymin) ymin=y0(ib)
            if(z0(ib).lt.zmin) zmin=z0(ib)
            if(x0(ib).gt.xmax) xmax=x0(ib)
            if(y0(ib).gt.ymax) ymax=y0(ib)
            if(z0(ib).gt.zmax) zmax=z0(ib)
            enddo
1000    continue
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
C THIS SUBROUTINE COMPUTES THE RADIUS OF GYRATION

        subroutine gyration(rg)
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/masscenter/xcm,ycm,zcm,xmcm,ymcm,zmcm
        common/bas/unit,men,lpdb
        common/mass/rmas(len),rsqmas(len)

        xcm=0.d0
        ycm=0.d0
        zcm=0.d0
        xmcm=0.d0
        ymcm=0.d0
        zmcm=0.d0
        do i=1,men
        xcm=xcm+x0(i)
        ycm=ycm+y0(i)
        zcm=zcm+z0(i)
        xmcm=xmcm+x0(i)*rmas(i)
        ymcm=ymcm+y0(i)*rmas(i)
        zmcm=zmcm+z0(i)*rmas(i)
        enddo
        xcm=xcm/men
        ycm=ycm/men
        zcm=zcm/men
        xmcm=xmcm/men-xdown
        ymcm=ymcm/men-ydown
        zmcm=zmcm/men-zdown
        rg=0.d0
        do i=1,men
        dx=x0(i)-xcm
        dy=y0(i)-ycm
        dz=z0(i)-zcm
        rg=rg+dx*dx+dy*dy+dz*dz
        enddo
        rg=rg/men
        rg=sqrt(rg)
        return
        end

C THIS SUBROUTINE COMPUTES W AND RADIUS OF GYRATION OF INDIVIDUAL CHAINS        
        subroutine cgyration()
        implicit double precision (a-h,o-z)
        parameter(len=10000)
        
        logical lconect(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/bas/unit,men,lpdb
        common/mass/rmas(len),rsqmas(len)
        common/gyr/cgyr(len),cend(len),xyzcm(3,len),w(len),knts(2,len)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4

        do ic=1,nchains
        chainlength=menchain(ic+1)-menchain(ic)

        knts(1,ic) = 0 ! no knots
        knts(2,ic) = 0 ! no knots
        
        cxcm=0.d0
        cycm=0.d0
        czcm=0.d0
        do i=menchain(ic)+1,menchain(ic+1)
        cxcm=cxcm+x0(i)
        cycm=cycm+y0(i)
        czcm=czcm+z0(i)
        enddo
        cxcm=cxcm/chainlength
        cycm=cycm/chainlength
        czcm=czcm/chainlength
        xyzcm(1,ic)=cxcm
        xyzcm(2,ic)=cycm
        xyzcm(3,ic)=czcm
        crg=0.d0
        do i=menchain(ic)+1,menchain(ic+1)
        dx=x0(i)-cxcm
        dy=y0(i)-cycm
        dz=z0(i)-czcm
        crg=crg+dx*dx+dy*dy+dz*dz
        enddo
        crg=crg/chainlength
        cgyr(ic)=sqrt(crg)
        w(ic)=0 ! dr/rb
        dx=x0(menchain(ic)+1)-x0(menchain(ic+1))
        dy=y0(menchain(ic)+1)-y0(menchain(ic+1))
        dz=z0(menchain(ic)+1)-z0(menchain(ic+1))
        cend(ic)=sqrt(dx*dx+dy*dy+dz*dz) !/chainlength
        enddo
        
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES THE NATIVE COUPLINGS WITHIN GO-MODEL
c   based on cut-off length dnat
        subroutine compute_cmap(dnat)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lpbcx,lpbcy,lpbcz,lconect(len)
        common/kier/lpbcx,lpbcy,lpbcz
        character aseq*3
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/sig/sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/cmap/klont,klist(3,len*50)
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/bas/unit,men,lpdb
        common/pid/Pi,c216
        klont=0

        do 10000 ic=1,nchains
            do 10000 icc=ic,nchains
                do 2000 ib1=menchain(ic)+1,menchain(ic+1)
                    if(ic.eq.icc) then
                        kdist=ib1+3
                    else
                        kdist=menchain(icc)+1
                    endif
                    do 2000 ib2=kdist,menchain(icc+1)
                        dx=xn(ib1)-xn(ib2)
                        dy=yn(ib1)-yn(ib2)
                        dz=zn(ib1)-zn(ib2)
                        if(lpbcx) dx = dx-xsep*nint(dx*xinv)
                        if(lpbcy) dy = dy-ysep*nint(dy*yinv)
                        if(lpbcz) dz = dz-zsep*nint(dz*zinv)
                        dal=dx*dx+dy*dy+dz*dz
                        dal=dsqrt(dal)
                        dcut=dnat
                        if(dal.le.dcut) then
                          klont=klont+1
                          klist(1,klont)=ib1
                          klist(2,klont)=ib2
                          if(ic.eq.icc) then
                            klist(3,klont)=1
                          else
                            klist(3,klont)=-1
                          endif
                        endif
2000            continue
10000   continue

        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   THIS SUBROUTINE COMPUTES THE NATIVE COUPLINGS WITHIN GO-MODEL
c   based on all-atom VdW spheres covering
        subroutine compute_contact_map(filename)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3,filename*32,aname*3
        logical lconect(len)
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/rall/rx_a(len,14),ry_a(len,14),rz_a(len,14)
        common/nall/na(len),aname(len,14)
        common/radi/vrad(len,14)
        common/cmap/klont,klist(3,len*50)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/bas/unit,men,lpdb

        alpha=(26/7.)**(1./6.)

        call load_allatom(filename)
        call assign_VdW_radius

        klont=0

        ! BUILD CONTACT MAP
        do 10000 ic=1,nchains
            do 10000 icc=ic,nchains
c               do 2000 ib1=1,men-2
                do 2000 ib1=menchain(ic)+1,menchain(ic+1)
                    na1=na(ib1)
                    if(ic.eq.icc) then
                        kdist=ib1+3
                    else
                        kdist=menchain(icc)+1
                    endif
c                   do 2000 ib2=ib1+2,men
                    do 2000 ib2=kdist,menchain(icc+1)
                        na2=na(ib2)
                        rmin=1.E+6
                        kcc=0
                        kccbb=0
                        kccbs=0
                        kccss=0
                        kccbsbs=0
                        kccbssb=0
                        do 200 ja1=1,na1
                        do 200 ja2=1,na2
                            rx=rx_a(ib1,ja1)-rx_a(ib2,ja2)
                            ry=ry_a(ib1,ja1)-ry_a(ib2,ja2)
                            rz=rz_a(ib1,ja1)-rz_a(ib2,ja2)
                            rr=sqrt(rx*rx+ry*ry+rz*rz)
                            vrsum=vrad(ib1,ja1)+vrad(ib2,ja2)
                            if(rr.le.vrsum*alpha) then
                                if(ja1.lt.5.and.ja2.lt.5) then
                                    kccbb=kccbb+1
                                else if(ja1.gt.4.and.ja2.gt.4) then
                                    kccss=kccss+1
                                else
                                    kccbs=kccbs+1
                                    if(ja1.lt.5.and.ja2.gt.4) then
                                        kccbsbs=kccbsbs+1
                                    endif
                                    if(ja1.gt.4.and.ja2.lt.5) then
                                        kccbssb=kccbssb+1
                                    endif
                                endif
                                kcc=kcc+1
                            endif
                            if(rr.lt.rmin) rmin=rr
200                     continue
                        if(kcc.gt.0) then
                        klont=klont+1
                        klist(1,klont)=ib1
                        klist(2,klont)=ib2
                        if(ic.eq.icc) then
                            klist(3,klont)=1
                        else
                            klist(3,klont)=-1
                        endif
                        if(kccbb.gt.0) klist(3,klont)=klist(3,klont)*2
                        if(kccbs.gt.0) klist(3,klont)=klist(3,klont)*3
                        if(kccss.gt.0) klist(3,klont)=klist(3,klont)*5
                        if(kccbsbs.gt.0) klist(3,klont)=klist(3,klont)*3
                        if(kccbssb.gt.0) klist(3,klont)=klist(3,klont)*7
                        endif
2000            continue
10000   continue
        
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function ran2(iseed)
      implicit double precision(a-h,o-z)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
c      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      common/rans/idum
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if(iseed.ne.0) idum=-iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,)B.n.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE LOADS THE NATIVE BACKBONE CONFORMATION OF A PROTEIN
C FROM ITS PDB FILE
        subroutine load_protein(filn)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lfrompdb(len),lpdb
        logical l3rdcn(len),lpbcx,lpbcy,lpbcz
        character aseq*3,ares*3,filn*32,bb*2,buffer*128,ch1*1,ch2*1,ch*1
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/kier/lpbcx,lpbcy,lpbcz
        common/bas/unit,men,lpdb
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/angnat/the0(len),phi0(len),lfrompdb
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/ssb/knct3rd(len),l3rdcn
        dimension ch(len)
        character bb2*2,bb4*4
        open(8,file=filn,status='old',iostat=ierr)
        if(ierr.ne.0) then
        write(*,'(a,2x,a)')'ERROR OPENING FILENAME',filn
        write(*,*)'PROGRAM STOPPED.'
        stop
        endif
        
        do ib=1,len
            lconect(ib)=.false.
        enddo
        nchains=0
        menchain(1)=0
        ib=0
        xdown=0.d0
        ydown=0.d0
        zdown=0.d0
15      read(8,'(a)',end=20) buffer
        if(buffer(1:4).ne.'ATOM') then
            if(buffer(1:3).eq.'END') goto 20
            if(buffer(1:3).eq.'TER') then
                if(ib.gt.menchain(nchains+1)) then !filter out DNA chain
                    lconect(ib)=.false.
                    nchains=nchains+1
                    menchain(nchains+1)=ib
                endif
            elseif(buffer(1:6).eq.'CRYST1') then
                read(buffer(1:33),'(6x,3f9.3)') xup,yup,zup
                xup=xup/unit
                yup=yup/unit
                zup=zup/unit
                xsep=xup-xdown
                ysep=yup-ydown
                zsep=zup-zdown
                xinv=1.d0/xsep
                yinv=1.d0/ysep
                zinv=1.d0/zsep
            endif
            goto 15
        endif
        if((buffer(17:17).ne.'A').and.(buffer(17:17).ne.' ')) goto 15
        read(buffer,'(13x,a4,a3,x,a1,i4,4x,3f8.3)')
     +  bb4,ares,ch1,ival,xval,yval,zval
        read(bb4,'(a2,2x)') bb
        read(bb4,'(x,a2,x)') bb2
        if(bb.eq.'CA'.or.bb2.eq.'CA') then
        ib=ib+1
        lconect(ib)=.true.
        if(lpdb) lfrompdb(ib)=.true.
        aseq(ib)=ares
        iseq(ib)=ival
        ch(ib)=ch1
        if(ib.gt.2) then
            if(lconect(ib-1).and.ch(ib-1).ne.ch(ib)) then
                lconect(ib-1)=.false.
                nchains=nchains+1
                menchain(nchains+1)=ib-1
            endif
        endif
        xn(ib)=xval/unit
        yn(ib)=yval/unit
        zn(ib)=zval/unit
        if(aseq(ib).eq.'GLY') then
        inameseq(ib)=1
        else if(aseq(ib).eq.'PRO') then
        inameseq(ib)=2
        else if(aseq(ib).eq.'GLN') then
        inameseq(ib)=3
        else if(aseq(ib).eq.'CYS') then
        inameseq(ib)=4
        else if(aseq(ib).eq.'ALA') then
        inameseq(ib)=5
        else if(aseq(ib).eq.'SER') then
        inameseq(ib)=6
        else if(aseq(ib).eq.'VAL') then
        inameseq(ib)=7
        else if(aseq(ib).eq.'THR') then
        inameseq(ib)=8
        else if(aseq(ib).eq.'ILE') then
        inameseq(ib)=9
        else if(aseq(ib).eq.'LEU') then
        inameseq(ib)=10
        else if(aseq(ib).eq.'ASN') then
        inameseq(ib)=11
        else if(aseq(ib).eq.'ASP') then
        inameseq(ib)=12
        else if(aseq(ib).eq.'LYS') then
        inameseq(ib)=13
        else if(aseq(ib).eq.'GLU') then
        inameseq(ib)=14
        else if(aseq(ib).eq.'MET') then
        inameseq(ib)=15
        else if(aseq(ib).eq.'HIS') then
        inameseq(ib)=16
        else if(aseq(ib).eq.'PHE') then
        inameseq(ib)=17
        else if(aseq(ib).eq.'ARG') then
        inameseq(ib)=18
        else if(aseq(ib).eq.'TYR') then
        inameseq(ib)=19
        else if(aseq(ib).eq.'TRP') then
        inameseq(ib)=20
        endif
        endif
        goto 15
20      continue
        close(8)
        men=ib
        
        open(8,file=filn,status='old',iostat=ierr)
16      read(8,'(a)',end=21) buffer
        if(buffer(1:6).eq.'SSBOND') then
            read (buffer(16:21), '(a1,i5)' ) ch1,icys1
            read (buffer(30:35), '(a1,i5)' ) ch2,icys2
            nssb=nssb+1
            do j=1,men
             if (ch(j).eq.ch1 .and. iseq(j).eq.icys1) then
              ksb(1,nssb)=j
             endif
             if (ch(j).eq.ch2 .and. iseq(j).eq.icys2) then
              ksb(2,nssb)=j
             endif
            enddo
        endif
        if(buffer(1:3).eq.'END') goto 21
        goto 16
21      continue
        close(8)
c        endif
        if(men.eq.0) then
            write(*,'(a,a)')'NO ATOMS IN THE FILE ',filn
            write(*,*)'PROGRAM STOPPED.'
            stop
        endif
        if(lconect(men)) then ! if pdb file does not have TER or CH tags
            lconect(men)=.false.
            nchains=nchains+1
            menchain(nchains+1)=men
        endif
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINES COMPUTE THE RMSD OF THE C ALPHA BACKBONE TO THE
C NATIVE BACKBONE TAKEN FROM PDB
        subroutine compute_rmsd(rms)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/bas/unit,men,lpdb
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        dimension r(len,3),rn(len,3)
        do ib=1,men
        r(ib,1)=x0(ib)
        r(ib,2)=y0(ib)
        r(ib,3)=z0(ib)
        rn(ib,1)=xn(ib)
        rn(ib,2)=yn(ib)
        rn(ib,3)=zn(ib)
        enddo
        rms = 0
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE PRINTS THE CHAIN COORDINATES IN PDB FORMAT
        subroutine print_conformation(iun,time,energy,rms)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character*3 aseq,ares
        character*2 chainid(len)
        logical lsim,lconect(len),l3rdcn(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/misc/ad,lconftm,lsim
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/bas/unit,men,lpdb
        common/gyr/cgyr(len),cend(len),xyzcm(3,len),w(len),knts(2,len)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/ssb/knct3rd(len),l3rdcn
        common/plates/zdown,zup,xup,xdown,yup,ydown
        common/neigh/nei(2,len),rnei,rneisq
        dimension nserial(len)
        
            xmin=x0(1)
            ymin=y0(1)
            zmin=z0(1)
            do ib=2,men
                if(x0(ib).lt.xmin) xmin=x0(ib)
                if(y0(ib).lt.ymin) ymin=y0(ib)
                if(z0(ib).lt.zmin) zmin=z0(ib)
            enddo
        
          write(iun,'(/,a,f10.1,2x,a,f12.2,2x,a,f9.4)')
     +   'REMARK   0 TIME =',time,'ENERGY =',energy,'RMSD =',rms*unit
          xmin=0.d0
          ymin=0.d0
          zmin=0.d0
        
        ic=1
        iatom=0
        cmin=-999/unit
        do 1000 ib=1,men
        ares=aseq(ib)
        iatom=iatom+1
        ic2=ic
1532    continue
        if(ic2 .gt. 153) then ! program supports up to 153 unique chains
            ic2=ic2-153
            goto 1532
        endif
        if(ic2 .le. 26) then ! A IS ASCII 65, Z IS 90 (26 letters)
            icchar=ic2+64 ! A,B,C,D...
        else
            if(ic2 .le. 52) then
                icchar=ic2+70 ! a,b,c,d...
            else
                icchar=ic2-5 ! 0,1,2,3...
            endif
        endif
        if(ic2 .le. 62) then
            write(chainid(ib),'(a,a)') ' ',char(icchar)
        else
            write(chainid(ib),'(i2)') ic2-53 ! 10,11,12,13...
        endif
        nserial(ib)=iatom
        if(x0(ib)-xmin.lt.cmin.or.y0(ib)-ymin.lt.cmin.or.
     +  z0(ib)-zmin.lt.cmin) then
        write(iun,'(a,i7,2x,a,a,a,i4,f12.2,2f8.2,f6.2)')
     +  'ATOM',iatom,'CA  ',ares,chainid(ib),iseq(ib),(x0(ib)-xmin)*unit
     +  ,(y0(ib)-ymin)*unit,(z0(ib)-zmin)*unit,float(nei(1,ib))
        else
        write(iun,'(a,i7,2x,a,a,a,i4,f12.3,2f8.3,f6.2)')
     +  'ATOM',iatom,'CA  ',ares,chainid(ib),iseq(ib),(x0(ib)-xmin)*unit
     +  ,(y0(ib)-ymin)*unit,(z0(ib)-zmin)*unit,float(nei(1,ib))
        endif
        if(.not. lconect(ib)) then
        iatom=iatom+1
        write(iun,'(a,i7,6x,a,a,i6,28x)') 'TER ',iatom,ares,
     +  chainid(ib),iseq(ib)
        iatom=iatom+1
        if(xyzcm(1,ic)-xmin.lt.cmin.or.xyzcm(2,ic)-ymin.lt.cmin.or.
     +  xyzcm(3,ic)-zmin.lt.cmin) then
        write(iun,'(a,i5,2x,3a,i4,f12.2,2f8.2)') 'HETATM',iatom,'C   ',
     +  'COG',chainid(ib),iseq(ib)+1,(xyzcm(1,ic)-xmin)*unit,
     +  (xyzcm(2,ic)-ymin)*unit,(xyzcm(3,ic)-zmin)*unit
        else
        write(iun,'(a,i5,2x,3a,i4,f12.3,2f8.3)') 'HETATM',iatom,'C   ',
     +  'COG',chainid(ib),iseq(ib)+1,(xyzcm(1,ic)-xmin)*unit,
     +  (xyzcm(2,ic)-ymin)*unit,(xyzcm(3,ic)-zmin)*unit
        endif
        write(iun,'(a,i4,a,f10.1,a,f8.2,a,f7.2,a,i4,a,i3,a,i3,a,f6.3)') 
     +  'REMARK',ic,' T= ',time,' RG= ',cgyr(ic)*unit,
     +  ' R_end_to_end= ',cend(ic)*unit,' N= ',menchain(ic+1)
     +  -menchain(ic),' K1= ',knts(1,ic),' K2= ',knts(2,ic),' W= ',w(ic)
c       write(iun,'(a,i4,a,f10.1,i6,3f8.3)')'COG   ',ic,' T= ',time,
c    +  menchain(ic+1)-menchain(ic),(xyzcm(1,ic)-xmin)*unit,
c    +  (xyzcm(2,ic)-ymin)*unit,(xyzcm(3,ic)-zmin)*unit
        ic=ic+1
        endif
1000    continue
        nssb2=0
        do 1100 ib=1,men
        if(l3rdcn(ib).and.(ib.lt.knct3rd(ib))) then
            nssb2=nssb2+1
            j=knct3rd(ib)
            write(iun,'(a,i4,2a,i5,2a,i5)') 'SSBOND',nssb2,' CYS',
     +      chainid(ib),iseq(ib),'    CYS',chainid(j),iseq(j)
            write(iun,'(a,2i5)') 'CONECT',nserial(ib),nserial(j)
        endif
1100    continue        
        write(iun,'(a)')'END'
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C THIS SUBROUTINE LOADS THE NATIVE BACKBONE CONFORMATION OF A PROTEIN 
C FROM ITS PDB FILE AND THEN CALCULATES THE TORSION ANGLES OF THE
C NATIVE STATE
        subroutine load_allatom(filename)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3,ares*3,filename*32,bb*3,buffer*128,bbo*3
        common/sequence/iseq(len),inameseq(len),aseq(len)
        character aname*3
        common/rall/rx_a(len,14),ry_a(len,14),rz_a(len,14)
        common/nall/na(len),aname(len,14)
        common/bas/unit,men,lpdb
        open(15,file=filename,status='old',iostat=ierr)
        if(ierr.ne.0) then
        write(*,'(a,2x,a)')'ERROR OPENING FILENAME',filename
        write(*,*)'PROGRAM STOPPED.'
        stop
        endif
        ivalold=-1
        ib=0
        do i=1,len
        na(i)=0
        enddo
        bbo='   '
15        read(15,'(a)',end=20) buffer
        if(buffer(1:3).eq.'END') goto 20
        if(buffer(1:4).ne.'ATOM') goto 15
        read(buffer,'(13x,a3,1x,a3,2x,i4,4x,3f8.3)')
     +  bb,ares,ival,xval,yval,zval
        if(bb.eq.bbo) goto 15
        if((buffer(17:17).ne.'A').and.(buffer(17:17).ne.' ')) goto 15
        bbo=bb
        if(ival.ne.ivalold) then !if(bb.eq.'N  ')then
        ivalold=ival
        ib=ib+1
        ja=0
        aseq(ib)=ares
        iseq(ib)=ival
        endif
        if(bb(1:1).eq.'N'.or.bb(1:1).eq.'C'.or.bb(1:1).eq.'O'.or.
     +     bb(1:1).eq.'S')then
        ja=ja+1
        rx_a(ib,ja)=xval
        ry_a(ib,ja)=yval
        rz_a(ib,ja)=zval
        na(ib)=ja
        aname(ib,ja)=bb
        endif
        goto 15
20        continue
        close(15)
        men=ib
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine assign_VdW_radius
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        character aseq*3,ares*3,ana*3
        common/sequence/iseq(len),inameseq(len),aseq(len)
        character aname*3
        common/rall/rx_a(len,14),ry_a(len,14),rz_a(len,14)
        common/nall/na(len),aname(len,14)
        common/radi/vrad(len,14)
        common/bas/unit,men,lpdb
        dimension nb(len)

        do 1000 ib=1,men
        vrad(ib,1)=1.64                ! BACKBONE NITROGEN
        vrad(ib,2)=1.88                ! BACKBONE C-ALPHA
        vrad(ib,3)=1.61                ! BACKBONE C'
        vrad(ib,4)=1.42                ! BACKBONE OXYGEN
        ares=aseq(ib)
        if(ares.eq.'GLY') then
          nb(ib)=4
          goto 1000
        else if(ares.eq.'PRO') then
          nb(ib)=7
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
        else if(ares.eq.'GLN') then
          nb(ib)=9
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.61
          vrad(ib,8)=1.42
          vrad(ib,9)=1.64
        else if(ares.eq.'CYS') then
          nb(ib)=6
          vrad(ib,5)=1.88
          vrad(ib,6)=1.77
        else if(ares.eq.'VAL') then
          nb(ib)=7
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
        else if(ares.eq.'PHE') then
          nb(ib)=11
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.61
          vrad(ib,8)=1.76
          vrad(ib,9)=1.76
          vrad(ib,10)=1.76
          vrad(ib,11)=1.76
        else if(ares.eq.'MET') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.77
          vrad(ib,8)=1.88
        else if(ares.eq.'ILE') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.88
        else if(ares.eq.'ASP') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.46
          vrad(ib,8)=1.42
        else if(ares.eq.'GLU') then
          nb(ib)=9
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.61
          vrad(ib,8)=1.46
          vrad(ib,9)=1.42
        else if(ares.eq.'LYS') then
          nb(ib)=9
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.88
          vrad(ib,9)=1.64
        else if(ares.eq.'ARG') then
          nb(ib)=11
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.64
          vrad(ib,9)=1.61
          vrad(ib,10)=1.64
          vrad(ib,11)=1.64
        else if(ares.eq.'SER') then
          nb(ib)=6
          vrad(ib,5)=1.88
          vrad(ib,6)=1.46
        else if(ares.eq.'THR') then
          nb(ib)=7
          vrad(ib,5)=1.88
          vrad(ib,6)=1.46
          vrad(ib,7)=1.88
        else if(ares.eq.'TYR') then
          nb(ib)=12
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.76
          vrad(ib,8)=1.76
          vrad(ib,9)=1.76
          vrad(ib,10)=1.76
          vrad(ib,11)=1.61
          vrad(ib,12)=1.46
        else if(ares.eq.'HIS') then
          nb(ib)=10
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.64
          vrad(ib,8)=1.76
          vrad(ib,9)=1.76
          vrad(ib,10)=1.64
        else if(ares.eq.'ASN') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.42
          vrad(ib,8)=1.64
        else if(ares.eq.'TRP') then
          nb(ib)=14
          vrad(ib,5)=1.88
          vrad(ib,6)=1.61
          vrad(ib,7)=1.76
          vrad(ib,8)=1.61
          vrad(ib,9)=1.64
          vrad(ib,10)=1.61
          vrad(ib,11)=1.76
          vrad(ib,12)=1.76
          vrad(ib,13)=1.76
          vrad(ib,14)=1.76
        else if(ares.eq.'ALA') then
          nb(ib)=5
          vrad(ib,5)=1.88
        else if(ares.eq.'LEU') then
          nb(ib)=8
          vrad(ib,5)=1.88
          vrad(ib,6)=1.88
          vrad(ib,7)=1.88
          vrad(ib,8)=1.88
        endif
1000        continue

        ! CHECK AND CORRECTION
        do ib=1,men
        nat=na(ib)
        if(nat.lt.nb(ib)) then
        write(1,*)'AMINO ACID HAS FEWER ATOMS THAN SHOULD BE'
        write(1,'(i5,2x,a3,2x,2i6)')iseq(ib),aseq(ib),nat,nb(ib)
        else if(nat.gt.nb(ib)) then
c       write(1,*)'AMINO ACID HAS MORE ATOMS THAN SHOULD BE'
c       write(1,'(i5,2x,a3,2x,2i6)')iseq(ib),aseq(ib),nat,nb(ib)
c       do j=1,nat
c       write(1,'(a3,2x,i3,2x,a3)')aseq(ib),iseq(ib),aname(ib,j)
c       enddo
        nat=nb(ib)
        na(ib)=nb(ib)
        endif
        do j=1,nat
        ares=aseq(ib)
        ana=aname(ib,j)
        rad=vrad(ib,j)
        if(ana(1:1).eq.'N'.and.rad.ne.1.64) then
           vrad(ib,j)=1.64
        else if(ana(1:1).eq.'S'.and.rad.ne.1.77) then
           vrad(ib,j)=1.77
        else if(ana(1:1).eq.'O'.and.rad.ne.1.42.and.rad.ne.1.46) then
           vrad(ib,j)=1.46
        else if(ana(1:1).eq.'C'.and.rad.ne.1.88.and.rad.ne.1.76.
     +          and.rad.ne.1.61) then
           vrad(ib,j)=1.88
        endif
        if(rad.eq.0.) then
         write(1,'(a)')'ATOM ERROR:'
         write(1,'(a3,2x,i3,2x,a3)')aseq(ib),iseq(ib),aname(ib,j)
         stop
        endif
        enddo
        enddo
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine gopotential(asigma)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lpbcx,lpbcy,lpbcz
        common/sig/sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
        common/cmap/klont,klist(3,len*50)
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/kier/lpbcx,lpbcy,lpbcz
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/vel/x1(len),y1(len),z1(len)
        common/bas/unit,men,lpdb

        icn=0
        delh=0.5
        mah=0

        bmin=1000
        bmax=-1000

        dsum1=0.d0
        do 2000 k=1,klont
           i=klist(1,k)
           j=klist(2,k)
           dx=xn(i)-xn(j)
           dy=yn(i)-yn(j)
           dz=zn(i)-zn(j)
           if(lpbcx) dx = dx-xsep*nint(dx*xinv)
           if(lpbcy) dy = dy-ysep*nint(dy*yinv)
           if(lpbcz) dz = dz-zsep*nint(dz*zinv)
           dal=dx*dx+dy*dy+dz*dz
           dal=dsqrt(dal)
           dsum1=dsum1+dal
           sont(k)=dal*0.5d0**(1.d0/6.d0) ! dal for 1012
 2000   continue

        if(klont.gt.0) then
           asigma=dsum1/klont
        else
           asigma=0.d0
        endif

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine compute_native_angles
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lfrompdb(len)
        common/nat/xn(len),yn(len),zn(len),ksb(2,len/10),nssb
        common/four/rx(4),ry(4),rz(4)
        common/bas/unit,men,lpdb
        common/angnat/the0(len),phi0(len),lfrompdb

        do 100 ib=2,men-1
        do ibb=1,3
        ib1=ib+ibb-2
        rx(ibb)=xn(ib1)
        ry(ibb)=yn(ib1)
        rz(ibb)=zn(ib1)
        enddo
        call bondangle(theta)
        the0(ib)=theta
100        continue

        do 200 ib=3,men-1
        do ibb=1,4
        ib1=ib+ibb-3
        rx(ibb)=xn(ib1)
        ry(ibb)=yn(ib1)
        rz(ibb)=zn(ib1)
        enddo
        call dihedral(phi)
        phi0(ib)=phi
200        continue

        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine evalangles(epot,lsldh,cofdih)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        logical lconect(len),lfrompdb(len),langle,ldi,ldisimp
        logical lsldh
        character*3 aseq
        common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4
        common/sequence/iseq(len),inameseq(len),aseq(len)
        common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
        common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
        common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
        common/chiral/CDH,langle,ldisimp
        common/four/rx(4),ry(4),rz(4)
        common/angnat/the0(len),phi0(len),lfrompdb
        common/angtemp/thetemp(len),phitemp(len),chir(len)
        common/bas/unit,men,lpdb
        dimension fi(3),fj(3),fk(3),fl(3),qij(3),qkj(3),qp(3)
        dimension ra(3),rb(3),rm(3),rn(3)
        ! CALCULATE ENERGY AND FORCE RELATED TO THE BOND ANGLE
!$omp parallel default(shared), private(ib,i1,i2,i3,rijrkj,
!$omp& d2ij,d2kj,k,dij,dkj,costh,theta,fi,fj,fk,rp,drp,qij,
!$omp& qkj,qp,ra,rb,ene,dvdp,j,theta2,theta3,theta4,theta5)
!$omp do reduction(+:epot)
        do 100 ib=2,men-1
        if(lconect(ib-1).and.lconect(ib)) then
        ! COMPUTING FORCE (THE SAME FOR EVERY TYPE OF POTENTIAL)
C      SWOPE, W. C. AND FERGUSON, D. M. (1992) ALTERNATIVE EXPRESSIONS 
C      FOR ENERGIES AND FORCES DUE TO ANGLE BENDING AND TORSIONAL ENERGY
C      J. COMPUT. CHEM. 13: 585-594. DOI: 10.1002/JCC.540130508
        ! v(i1)=-rij, v(i2)=rkj
        i1=ib-1
        i2=ib
        i3=ib+1
        rijrkj=0.d0
        d2ij=0.d0
        d2kj=0.d0
        do k=1,3
         rijrkj=rijrkj-v(k,i1)*v(k,i2)
         d2ij=d2ij+v(k,i1)*v(k,i1)
         d2kj=d2kj+v(k,i2)*v(k,i2)
        enddo
        dij=sqrt(d2ij)
        dkj=sqrt(d2kj)
        costh=rijrkj/(dij*dkj)
        costh=min(costh,1.d0)
        costh=max(costh,-1.d0)
        theta=acos(costh)
        if (costh.eq.0.d0) then
         do k=1,3
          fi(k)=0.d0
          fj(k)=0.d0
          fk(k)=0.d0
         enddo
        else
         do k=1,3
          qij(k)=-v(k,i1)/dij
          qkj(k)=v(k,i2)/dkj
          qp(k)=vxv(k,i2)
         enddo
         ra(1)=qij(2)*qp(3)-qij(3)*qp(2)
         ra(2)=qij(3)*qp(1)-qij(1)*qp(3)
         ra(3)=qij(1)*qp(2)-qij(2)*qp(1)
         rb(1)=qkj(2)*qp(3)-qkj(3)*qp(2)
         rb(2)=qkj(3)*qp(1)-qkj(1)*qp(3)
         rb(3)=qkj(1)*qp(2)-qkj(2)*qp(1)
         do k=1,3
          fi(k)=ra(k)/dij
          fk(k)=-rb(k)/dkj
          fj(k)=-fi(k)-fk(k)
         enddo
        endif
        ! COMPUTING POTENTIAL-DEPENDENT PART
        thetemp(ib)=theta
        if(langle) then
        if(lfrompdb(ib-1).and.lfrompdb(ib+1)) then
            theta=theta-the0(ib)
            ene=CBA*theta*theta
            dvdp=2.d0*CBA*theta
        else
            j=min(inameseq(ib),3)
            k=min(inameseq(ib+1),3)
            theta2=theta*theta
            theta3=theta2*theta
            theta4=theta3*theta
            theta5=theta4*theta
            ene=angpot(1,j,k)+angpot(2,j,k)*theta+angpot(3,j,k)*theta2
     +           +angpot(4,j,k)*theta3+angpot(5,j,k)*theta4
     +           +angpot(6,j,k)*theta5+angpot(7,j,k)*theta3**2
            dvdp=angpot(9,j,k)+angpot(10,j,k)*theta
     +           +angpot(11,j,k)*theta2+angpot(12,j,k)*theta3
     +           +angpot(13,j,k)*theta4+angpot(14,j,k)*theta5
        endif
        epot=epot+ene
        fx(i1)=fx(i1)-dvdp*fi(1)
        fy(i1)=fy(i1)-dvdp*fi(2)
        fz(i1)=fz(i1)-dvdp*fi(3)
        fx(i2)=fx(i2)-dvdp*fj(1)
        fy(i2)=fy(i2)-dvdp*fj(2)
        fz(i2)=fz(i2)-dvdp*fj(3)
        fx(i3)=fx(i3)-dvdp*fk(1)
        fy(i3)=fy(i3)-dvdp*fk(2)
        fz(i3)=fz(i3)-dvdp*fk(3)
        endif
        endif
100     continue
!$omp enddo nowait
!$omp end parallel 
        ! CALCULATE ENERGY AND FORCE RELATED TO THE DIHEDRAL ANGLE
!$omp parallel default(shared), private(ib,i1,i2,i3,i4,rm,rn,d2n,d2m,
!$omp& d2rkj,k,phi,fi,fj,fk,fl,dmn,cosphi,rijn,drkj,rijrkj,rklrkj,df,
!$omp& ene,dvdp,j,sinfi,cosfi,sin2fi,cos2fi,sincosfi)
!$omp do reduction(+:epot)
        do 200 ib=3,men-1
        if(lconect(ib-2) .and. lconect(ib-1) .and. lconect(ib)) then
        ! COMPUTING FORCE (THE SAME FOR EVERY TYPE OF POTENTIAL)
C     B EKKER, H., BERENDSEN, H. J. C. AND VAN GUNSTEREN, W. F. (1995)
C     FORCE AND VIRIAL OF TORSIONAL-ANGLE-DEPENDENT POTENTIALS
C     J. COMPUT. CHEM., 16: 527-533. DOI: 10.1002/JCC.540160502
        ! v(i1)=-rij, v(i2)=rkj, v(i3)=-rkl
        i1=ib-2
        i2=ib-1
        i3=ib
        i4=ib+1
        rm(1)=vxv(1,i2)
        rm(2)=vxv(2,i2)
        rm(3)=vxv(3,i2)
        rmnrm=vnrm(i2)
        rn(1)=vxv(1,i3)
        rn(2)=vxv(2,i3)
        rn(3)=vxv(3,i3)
        rnnrm=vnrm(i3)
        !d2n=0.d0
        !d2m=0.d0
        d2rkj=0.d0
        do k=1,3
        ! d2n=d2n+rn(k)*rn(k)
        ! d2m=d2m+rm(k)*rm(k)
         d2rkj=d2rkj+v(k,i2)*v(k,i2)
        enddo
!        if (d2n.eq.0.d0 .or. d2m.eq.0.d0) then
!         phi=0.d0
!         do k=1,3
!          fi(k)=0.d0
!          fj(k)=0.d0
!          fk(k)=0.d0
!          fl(k)=0.d0
!         enddo
!       else
         dmn=rn(1)*rm(1)+rn(2)*rm(2)+rn(3)*rm(3)
         cosphi=max(min(dmn,1.d0),-1.d0) !/sqrt(d2n*d2m)
         rijn=0.d0
         do k=1,3
          rijn=rijn-v(k,i1)*rn(k)
         enddo
         if (rijn.lt.0.d0) then
          phi=-1.d0*acos(cosphi)
         else
          phi=acos(cosphi)
         endif
         drkj=sqrt(d2rkj)
         do k=1,3
          fi(k)=rm(k)*drkj/rmnrm !/d2m
          fl(k)=-rn(k)*drkj/rnnrm !/d2n
         enddo
         rijrkj=0.d0
         rklrkj=0.d0
         do k=1,3
          rijrkj=rijrkj-v(k,i1)*v(k,i2)
          rklrkj=rklrkj-v(k,i3)*v(k,i2)
         enddo
         do k=1,3
          df=(fi(k)*rijrkj-fl(k)*rklrkj)/d2rkj
          fj(k)=-fi(k)+df
          fk(k)=-fl(k)-df
         enddo
!        endif
        ! CALCULATING PART DEPENDENT ON POTENTIAL
        phitemp(ib)=phi
        if(langle.and.ldi) then
        if(lfrompdb(ib-2).and.lfrompdb(ib+1)) then
            phi=phi-phi0(ib)
            if(ldisimp) then
                ene=0.5*CDH*phi*phi
                dvdp=-CDH*phi
            else
                ene=CDA*(1.d0-dcos(phi))+CDB*(1.d0-dcos(3.d0*phi))
                dvdp=CDA*sin(phi)+3.d0*CDB*sin(3.d0*phi)
            endif
        else
            j=min(inameseq(ib-1),3)
            k=min(inameseq(ib),3)
            sinfi=dsin(phi)
            cosfi=cosphi
            sin2fi=sinfi**2
            cos2fi=cosfi**2
            sincosfi=sinfi*cosfi
            ene=dihpot(1,j,k)+dihpot(2,j,k)*sinfi+dihpot(3,j,k)*cosfi
     +      +dihpot(4,j,k)*sin2fi+dihpot(5,j,k)*cos2fi
     +      +dihpot(6,j,k)*sincosfi
            dvdp=dihpot(2,j,k)*cosfi-dihpot(3,j,k)*sinfi
     +      +2.d0*(dihpot(4,j,k)-dihpot(5,j,k))*sincosfi
     +      +dihpot(6,j,k)*(cos2fi-sin2fi)
        endif
        if(lsldh) then
            dvdp=dvdp*cofdih
            ene=ene*cofdih
        endif
        epot=epot+ene    
        fx(i1)=fx(i1)-dvdp*fi(1)
        fy(i1)=fy(i1)-dvdp*fi(2)
        fz(i1)=fz(i1)-dvdp*fi(3)
        fx(i2)=fx(i2)-dvdp*fj(1)
        fy(i2)=fy(i2)-dvdp*fj(2)
        fz(i2)=fz(i2)-dvdp*fj(3)
        fx(i3)=fx(i3)-dvdp*fk(1)
        fy(i3)=fy(i3)-dvdp*fk(2)
        fz(i3)=fz(i3)-dvdp*fk(3)
        fx(i4)=fx(i4)-dvdp*fl(1)
        fy(i4)=fy(i4)-dvdp*fl(2)
        fz(i4)=fz(i4)-dvdp*fl(3)
        endif
        endif
200     continue
!$omp enddo nowait
!$omp end parallel 
        return
        end

C THIS SUBROUTINE RETURN THE BOND ANGLE AT THE SECOND SITE
c ux1=x0(ib)-x0(ib-1)
c ux2=x0(ib)-x0(ib+1)
        subroutine bondangle(theta)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/four/rx(4),ry(4),rz(4)
        common/bas/unit,men,lpdb

        ux1=rx(2)-rx(1)
        uy1=ry(2)-ry(1)
        uz1=rz(2)-rz(1)
        uu1=sqrt(ux1*ux1+uy1*uy1+uz1*uz1)

        ux2=rx(2)-rx(3)
        uy2=ry(2)-ry(3)
        uz2=rz(2)-rz(3)
        uu2=sqrt(ux2*ux2+uy2*uy2+uz2*uz2)

        u12=ux1*ux2+uy1*uy2+uz1*uz2
        u12=u12/(uu1*uu2)

        u12=min(u12,1.d0)
        u12=max(u12,-1.d0)
        theta=dacos(u12)

        return
        end

C THIS SUBROUTINE RETURN THE DIHEDRAL ANGLE AT THE THIRD SITE

        subroutine dihedral(phi)
        implicit double precision(a-h,o-z)
        parameter(len=10000)
        common/four/rx(4),ry(4),rz(4)
        common/bas/unit,men,lpdb

        ux1=rx(2)-rx(1)
        uy1=ry(2)-ry(1)
        uz1=rz(2)-rz(1)

        ux2=rx(3)-rx(2)
        uy2=ry(3)-ry(2)
        uz2=rz(3)-rz(2)

        ux3=rx(4)-rx(3)
        uy3=ry(4)-ry(3)
        uz3=rz(4)-rz(3)

        vx1=uy1*uz2-uz1*uy2
        vy1=uz1*ux2-ux1*uz2
        vz1=ux1*uy2-uy1*ux2
        vv1=vx1*vx1+vy1*vy1+vz1*vz1

        vx2=uy2*uz3-uz2*uy3
        vy2=uz2*ux3-ux2*uz3
        vz2=ux2*uy3-uy2*ux3
        vv2=vx2*vx2+vy2*vy2+vz2*vz2

        v12=vx1*vx2+vy1*vy2+vz1*vz2
        v12=v12/sqrt(vv1*vv2)
        v12=min(v12,1.d0)
        v12=max(v12,-1.d0)
        phi=dacos(v12)

        di=vx1*ux3+vy1*uy3+vz1*uz3
        if(di.lt.0.d0) phi=-phi

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE sort2(n,arr,ibarr)
        implicit double precision(a-h,o-z)
      INTEGER n,M,NSTACK           !    from "Numerical Recipes"
c      REAL arr(n)
       dimension arr(n)
                integer ibarr(n)
      PARAMETER (M=7,NSTACK=500)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
c      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
       a=arr(j)
       ib=ibarr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
       ibarr(i+1)=ibarr(i)
11        continue
          i=0
2         arr(i+1)=a
      ibarr(i+1)=ib
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
      itempb=ibarr(k)
        arr(k)=arr(l+1)
      ibarr(k)=ibarr(l+1)
        arr(l+1)=temp
      ibarr(l+1)=itempb
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
      itempb=ibarr(l+1)
          arr(l+1)=arr(ir)
       ibarr(l+1)=ibarr(ir)
          arr(ir)=temp
       ibarr(ir)=itempb
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
       itempb=ibarr(l)
          arr(l)=arr(ir)
       ibarr(l)=ibarr(ir)
          arr(ir)=temp
       ibarr(ir)=itempb
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
      itempb=ibarr(l+1)
          arr(l+1)=arr(l)
       ibarr(l+1)=ibarr(l)
          arr(l)=temp
       ibarr(l)=itempb
        endif
        i=l+1
        j=ir
        a=arr(l)
      ib=ibarr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
      itempb=ibarr(i)
        arr(i)=arr(j)
       ibarr(i)=ibarr(j)
        arr(j)=temp
       ibarr(j)=itempb
        goto 3
5       arr(l)=arr(j)
       ibarr(l)=ibarr(j)
        arr(j)=a
       ibarr(j)=ib
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          write(*,*)'NSTACK too small in sort'
          stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
