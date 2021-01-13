c       The i,i+2 contacts purely repulsive

      program cg
      implicit double precision(a-h,o-z)

      parameter(len=10000)      !maximum number of all residues together
      character aseq*3,pdbfile*32,outfile*32,savfile*32,state*2,buf1*4
      character rstfile*64,seqfile*32,arg*32,stafile*32,filname*32
      character buffer*128,buf2*4,paramfile*32,mapfile*32,cmapf*32

      logical lconftm,lmedian,lthermo,lforce,lunfold,lvelo,lmass,ldisp
      logical lchiral,lpullrel,lallatom,langle,lj1012,lconstvol,lradii
      logical loscillate,l3rdcn(len),lfrmscrtch,lwarmup,lshear,ldelrst
      logical lfrompdb(len),lwritemap,lsimpang,ldi,lsim,lstartpdb,lpbc
      logical lgln,lmrs,lparam,lchargend,ldet,lpbcx,lpbcy,lpbcz,lcintr
      logical lcoilang,lcoildih,lsawconftm,lpid,lsqpbc,ldynss,lpullfin
      logical lteql,lposcrd,lsslj,lnatend,ldisimp,lkmt,lsselj,lminfbox
      logical lunwrap,lampstrict,ljwal,lenetab,lfcc,lcleanrst,lrepcoul
      logical lnowal,lmj,lbar,lcmap,lwals,lsink,lwrtang,lwal,lcpb,lii4
      logical lobo,lwritego,lcdnat,lcospid,lepid,lrmsmax,lecperm,lsldh
      logical lrst,lcontin,lconect(len),lpdb,ldens,lmaxforce,lwritexyz
      logical lcpot

      common/sequence/iseq(len),inameseq(len),aseq(len)
      common/pos/x0(len),y0(len),z0(len),v(6,len),vxv(6,len),vnrm(len)
      common/vel/x1(len),y1(len),z1(len)
      common/der/x2(len),y2(len),z2(len),x3(len),y3(len),z3(len)
      common/der2/x4(len),y4(len),z4(len),x5(len),y5(len),z5(len)
      common/for/fx(len),fy(len),fz(len),xsep,ysep,zsep,xinv,yinv,zinv
      common/parm/f02,f12,f32,f42,f52
      common/sig/sigma1(500),sigma0,sont(len*99),cut,rcut,rcutsq,cutsq
      common/verl/oxv(3,len),vrcut2sq,af,kfcc(3,len*500,2),kfccw,menw
      common/cmap/kront,krist(3,len*500),klont,klist(3,len*50),lcintr
      common/cmp2/kcont,kcist(3,len*500),kqont,kqist(4,len*1500,2),jq
      common/sdch/icnt(6),ksdchn(21,4),ksdchns(len),khbful(4,len)
      common/cmapi/cntfct,imap(len*50),icn,intrhc,intrsc,intehc,intesc
      common/equil/kteql,ktrest,nratvel,nratveld,kconnecttime
      common/bon/bond,b(len-1),lconect,menchain(len),nchains,lii4,lcpb
      common/nat/xn(len),yn(len),zn(len),enative,ksb(2,len/10),nssb
      common/hhar/HH1,HH2,H1,H2,potcoeff,tolerance,screend,coul,factor
      common/misc/adia(len),ad,ethi,mcmr,mchi,ncord,lconftm,lcpot,lsim
      common/wal/walmindst,kbwal(9),icw(2),lobo,ljwal,lwal,lwals,ldens
      common/bas/unit,men,lsqpbc,lpdb,lwritemap,lradii,lsink,lkmt,lfcc
      common/kier/afx,afy,afz,shear,lshear,lpbcx,lpbcy,lpbcz,kbperiod
      common/pull/xpul(len),ypul(len),zpul(len),vpulx,vpuly,vpulz,cofp
      common/mass/rmas(len),rsqmas(len)
      common/pid/alphacos(3),rbb(2),psi0ss,psi0bb(2),Pi,c216,vmp,lbar
      common/masscenter/xcm,ycm,zcm,xmcm,ymcm,zmcm
      common/gyr/cgyr(len),cend(len),xyzcm(3,len),w(len),knts(2,len)
      common/chiral/chirn(len),echi,CDH,lchiral,langle,ldisimp,lwrtang
      common/ang/CBA,CDA,CDB,dar,angpot(16,20,20),dihpot(16,20,20),ldi
      common/angnat/the0(len),phi0(len),lfrompdb,lsimpang,lcoilang
      common/angtemp/thetemp(len),phitemp(len),chir(len),lcoildih
      common/restr/bckbmin,bckb2min,sdchnmax,sfact,kmaxbckb,lecperm
      common/restart/delta,work,sep0,rstfile,stafile,filname,klenstr
      common/ssb/knct3rd(len),l3rdcn,disul,dmrs,amrs,rmrs,lmrs,ldynss
      common/ssb2/dislj,icnss,icdss,sigss,lsselj,lrepcoul,lepid
      common/respul/z0temp(len),ksorted(len),ip1,ip2,ipwn,ipw(2,len)
      common/plates/zdown,zup,zuforce,zdforce,fwal,xup,xdown,yup,ydown
      common/xyforces/xuforce,xdforce,yuforce,ydforce
      common/nmapi/rangli(8),emj(500),ncnt,nli(5,len*50),ldet,lpid,lmj
      common/neigh/nei(2,len),rnei,rneisq,lposcrd,neimin,neimaxdisul
      common/mr/tene(18001,3,3,2),pp,lenetab
      dimension ttab(150)       ! table of temperatures
      dimension kbt(len*50),fbt(len*50),fbt2(len*50)
      dimension kut(len*50),fut(len*50),fut2(len*50)
      dimension mtraj(len*20),ncorders(len)
      dimension tfold(501)      ! folding time for each trajectory
      dimension works(501),youngmod(4,501) ! work and E per oscillation
      dimension aufres(100000),aufres2(100000),aree(100000)
      dimension adfres(100000),adfres2(100000)
      dimension caac(len*50),kaak(len*50)

        pdbfile='1ubq.pdb'         ! PDB INPUT FILEAME
        seqfile='glut.txt'         ! SEQUENCE INPUT FILE
        outfile='saw.out'          ! OUTPUT FILENAME
        savfile='saw.pdb'          ! FILENAME FOR SAVING CONFORMATIONS
        rstfile='saw.rst'          ! RESTART FILE
        paramfile='parameters.txt' ! PARAMETER FILE
        stafile='(a,0,a)'          ! FORMAT STRING FOR OTHER FILES
        ! SIMULATION MODE ----------------------------------------------
        lpdb=.false.      ! USING A PDB FILE FOR GO-MODEL SIMULATIONS
        lforce=.false.    ! USING EXTERNAL FORCE FOR AFM-LIKE STRETCHING
        lvelo=.false.     ! CONST VELO. WITH AFM, MUST BE FALSE IF LWALL
        lwal=.false.      ! EXISTENCE OF THE SIMULATION BOX
        lwals=.false.     ! WALLS IN ALL 3 DIRECTIONS, NOT ONLY Z
        ldens=.false.     ! SYSTEM IS SQUEEZED TO REACH TARGET DENSITY
        lmedian=.false.   ! CALCULATE THE MEDIAN FOLDING TIME
        lthermo=.false.   ! CALCULATE THERMODYNAMIC PROPERTIES
        lunfold=.false.   ! .TRUE. FOR STUDY OF UNFOLDING
        loscillate=.false.! WALL OSCILLATION MEASUREMENTS
        lshear=.false.    ! TRUE FOR SHEARING SIMULATIONS
        ! SIMULATION OPTIONS -------------------------------------------
        lfcc=.false.      ! FCC WALLS
        lstartpdb=.false. ! START FROM PDB FILE
        lnatend=.false.   ! END WHEN ALL NATIVE CONTACTS PRESENT
        lconftm=.false.   ! CALCULATE TIME FOR CONTACT FORMATION
        lsawconftm=.true. ! USE SELF-AVOIDING WALK FOR THAT TASK
        lallatom=.false.  ! GO MODEL WITH ALL-ATOM CONTACT MAP
        lnowal=.false.    ! PBC EVEN FOR WALL DIRECTIONS
        lwarmup=.false.   ! WARM UP THE SYSTEM STARTING FROM NAT
        lteql=.false.     ! USE Teql FOR EQUILIBRATION, THEN USE T
        lrst=.false.      ! RESTART FROM SAVED STATE
        ldelrst=.true.    ! DELETE OLD RESTART FILES DURING SIMULATION
        lfrmscrtch=.false.! RESTART FROM PDB FILE (NO VELOCITIES)
        lminfbox=.false.  ! EXTEND SQUEEZED BOX TO FIND FORCE MINIMUM
        lpullfin=.true.   ! PULL THE BOX WHEN SIMULATION FINISHES
        lobo=.false.      ! ATTACH RESIDUES TO WALL ONE BY ONE
        lwritemap=.false. ! WRITE MAPFILE WITH CONTACTS
        lwritego=.false.  ! WRITE GO CONTACTS IN MAPFILE FORMAT
        lwrtang=.true.    ! IF LWRITEMAP IS TRUE WRITE ALSO ANGLES
        lwritexyz=.false. ! WRITE COORDINATES IN XYZ FILE INSTEAD OF PDB
        lposcrd=.false.   ! SHIFT COORDINATES TO STAY POSITIVE
        lcdnat=.false.    ! USE CUSTOM CUTOFF PARAMETERS FOR GO MAP
        lconstvol=.false. ! KEEP CONSTANT VOLUME DURING STRETCHING
        lcmap=.false.     ! FOR CONTACT MAP LOADED FROM A FILE
        ldisp=.false.     ! FOR DISPLACING ONE PROTEIN FROM ANOTHER
        lpbc=.false.      ! PERIODIC BOUNDARY COND-S IN EVERY DIRECTION
        lpbcx=.false.     ! PERIODIC BOUNDARY CONDITIONS IN X DIRECTION
        lpbcy=.false.     ! PERIODIC BOUNDARY CONDITIONS IN Y DIRECTION
        lpbcz=.false.     ! PERIODIC BOUNDARY CONDITIONS IN Z DIRECTION
        lcpb=.false.      ! PBC ALSO DURING GENERATING CHAINS BY SAW
        lunwrap=.false.   ! UNWRAP INPUT PDB FILE FROM PBC
        lsqpbc=.true.     ! PBC ALSO DURING SQUEEZING
        lampstrict=.true. ! AMPLITUDE OF OSCIL. IS SET TO ampstrict
        lmaxforce=.false. ! GATHER DATA ABOUT FMAX BETWEEN 2 THRESHOLDS
        lpullrel=.false.  ! PULL WITH CONSTANT VEL THEN RELEASE (mpull)
        ldet=.false.      ! FOR WRITING DETAILS IN THE CONTACT MAP
        lkmt=.true.       ! USE KMT ALGORITHM TO COMPUTE KNOT POSITIONS
        lrmsmax=.false.   ! STOP SIMULATION IF RMS IS BIGGER THAN RMSMAX
        ! SIMULATION ITSELF --------------------------------------------
        lparam=.true.     ! IF NO PARAMFILE IS SPECIFIED, LSIM=.TRUE.
        lcleanrst=.true.  ! IF CONTACTS SHOULD BE PURGED EVERY KRST TAU
        lsim=.false.      ! FOR MODEL WITH SIMPLIFIED CONTACT POTENTIAL
        ldisimp=.false.   ! HARMONIC DIHEDRAL POTENTIAL
        lmass=.false.     ! FOR TAKING INTO ACCOUNT DIFF. MASSES
        lchargend=.false. ! SHOULD PROTEIN ENDS BE CHARGED
        lsimpang=.false.  ! FOR MODEL WITH SIMPLIFIED ANGLE POT.
        lcoilang=.false.  ! RAND. COIL ANGLE POTENTIAL FOR GOMODEL
        lcoildih=.false.  ! RAND. COIL ANG. POT. FOR DIHEDRAL POTENTIAL
        lrepcoul=.false.  ! ONLY REPULSIVE PART OF COULOMB POT, NO ADIAB
        lecperm=.true.    ! USE CONSTANT REL. PERMITTIVITY FOR COULOMB
        lcpot=.true.   ! CUSTOM POTENTIAL NOT BASED ON NATIVE STRUCTURE
        lcintr=.true.  ! CUSTOM POTENTIAL ALSO FOR INTRACHAIN CONTACTS
        lpid=.false.   ! PSEUDO IMPROPER DIHEDRAL POTENTIAL
        lcospid=.false.! COSINE FUNCTION FOR THE PID POTENTIAL
        lepid=.false.  ! USE PID POTENTIAL TO CALCULATE ELECTROSTATICS
        lbar=.false.   ! SOLVATION BARRIER FOR PID POTENTIAL
        lmj=.false.    ! CONTACT DEPTH SCALED BY MIYAZAWA-J-LIKE MATRIX
        lchiral=.false.! POTENTIAL OF CHIRALITY
        langle=.true.  ! FOR MODEL WITH ANGLE POTENTIALS
        ldi=.true.     ! FOR DIHEDRAL TERM IN POTENTIALS, IF LANGLE IS T
        lii4=.true.    ! i+4 CONTACTS ARE INCLUDED BY DEFAULT, EXCEPT BB
        lsink=.false.  ! FOR SINK-LIKE POTENTIAL FOR NON-NATIVE CONTACTS
        lradii=.false. ! IF AMINO ACIDS SHOULD HAVE RADII FROM PARAMFILE
        lgln=.false.   ! FOR COUNTING ONLY GLNs IN CONTACTS
        lmrs=.false.   ! MORSE POTENTIAL FOR SSBONDS (DEFAULT HARMONIC)
        ldynss=.false. ! DYNAMIC FORMING OF DISULFIDE BRIDGES
        lsslj=.false.  ! TREAT SS BONDS AS NORMAL L-J CONTACTS
        lsselj=.false. ! TREAT SS BONDS AS EXCLUSIVE L-J CONTACTS
        lj1012=.false. ! FOR L-J 10-12 POTENTIALS
        ljwal=.true.   ! LENARD-JONES POTENTIAL FOR WALL INTERACTIONS
        lenetab=.false.! TABULARIZED VALUES FOR LOCAL POTENTIALS
        lsldh=.false.  ! SLOWLY TURN ON DIHEDRAL POTENTIAL IN THE START
        ! SIMULATION PARAMETERS - GENERAL ------------------------------
        c216=2.0d0**(1.d0/6.d0)
        Pi=dacos(-1.d0)! PI
        twopi=2.d0*Pi  ! 2PI
        unit=5.d0      ! LENGTH UNIT
        iseed=448      ! RANDOM SEED
        ntraj=1        ! NUMBER OF TRAJECTORIES
        mskip=0        ! SKIPPING STEPS [tau]
        af=4.d0        ! UNIT CELL SIZE FOR FCC WALLS [Angstrem]
        mstep=3000000  ! TOTAL SIMULATION TIME [tau]
        ktrest=10000   ! TIME AFTER SQUEEZING, BEFORE PULLING [TAU]
        kteql=0        ! EQUILIBRATION TIME BEFORE STRETCHING ETC [TAU]
        kwrite=100     ! HOW OFTEN TO PRINT ENERGY, 0: DON'T PRINT
        ksave= 1000    ! HOW OFTEN TO SAVE CONFORMATION, 0: DON'T SAVE
        kksave=1000    ! KSAVE BEFORE COOLDOWN (for kteql equilibration)
        krst=0         ! RESTART
        tstart=0.35    ! STARTING TEMPERATURE
        tend=0.35      ! ENDING TEMPERATURE
        tstep=0.05     ! TEMPERATURE STEP
        nen1=1         ! FIRST SIMULATED INDEX (INDICES<nen1 ARE FROZEN)
        rmsmax=10.0    ! STOP THE SIMULATION IF RMSD IS BIGGER THAN THIS
        ! SIMULATION PARAMETERS - DISULFIDE BRIDGES
        nssb=0         ! NUMBER OF NATIVE DISULFIDE BRIDGE, 0: NO BRIDGE
c        ksb(1,1)=56   ! FOR SSBONDS PUT RES NUMBER AS IS IN THE PROGRAM
c        ksb(2,1)=268  ! THIS ADDS 1 SSBOND BETWEEN RESIDUES 56 AND 268
        disul=1.d0     ! 5.d0 STRENGTH OF A DISULFIDE BRIDGE
        dislj=4.d0     ! HOW MANY TIMES L-J SSBOND IS STRONGER THAN L-J
        amrs=0.5       ! REVERSE WIDTH OF MORSE POTENTIAL [Angstrem^-1]
        smorse=15.0    ! SPRING CONSTANT EQUIVALENT TO MORSE [eps/A^2]
        rmrs=5.9       ! MORSE POT. MINIMUM FOR C-C PAIR [Angstrems]
        neimin=0       ! MAX. NUM. OF NEIGHBORS TO INCREASE COORD NUMBER
        neimaxdisul=9  ! MIN. NUM. OF NEIGHBORS TO PREVENT SS OXIDATION
        ! SIMULATION PARAMETERS - ANGLE POTENTIALS
        echi=1.0       ! COEFF. FOR CHIRALITY POTENTIAL [epsilon]
        CBA=30.d0!20.0COEFF FOR BOND ANGLE POTENTIALS,   20 for Clementi
        CDA=0.66 !1.0 COEFF FOR DIHED. ANGLE POT. (K1),   1 for Clementi
        CDB=0.66 !0.5 COEFF FOR DIHED. ANGLE POT. (K3), 0.5 for Clementi
        CDH=3.33 !1.0 COEFF FOR DIHED. ANGLE POT. (HARMONIC), 1 for Clem
        ! SIMULATION PARAMETERS - PULLING AND STRETCHING
        dar=1.d-10     ! DISPLACEMENT FOR COMPUTING THE FORCES
        coef=0.01      ! COEFFICIENT FOR CONST FORCE  [epsilon/(a*tau)]
        velo=0.005     ! PULLING VELOCITY IN ANGSTROMS/TAU
        veldist=12.d0  ! DISTANCE BEFORE ACQUIRING FULL VELOCITY
        mpull=100      ! TIME OF PULLING BEFORE RELEASE [tau]
        HH1=30.0       ! THE PULLING SPRING [epsilon/A^2]
        HH2=0.         ! THE PULLING SPRING [epsilon/A^2]
        cofp=1.0       ! COEFFICIENT OF THE PULLING FORCE
        dnaver = 0.0   ! LENGTH OVER WHICH FORCE IS AVERAGED (Angstrem)
        naver= 100     ! AVERAGING INSTATANEOUS F OVER MD STEPS (tau)
        kwforce=naver  ! 100     HOW OFTEN TO RECORD THE FORCE [tau]
        !       100000 is maximum simulation time divided by kwforce
        ! SIMULATION PARAMETERS - SIMULATION BOX
        tdens=0.001    ! TARGET DENSITY (IN RESIDUES/ANSTREM^3)
        sdens=0.0001   ! INITIAL DENSITY (THE BOX WILL SQUEEZE TO TDENS)
        densvelo=0.02  ! VELOCITY OF SQUEEZING FROM EVERY DIRECTION
        sepmin=10.0    ! FINAL DISTANCE BETWEEN WALLS [Angstrems]
        kbperiodmax=6  ! NUMBER OF OSCILLATIONS AFTER WHICH IS PULLING
        amp=0.1        ! AMPLITUDE OF OSCILLATIONS, RELATIVE TO WALLDIST
        ampstrict=10.0 ! AMPLITUDE OF OSCILLATIONS, STRICT [Angstrem]
        omega=0.0001   ! ANGULAR FREQUENCY OF OSCILLATIONS [1/tau]
c       period around 18000 tau means vmax=0.01 A/tau for large systems
        period=twopi/omega   ! TIME PERIOD OF OSCILLATIONS [tau]
        walmindst=5.d0 ! MINIMAL DISTANCE FOR STICKING TO WALL [A]
        kconnecttime=9 ! WHEN TO CONNECT BEADS TO WALL (1,3,5 OR 7)
        ipwn=-1        ! NR OF RES. STICKED TO WALL (NEGATIVE-AUTO)
        fwal = 4.0     ! WALL POTENTIAL COEFFICIENT
        ! SIMULATION PARAMETERS - TECHNICAL
        alphacos(1)=6.4! INVERSE WIDTH OF THE 1ST BB PID POTENTIAL
        alphacos(2)=6.0!3.6! INVERSE WIDTH OF THE 2ND BB PID POTENTIAL
        alphacos(3)=1.2! INVERSE WIDTH OF THE SS PID POTENTIAL
        rbb(1)=5.6     ! R_MIN OF THE 1ST BB PID POTENTIAL
        rbb(2)=6.2     ! R_MIN OF THE 2ND BB PID POTENTIAL
        psi0ss=-0.23   ! ANGLE MINIMUM OF THE 1ST BB PID POTENTIAL
        psi0bb(1)=1.05 ! ANGLE MINIMUM OF THE 2ND BB PID POTENTIAL
        psi0bb(2)=-1.44! ANGLE MINIMUM OF THE SS PID POTENTIAL
        epsbb=0.2      ! AMPLITUDE OF BB INTERACTIONS in PID [epsilon]
        rnei=7.5       ! CUT-OFF DISTANCE FOR NEIGHBOUR RES. [Angstrem]
        ad=2000        ! TIME FOR ADIABATIC POTENTIAL TURNOFF [steps]
        potcoeff=1.d0  ! SCALE OF LCPOT FORCES
        cntfct=1.3d0   ! MULTIPLIES SIGMA TO CHECK IF CONTACT PRESENT
c        factor=2.0    ! HOW MANY TIMES BACKB-BACKB. TOLERANCE IS BIGGER
        confcut=4.56   ! MIN. DISTANCE FOR CHAINS GENERATED BY SAW [A]
        bond=3.8       ! BOND LENGTH [Angstrem]
        ethi=1.0       ! THIRUMALAI OVERLAP FUNCTION [Angstrem] unused
        dnat=0.0       ! CUT-OFF DISTANCE FOR NATIVE CONTACTS [Angstrem]
        delta=0.005    ! INTEGRATION TIME STEP [tau]
        gamma=2.0      ! LANGEVIN PARAMETER [m/tau]
        H1 = 50.0      ! HARMONIC COEFF. [epsilon/A^2]
        H2 = 0.0       ! ANHARMONIC COEFF. [epsilon/A^4]
        verlcut=10.0   ! CUT-OFF DISTANCE FOR VERLET LIST[Angstrem]
        tolerance=0.d0 ! HOW FAR FROM THE MINIMUM CONTACT TURNS ON [%]
        sigma1(1)=1.5  ! FOR LEGACY (UNLESS LSIM IS T SHOULD NOT MATTER)
        sigma1(2)=1.5  ! FOR LEGACY (UNLESS LSIM IS T SHOULD NOT MATTER)
        sigma1(3)=1.5  ! FOR LEGACY (UNLESS LSIM IS T SHOULD NOT MATTER)
        sigma1(4)=5.0  ! L-J MINIMUM FOR BACKBONE CONTACTS [Angstrems]
        sigma1(5)=7.5  ! L-J MINIMUM FOR SIDECHAIN CONTACTS [Angstrems]
        sigma1(6)=6.6  ! L-J MINIMUM FOR BCKB-SDCH CONTACTS [Angstrems]
        sigma1(7)=6.6  ! L-J MINIMUM FOR SDCH-BCKB CONTACTS [Angstrems]
        sigma1(8)=6.1  ! L-J MINIMUM FOR i,i+4 CONTACTS [Angstrems]
        screend=10.0   ! ELECTROSTATIC SCREENING LENGTH [Angstrems]
        coul=85.0      ! CONSTANT FOR COULOMBIC INTERACTION [eps*A*A]
        if(lecperm) coul=2.63 ! 210 IF REL. PERMITTIVITY=1 [eps*A]
        cut=5.0        ! CUT-OFF DISTANCE FOR REPULSIVE TERM [Angstrem]
        rcut=18.d0     ! CUT-OFF DISTANCE FOR THE REST [Angstrem]
        kmaxbckb=2     ! MAXIMAL NUMBER OF BACKBONE CONTACTS
        bckbmin=0.75   ! MINIMUM VALUE OF BACKBONE ANGLE
        bckb2min=0.92  ! MINIMUM VALUE OF BACKBONE ANGLE
        sdchnmax=0.5   ! MAXMUM VALUE OF SIDECHAIN ANGLE
        ndomain=1      ! NUMBER OF DOMAIN (FOR TITIN ONLY)
      ! READING VARIABLES --------------------------------------------
      if(iargc().gt.0) then     ! reading parameters from external file
         call getarg(1,arg)     ! for more than 1 args, do i_arg=1,iargc()
         open(7,file=arg,status='old') ! no spaces in filenames allowed
         write(*,*) 'RUNNING INPUTFILE ',arg
 11      read(7,'(a)',end=12, err=12) buffer
         if(buffer(1:4).eq.'velo') then
            read(buffer(5:),*) velo
         elseif(buffer(1:7).eq.'pdbfile') then
            read(buffer(8:),*) pdbfile
         elseif(buffer(1:7).eq.'seqfile') then
            read(buffer(8:),*) seqfile
         elseif(buffer(1:7).eq.'outfile') then
            read(buffer(8:),*) outfile
         elseif(buffer(1:7).eq.'savfile') then
            read(buffer(8:),*) savfile
         elseif(buffer(1:7).eq.'rstfile') then
            read(buffer(8:),*) rstfile
         elseif(buffer(1:7).eq.'mapfile') then
            read(buffer(8:),*) mapfile
            lwritemap=.true.
         elseif(buffer(1:7).eq.'verlcut') then
            read(buffer(8:),*) verlcut
         elseif(buffer(1:6).eq.'dnaver') then
            read(buffer(7:),*) dnaver
         elseif(buffer(1:6).eq.'factor') then
            read(buffer(7:),*) factor
         elseif(buffer(1:9).eq.'lconstvol') then
            read(buffer(10:),*) lconstvol
         elseif(buffer(1:9).eq.'lstartpdb') then
            read(buffer(10:),*) lstartpdb
         elseif(buffer(1:8).eq.'lslowdih') then
            read(buffer(9:),*) lsldh
         elseif(buffer(1:5).eq.'lsldh') then
            read(buffer(6:),*) lsldh
         elseif(buffer(1:5).eq.'lcmap') then
            read(buffer(6:),*) lcmap
         elseif(buffer(1:5).eq.'lsink') then
            read(buffer(6:),*) lsink
         elseif(buffer(1:4).eq.'lpid') then
            read(buffer(5:),*) lpid
         elseif(buffer(1:4).eq.'lfcc') then
            read(buffer(5:),*) lfcc
         elseif(buffer(1:4).eq.'lii4') then
            read(buffer(5:),*) lii4
         elseif(buffer(1:4).eq.'lbar') then
            read(buffer(5:),*) lbar
         elseif(buffer(1:4).eq.'lrst') then
            read(buffer(5:),*) lrst
         elseif(buffer(1:6).eq.'lradii') then
            read(buffer(7:),*) lradii
         elseif(buffer(1:6).eq.'lcdnat') then
            read(buffer(7:),*) lcdnat
         elseif(buffer(1:5).eq.'ljwal') then
            read(buffer(6:),*) ljwal
         elseif(buffer(1:5).eq.'lepid') then
            read(buffer(6:),*) lepid
         elseif(buffer(1:5).eq.'disul') then
            read(buffer(6:),*) disul
         elseif(buffer(1:5).eq.'dislj') then
            read(buffer(6:),*) dislj
         elseif(buffer(1:5).eq.'epsbb') then
            read(buffer(6:),*) epsbb
         elseif(buffer(1:5).eq.'gamma') then
            read(buffer(6:),*) gamma
         elseif(buffer(1:5).eq.'cmapf') then
            read(buffer(6:),*) cmapf
            lcmap=.true.
         elseif(buffer(1:9).eq.'tolerance') then
            read(buffer(10:),*) tolerance
         elseif(buffer(1:9).eq.'lcleanrst') then
            read(buffer(10:),*) lcleanrst
         elseif(buffer(1:9).eq.'lwritemap') then
            read(buffer(10:),*) lwritemap
         elseif(buffer(1:9).eq.'lwritexyz') then
            read(buffer(10:),*) lwritexyz
         elseif(buffer(1:7).eq.'lwrtang') then
            read(buffer(8:),*) lwrtang
         elseif(buffer(1:7).eq.'lcospid') then
            read(buffer(8:),*) lcospid
         elseif(buffer(1:7).eq.'lecperm') then
            read(buffer(8:),*) lecperm
            coul=210.0
         elseif(buffer(1:7).eq.'lrmsmax') then
            read(buffer(8:),*) lrmsmax
         elseif(buffer(1:11).eq.'wallmindist') then
            read(buffer(12:),*) walmindst
         elseif(buffer(1:12).eq.'lfromscratch') then
            read(buffer(13:),*) lfrmscrtch
         elseif(buffer(1:11).eq.'lwritegomap') then
            read(buffer(12:),*) lwritego
         elseif(buffer(1:8).eq.'lwritego') then
            read(buffer(9:),*) lwritego
         elseif(buffer(1:7).eq.'lunwrap') then
            read(buffer(8:),*) lunwrap
         elseif(buffer(1:10).eq.'lfrmscrtch') then
            read(buffer(11:),*) lfrmscrtch
         elseif(buffer(1:10).eq.'lampstrict') then
            read(buffer(11:),*) lampstrict
         elseif(buffer(1:12).eq.'kconnecttime') then
            read(buffer(13:),*) kconnecttime
         elseif(buffer(1:8).eq.'displace') then
            read(buffer(9:),*) iprota,jprotb,away
            ldisp=.true.
         elseif(buffer(1:8).eq.'lpullfin') then
            read(buffer(9:),*) lpullfin
         elseif(buffer(1:8).eq.'lrepcoul') then
            read(buffer(9:),*) lrepcoul
         elseif(buffer(1:8).eq.'kmaxbckb') then
            read(buffer(9:),*) kmaxbckb
         elseif(buffer(1:7).eq.'screend') then
            read(buffer(8:),*) screend
         elseif(buffer(1:7).eq.'bckbmin') then
            read(buffer(8:),*) bckbmin
         elseif(buffer(1:7).eq.'bckbmax') then
            read(buffer(8:),*) bckbmax
         elseif(buffer(1:8).eq.'sdchnmin') then
            read(buffer(9:),*) sdchnmin
         elseif(buffer(1:8).eq.'sdchnmax') then
            read(buffer(9:),*) sdchnmax
         elseif(buffer(1:8).eq.'bckb2min') then
            read(buffer(9:),*) bckb2min
         elseif(buffer(1:8).eq.'bckb2max') then
            read(buffer(9:),*) bckb2max
         elseif(buffer(1:8).eq.'lrestart') then
            read(buffer(9:),*) lrst
         elseif(buffer(1:8).eq.'lsimpang') then
            read(buffer(9:),*) lsimpang
         elseif(buffer(1:8).eq.'lcoilang') then
            read(buffer(9:),*) lcoilang
         elseif(buffer(1:8).eq.'lcoildih') then
            read(buffer(9:),*) lcoildih
         elseif(buffer(1:8).eq.'lallatom') then
            read(buffer(9:),*) lallatom
         elseif(buffer(1:9).eq.'lchargend') then
            read(buffer(10:),*) lchargend
         elseif(buffer(1:9).eq.'ampstrict') then
            read(buffer(10:),*) ampstrict
         elseif(buffer(1:8).eq.'potcoeff') then
            read(buffer(9:),*) potcoeff
         elseif(buffer(1:8).eq.'densvelo') then
            read(buffer(9:),*) densvelo
         elseif(buffer(1:6).eq.'cntfct') then
            read(buffer(7:),*) cntfct
         elseif(buffer(1:6).eq.'lsqpbc') then
            read(buffer(7:),*) lsqpbc
         elseif(buffer(1:5).eq.'lpbcx') then
            read(buffer(6:),*) lpbcx
         elseif(buffer(1:5).eq.'lpbcy') then
            read(buffer(6:),*) lpbcy
         elseif(buffer(1:5).eq.'lpbcz') then
            read(buffer(6:),*) lpbcz
         elseif(buffer(1:4).eq.'lobo') then
            read(buffer(5:),*) lobo
         elseif(buffer(1:4).eq.'lsim') then
            read(buffer(5:),*) lsim
         elseif(buffer(1:4).eq.'ldet') then
            read(buffer(5:),*) ldet
         elseif(buffer(1:4).eq.'lpbc') then
            read(buffer(5:),*) lpbc
         elseif(buffer(1:4).eq.'lcpb') then
            read(buffer(5:),*) lcpb
         elseif(buffer(1:4).eq.'lkmt') then
            read(buffer(5:),*) lkmt
         elseif(buffer(1:4).eq.'ldih') then
            read(buffer(5:),*) ldi
         elseif(buffer(1:7).eq.'ldisimp') then
            read(buffer(8:),*) ldisimp
         elseif(buffer(1:7).eq.'lmedian') then
            read(buffer(8:),*) lmedian
         elseif(buffer(1:7).eq.'lenetab') then
            read(buffer(8:),*) lenetab
         elseif(buffer(1:3).eq.'CBA') then
            read(buffer(4:),*) CBA
         elseif(buffer(1:3).eq.'CDA') then
            read(buffer(4:),*) CDA
         elseif(buffer(1:3).eq.'CDB') then
            read(buffer(4:),*) CDB
         elseif(buffer(1:3).eq.'lmj') then
            read(buffer(4:),*) lmj
         elseif(buffer(1:4).eq.'coul') then
            read(buffer(5:),*) coul
         elseif(buffer(1:4).eq.'nen1') then
            read(buffer(5:),*) nen1
         elseif(buffer(1:4).eq.'fwal') then
            read(buffer(5:),*) fwal
         elseif(buffer(1:5).eq.'adiab') then
            read(buffer(6:),*) ad
         elseif(buffer(1:5).eq.'lvelo') then
            read(buffer(6:),*) lvelo
         elseif(buffer(1:5).eq.'lmass') then
            read(buffer(6:),*) lmass
         elseif(buffer(1:6).eq.'langle') then
            read(buffer(7:),*) langle
         elseif(buffer(1:7).eq.'lchiral') then
            read(buffer(8:),*) lchiral
         elseif(buffer(1:7).eq.'lconftm') then
            read(buffer(8:),*) lconftm
         elseif(buffer(1:7).eq.'confcut') then
            read(buffer(8:),*) confcut
         elseif(buffer(1:7).eq.'lnatend') then
            read(buffer(8:),*) lnatend
         elseif(buffer(1:7).eq.'lposcrd') then
            read(buffer(8:),*) lposcrd
         elseif(buffer(1:7).eq.'lthermo') then
            read(buffer(8:),*) lthermo
         elseif(buffer(1:7).eq.'ldelrst') then
            read(buffer(8:),*) ldelrst
         elseif(buffer(1:6).eq.'lnowal') then
            read(buffer(7:),*) lnowal
         elseif(buffer(1:10).eq.'lsawconftm') then
            read(buffer(11:),*) lsawconftm
         elseif(buffer(1:4).eq.'lmrs') then
            read(buffer(5:),*) lmrs
         elseif(buffer(1:6).eq.'lmorse') then
            read(buffer(7:),*) lmrs
         elseif(buffer(1:6).eq.'lparam') then
            read(buffer(7:),*) lparam
         elseif(buffer(1:6).eq.'ldynss') then
            read(buffer(7:),*) ldynss
         elseif(buffer(1:6).eq.'lsselj') then
            read(buffer(7:),*) lsselj
         elseif(buffer(1:6).eq.'lcintr') then
            read(buffer(7:),*) lcintr
         elseif(buffer(1:5).eq.'lsslj') then
            read(buffer(6:),*) lsslj
         elseif(buffer(1:5).eq.'ldens') then
            read(buffer(6:),*) ldens
         elseif(buffer(1:5).eq.'tdens') then
            read(buffer(6:),*) tdens
         elseif(buffer(1:5).eq.'sdens') then
            read(buffer(6:),*) sdens
         elseif(buffer(1:4).eq.'krst') then
            read(buffer(5:),*) krst
         elseif(buffer(1:5).eq.'kteql') then
            read(buffer(6:),*) kteql
         elseif(buffer(1:5).eq.'ksave') then
            read(buffer(6:),*) ksave
         elseif(buffer(1:6).eq.'kksave') then
            read(buffer(7:),*) kksave
         elseif(buffer(1:6).eq.'ktrest') then
            read(buffer(7:),*) ktrest
         elseif(buffer(1:6).eq.'kwrite') then
            read(buffer(7:),*) kwrite
         elseif(buffer(1:6).eq.'sepmin') then
            read(buffer(7:),*) sepmin
         elseif(buffer(1:6).eq.'lwall ') then
            read(buffer(7:),*) lwal
         elseif(buffer(1:6).eq.'lwalls') then
            read(buffer(7:),*) lwals
         elseif(buffer(1:6).eq.'lshear') then
            read(buffer(7:),*) lshear
         elseif(buffer(1:9).eq.'paramfile') then
            read(buffer(10:),*) paramfile
            lparam=.true.
         elseif(buffer(1:10).eq.'loscillate') then
            read(buffer(11:),*) loscillate
         elseif(buffer(1:11).eq.'kbperiodmax') then
            read(buffer(12:),*) kbperiodmax
         elseif(buffer(1:5).eq.'omega') then
            read(buffer(6:),*) omega
            period=twopi/omega
         elseif(buffer(1:6).eq.'period') then
            read(buffer(7:),*) period
            omega=twopi/period
         elseif(buffer(1:5).eq.'iseed') then
            read(buffer(6:),*) iseed
         elseif(buffer(1:5).eq.'ntraj') then
            read(buffer(6:),*) ntraj
         elseif(buffer(1:5).eq.'mstep') then
            read(buffer(6:),*) mstep
         elseif(buffer(1:5).eq.'lcpot') then
            read(buffer(6:),*) lcpot
         elseif(buffer(1:4).eq.'lpdb') then
            read(buffer(5:),*) lpdb
         elseif(buffer(1:4).eq.'ipwn') then
            read(buffer(5:),*) ipwn
         elseif(buffer(1:4).eq.'bond') then
            read(buffer(5:),*) bond
         elseif(buffer(1:4).eq.'dnat') then
            read(buffer(5:),*) dnat
         elseif(buffer(1:4).eq.'rcut') then
            read(buffer(5:),*) rcut
         elseif(buffer(1:4).eq.'cofp') then
            read(buffer(5:),*) cofp
         elseif(buffer(1:4).eq.'rnei') then
            read(buffer(5:),*) rnei
         elseif(buffer(1:6).eq.'neimin') then
            read(buffer(7:),*) neimin
         elseif(buffer(1:11).eq.'neimaxdisul') then
            read(buffer(7:),*) neimaxdisul
         elseif(buffer(1:3).eq.'HH1') then
            read(buffer(4:),*) HH1
         elseif(buffer(1:3).eq.'amp') then
            read(buffer(4:),*) amp
         elseif(buffer(1:3).eq.'cut') then
            read(buffer(4:),*) cut
         elseif(buffer(1:5).eq.'acos1') then
            read(buffer(6:),*) alphacos(1)
         elseif(buffer(1:5).eq.'acos2') then
            read(buffer(6:),*) alphacos(2)
         elseif(buffer(1:5).eq.'acos3') then
            read(buffer(6:),*) alphacos(3)
         elseif(buffer(1:6).eq.'psi0ss') then
            read(buffer(7:),*) psi0ss
         elseif(buffer(1:7).eq.'psi0bb1') then
            read(buffer(8:),*) psi0bb(1)
         elseif(buffer(1:7).eq.'psi0bb2') then
            read(buffer(8:),*) psi0bb(2)
         elseif(buffer(1:4).eq.'rbb1') then
            read(buffer(5:),*) rbb(1)
         elseif(buffer(1:4).eq.'rbb2') then
            read(buffer(5:),*) rbb(2)
         elseif(buffer(1:4).eq.'bbrm') then
            read(buffer(5:),*) sigma1(4)
         elseif(buffer(1:4).eq.'ssrm') then
            read(buffer(5:),*) sigma1(5)
         elseif(buffer(1:4).eq.'bsrm') then
            read(buffer(5:),*) sigma1(6)
            sigma1(7)=sigma1(6)
         elseif(buffer(1:4).eq.'i4rm') then
            read(buffer(5:),*) sigma1(8)
         elseif(buffer(1:4).eq.'temp') then
            read(buffer(5:),*) tstart
            tend=tstart
         elseif(buffer(1:4).eq.'teql') then
            read(buffer(5:),*) teql
            lteql=.true.
         elseif(buffer(1:4).eq.'tend') then
            read(buffer(5:),*) tend
         elseif(buffer(1:5).eq.'tstep') then
            read(buffer(6:),*) tstep
         elseif(buffer(1:6).eq.'tstart') then
            read(buffer(7:),*) tstart
         elseif(buffer(1:7).eq.'klenstr') then
            read(buffer(8:),*) klenstr
         elseif(buffer(1:4).eq.'file') then
            read(buffer(5:),*) filname
            write(stafile,*) '(a',klenstr,',a)'
            write(outfile,stafile) filname,'.out'
            write(mapfile,stafile) filname,'.map'
            write(savfile,stafile) filname,'.pdb'
         else                   ! writing to console, unless file indexes 5 or 6 are in use
            write(*,*) 'UNRECOGNIZED OPTION: ',buffer
         endif
         goto 11
 12      close(7)
      endif
      if(lparam) call load_paramfile(paramfile)
      do i=1,len                ! LFROMPDB MUST BE ZEROED BEFORE LOADING CMAPS
         the0(i)=-1.0
         phi0(i)=0.0
         lfrompdb(i)=lsimpang
      enddo
      if(lpbc) then
         lpbcx=.true.
         lpbcy=.true.
         lpbcz=.true.
      endif
      klont=0
      if(lwal.and..not.lnowal) lpbcz=.false. ! no PBC in Z
      if(lsselj) lmrs=.true.    ! lmorse is overriden by lsselj
      if(dnaver.gt.0.0) naver=nint(dnaver/velo)
      nratvel=nint(veldist/velo) ! NUMBER OF STEPS TO PULLING VELOCITY
      nratveld=nint(veldist/densvelo) ! N OF STEPS TO SQUEEZE VELOCITY
      open(1,file=outfile,status='unknown')
      if(lwritemap) open(22,file=mapfile,status='unknown')
      if(ksave.ne.0) open(2,file=savfile,status='unknown')

      write(1,*)'#I,I+2 CONTACTS PURELY REPULSIVE'

!     SCALE LENGTHS
      cut=cut/unit              ! REPULSIVE INTERACTIONS CUT-OFF
      sfact=(1.d0+tolerance)*c216 !2.0d0**(1.d0/6.d0)
      sigma0=cut/c216           !*0.5d0**(1.d/6.d)! REPULSIVE INTERACTIONS SIGMA
c     dfold=dfold/unit
      rcut = rcut/unit          ! OTHER POTENTIALS CUT-OFF
      rnei=rnei/unit
      rcutsq=rcut*rcut
      cutsq=cut*cut
      rneisq=rnei*rnei
c     bckbmin=bckbmin*bckbmin
c     bckb2min=bckb2min*bckb2min
      sdchnmax=sdchnmax*sdchnmax
      verlcut=verlcut/unit
      confcut=confcut/unit
      vrcut2sq=verlcut*verlcut/4.d0
      dnat=dnat/unit
      H1=H1*unit*unit
      H2=H2*unit**4
      HH1=HH1*unit*unit
      HH2=HH2*unit**4
      ethi=ethi/unit
      bond=bond/unit
      af=af/unit
      rbb(1)=rbb(1)/unit*0.5d0**(1.d0/6.d0)
      rbb(2)=rbb(2)/unit*0.5d0**(1.d0/6.d0)
      if(.not.lcospid) then
         alphacos(1)=alphacos(1)/Pi
         alphacos(2)=alphacos(2)/Pi
         alphacos(3)=alphacos(3)/Pi
         vmp=1.d0
      else
         vmp=Pi
      endif
      ampstrict=ampstrict/unit
      sepmin=sepmin/unit
      sigma1(9)=walmindst       ! for wall interactions via L-J potential
      walmindst=walmindst/unit
      tdens=tdens*unit**3
      sdens=sdens*unit**3
      vpull=velo/unit
      vpull=vpull*delta         ! corresponding infinitesimal displacement
      vpulld=densvelo/unit*delta
      vtarget=vpull             ! target velocity acquired after nratvel steps
!     SETUP TABLE OF TEMPERATURE
      nt=nint(abs(tstart-tend)/tstep)+1
      if(nt.gt.150) nt=150
      ttstep=tstep
      if(tstart.gt.tend) ttstep=-tstep
      do i=1,nt
         ttab(i)=tstart+(i-1)*ttstep
      enddo
c     bckbmax=bckbmax-bckbmin
c     bckb2max=bckb2max-bckb2min
c     sdchnmin=sdchnmax-sdchnmin

!     LOAD PROTEIN CONFORMATION
      if(lpdb.or.lstartpdb) then
         write(1,'(/,a,2x,a,/)')'#PDB FILE =',pdbfile
         call load_protein(pdbfile,lunwrap)
         do i=1,men
            x0(i)=xn(i)
            y0(i)=yn(i)
            z0(i)=zn(i)
         enddo
      else
         write(1,'(/,a,2x,a,/)')'#SEQ FILE =',seqfile
         call load_sequence(seqfile)
         call confstart(sdens,confcut)
         do i=1,men
            xn(i)=x0(i)
            yn(i)=y0(i)
            zn(i)=z0(i)
         enddo
      endif

      if(lallatom) then
         write(1,'(a)')'#CONSTRUCT CONTACT-MAP BASED ON ALL-ATOM'
         call compute_contact_map(pdbfile)
      elseif(klont.eq.0) then
         write(1,'(a,f6.2,a)')
     $        '#CONSTRUCT CONTACTMAP BASED ON CUT-OFF =',
     +        dnat*unit,' ANGSTROM'
         call compute_cmap(dnat,lcdnat)
      else
         write(1,'(a)')'#CONTACT-MAP BASED ON THE SEQUENCE FILE'
      endif

      if(lcmap) then
         write(1,'(a,a)')'#CONSTRUCT CONTACT-MAP BASED ON FILE ',cmapf
         call load_cmap(cmapf)
      endif
      if(ipwn.lt.0) then
         ipwn=(men-nint((men**(1./3)-walmindst*tdens**(1./3))**3))/3
      endif
      ipwn=2*(ipwn/2)
!     BUILD TITIN
      if(ndomain.gt.1) then
         write(1,'(a,i10)')'#NUMBER OF DOMAINS',ndomain
         call build_titin(ndomain)
         write(1,'(a)')'#DO NOT ALLOW CONTACTS BETWEEN DOMAINS'
         call interdomain(ndomain)
      endif

      ip1=1
      ip2=men
      if(lmass) then
         write(1,'(a)')'#CONSIDERING AMINO ACID MASSES'
         call amino_acid_mass
      else
         do i=1,men
            rmas(i)=1.d0
         enddo
      endif
      jq=1
      part=men-nen1+1
      mcmr=men*(men-1)/2-(men-1)
      mchi=(men*men-5*men+6)/2
      targetvolume=men/tdens
      reemax=3.0*targetvolume**(1./3)
      xmin=x0(1)
      ymin=y0(1)
      zmin=z0(1)
      xmax=x0(1)
      ymax=y0(1)
      zmax=z0(1)
      if(.not.lcpb) then
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
      endif
      if(lobo) kconnecttime=3
      call update_verlet_list(verlcut,nen1)
      corder=0
      icor=0
      do k=1,kront
C     if(krist(3,k).ne.0) then
         corder=corder + abs(krist(1,k)-krist(2,k))
         icor=icor+1
C     endif
      enddo
      if(icor.gt.0) then
         corder=corder/icor
         corder=corder/men
         write(1,'(a,f8.4)') '#RELATIVE TOTAL CONTACT ORDER', corder
      else
         write(1,'(a)') '#NO NON-REPULSIVE CONTACTS'
      endif
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
            write(1,*)'#RESIDUE ',ksb(ks,1),
     $           ' IS NOT A CYSTEINE. PLS CHECK!'
            stop
         endif
         if(aseq(i2).ne.'CYS') then
            write(1,*)'#RESIDUE ',ksb(ks,2),
     $           ' IS NOT A CYSTEINE. PLS CHECK!'
            stop
         endif
         icheck=0
         if(.not.(lsslj.or.lsselj)) then
            do k=1,klont
               ki1=klist(1,k)
               ki2=klist(2,k)
               if((ki1.eq.i1.and.ki2.eq.i2)
     $              .or.(ki1.eq.i2.and.ki2.eq.i1)) then
               klist(3,k)=sign(631,klist(3,k)) ! SS BONDS HAVE VALUE +-631
               icheck=1
               write(1,'(a,i4,5x,2(a3,i4,3x))')'#NATIVE SS BOND',ks,
     +              aseq(i1),iseq(i1),aseq(i2),iseq(i2)
            endif
         enddo
         if(icheck.eq.0) then
            write(1,'(a,i4,5x,2(a3,i4,3x))')
     $           '#SS BOND COULD NOT BE MADE',ks,
     +           aseq(i1),iseq(i1),aseq(i2),iseq(i2)
         endif
      endif
 47   continue
c     write(1,'(a,f7.1)')'SS BOND STRENGTH',disul
      if(lcpot) then
         write(1,'(a)')
     $        '#USING CUSTOM ATTRACTIVE L-J CONTACT POTENTIAL!'
         if(lsselj) then
            write(1,'(a)')
     $           '#USING L-J POTENTIAL FOR NON-NATIVE SS BONDS!'
            write(1,'(a,f6.2,a,f6.2)')'#R_LJ=',rmrs,' DEPTH=',dislj
         else
            write(1,'(a)')
     $           '#USING MORSE POTENTIAL FOR NON-NATIVE SS BONDS!'
            write(1,'(a,f6.2,a,f6.2)')'#R_MORSE=',rmrs,' A_MORSE=',amrs
         endif
         if(lsim.or.(lcpot.and..not.lparam)) then
            lsim=.true.
            write(1,'(a,f8.2)')
     $           '#SS CONTACT EQUILIBRIUM DISTANCE=',sigma1(5)
         else
            write(1,'(a,a)')
     $           '#CONTACTS BASED ON DATA FROM FILE ',paramfile
         endif
         write(1,'(a,f6.2)') '#USING DEBYE SCREENING LENGTH',screend
      endif

!     PUT CUSTOM POTENTIALS HERE

      smorse=smorse*unit*unit
      amrs=amrs*unit
      screend=screend/unit
      coul=coul/unit
      if(.not.lecperm) coul=coul/unit
      dmrs=smorse/amrs**2
      rmrs=rmrs/unit
      sigma1(88)=rmrs*0.5d0**(1.d0/6.d0)
      sigss=sigma1(88)
      do i=4,9
         sigma1(i)=sigma1(i)/unit*0.5d0**(1.d0/6.d0)
      enddo
c     sigma1(2)=dalQQ*0.5d0**(1.d0/6.d0)

      write(1,'(a)')'#USING HARMONIC POTENTIALS FOR NATIVE SS BONDS!'
      olno=ran2(iseed)          ! FOR LEGACY (COMPARE OLD VERSION, SAME SEED)
      write(1,'(a,2(a,f6.2))')'#USING ANHARMONIC POTENTIAL',
     +     '   H1 =',H1/unit/unit,'   H2 =',H2/unit**4

      if(lchiral) then
         write(1,'(a)')'#USING CHIRALITY POTENTIALS'
         call model_chirality
      endif

      if(langle) then
         if(lsimpang) ldi=.false.
         if(ldi) then
            write(1,'(a)')
     $           '#USING POTENTIALS FOR BOND AND DIHEDRAL ANGLES'
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
         if(ldynss .and. inameseq(i).eq.4 .and. inameseq(j).eq.4
     +        .and. abs(klist(3,k)).eq.631) then
            ijdiff=0
            write(1,'(a,2i4,5x,2(a3,i4,3x))')'#DYNAMIC SS BOND',i,j,
     +           aseq(i),iseq(i),aseq(j),iseq(j)
         endif
         if(.not.lconect(i)) ijdiff=5
         if(ijdiff.ge.3) then   ! 4) then
            km=km+1
            klist(1,km)=i
            klist(2,km)=j
            klist(3,km)=klist(3,k)
         endif
      enddo
      klont=km

      ngln=0
      do ib=1,men
         x0(ib)=xn(ib)          ! assign native coordinates to actual ones
         y0(ib)=yn(ib)
         z0(ib)=zn(ib)
         adia(ib)=0             ! zero the table for adiabatic turning off
         z0temp(ib)=zn(ib)      ! initialize the table to be sorted by z
         if(ljwal) z0temp(ib) = 0
         ksorted(ib)=ib         ! indexes to be sorted
         xpul(ib)=0.0           ! initialize ref. values of pulled resid.
         ypul(ib)=0.0
         zpul(ib)=0.0
         nei(1,ib)=0            ! set neighbour counter to zero
         nei(2,ib)=0
         if(aseq(ib).eq.'GLN') ngln=ngln+1 ! count glutamines
         ksdchns(ib)=ksdchn(inameseq(ib),1) ! type of sidechain
!     khbful(3,ib)=ksdchn(inameseq(ib),2) !nr of potential hbonds
      enddo
      if(lchargend) then
         do ic=1,nchains
            if(ksdchns(menchain(ic)+1).eq.4) then
               ksdchns(menchain(ic)+1)=2 !N-terminal is zwitterion
            else
               ksdchns(menchain(ic)+1)=5 !N-terminal of protein is positive
            endif
            if(ksdchns(menchain(ic+1)).eq.5) then
               ksdchns(menchain(ic+1))=2 ! C-terminal is zwitterion
            else
               ksdchns(menchain(ic+1))=4 ! C-terminal is negative
            endif
         enddo
      endif

      bond=0.d0                 ! length of the Ca-Ca bond is averaged over all bonds
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
      bond=bond/(men-1)         ! this averaging previously was in gopotential

      write(1,'(a)')'#USING THE GO-LIKE 6-12 LJ POTENTIALS'
      if(lpdb) call gopotential(asigma)
      if(lwritego) then
         call print_cmap(22,0)
         close(1)
         close(2)
         close(22)
         stop
      endif

      call prepare(edsg)
      call evalgo(edsg,chi)
      if(lwal.and..not.lnowal) call evalwall(edsg)
      if(lpid) then
         call evalimproper(edsg,lcospid,epsbb)
      else
         call evalcpot(edsg)
      endif
      if(langle.or.lwritemap) call evalangles(edsg,lsldh,1.d0)
      if(lchiral) then
         call eval_chirality(enechi)
         edsg=edsg+enechi
      endif
      if(lwritemap) call print_map(22,0)
      if(lpullrel) lvelo=.TRUE.

      write(1,'(/,a,i10)') '#TOTAL PROTEIN LENGTH      ',men
      write(1,'(a,i10)')  '#NUMBER OF NATIVE CONTACTS ',klont
      if(lpdb) then
         write(1,'(a,f10.4)')'#AVERAGE LENGTH OF CONTACTS',asigma*unit
      else
         write(1,'(a)')'#NO PDB FILE USED FOR GO MODEL CONSTRUCTION'
      endif
      write(1,'(a,f10.2)')'#ENERGY OF NATIVE STATE    ',edsg
      write(1,*)
      if(lforce) write(1,'(a,f7.4)')'#USING AFM FORCE  ',coef
      if(lvelo.or.lforce) write(1,'(a,2f10.2)')
     +     '#PULLING SPRING CONSTANTS  ',HH1/unit/unit,HH2/unit**4
      if(naver.ne.0)
     +     write(1,'(a,i8,a)') '#FORCE AVERAGED OVER ',naver,' TAU'
      write(1,'(a,f7.2)') '#VERLET LIST CUTOFF ',verlcut*unit
      if(loscillate) then
         write(1,'(a,f8.6,a,f9.1)')
     +        '#ANGULAR FREQUENCY ',omega,' PERIOD ',period
      else if(lvelo .OR. lwal) then
         write(1,'(a,f7.4)') '#CONSTANT VELOCITY',velo
      endif
      if(ldens) then
         write(1,'(a,f8.5)')'#SQEEZING VELOCITY',densvelo
         write(1,'(a,f8.5,a)')'#TARGET DENSITY ',
     +        tdens/(unit**3),' RESIDUES/A^3'
      endif
      if(lvelo) write(1,'(a,i10,a)') '#KWFORCE ',kwforce,' TAU'
      if(lwal) write(1,'(a,f7.4)')'#WALL POTENTIAL COEFFICIENT ',fwal
      if(lpullrel) then
         write(1,'(a)')'#STUDY PULLING AT CONSTANT VELOCITY THEN STOP'
         write(1,'(a)')'#PULLING AT CERTAIN DISTANCE'
         write(1,'(a,i10)')'#PULLING TIME [tau]        ',mpull
      else
         if(lmedian) write(1,'(a)')'#COMPUTING MEDIAN FOLDING TIMES'
         if(lthermo) write(1,'(a)')
     $        '#COMPUTING THERMODYNAMIC PROPERTIES'
         if(lunfold) write(1,'(a)')'#STUDYING UNFOLDING'
      endif
      if(lconftm) then
         if(lunfold) then
            write(1,'(a)')
     $           '#COMPUTING AVERAGED TIMES FOR CONTACT BREAKING'
            write(1,'(a)')
     $           '#AND AVERAGED LIFE TIMES OF CONTACTS'
         else
            write(1,'(a)')
     $           '#COMPUTING AVERAGED TIMES FOR CONTACT FORMATION'
         endif
      endif
      write(1,'(/,a,f10.3)')'#DELTA    =',delta
      write(1,'(a,f10.3)')'#GAMMA    =',gamma
      write(1,'(/,a,i10)') '#NUMBER OF TRAJECTORIES ',ntraj
      write(1,'(a,i10)') '#SIMULATION TIME        ',mstep
      write(1,'(a,i10)') '#SKIPPING STEPS         ',mskip
      write(1,'(a,i10)') '#EQUILIBRATION TIME     ',kteql
      write(1,'(a,i10)')'#RANDOM SEED            ',iseed
      write(1,'(/,a,7x,3f7.3)')'#TSTART TEND TSTEP',tstart,tend,tstep
      write(1,'(a,i6)')'#NUMBER OF TEMPERATURE STEPS',nt

!     RESCALE TIME PARAMETERS
      kunit=nint(1.d0/delta)    ! if delta = 0.005, this is 200
      naver=naver*kunit
      krst=krst*kunit
      mskip=mskip*kunit
      mstep=mstep*kunit
      ktrest=ktrest*kunit
      kteql=kteql*kunit
      nratvel=nratvel*kunit
      nratveld=nratveld*kunit
      kwrite=kwrite*kunit
      ksave=ksave*kunit
      kksave=kksave*kunit
      mpull=mpull*kunit
      kwforce=kwforce*kunit
      kwquarterperiod=nint(0.25*period)*kunit ! 1/4 of period
      omega=omega*delta
!     SCALE FACTORS FOR VELOCITIES DURING EQUILIBRATION
      delsq=delta*delta
      deltsq=0.5d0*delsq
!     SET PARAMETERS IN PREDICTOR-CORRECTOR METHOD
      f02=dble(3./16.)
      f12=dble(251./360.)
      f32=dble(11./18.)
      f42=dble(1./6.)
      f52=dble(1./60.)
      dfr=ran2(iseed)           ! dfr is not used anywhere, just here (legacy)

C     ===============================================
!     LOOP OVER TEMPERATURES
      do 2000 it=1,nt
         if(lteql  .and. kteql.gt.0) then
            tr = teql
            write(1,'(/,a,f7.3,a,i7,a,f7.3)')
     +           '#Teql',tr,' for the first ',
     $           kteql/kunit,' Tau, then ',ttab(it)
         else
            tr = ttab(it)
            write(1,'(/,a,f7.3)')'#TEMPERATURE ',tr
         endif

!     LANGEVIN PARAMETERS
         gamma2=gamma/delta
         const2=2.d0*tr*gamma*delta ! assume xmas=1
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

         if(lthermo) then
            acv=0.d0
            acv2=0.d0
            apnat=0.d0
            apnat2=0.d0
            achi=0.d0
            achi2=0.d0
         endif

         nfm=mstep/kwforce
         if(lvelo.or.lwal) then
            do i=1,nfm          ! averages over multiple trajectories
               aufres(i)=0.d0
               aufres2(i)=0.d0
               aree(i)=0.d0
               adfres(i)=0.d0
               adfres2(i)=0.d0
            enddo
            afmax1=0
            afmax2=0
            atmax1=0
            atmax2=0
         endif

C     =========================================
!     LOOP OVER STARTING CONFIGURATIONS
         iterate=0
         do 1000 itraj=1,ntraj

            if(lteql) then
               tr = teql
               const2=2.d0*tr*gamma*delta ! assume xmas=1
               const2=dsqrt(const2)*delta
               aheat=delsq*part*3*tr
            endif

!     LANGEVIN PARAMETERS

            kbwal(1)=0          ! TIME COUNTER AFTER SQUEEZING
            kbwal(2)=0          ! T COUNTER BEFORE REACHING VELOCITY TO FIND MINIMUM
            kbwal(3)=0          ! T COUNTER AFTER REACHING MINIMUM
            kbwal(4)=0          ! T COUNTER BEFORE REACHING SQUEEZING VELOCITY
            kbwal(5)=0          ! T COUNTER AFTER REACHING MAXIMUM AMPLITUDE
            kbwal(6)=-2         ! T COUNTER OF OSCILLATIONS
            kbwal(7)=0          ! T COUNTER AFTER OSCILLATIONS
            kbwal(8)=0          ! T COUNTER BEFORE REACHING PULLING VELOCITY
            kbwal(9)=0          ! ALWAYS 0
            lcontin=.true.      !true: simulation proceeds; false: it stops
            vtarget=abs(vtarget)
            fresist = 0.d0
            afresist = 0.d0
            bfresist = 0.d0
            fresistperp = 0.d0
            afresistperp = 0.d0
            bfresistperp = 0.d0
            aufresist = 0.d0
            adfresist = 0.d0
            axufresist = 0.d0
            axdfresist = 0.d0
            ayufresist = 0.d0
            aydfresist = 0.d0
            bufresist=0.0
            bdfresist=0.0
            bxufresist = 0.d0
            bxdfresist = 0.d0
            byufresist = 0.d0
            bydfresist = 0.d0
            work=0.d0           !    WORK DONE IN ONE OSCILLATION CYCLE
            kbperiod=1          !    NUMBER OF OSCILLATION CYCLE
            intrsc=nssb         !    NUMBER OF INTRACHAIN DISULFIDE BONDS
            intesc=0            !    NUMBER OF INTERCHAIN DISULFIDE BONDS
            icnss=0             !    NUMBER OF NATIVE DISULFIDE BONDS
            icdss=0             !    NUMBER OF NON-NATIVE DISULFIDE BONDS
            cofdih=0.0          !    COEFFICIENT FOR SLOWLY TURNING ON DIHEDRAL POT.
            shear=0.0           !    THE SYSTEM IS NOT SHEARED AT THE BEGINNING
            menw=0              !    NUMBER OF FCC WALL BEADS IS 0 AT THE BEGINNING
            kfccw=0             !    LENGTH OF THE FCC WALL CONTACT LIST IS ALSO 0
            kqont=0             !    NUMBER OF NON-NATIVE CONTACTS IS ALS0 0
            jq=1                !    VERLET LIST HAS 2 COPIES, JQ=1 AND JQ=2
            if(lpullrel) then
               lvelo=.TRUE.
               lunfold=.TRUE.
               lmedian=.FALSE.
            endif

            if(lmedian.and.(.not.lconftm).and.(.not.lthermo)) then
               if(inot.gt.ntraj/2+1) goto 1001
            endif
            iterate=iterate+1

            do i=1,klont
               kbt(i)=0
               kut(i)=0
               imap(i)=0
            enddo

            do i=1,4            ! first row of E is preparation
               youngmod(i,1)=0.d0
            enddo               ! work is the work required to place polymers together

            if(lthermo) then
c     initialize averages over the trajectories
               ave=0.d0
               ave2=0.d0
               pnat=0.d0
               bchi=0.d0
               bchi2=0.d0
            endif

!     STARTING CONFORMATION
            if(lconftm.or.lmedian) then ! A STRAIGHT-LINE
               do i=1,men
                  x0(i)=0.d0
                  y0(i)=0.d0
                  z0(i)=(i-1)*bond
               enddo
               if(lsawconftm) call confstart(sdens,confcut)
            else if(lwarmup) then
               if(it.eq.1.and.itraj.eq.1) then
                  do i=1,men
                     x0(i)=xn(i)
                     y0(i)=yn(i)
                     z0(i)=zn(i)
                  enddo
               endif
            else if(lpdb.or.lstartpdb) then ! THE NATIVE STATE
               do i=1,men
                  x0(i)=xn(i)
                  y0(i)=yn(i)
                  z0(i)=zn(i)
               enddo
            else
               call confstart(sdens,confcut)
               do i=1,men
                  xn(i)=x0(i)
                  yn(i)=y0(i)
                  zn(i)=z0(i)
               enddo
            endif

            if(ldisp) call displace(iprota,jprotb,away)

            if(.not.lcpb) then
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
            endif
            if(lwal) then       ! CALCULATE WALL ENDS
               ip1=0
               ip2=0
               vpull=0.d0
               do ib=1,men
                  z0temp(ib)=z0(ib)
                  if(ljwal) z0temp(ib)=0
                  ksorted(ib)=ib
                  ipw(1,ib)=0
                  ipw(2,ib)=ib
               enddo
               if(ldens.and..not.lcpb) then
                  startvolume=(men)/(sdens) !*unit**3)
                  startboxsize=0.5*startvolume**(1.0/3.0)
                  xyzmin=-startboxsize
                  xyzmax=startboxsize
                  do ib=1,men   ! to form a cube, wals must have same distance
                     if(x0(ib)-2*bond.lt.xyzmin) xyzmin=x0(ib)-2*bond
                     if(y0(ib)-2*bond.lt.xyzmin) xyzmin=y0(ib)-2*bond
                     if(z0(ib)-2*bond.lt.xyzmin) xyzmin=z0(ib)-2*bond
                     if(x0(ib)+2*bond.gt.xyzmax) xyzmax=x0(ib)+2*bond
                     if(y0(ib)+2*bond.gt.xyzmax) xyzmax=y0(ib)+2*bond
                     if(z0(ib)+2*bond.gt.xyzmax) xyzmax=z0(ib)+2*bond
                  enddo
                  xdown=xyzmin
                  xup=xyzmax
                  ydown=xyzmin
                  yup=xyzmax
                  zdown=xyzmin
                  zup=xyzmax
                  xsep=xup-xdown
                  ysep=yup-ydown
                  zsep=zup-zdown
                  xinv=1.d0/xsep
                  yinv=1.d0/ysep
                  zinv=1.d0/zsep
               endif
            endif
            oldxup=xup
            oldxdown=xdown
            oldyup=yup
            oldydown=ydown
            oldzup=zup
            oldzdown=zdown

            if(lforce.or.lvelo) then ! COMPUTE THE DIRECTION OF FORCE
               afx=xn(ip2)-xn(ip1)
               afy=yn(ip2)-yn(ip1)
               afz=zn(ip2)-zn(ip1)
               aff=sqrt(afx*afx+afy*afy+afz*afz)
               afx=afx/aff
               afy=afy/aff
               afz=afz/aff
               vpulx=vpull*afx
               vpuly=vpull*afy
               vpulz=vpull*afz
            endif

            if(lvelo) then
               reemax=0.975*bond*men ! maximum length, after it trajectory ends
               xpul(1)=xn(ip2)
               ypul(1)=yn(ip2)
               zpul(1)=zn(ip2)
               if(lmaxforce) then
                  fmax1=0
                  fmax2=0
                  tresh1=70     ! threshold for measuring maximum force for 1tit hard
                  tresh2=160    ! numbers are in angstrems?
                  tresh1=tresh1/velo ! thresh1 and 2 are in taus
                  tresh2=tresh2/velo
               endif
            endif

!     LOAD INITIAL VELOCITIES OF PARTICLES
            call intvel3d(aheat,part,nen1)
            kb0=0
            sep0=0.0
            do i=1,men          ! ZERO THE CONNECTION TABLES
               l3rdcn(i)=.false. ! true if residue i forms a disulfide bond
               knct3rd(i)=0     ! index of residue that is bonded to i
               khbful(1,i)=0
               khbful(2,i)=0
               khbful(3,i)=0
               khbful(4,i)=0
               nei(1,ib)=0
               nei(2,ib)=0
            enddo
            if(.not.ldynss) then
               do ks=1,nssb
                  i1=ksb(1,ks)
                  i2=ksb(2,ks)
                  l3rdcn(i1)=.true.
                  l3rdcn(i2)=.true.
                  knct3rd(i1)=i2
                  knct3rd(i2)=i1
               enddo
            endif
            if(itraj.gt.1) lrst=.false.
            if(lrst) then
               open(21,file=rstfile,status='unknown')
               read(21,*)time,work
               kb0=int(time/delta)
               mstep=mstep+kb0
               read(21,*)zup,zdown,yup,ydown,xup,xdown,sep0,shear
               if(lampstrict) amp=ampstrict/sep0
               read(21,*)(kbwal(i),i=1,8)
               state='E '       ! EQUILIBRATION
               if(kbwal(1).gt.0) vtarget=-1.0*vtarget
               if(loscillate.and.kbwal(3).gt.0) vtarget=-1.0*vtarget
               if(kbwal(kconnecttime).gt.0) then
                  if(lfcc) call make_fcc()
                  if(kbwal(kconnecttime+1).eq.0) state='B ' ! AFTER BEADS
               endif
               if(loscillate.and.kbwal(6).gt.-2) then
                  state='O '
                  vtarget=-1.0*vtarget
               endif
               if(kbwal(8).gt.0) state='P '
               if(lfrmscrtch) then
                  kteql=kb0+ktrest ! EQUILIBRATION AFTER RESTART
c     call update_verlet_list(verlcut,nen1) ! only for ssbond
c     call compute_ssbonds()                ! now obsolete
                  do i=1,men    ! this is OK, positions are from the PDB file
                     x0(i)=x0(i)+xdown
                     y0(i)=y0(i)+ydown
                     z0(i)=z0(i)+zdown
                  enddo
               else
                  read(21,*)intrsc,intesc
                  do k=1,2*(intrsc+intesc)
                     read(21,*)i,j
                     l3rdcn(i)=.true.
                     l3rdcn(j)=.true.
                     knct3rd(i)=j
                     knct3rd(j)=i
                  enddo
                  read(21,*)(ksorted(i),i=1,men)
                  j=1
c     do k=1,ipwn/2
c     i=ksorted(k)
c     ipw(1,j)=0
c     ipw(2,j)=i
c     j=j+1
c     i=ksorted(men+1-k)
c     ipw(1,j)=1
c     ipw(2,j)=i
c     j=j+1
c     enddo
                  read(21,*)(xpul(i),ypul(i),zpul(i),i=1,men)
                  read(21,*)(x0(i),y0(i),z0(i),i=1,men)
                  read(21,*)(x1(i),y1(i),z1(i),i=1,men)
                  read(21,*)icnss,icdss,ip1,ip2
                  if(lwal) then
                     read(21,*)(ipw(1,i),i=1,men)
                     read(21,*)(ipw(2,i),i=1,men)
                  endif
                  read(21,*)kbperiod
               endif
               close(21)
            endif

            xsep=xup-xdown
            ysep=yup-ydown
            zsep=zup-zdown
            xinv=1.d0/xsep
            yinv=1.d0/ysep
            zinv=1.d0/zsep
            call update_verlet_list(verlcut,nen1) ! set up Verlet list
!     ASSIGN INITIAL ACCELERATION BASED ON INITIAL POSITIONS
!     if(lj1012) then
!     call evalgo_1012(epot)
!     else
            call prepare(epot)
            call evalgo(epot,chi)
            if(lwal.and..not.lnowal) call evalwall(epot)
            if(lpid) then
               call evalimproper(epot,lcospid,epsbb)
            else
               call evalcpot(epot)
            endif
!     endif
            if(lchiral) then
               call eval_chirality(enechi)
               epot=epot+enechi
            endif
            if(langle.or.lwritemap) call evalangles(epot,lsldh,0.d0)

!     SCALE ACCELERATIONS
            do 530 i=1,men
               x2(i)=fx(i)*deltsq
               y2(i)=fy(i)*deltsq
               z2(i)=fz(i)*deltsq
 530        continue

!     TABLE HEADING
            if(ksave.ne.0) write(2,'(a,i9)')'MODEL',iterate
            if(kwrite.ne.0) then
               write(1,'(//,a,i4)')'#TRAJECTORY',iterate
               if(lforce) then
                  write(1,'(a,a)')
     $                 '#    TIME   AFORCE      EPOT      ETOT',
     +'   ICN     RG    RMSD  D(1,N)   VEL'
               else if(lvelo) then
                  write(1,'(a,a)')'#    TIME      EPOT      ETOT',
     +                 '   ICN     RG    RMSD   D(1,N)   FORCE'
               else if(lwal) then
                  if(lwals) then
                     write(1,'(a,a,a,a)')
     $                    '#S     TIME        EPOT        ETOT',
     +' INTRHC INTEHC INTRSC INTESC    ICN     RMSD     FZ_UP',
     +'    FZDOWN       |F|       SEP    FPULLZ    FPULLX',
     +'         W NCORD     FX_UP    FXDOWN     FY_UP    FYDOWN'
c     else if(lwals) then
c     write(1,'(a,a,a)')'#S    TIME      EPOT      ETOT   ICN     RG',
c     +  '    RMSD   FZ_UP  FZDOWN     |F|      SEP',
c     +  '   FX_UP  FXDOWN   FY_UP  FYDOWN   FPULL'
                  else
                     write(1,'(a,a,a,a)')
     $                    '#S     TIME      EPOT        ETOT',
     +' INTRHC INTEHC INTRSC INTESC    ICN     RMSD     FZ_UP',
     +'    FZDOWN       |F|       SEP    FPULLZ    FPULLX',
     +'         W NCORD     XCM     YCM     ZCM   ICW'
                  endif
               else if(ldynss) then
                  write(1,'(a,a,a)')
     $                 '#    TIME          EPOT          ETOT',
     +'   ICN ICNss ICDss',
     +'      RG       L    RMSD NCORD     W CORDR KNOTS KNOTE'
               else
                  write(1,'(a,a,a)')
     $                 '#    TIME          EPOT          ETOT',
     +'   ICN B1-B2 S1-S2 B1-S2 B1-B1 S1-S1 B1-S1',
     +'      RG       L    RMSD NCORD     W CORDR KNOTS KNOTE'
               endif
            endif

            if(lpullrel) then
               write(1,'(a)')'#START PULLING WITH CONSTANT VELOCITY'
               call flush(1)
               do i=1,mpull
                  call lang(twopi,gamma2,const2,nen1)
                  call corr(deltsq,nen1)
                  call predct(nen1)
                  call prepare(epot)
                  call evalgo(epot,chi)
                  if(lwal.and..not.lnowal) call evalwall(epot)
                  if(lpid) then
                     call evalimproper(epot,lcospid,epsbb)
                  else
                     call evalcpot(epot)
                  endif
                  if(lchiral) then
                     call eval_chirality(enechi)
                     epot=epot+enechi
                  endif
                  if(langle.or.lwritemap) then
                     call evalangles(epot,lsldh,min(i*1.d0/ad,1.d0))
                  endif
                  call vafm(fresist,fresistperp,epot)
               enddo
               write(1,'(a,f12.2)')
     $              '#STOP PULLING AT D= ',mpull*delta*velo*unit
               lmedian=.TRUE.
               lunfold=.FALSE.
               lvelo=.FALSE.
               lforce=.FALSE.
               kb0=mpull
            endif

c     -----------------------------------------
c     ENTER MAIN LOOP OF SIMULATION
c     -----------------------------------------

            kb=kb0
 533        continue
            kb=kb+1

            if(lmass) then
               call lang_mass(twopi,gamma2,const2,nen1)
            else
               call lang(twopi,gamma2,const2,nen1)
            endif
            call corr(deltsq,nen1)
            call predct(nen1)
!     if(lj1012) then
!     call evalgo_1012(epot)
!     else
            call prepare(epot)
            call evalgo(epot,chi)
            if(lwal.and..not.lnowal) call evalwall(epot)
            if(lpid) then
               call evalimproper(epot,lcospid,epsbb)
            else
               call evalcpot(epot)
            endif
!     endif
            if(lchiral) then
               call eval_chirality(enechi)
               epot=epot+enechi
            endif
            if(langle.or.lwritemap) then
               if(lsldh.and.kb.le.ad) cofdih=kb*1.d0/ad
               call evalangles(epot,lsldh,cofdih)
            endif
            if(lmass) then
               do ib=1,men
                  rma=rmas(ib)
                  fx(ib)=fx(ib)/rma
                  fy(ib)=fy(ib)/rma
                  fz(ib)=fz(ib)/rma
               enddo
            endif
            if(kb.gt.kteql) then ! START OF A GLOBAL EQUILIBRATION CONDITION

               if(lwal) then
                  totforce=zuforce+zdforce+fresist ! fresist<0, zforce>0
                  btotforce=bufresist+bdfresist+bfresist !averaged total force
                  if(ldens) then
!---------PHASE 1: SQUEEZING PROTEINS TOGETHER----------
                     if(kbwal(1).eq.0) then
                        if(xsep*ysep*zsep.le.2.0*targetvolume) then
                           if((kb-kteql).lt.nratveld) then !REACHING VPULL
                              vpull=vtarget*float(kb-kteql)/nratveld
                              state='S ' ! SQUEEZING
                           else
                              vpull=vtarget !SQUEEZING AND CHECKING VOLUME
                           endif
                           if(xsep*ysep*zsep.le.targetvolume) then
                              call print_restart(kb,itraj)
                              kbwal(1)=1
                              vtarget=-1.0*vtarget
                              vpull=0.d0
                              if(kconnecttime.eq.1) then
                                 if(.not.ljwal) call connect_to_wal()
                                 if(lfcc) call make_fcc()
                                 state='B ' ! AFTER ATTACHING THE BEADS
                              else
                                 state='E ' ! EQUILIBRATION
                              endif
                           endif
                        else
                           if((kb-kteql).lt.nratveld) then !REACHING VPULLD
                              vpull=vpulld*float(kb-kteql)/nratveld
                              state='S ' ! SQUEEZING
                           else
                              vpull=vpulld ! SQUEEZING AND CHECKING VOLUME
                           endif
                        endif
!---------PHASE 2: RESTING AFTER SQUEEZING PROTEINS-----
                     else if(kbwal(1).lt.ktrest) then
                        kbwal(1)=kbwal(1)+1
!---------PHASE 3: EXTENDING THE BOX TO FIND MINIMUM----
                     else if(kbwal(3).eq.0) then
                        bceilf=sign(min(
     $                       abs(0.5*btotforce),1.d0),
     $                       btotforce)
                        if(kbwal(2).lt.nratvel) then
                           vpull=vtarget*bceilf*
     $                          float(kbwal(2))/nratvel
                           kbwal(2)=kbwal(2)+1
                           state='M ' ! SEARCHING FOR MINIMUM (F=0)
                        else    ! RELAXATION SPEED IS PROPORTIONAL TO FORCE
                           vpull=vtarget*bceilf
                        endif
                        if(abs(btotforce).lt.(0.05+0.00005*men).or.
     +                       .not.lminfbox) then
                           vpull=0.0
                           icheck=0
                           if(lobo) then
                              state='B ' ! AFTER ATTACHING THE BEADS
                              if(.not.ljwal)
     $                             call connect_to_wal_one_by_one()
                              if(ip1+ip2.ge.ipwn) icheck=1
                           else
                              icheck=1
                           endif
                           if(icheck.eq.1) then
                              call print_restart(kb,itraj)
                              kbwal(3)=1
                              if(lshear) then
                                 sep0=xsep
                              else
                                 sep0=zsep
                              endif
                              if(lampstrict) amp=ampstrict/sep0
                              vpull=0.d0
                              if(loscillate) vtarget=-1.0*vtarget
                              if(.not.lobo) then
                                 if(kconnecttime.eq.3) then
                                    if(.not.ljwal)
     $                                   call connect_to_wal()
                                    if(lfcc) call make_fcc()
                                    state='B ' ! AFTER ATTACHING THE BEADS
                                 else
                                    if(.not.lminfbox) kbwal(3)=ktrest
                                    state='E ' ! EQUILIBRATION
                                 endif
                              endif
                           endif
                        endif
!---------PHASE 4: RESTING AFTER FINDING MINIMUM--------
                     else if(kbwal(3).lt.ktrest) then
                        kbwal(3)=kbwal(3)+1
!---------PHASE 5: SQUEEZING BOX TO FULL AMPLITUDE------
                     else if(loscillate.and.kbwal(5).eq.0) then
                        if(((.not.lshear).and.zsep.gt.sep0*(1.0-amp))
     +                       .or.(lshear.and.shear.lt.sep0*amp*0.5))
     $                       then
                        if(kbwal(4).lt.nratvel) then
                           kbwal(4)=kbwal(4)+1
                           vpull=vtarget*float(kbwal(4))/nratvel
                           state='A ' ! REACHING MAXIMAL AMPLITUDE
                        else
                           vpull=vtarget
                        endif
                     else
                        call print_restart(kb,itraj)
                        kbwal(5)=1
                        vpull=0.0
                        if(kconnecttime.eq.5) then
                           if(.not.ljwal) call connect_to_wal()
                           if(lfcc) call make_fcc()
                           state='B ' ! AFTER ATTACHING THE BEADS
                        else
                           state='E ' ! EQUILIBRATION
                        endif
                     endif
!---------PHASE 6: RESTING AFTER FINDING MAX AMPLITUDE--
                  else if(loscillate.and.kbwal(5).lt.ktrest) then
                     kbwal(5)=kbwal(5)+1
!---------PHASE 7: OSCILLATIONS-------------------------
                  else if(loscillate.and.kbwal(6).eq.-2) then
                     kbwal(6)=-1
                     state='O ' ! OSCILLATIONS
                     vtarget=-1.0*vtarget
                  else if(loscillate.and.kbperiod.le.kbperiodmax) then
                     kbwal(6)=kbwal(6)+1
                     vpull=-1.d0*sep0*amp*omega*sin(omega*kbwal(6))
                     if(mod(kbwal(6),kwquarterperiod).eq.0) then
                        kremainder=mod(kbwal(6)/kwquarterperiod,4)
                        if(kremainder.eq.0) then
                           works(kbperiod)=work
                           k=kbperiod
                           write(1,'(a,4f13.9,f10.2)') '#',
     +                          (youngmod(j,k)/unit**3,j=1,4),works(k)
                           work=0.d0
                           kbperiod=kbperiod+1
                           call print_restart(kb,itraj)
                           kbwal(6)=kbwal(6)+1 ! for the last check
                        endif
                        emd=btotforce/(amp*xsep*ysep)
                        youngmod(kremainder+1,kbperiod)=emd
                     endif
!---------PHASE 8: RESTING AFTER OSCILLATIONS-----------
                  else if(kbwal(7).lt.ktrest) then
                     klastquarter=mod(kbwal(6),kwquarterperiod)
                     if(loscillate.and.klastquarter.ne.0) then
                        kbwal(6)=kbwal(6)+1 !last check return to sep0
                        vpull=-1.d0*sep0*amp*omega*sin(omega*kbwal(6))
                        if(klastquarter.eq.0) then
                           emd=btotforce/(amp*xsep*ysep)
                           write(1,'(a,4f13.9,f10.2)') '#',emd/unit**3,
     +                          0.0,0.0,0.0,work
                           work=0.d0
                        endif
                     else
                        if(kconnecttime.eq.8 .and. kbwal(7).eq.0) then
                           if(.not.ljwal) call connect_to_wal()
                           if(lfcc) call make_fcc()
                           state='B ' ! AFTER ATTACHING THE BEADS
                        else
                           state='E ' ! EQUILIBRATION
                        endif
                        kbwal(7)=kbwal(7)+1
                        vpull=0.d0
                     endif
!---------PHASE 9: PULLING------------------------------
                  else if(vpull.gt.vtarget.and.lpullfin) then
                     vpull=vtarget*float(kbwal(8))/nratvel
                     kbwal(8)=kbwal(8)+1
                     state='P ' ! PULLING
                     if(lshear) sep0=zsep
                     kbwal(4)=0
                  else if(lpullfin) then
                     vpull=vtarget
                  else
                     vpull=0.0
                  endif
               else             !TODO
                  if(vpull.lt.abs(vtarget)) then ! REACH VTARGET
                     if(vtarget.gt.0) then
                        state='S ' ! SQUEEZING
                        vpull=vtarget*float(kb-kteql)/nratvel
                     else if(kbwal(3).le.ktrest) then !EQUILIBRATION
                        state='R ' ! RESTING
                        kbwal(3)=kbwal(3)+1
                        vpull=0.d0
                     else
                        state='P ' ! PULLING
                        vpull=vtarget*float(kbwal(8))/nratvel
                        kbwal(8)=kbwal(8)+1
                     endif
                  else
                     vpull=vtarget
                  endif
               endif

               if(lshear.and.kbwal(4).gt.0) then
                  shear=shear+vpull*0.5
                  work=work+fresistperp*vpull
               else
                  zup=zup-vpull*0.5
                  zdown=zdown+vpull*0.5
                  work=work+totforce*vpull ! if box shrinks, vpull > 0
                  if(ldens.and.kbwal(1).eq.0) then
                     xup=xup-vpull*0.5
                     xdown=xdown+vpull*0.5
                     work=work+(yuforce+ydforce+xuforce+xdforce)*vpull
                     yup=yup-vpull*0.5
                     ydown=ydown+vpull*0.5 ! bad indent in the next line
                  else if(lconstvol.and.kbwal(8).eq.0
     $                    .and.kbwal(3).gt.0) then
                     vpullcv=-xsep*(1.d0-1.d0/sqrt(1.d0+vpull/zsep))
                     xup=xup-vpullcv*0.5
                     xdown=xdown+vpullcv*0.5
                     yup=yup-vpullcv*0.5
                     ydown=ydown+vpullcv*0.5
                     work=work+
     $                    (yuforce+ydforce+xuforce+xdforce)*vpullcv
                  endif
                  xsep=xup-xdown
                  ysep=yup-ydown
                  zsep=zup-zdown
                  xinv=1.d0/xsep
                  yinv=1.d0/ysep
                  zinv=1.d0/zsep
               endif
               lcontin=zsep.gt.sepmin
               if(kbwal(8).gt.0) then
                  lcontin=lcontin.and.zsep.lt.reemax
               endif
               if((.not.ldens) .AND. (vtarget.gt.0)) then
                  if((icw(1).ge.ipwn/2).and.(icw(2).ge.ipwn/2)) then
                     vtarget=-1.0*vtarget
                     if(.not.ljwal) call connect_to_wal()
                     if(lfcc) call make_fcc()
                  endif
               endif
            endif

            do kwal=1,menw/2
               z0(men+kwal)=zdown
               z0(men+menw/2+kwal)=zup
            enddo

            if(lforce) then     ! APPLY THE AFM FORCE
               aforce=coef*kb*delta
               call afm(aforce,epot)
            else if(lvelo .or. kbwal(kconnecttime).gt.0) then ! APPLY FORCE
               call vafm(fresist,fresistperp,epot)
               if(naver.ne.0) then
                  afresist=afresist+fresist
                  afresistperp=afresistperp+fresistperp
                  if(mod(kb,naver).eq.0) then
                     bfresist=afresist/naver
                     afresist=0.
                     bfresistperp=afresistperp/naver
                     afresistperp=0.
                  endif
               endif
            endif

         else                   ! Continuation of equlibration condition
            state='E '
            if(kb.eq.kteql .and. kteql.gt.0) then ! last step of equil.
               tr = ttab(it)    ! matters only if lteql is true
               const2=2.d0*tr*gamma*delta ! assume xmas=1
               const2=dsqrt(const2)*delta
               aheat=delsq*part*3*tr
               if(lforce.or.lvelo) then
                  afx=x0(ip2)-x0(ip1)
                  afy=y0(ip2)-y0(ip1)
                  afz=z0(ip2)-z0(ip1)
                  aff=sqrt(afx*afx+afy*afy+afz*afz)
                  afx=afx/aff
                  afy=afy/aff
                  afz=afz/aff
                  vpulx=vpull*afx
                  vpuly=vpull*afy
                  vpulz=vpull*afz
                  if(lvelo) then
                     xpul(1)=x0(ip2)
                     ypul(1)=y0(ip2)
                     zpul(1)=z0(ip2)
                  endif
                  do i=1,men
                     xn(i)=x0(i)
                     yn(i)=y0(i)
                     zn(i)=z0(i)
                  enddo
               endif
            endif
         endif                  ! end of equilibration condition

         if(lwal) then          ! MEASURE THE FORCE TO WALLS (A AND D)
            if(naver.ne.0) then
               aufresist=aufresist+zuforce
               adfresist=adfresist+zdforce
               if(mod(kb,naver).eq.0) then
                  bufresist=aufresist/naver
                  aufresist=0.
                  bdfresist=adfresist/naver
                  adfresist=0.
               endif
               if(lwals.or.
     $              (ldens.and.kbwal(1).eq.0.and..not.lsqpbc))then
               axufresist=axufresist+xuforce
               axdfresist=axdfresist+xdforce
               ayufresist=ayufresist+yuforce
               aydfresist=aydfresist+ydforce
               if(mod(kb,naver).eq.0) then
                  bxufresist=axufresist/naver
                  axufresist=0.
                  bxdfresist=axdfresist/naver
                  axdfresist=0.
                  byufresist=ayufresist/naver
                  ayufresist=0.
                  bydfresist=aydfresist/naver
                  aydfresist=0.
               endif
            endif
         endif
      endif

      if(lnatend.and. .not. lmedian) then
         lcontin=(icn.lt.klont)
         if(ldynss) lcontin=((icn+icnss).lt.(klont+nssb))
      endif
      if(kb.lt.mskip) goto 533  ! SKIPPING STEPS

      if(lthermo.and.(kb.eq.mskip)) then
c     initialize things for one trajectory after skipping
      endif

!     CALCULATE KINETIC ENERGY AND MEAN COORDINATION NUMBER
      sumvel=0.d0
      do 540 i=1,men
         sumvel=sumvel+(x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i))*rmas(i)
 540  continue
      ekin=sumvel/(2*delsq)     ! kinetic energy
      etot=epot+ekin            ! total energy
      if(lgln) then
         ancord=(1.0*ncord)/ngln+2.0 ! (mean coord. num.) for GLN
      else
         ancord=(1.0*ncord)/men+2.0 ! (mean coord. num.) for all AAs
      endif
!     FIND CONTACTS WHICH APPEAR FOR THE FIRST TIME DURING FOLDING
      if(lconftm) then
         do k=1,klont
            i=klist(1,k)
            j=klist(2,k)
            if(lunfold.and.imap(k).eq.1) kut(k)=kb ! CONTACT LAST APPERANCE
            if(kbt(k).eq.0) then
               if(lunfold) then
                  if(imap(k).eq.0) kbt(k)=kb ! CONTACT DISAPPEARS
               else
                  if(imap(k).eq.1) kbt(k)=kb ! CONTACT APPEARS
               endif
            endif
         enddo
      endif

!     UPDATE VERLET LIST IF CERTAIN DISTANCE IS CROSSED
      vrcut2sqlim = vrcut2sq
      if(lpbcx) vrcut2sqlim = vrcut2sqlim-max(xdown-oldxdown,0.0)
      if(lpbcx) vrcut2sqlim = vrcut2sqlim-max(oldxup-xup,0.0)
      if(lpbcy) vrcut2sqlim = vrcut2sqlim-max(ydown-oldydown,0.0)
      if(lpbcy) vrcut2sqlim = vrcut2sqlim-max(oldyup-yup,0.0)
      if(lpbcz) vrcut2sqlim = vrcut2sqlim-max(zdown-oldzdown,0.0)
      if(lpbcz) vrcut2sqlim = vrcut2sqlim-max(oldzup-zup,0.0)
      do 1525 i=1,men
         if((x0(i)-oxv(1,i))**2+(y0(i)-oxv(2,i))**2+(z0(i)-oxv(3,i))**2
     +        .gt.vrcut2sqlim) then
            call update_verlet_list(verlcut,nen1)
            oldxup=xup
            oldxdown=xdown
            oldyup=yup
            oldydown=ydown
            oldzup=zup
            oldzdown=zdown
            goto 1526
         endif
 1525 continue
 1526 continue

      if(lvelo.and.mod(kb,kwforce).eq.0) then
         ifn=kb/kwforce
         dx=x0(ip1)-x0(ip2)     ! END TO END DISTANCE
         dy=y0(ip1)-y0(ip2)
         dz=z0(ip1)-z0(ip2)
         ree=sqrt(dx*dx+dy*dy+dz*dz)
         aree(ifn)=aree(ifn)+ree
         if(naver.eq.0) then
            aufres(ifn)=aufres(ifn)+fresist
            aufres2(ifn)=aufres2(ifn)+fresist*fresist
         else
            aufres(ifn)=aufres(ifn)+bfresist
            aufres2(ifn)=aufres2(ifn)+bfresist*bfresist
         endif
      endif

      if(lwal.and.mod(kb,kwforce).eq.0) then
         ifn=kb/kwforce
!     END TO END DISTANCE
         aree(ifn)=aree(ifn)+zsep ! aree measures avg distance of ends
         if(naver.eq.0) then    ! in this context avg are not relevant?
            aufres(ifn)=aufres(ifn)+zuforce
            aufres2(ifn)=aufres2(ifn)+zuforce*zuforce
            adfres(ifn)=adfres(ifn)+zdforce
            adfres2(ifn)=adfres2(ifn)+zdforce*zdforce
         else                   ! the averages for aufres and adfres are over kwforce
            aufres(ifn)=aufres(ifn)+bufresist
            aufres2(ifn)=aufres2(ifn)+bufresist*bufresist
            adfres(ifn)=adfres(ifn)+bdfresist
            adfres2(ifn)=adfres2(ifn)+bdfresist*bdfresist
         endif                  ! force for output file (resist) is averaged over naver
      endif                     ! force for after-simulation statistis (res) over kwforce

!     PRINTING TO FILE EVERY KWRITE STEPS AND AT THE END
      if(kwrite.ne.0) then
         if(kb.eq.1.or.mod(kb,kwrite).eq.0 .or. .not. lcontin) then
            ktime=nint(kb*delta)
            call gyration(rg)   ! RADIUS OF GYRATION
            call cgyration()    ! Rg and other parameters like W
            call compute_rmsd(rms) ! RMSD

            if(.not.lwal) then
               dx=x0(ip1)-x0(ip2) ! END TO END DISTANCE
               dy=y0(ip1)-y0(ip2)
               dz=z0(ip1)-z0(ip2)
               ree=sqrt(dx*dx+dy*dy+dz*dz)
            endif

            icor=0
            corder=0.d0
            do k=1,kqont
               if(abs(kqist(4,k,jq)).gt.1) then
                  corder=corder + abs(kqist(1,k,jq)-kqist(2,k,jq))
                  icor=icor+1
               endif
            enddo
            if(icor.gt.0) corder=corder/icor
            corder=corder/men

            if(lvelo) then
               lcontin=ree.lt.reemax
c     projection of the restoring force on the direction of the force
               if(naver.eq.0) then
                  write(1,'(i9,2f11.2,i6,f7.2,f8.3,f8.2,f9.2)')
     +                 ktime,epot,etot,icn,rg*unit,
     $                 rms*unit,ree*unit,fresist/unit
               else
                  write(1,'(i9,2f11.2,i6,f7.2,f8.3,f8.2,f9.2)')
     +                 ktime,epot,etot,icn,rg*unit,
     $                 rms*unit,ree*unit,bfresist/unit
                  if(lmaxforce) then
                     if(ktime.le.nint(tresh1)
     $                    .and.bfresist.gt.fmax1) then
                     fmax1=bfresist
                     tmax1=dble(ktime)
                  endif
                  if(ktime.gt.nint(tresh1)
     $                 .and.ktime.le.nint(tresh2).and.
     +                 bfresist.gt.fmax2) then
                     fmax2=bfresist
                     tmax2=dble(ktime)
                  endif
               endif
            endif
            go to 9997
         endif

         if(lwal) then
            ubr=bufresist/unit
            dbr=bdfresist/unit
            fbr=-1.d0*bfresist/unit
            fbrp=-1.d0*bfresistperp/unit
            if(lshear.and.kbwal(4).gt.0) then
               zse=shear*2.0*unit ! shear is applied to both wals, so 2x
            else
               zse=zsep*unit    ! 0.5vpull is applied to both wals, so 1x
            endif
            if(lwals) then
               xubr=bxufresist/unit
               xdbr=bxdfresist/unit
               yubr=byufresist/unit
               ydbr=bydfresist/unit
               abr=sqrt((xubr-xdbr)**2+(yubr-ydbr)**2+(ubr-dbr)**2)
c     '#S     TIME       EPOT       ETOT INTRHC INTEHC INTRSC INTESC       RG
c     RMSD     FZ_UP    FZDOWN       |F|       SEP    FPULLZ    FPULLX
c     W NCORD     FX_UP    FXDOWN     FY_UP    FYDOWN'
               write(1,
     $'(a2,i9,2f11.2,5i7,f9.2,3f10.3,4f10.2,f6.2,4f10.2)')
     +state,ktime,epot,etot,intrhc,intehc,intrsc,intesc,icn,
     +rms*unit,ubr,dbr,abr,zse,fbr,fbrp,work,ancord,
     +xubr,xdbr,yubr,ydbr
            else
               abr=(ubr+dbr)/2
               write(1,
     $'(a2,i9,2f11.2,5i7,f9.2,3f10.3,4f10.2,f6.2,3f8.2,i6)')
     +state,ktime,epot,etot,intrhc,intehc,intrsc,intesc,icn,
     +rms*unit,ubr,dbr,abr,zse,fbr,fbrp,work,ancord,
     +xmcm*unit,ymcm*unit,zmcm*unit,icw(1)+icw(2)
            endif
            go to 9997
         endif

         if(lforce) then
            dx=x0(ip1)-x0(ip2)  ! END TO END DISTANCE
            dy=y0(ip1)-y0(ip2)
            dz=z0(ip1)-z0(ip2)
            ree=sqrt(dx*dx+dy*dy+dz*dz)
c     projection of the velocity on the direction of the force
            vl=x1(men)*afx + y1(men)*afy + z1(men)*afz
            vl=vl*unit/delta
            write(1,'(i9,f9.3,2f10.3,i6,f7.2,f8.3,f8.2,f8.2)')
     +           ktime,aforce,epot,etot,icn,
     $           rg*unit,rms*unit,ree*unit,vl
            go to 9997
         endif

         if(ldynss) then
            write(1,'(i9,2f14.3,3i6,2f8.2,f8.3,3f6.2,2i6)')
     $           ktime,epot,etot,
     +           icn,icnss,icdss,
     +           rg*unit,ree*unit,rms*unit,
     $           ancord,w(1),corder,knts(1,1),knts(2,1)
            go to 9997
         endif

         write(1,'(i9,2f14.3,7i6,2f8.2,f8.3,3f6.2,2i6)')
     $        ktime,epot,etot,
     +        icn,icnt(1),icnt(2),icnt(3),icnt(4),icnt(5),icnt(6),
     +        rg*unit,ree*unit,rms*unit,
     $        ancord,w(1),corder,knts(1,1),knts(2,1)

 9997    continue
         call flush(1)

      endif
      endif
      if(ksave.ne.0.and.((kb.ge.(kteql).and.(mod(kb,ksave).eq.0)).or.
     +     (kb.eq.1.or.(kb.lt.(kteql).and.mod(kb,kksave).eq.0)))) then
!     at first kteql time steps, sequence can be saved more often
         time=kb*delta
         if(mod(ksave,kwrite).ne.0) then
            call compute_rmsd(rms) ! RMSD
            call cgyration()    ! RG FOR INDIVIDUAL CHAINS
         endif
         if(lwritexyz) then
            call print_conf_xyz(2,time,epot,rms)
         else
            call print_conformation(2,time,epot,rms)
         endif
         if(lwritemap) call print_map(22,nint(kb*delta))
         call flush(2)
      endif

      if(krst.ne.0.and.mod(kb,krst).eq.0) then
         if(lcleanrst) then
            intrsc=0
            intesc=0
            icnss=0
            icdss=0
            do ib=1,men
               khbful(1,ib)=0
               khbful(2,ib)=0
               khbful(3,ib)=0
               khbful(4,ib)=0
               l3rdcn(ib)=.false.
               knct3rd(ib)=0
            enddo
            jq2=3-jq
            do k=1,kqont
               kqist(3,k,jq)=0
               kqist(4,k,jq)=sign(1,kqist(4,k,jq))
               kqist(3,k,jq2)=0
               kqist(4,k,jq2)=sign(1,kqist(4,k,jq2))
            enddo
         endif
         call print_restart(kb,itraj)
         if(ldelrst) then
            time=(kb-2*krst)*delta
            if(time.gt.0.d0) then
               write(stafile,*) '(a',klenstr,',i2.2,i10.10,a)'
               write(rstfile,stafile) filname,itraj,nint(time),'.rst'
               open(38,iostat=irstat,file=rstfile,status='unknown')
               if (irstat.eq.0) close(38, status='DELETE')
            endif
         endif
      endif

      if(lthermo.and.kb.gt.kteql) then
         ave=ave+etot
         ave2=ave2+etot*etot
         if(icn.eq.klont) pnat=pnat+1.0
         bchi=bchi+chi
         bchi2=bchi2+chi*chi
      endif

      if(abs(etot).gt.99999999.9) lcontin=.false. ! stop if it blew up
      if(lrmsmax.and.rms*unit.gt.rmsmax) lcontin=.false. ! stop if rms

      if(lmedian) then
         if(icn.lt.klont.and.kb.lt.mstep+mskip) goto 533
      else
         if(kb.lt.mstep+mskip .AND. lcontin) goto 533
      endif

c     -----------------------------------------
!     END LOOP OF SIMULATION
c     -----------------------------------------
 544  continue

      if(ksave.ne.0.and.mod(kb,ksave).ne.0) then
         time=kb*delta
         call compute_rmsd(rms) ! RMSD
         call cgyration()
         if(lwritexyz) then
            call print_conf_xyz(2,time,epot,rms)
         else
            call print_conformation(2,time,epot,rms)
         endif
         call flush(2)
      endif

!     ACCUMULATE CONTACT BREAKING TIMES
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

!     THE FOLDING TIME
      if(lmedian) then
         time=kb*delta
         tfold(iterate)=time
         if(icn.ne.klont) inot=inot+1
         if(kwrite.ne.0.and.kb.lt.mstep+mskip) then
            etot=epot+ekin      ! total energy
            call gyration(rg)   ! RADIUS OF GYRATION
            call compute_rmsd(rms) ! RMSD
            if(lforce) then
c     projection of the velocity on the direction of the force
               vl=x1(men)*afx + y1(men)*afy + z1(men)*afz
               vl=vl*unit/delta
               write(1,'(f9.1,f9.3,2f10.3,i6,f7.2,f8.3,f8.2,f8.2)')
     +              time,aforce,epot,etot,
     $              icn,rg*unit,rms*unit,ree*unit,vl
            else
               write(1,'(f9.1,2f10.3,i6,f7.2,f8.3)')
     +              time,epot,etot,icn,rg*unit,rms*unit
            endif
            call flush(1)
         endif
      endif

      if(loscillate) then
         write(1,'(//,a,x,5(9x,a))') '#','E1','E2','E3','E4','W'
         do k=1,kbperiod-1
            write(1,'(a,4f11.5,f10.2)')
     +           '# ',(youngmod(j,k)/unit**3,j=1,4),works(k)
         enddo
         write(1,'(a,f11.2)') '# ',work
      endif

      if(lthermo) then
         tav=kb-kteql-mskip
         ave=ave/tav
         ave2=ave2/tav
         pnat=pnat/tav
         cv=(ave2-ave*ave)/(tr*tr) ! SPECIFIC HEAT
         apnat=apnat+pnat
         apnat2=apnat2+pnat*pnat
         acv=acv+cv
         acv2=acv2+cv*cv
         bchi=bchi/tav
         bchi2=bchi2/tav
         dchi=bchi2-bchi*bchi   ! STRUCTURAL SUSCEPTIBILITY
         achi=achi+dchi
         achi2=achi2+dchi*dchi
      endif
!     nfm is maximal 'time' for which all trajectories are recorded
      if(lvelo.and.(kb/kwforce).lt.nfm) nfm = kb/kwforce
      if(lvelo.and.lmaxforce) then
         write(1,*)'#the first peak '
         write(1,*)'#max force ',fmax1/unit,
     $        '  displacement ',tmax1*velo
         write(1,*)'#the second peak '
         write(1,*)'#max force ',fmax2/unit,
     $        '  displacement ',tmax2*velo
         afmax1=afmax1+fmax1/unit
         afmax2=afmax2+fmax2/unit
         atmax1=atmax1+tmax1*velo
         atmax2=atmax2+tmax2*velo
      endif
c     print distances in the contacts

 1000 continue
!     END LOOP OVER CONFIGURATIONS

 1001 continue

      if(lconftm) then
         if(lunfold) then
            write(1,'(/,a)')
     $           '#AVERAGE TIME NEEDED FOR BREAKING EACH CONTACT'
            write(1,'(2a)')'# ICN    I    J  J-I       t0    DISP.',
     +           '       t1    DISP.'
         else
            write(1,'(/,a)')
     $           '#AVERAGE TIME NEEDED FOR FORMING EACH CONTACT'
            write(1,'(a)')'# ICN    I    J  J-I       t0    DISP.'
         endif
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
            if(lunfold) then
               cc=fut(k)
               dd=sqrt(abs(fut2(k)-cc*cc))/2
               write(1,'(4i5,4f11.2)')k,iseq(i),iseq(j),j-i,aa,bb,cc,dd
            else
c     write(1,'(4i5,2f11.2)')k,iseq(i),iseq(j),j-i,aa,bb
               caac(k)=aa       ! time of contact break
               kaak(k)=k        ! index of contact
            endif
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
               if(fcaa.gt.0.01)
     $              write(1,'(4i5,2f11.2)')ll,ise,jse,j-i,fcaa
            endif
         enddo
      endif

!     THE MEDIAN FOLDING TIME
      if(lmedian) then
         call sort(iterate,tfold)
         tmed=tfold(ntraj/2+1)
         write(1,'(/,a)')'#  TEMP     TMEDIAN  NTRAJ  INOT'
         write(1,'(a,f6.2,f12.2,i7,i6)')'#',tr,tmed,iterate,inot
      endif

      if(lthermo) then
         pnat=apnat/iterate
         pnat2=apnat2/iterate
         pnat2=sqrt(abs(pnat2-pnat*pnat))/2.d0
         cv=acv/iterate
         cv2=acv2/iterate
         cv2=sqrt(abs(cv2-cv*cv))/2.d0
         dchi=achi/iterate
         dchi2=achi2/iterate
         dchi2=sqrt(abs(dchi2-dchi*dchi))/2.d0
         write(1,'(/,a,a)')
     +        '#  TEMP      PNAT    DISP/2           CV       DISP/2',
     +        '         CHI      DISP/2'
         write(1,'(a,f6.2,2f10.5,2f13.3,2f12.4)')
     +        '#',tr,pnat,pnat2,cv,cv2,dchi,dchi2
      endif


      if(lvelo.and.iterate.gt.1) then
         write(1,'(//,a)')
     +        '#AVERAGED FORCE ON STRETCHING WITH CONSTANT VELOCITY'
         write(1,'(a)') '#    TIME  <D(1,N)>  <FORCE>  DISP./2'
         do i=1,nfm
            af1=aufres(i)/iterate
            af2=aufres2(i)/iterate-af1*af1
            af2=sqrt(abs(af2))*0.5d0
            ree=aree(i)/iterate*unit
            write(1,'(f9.1,f10.2,2f9.2)')i*kwforce*delta,ree,af1,af2
         enddo
         if(lmaxforce) then
            write(1,*)' '
            write(1,*)'#average over trajectories '
            at=1./ntraj
            write(1,*)'#the first peak '
            write(1,*)
     $           '#max force ',afmax1*at,'  displacement ',atmax1*at
            write(1,*)'#the second peak '
            write(1,*)
     $           '#max force ',afmax2*at,'  displacement ',atmax2*at
         endif
      endif

 2000 continue
!     END LOOP OVER TEMPERATURES
C     ===============================================

      close(1)
c     write(2,'(a)') 'ENDMDL'
      close(2)
      close(22)
c     ----------------------------
c     END OF EXECUTION
c     ----------------------------
      end