      subroutine main_loop_th(nres,kcont,xn,yn,zn,seq,iseq,unit,
     $ ch,icmap,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta,isi,
     $ tsim,tsav,tpdb,teql,kupdate,ntrj,nch,ichain,kk_out,kk_pdb,
     $ eij,lchir,langl,l612,l612s,l1012,l61012,lpc,llf,lbd,
     $ tstart,tend,tstep,cutparam,s612,nr_l,rcut3)
     
      implicit none

      integer kk_out,kk_pdb
      integer nres,kcont,kront,nch,isi
      integer icmap(nres,nres),imap(nres,nres)
      integer iseq(nres)
      character*3 seq(nres)
      character*1 ch(nres)
      integer conttime(kcont)
      real*8 xn(nres),yn(nres),zn(nres)
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 x1(nres),y1(nres),z1(nres)
      real*8 x2(nres),y2(nres),z2(nres)
      real*8 x3(nres),y3(nres),z3(nres)
      real*8 x4(nres),y4(nres),z4(nres)
      real*8 x5(nres),y5(nres),z5(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 b(nres),chirn(nres)
      real*8 phin(nres),thetan(nres)
      real*8 dij(kcont),sij(kcont)
      real*8 eij(nres,nres)
      real*8 gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta
      real*8 epot,ekin,etot,enechi,enephi,eneth,bond,s612
      real*8 tstart,tend,tstep,rcut3
      integer klist(kcont,2)
      integer krist(nres*nres,2)
      integer ichain(nch,2)
      real*8 rcut,sigma0,unit,cutparam
      real*8 tsim,tsav,tpdb,teql,tim
      real*8 rms,rg,pnat,avene,avene2,cv,avq,w
      real*8 kw(nres,3),kw2(nres,3)
      real*8 rmsf(nres)
      real*8 drij(kcont),dr2ij(kcont),dfij(kcont)
      integer kupdate,ntrj
      integer*4 nsim,nsav,npdb,neql
      integer i,j,k,it1,it2,nr_l,inc,imodel
      logical lpc,llf,lbd,lchir,langl,l612,l612s,l1012,l61012
C GET READY
      call gopotential612(nres,xn,yn,zn,icmap,b,bond,unit,
     $     kcont,klist,sij,dij,rcut,sigma0,nch,ichain,rcut3)
      call model_chirality(nres,xn,yn,zn,bond,chirn,nch,ichain)
      call get_nat_angles(nres,xn,yn,zn,nch,ichain,phin,thetan)
      do i=1,nres
        x0(i)=xn(i)
        y0(i)=yn(i)
        z0(i)=zn(i)
        x1(i)=0.d0
        y1(i)=0.d0
        z1(i)=0.d0
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
      enddo
      call center(nres,x0,y0,z0)
      if (tstep.gt.0.d0) temp=tstart
      call intvel3d(nres,x1,y1,z1,temp,delta,unit)
      call update_vlist(nres,x0,y0,z0,icmap,krist,kront,isi,
     $     temp,delta,gamma,kupdate,rcut,unit)
      imodel=0
      call write_pdb(kk_pdb,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
      nsim=tsim/delta
      nsav=tsav/delta
      npdb=tpdb/delta
      neql=teql/delta
 200  continue
C LOOP OVER TRAJECTORIES
      do k=1,ntrj
       pnat=0.d0
       avene=0.d0
       avene2=0.d0
       avq=0.d0
       do i=1,nres
        kw(i,1)=0.d0
        kw(i,2)=0.d0
        kw(i,3)=0.d0
        kw2(i,1)=0.d0
        kw2(i,2)=0.d0
        kw2(i,3)=0.d0
       enddo
       do i=1,kcont
        drij(i)=0.d0
        dr2ij(i)=0.d0
        conttime(i)=0
       enddo
       it1=0
       it2=0
       write(kk_out,'(/,a12,i4)') '# TRAJECTORY',k
 1222  format(2x,a8,8x,a4,8x,a4,3x,a3,6x,a2,4x,a4,7x,a1,/)
       write(kk_out,1222) '#   TIME','EPOT','ETOT','INC','RG','RMSD','W'
C THE MAIN LOOP STARTS HERE
       do j=1,nsim
        do i=1,nres
         fx(i)=0.d0
         fy(i)=0.d0
         fz(i)=0.d0
        enddo
        tim=j*delta
C DYNAMICS 
        if (mod(j,kupdate).eq.0 ) then
          call update_vlist(nres,x0,y0,z0,icmap,krist,kront,isi,
     $         temp,delta,gamma,kupdate,rcut,unit)
        endif
        if (lpc) then
         call predct(nres,x0,y0,z0,x1,y1,z1,x2,y2,z2,
     $        x3,y3,z3,x4,y4,z4,x5,y5,z5)
        endif
        if (l612) then
         call evalgo612(nres,x0,y0,z0,fx,fy,fz,b,icmap,eij,
     $        kcont,klist,sij,dij,krist,kront,epot,
     $        H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,
     $        j,cutparam,nr_l)
        endif
        if (l612s) then
         call evalgo612shift(nres,x0,y0,z0,fx,fy,fz,b,icmap,eij,
     $        kcont,klist,sij,dij,krist,kront,epot,s612,
     $        H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,
     $        j,cutparam,nr_l)
        endif
        if (l1012) then
         call evalgo1012(nres,x0,y0,z0,fx,fy,fz,b,icmap,eij,
     $        kcont,klist,dij,krist,kront,epot,
     $        H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,
     $        j,cutparam,nr_l)
        endif
        if (l61012) then
        call evalgo61012(nres,x0,y0,z0,fx,fy,fz,b,icmap,eij,
     $       kcont,klist,dij,krist,kront,epot,
     $       H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,
     $       j,cutparam,sij,nr_l)
        endif
        if (lchir) then
         call eval_chirality(nres,x0,y0,z0,fx,fy,fz,
     $        bond,chirn,echi,enechi,nch,ichain)
         epot=epot+enechi
        endif
        if (langl) then
         call eval_angles(nres,x0,y0,z0,fx,fy,fz,nch,ichain,
     $        kphi1,kphi2,phin,enephi,ktheta,thetan,eneth)
         epot=epot+enephi+eneth
        endif
        if (lpc) then
         call langvin(nres,x1,y1,z1,fx,fy,fz,gamma,delta,temp)
         call corr(nres,x0,y0,z0,x1,y1,z1,x2,y2,z2,
     $        x3,y3,z3,x4,y4,z4,x5,y5,z5,fx,fy,fz,delta)
        endif
        if (llf) then
         call leapfrog_dyn(nres,x0,y0,z0,x1,y1,z1,fx,fy,fz,
     $        gamma,delta,temp)
        endif
        if (lbd) then
         call brown_dyn(nres,x0,y0,z0,fx,fy,fz,gamma,delta,temp)
        endif
        ekin=0.d0
        if (lpc .or. llf) then
         do i=1,nres
          ekin=ekin+x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i)
         enddo
         ekin=ekin/(2.d0*delta*delta)
        endif
        etot=epot+ekin
C SAVE RESULTS
        if (j.gt.neql) then
         it1=it1+1
         epot=epot/nres
         etot=etot/nres
         avene=avene+etot
         avene2=avene2+etot*etot
         if (inc.eq.kcont) pnat=pnat+1.d0
         avq=avq+real(inc)/kcont         
         if (mod(j,nsav).eq.0) then
          it2=it2+1
          call compute_rmsd(nres,x0,y0,z0,xn,yn,zn,rms,kw,kw2)
          call calc_w_rg(nres,x0,y0,z0,rg,w)
          write(kk_out,'(a1,1x,f8.1,2f12.4,i6,3f8.3)') 
     $          'H',tim,epot,etot,inc,rg*unit,rms*unit,w
          call contact_fluct(nres,kcont,klist,dij,x0,y0,z0,drij,dr2ij)
         endif
         if (mod(j,npdb).eq.0) then
          imodel=imodel+1
          call write_pdb(kk_pdb,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
         endif
        endif
       enddo
C WRITE OUT RESULTS
       pnat=pnat/it1
       avq=avq/it1
       avene=avene/it1
       avene2=avene2/it1
       cv=(avene2-avene**2)/(temp*temp)
       do i=1,nres
        kw(i,1)=kw(i,1)/it2
        kw(i,2)=kw(i,2)/it2
        kw(i,3)=kw(i,3)/it2
        kw2(i,1)=kw2(i,1)/it2
        kw2(i,2)=kw2(i,2)/it2
        kw2(i,3)=kw2(i,3)/it2
        rmsf(i)=sqrt(kw2(i,1)-kw(i,1)*kw(i,1)
     $              +kw2(i,2)-kw(i,2)*kw(i,2)
     $              +kw2(i,3)-kw(i,3)*kw(i,3))
       enddo
       do i=1,kcont
        drij(i)=drij(i)/it2
        dr2ij(i)=dr2ij(i)/it2
        dfij(i)=sqrt(dr2ij(i)-drij(i)*drij(i))
       enddo
 1223  format(/,a6,9x,a3,8x,a4,11x,a1,10x,a2,/)
       write(kk_out,1223) '# TEMP','<E>','PNAT','Q','CV'
       write(kk_out,'(a1,1x,f6.3,4f12.6)')'A',temp,avene,pnat,avq,cv
       write(kk_out,'(/,a6,12x,a4,/)') '#    I','RMSF'
       do i=1,nres
        write(kk_out,'(a1,1x,i6,f16.8)') 'R',i, rmsf(i)*unit
       enddo
       write(kk_out,'(/,a)') '# FLUCTUATIONS IN THE LENGTH OF CONTACTS'
 1224  format(a6,4x,a2,4x,a2,14x,a2,10x,a6,4x,a16,/)
       write(kk_out,1224) '#    I','I1','I2','d0','<d-d0>',
     $   'SQRT(<(d-d0)^2>)'
       do i=1,kcont
        write(kk_out,'(a1,1x,3i6,3f16.8)')'F',i,klist(i,1),klist(i,2),
     $   dij(i)*unit,drij(i)*unit,dfij(i)*unit
       enddo
      enddo
      if (tstep.gt.0.d0) then
       temp=temp+tstep
       if (temp.le.tend) goto 200
      endif
      return
      end
