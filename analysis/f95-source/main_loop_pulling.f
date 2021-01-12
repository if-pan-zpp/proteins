      subroutine main_loop_pulling(nres,kcont,xn,yn,zn,seq,iseq,unit,
     $ ch,icmap,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta,isi,
     $ tsim,tsav,tpdb,kupdate,ntrj,nch,ichain,kk_out,kk_pdb,
     $ eij,lchir,langl,l612,l612s,l1012,l61012,lpc,llf,lbd,
     $ vpull,HH1,HH2,ip1,ip2,tfave,max_pull_f,stop_pulling,cutparam,
     $ s612,nr_l,knot_ask,kk_knot,tksav,knumbers,rr,arg,c_velo,c_for,
     $ forcon,c_param,a_param,stop_pulling2,tstop,ktight,rcut3)
     
      implicit none

      integer kk_out,kk_pdb,kk_knot,knumbers,rr
      integer nres,kcont,kront,nch,isi,ip1,ip2,klewy,kprawy
      integer kupdate,ntrj,inc,inn,imodel,nfave,ki,ktight,kmax
      integer*4 nsim,nsav,npdb,ksav,nstop
      integer icmap(nres,nres),imap(nres,nres)
      integer iseq(nres)
      integer conttime(kcont)
      character*3 seq(nres)
      character*1 ch(nres)
      character*1 knot_ask
      integer, dimension(:,:), allocatable :: kb_t
      integer, dimension(:,:), allocatable :: ke_t 
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
      real*8 eij(nres,nres)
      real*8 dij(kcont),sij(kcont)
      real*8 dtab(kcont)
      real*8 gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta
      real*8 epot,ekin,etot,enechi,enephi,eneth,bond,s612
      integer klist(kcont,2)
      integer krist(nres*nres,2)
      integer ichain(nch,2)
      real*8 rcut,sigma0,unit,max_pull_f,forcon,vl,vl1,vl2,rcut3
      real*8 tsim,tsav,tpdb,tim,cutparam,tksav,tstop
      real*8 rg,vpull,HH1,HH2,fresist,tfave,stop_cond
      real*8 afx,afy,afz,ree,dx,dy,dz,aforce,c_param,a_param
      real*8 vpullx,vpully,vpullz,xpull,ypull,zpull,afres,bfres
      character*80 arg
      integer i,j,k,nr_l
      logical lpc,llf,lbd,lchir,langl,l612,l612s,l1012,l61012
      logical stop_pulling
      logical c_velo,c_for,stop_pulling2
C GET READY

      if (knumbers.ne.0) then
      allocate(kb_t(knumbers,2))
      allocate(ke_t(knumbers,2))
      call get_knot_limits(nres,iseq,seq,ch,icmap,kcont,
     $ rr,arg,kb_t,ke_t,knumbers,kk_out)
      else
      allocate(kb_t(1,2))
      allocate(ke_t(1,2))
      kb_t(1,1)=1
      ke_t(1,1)=nres
      endif 
              
      call gopotential612(nres,xn,yn,zn,icmap,b,bond,unit,
     $     kcont,klist,sij,dij,rcut,sigma0,nch,ichain,rcut3)
        
      call model_chirality(nres,xn,yn,zn,bond,chirn,nch,ichain)
      call get_nat_angles(nres,xn,yn,zn,nch,ichain,phin,thetan)
      nsim=tsim/delta
      nstop=tstop/delta
      nsav=tsav/delta
      npdb=tpdb/delta
      ksav=tksav/delta
      nfave=tfave/delta
      imodel=0
      stop_cond=c_param*a_param*(nres-1)
      
C LOOP OVER TRAJECTORIES
      do k=1,ntrj
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
       do i=1,kcont
        conttime(i)=0
        dtab(i)=0
       enddo
       call intvel3d(nres,x1,y1,z1,temp,delta,unit)
C C-TERMINUS ATTACHED TO THE PULLING SPRING
       xpull=x0(ip2)
       ypull=y0(ip2)
       zpull=z0(ip2)
       afx=x0(ip2)-x0(ip1)
       afy=y0(ip2)-y0(ip1)
       afz=z0(ip2)-z0(ip1)
       ree=sqrt(afx*afx+afy*afy+afz*afz)
       afx=afx/ree
       afy=afy/ree
       afz=afz/ree
       vpullx=vpull*afx*delta
       vpully=vpull*afy*delta
       vpullz=vpull*afz*delta       
       afres=0.d0
       bfres=0.d0
       kmax=0
       
       call update_vlist(nres,x0,y0,z0,icmap,krist,kront,isi,
     $     temp,delta,gamma,kupdate,rcut,unit)
       call write_pdb(kk_pdb,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
       write(kk_out,'(/,a12,i4)') '# TRAJECTORY',k
       
       if (c_velo) then
 1222  format(2x,a8,8x,a4,8x,a4,3x,a3,6x,a2,2x,a6,3x,a5/)
       write(kk_out,1222) '#   TIME','EPOT','ETOT',
     $                    'INC','RG','D(1,N)','FORCE'
       endif
                      
       if (c_for) then 
 1223  format(a12,a9,a12,a12,a6,a8,1x,a8/)
       write(kk_out,1223) '#       TIME','FORCE','EPOT','ETOT',
     $                    'INC','RG','D(1,N)'
       endif
C THE MAIN LOOP STARTS HERE
       do j=1,nsim
        
        do i=1,nres
         fx(i)=0.d0
         fy(i)=0.d0
         fz(i)=0.d0
        enddo
        tim=j*delta
        aforce=forcon
        xpull=xpull+vpullx
        ypull=ypull+vpully
        zpull=zpull+vpullz        
C DYNAMICS 
        if (mod(j,kupdate).eq.0) then
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
     $        H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,j,
     $        cutparam,nr_l)
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
     $        H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,j,
     $        cutparam,nr_l)
        endif
        if (l61012) then
        call evalgo61012(nres,x0,y0,z0,fx,fy,fz,b,icmap,eij,
     $       kcont,klist,dij,krist,kront,epot,
     $       H1,H2,rcut,sigma0,nch,ichain,imap,inc,conttime,j,
     $       cutparam,nr_l)
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
 
        if (c_velo) then
        call vafm(nres,x0,y0,z0,xn,yn,zn,ip1,ip2,fx,fy,fz,fresist,
     $      HH1,HH2,afx,afy,afz,xpull,ypull,zpull)
        endif
        
        if (c_for) then
        call afm(nres,x0,y0,z0,xn,yn,zn,ip1,ip2,fx,fy,fz,fresist,
     $      HH1,HH2,afx,afy,afz,xpull,ypull,zpull,aforce)
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
        afres=afres+fresist
        if (mod(j,nfave).eq.0) then
         bfres=afres/nfave
         afres=0.d0
        endif

        dx=x0(ip2)-x0(ip1)
        dy=y0(ip2)-y0(ip1)
        dz=z0(ip2)-z0(ip1)
        ree=sqrt(dx*dx+dy*dy+dz*dz)
        do i=1,kcont
         if (conttime(i).eq.j) dtab(i)=ree
        enddo
        
        if ((knot_ask.eq.'Y').and.(mod(j,ksav).eq.0).and.
     &   (knumbers.ne.0)) then
         kmax=0
         do ki=1,knumbers
          call kmt(nres,x0,y0,z0,kb_t(ki,2),ke_t(ki,2),
     &     kb_t(ki,1),ke_t(ki,1))  
         write(kk_knot,'(a1,1x,i2,1x,f10.1,8x,i7,8x,i7)')'K',ki,tim,
     &    kb_t(ki,2),ke_t(ki,2)
         if (kmax.lt.(abs(kb_t(ki,2)-ke_t(ki,2)))) 
     &    kmax=abs(kb_t(ki,2)-ke_t(ki,2))
         enddo
        write(kk_knot,'(a43)')
     &   '-------------------------------------------'

         call flush(kk_knot)
         else if ((knot_ask.eq.'Y').and.(mod(j,ksav).eq.0).and.
     &   (knumbers.eq.0)) then
         do ki=1,1
          call kmt(nres,x0,y0,z0,kb_t(ki,2),ke_t(ki,2),
     &     kb_t(ki,1),ke_t(ki,1))  
         write(kk_knot,'(a1,1x,i2,1x,f10.1,8x,i7,8x,i7)')'K',ki,tim,
     &    kb_t(ki,2),ke_t(ki,2)
         kmax=abs(kb_t(ki,2)-ke_t(ki,2))
         enddo
        write(kk_knot,'(a43)')
     &   '-------------------------------------------'
         call flush(kk_knot)
         endif
        
        if (mod(j,nsav).eq.0) then
         epot=epot/nres
         etot=etot/nres
         call gyration(nres,x0,y0,z0,rg)
         call get_nc(nres,imap,kcont,klist,ch,inn)

         if (c_velo) then
          write(kk_out,'(a1,1x,f10.1,2f12.4,2i6,f8.3,1x,f8.3,1x,f8.3)') 
     $        'F',tim,epot,etot,inc,inn,rg*unit,ree*unit,bfres/unit
         endif
         
         if (c_for) then
c         vl=x1(ip1)*afx+y1(ip1)*afy+z1(ip1)*afz
c         vl1=vl*unit/delta
c         vl=x1(ip2)*afx+y1(ip2)*afy+z1(ip2)*afz
c         vl2=vl*unit/delta
          write(kk_out,
     $     '(a1,1x,f10.1,f9.3,2f12.4,2i6,f8.3,1x,f8.3)') 
     $     'F',tim,aforce/unit,epot,etot,inc,inn,rg*unit,ree*unit
         endif
                
         call flush(kk_out)
            
            if ((inc.eq.0 .or. inn.eq.0).and.(stop_pulling)) then 
               write(kk_out,'(a)')'ICN=0 -- SIMULATION STOPPED'
               goto 2000
            endif
            if ((ree*unit.gt.stop_cond).and.(stop_pulling2)) 
     $       then 
               write(kk_out,'(a7,f4.2,a1,f3.1,a28)')'D(1,N)>',c_param,
     $         '*',a_param,'*(N-1) -- SIMULATION STOPPED'
               goto 2000
            endif
            if ((tstop.gt.1d-1).and.(nstop.le.j).and.(ktight.ge.kmax)) 
     $       then 
               write(kk_out,'(a5,f8.1,a53)')'TIME ',tstop,
     $         ' ELAPSED! THE KNOT IS TIGHTENED -- SIMULATION STOPPED'
               goto 2000
            endif
           if ((bfres/unit.ge.max_pull_f)) then 
               write(kk_out,'(a,f12.6)')'FORCE GREATER THAN ',
     $         max_pull_f
               write(kk_out,'(a)') 'SIMULATION STOPPED'          
               goto 2000
            endif 
        endif
        if (mod(j,npdb).eq.0) then
         imodel=imodel+1
         call write_pdb(kk_pdb,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
        endif
        call flush(0)
       enddo
 2000  continue
       write(kk_out,'(/,a,/)') 'CONTACT LAST APPEARANCE'
       do j=1,kcont
       write(kk_out,'(a1,1x,i4,1x,i3,1x,i3,1x,i3,1x,2f10.2)')
     $      'T',j,klist(j,1),klist(j,2),abs(klist(j,1)-klist(j,2)),
     $      conttime(j)*delta,dtab(j)*unit
       enddo 
       imodel=imodel+1
       call write_pdb(kk_pdb,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
      enddo      
      return
      end
