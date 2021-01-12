      subroutine main_loop_folding(nres,kcont,xn,yn,zn,seq,iseq,unit,
     $ ch,icmap,gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta,isi,
     $ tsim,tsav,tpdb,kupdate,ntrj,nch,ichain,kk_out,kk_pdb,
     $ eij,lchir,langl,l612,l612s,l1012,l61012,lpc,llf,lbd,xu,yu,zu,
     $ cutparam,s612,nr_l,rcut3)
     
      implicit none

      integer kk_out,kk_pdb
      integer nres,kcont,kront,nch,isi
      integer icmap(nres,nres),imap(nres,nres)
      integer iseq(nres)
      integer conttime(kcont)
      character*3 seq(nres)
      character*1 ch(nres)
      real*8 xn(nres),yn(nres),zn(nres)
      real*8 x0(nres),y0(nres),z0(nres)
      real*8 x1(nres),y1(nres),z1(nres)
      real*8 x2(nres),y2(nres),z2(nres)
      real*8 x3(nres),y3(nres),z3(nres)
      real*8 x4(nres),y4(nres),z4(nres)
      real*8 x5(nres),y5(nres),z5(nres)
      real*8 fx(nres),fy(nres),fz(nres)
      real*8 xu(nres,ntrj),yu(nres,ntrj),zu(nres,ntrj)
      real*8 b(nres),chirn(nres)
      real*8 phin(nres),thetan(nres)
      real*8 dij(kcont),sij(kcont)
      real*8 eij(nres,nres)
      real*8 gamma,delta,temp,H1,H2,echi,kphi1,kphi2,ktheta
      real*8 epot,ekin,etot,enechi,enephi,eneth,bond,s612
      integer klist(kcont,2)
      integer krist(nres*nres,2)
      integer ichain(nch,2)
      real*8 rcut,sigma0,unit,cutparam
      real*8 tsim,tsav,tpdb,tim
      real*8 tfold(ntrj)
      real*8 rg,w,rcut3
      integer kupdate,ntrj
      integer i,j,k,nr_l,inc,imodel
      integer*4 nsim,nsav,npdb
      logical lpc,llf,lbd,lchir,langl,l612,l612s,l1012,l61012,lunfold
C GET READY
      call gopotential612(nres,xn,yn,zn,icmap,b,bond,unit,
     $     kcont,klist,sij,dij,rcut,sigma0,nch,ichain,rcut3)
      call model_chirality(nres,xn,yn,zn,bond,chirn,nch,ichain)
      call get_nat_angles(nres,xn,yn,zn,nch,ichain,phin,thetan)
      nsim=tsim/delta
      nsav=tsav/delta
      npdb=tpdb/delta
      imodel=0      
C LOOP OVER TRAJECTORIES
      do k=1,ntrj
       do i=1,nres
        x0(i)=xu(i,k)
        y0(i)=yu(i,k)
        z0(i)=zu(i,k)
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
       enddo
       call center(nres,x0,y0,z0)
       call intvel3d(nres,x1,y1,z1,temp,delta,unit)
       lunfold=.true.
       tfold(k)=0.d0
       call update_vlist(nres,x0,y0,z0,icmap,krist,kront,isi,
     $     temp,delta,gamma,kupdate,rcut,unit)
       call write_pdb(kk_pdb,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
       write(kk_out,'(/,a12,i4)') '# TRAJECTORY',k
 1222  format(a8,8x,a4,8x,a4,3x,a3,6x,a2,7x,a1,/)
       write(kk_out,1222) '#   TIME','EPOT','ETOT','INC','RG','W'
C THE MAIN LOOP STARTS HERE
       do j=1,nsim
        do i=1,nres
         fx(i)=0.d0
         fy(i)=0.d0
         fz(i)=0.d0
        enddo
        tim=j*delta
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
        if (inc.eq.kcont.and.lunfold) then
         lunfold=.false.
         tfold(k)=tim
        endif
C SAVE RESULTS
        if (mod(j,nsav).eq.0) then
         epot=epot/nres
         etot=etot/nres
         call calc_w_rg(nres,x0,y0,z0,rg,w)
         write(kk_out,'(a1,1x,f8.1,2f12.4,i6,2f8.3)') 
     &        'Z',tim,epot,etot,inc,rg*unit,w
         if (.not.lunfold) goto 2000
        endif
        if (mod(j,npdb).eq.0) then
         imodel=imodel+1
         call write_pdb(kk_pdb,nres,x0,y0,z0,seq,iseq,ch,imodel,unit)
        endif
       enddo
 2000  continue
       if (lunfold) then
        write(kk_out,'(/,a)') '# NO FOLDING EVENT'
        tfold(k)=tsim
       else
        write(kk_out,'(/,a16,f16.4)') '# FOLDING TIME: ',tfold(k)
       endif
       write(kk_out,'(/,a,/)') 'CONTACT FIRST APPEARANCE'
       do j=1,kcont
        write(kk_out,'(a1,1x,i4,1x,i3,1x,i3,1x,i3,1x,f10.2)')
     $     'T',j,klist(j,1),klist(j,2),abs(klist(j,1)-klist(j,2)),
     $      conttime(j)*delta
       enddo        
      enddo
      call sort_set(ntrj,tfold)
      write(kk_out,'(/,a,/)') '# FOLDING TIME TABLE'
      do k=1,ntrj
       write(kk_out,'(a1,1x,i6,f12.4)') 'M',k,tfold(k)
      enddo      
      
      if (mod(ntrj,2).eq.1) then
       write(kk_out,'(/,a13,1x,f12.4,/)') 'MEDIAN TIME: ',
     $       tfold((ntrj-1)/2+1)
      endif
      
      if (mod(ntrj,2).eq.0) then
       write(kk_out,'(/,a13,1x,f12.4,/)') 'MEDIAN TIME: ',
     $       (tfold(ntrj/2)+tfold(ntrj/2+1))/2
      endif       
      
      return
      end
