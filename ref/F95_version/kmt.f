      subroutine kmt(nres,x0,y0,z0,klewy,kprawy,kbegin,kend)
      implicit none
      integer nres,i,ii,j,n,kbegin,kend,nres1
      integer klewy,kprawy,ifun,ifun3,knot,knot1
      integer kbinl,kbinr,kendl,kendr,koba,modify
      real*8 x0(nres),y0(nres),z0(nres)
      integer, dimension(:), allocatable :: na
      integer, dimension(:), allocatable :: nb
      real*8, dimension(:), allocatable :: xk
      real*8, dimension(:), allocatable :: yk
      real*8, dimension(:), allocatable :: zk
      character aseq*3,pdbfile*20,outfile*20
      double precision det0,det1,det2,det3,xji,yji,zji
      double precision ro, xro,yro,zro,f12xy,f23xy,f31xy
    
      knot=0
      knot1=0
      kbinl=0
      kbinr=0
      kendl=0
      kendr=0
      koba=0
      nres1=nres

      !kbegin=50
      !kend=120
      nres=kend-kbegin+1

      allocate(xk(kend-kbegin+1))
      allocate(yk(kend-kbegin+1))
      allocate(zk(kend-kbegin+1))
      allocate(na(kend-kbegin+1))
      allocate(nb(kend-kbegin+1))

      do i=1,nres!kbegin,kend
         xk(i)=x0(i+kbegin-1)
         yk(i)=y0(i+kbegin-1)
         zk(i)=z0(i+kbegin-1)
         na(i)=i+kbegin-1
         nb(i)=i+kbegin-1
      enddo

      ii=0
      j=0
      i=2
      n=nres
      modify=0

566   continue   
      ifun3=0
         det0=0d0
         det1=0d0
         det2=0d0
         ifun=0
         ro=0

!#########################################################################

         do j=1,n-1
       if (((j.lt.(i-2)).or.(j.gt.(i+1))).and.(ifun.ne.1)) then
    
          ro=0.d0
          xro=0.d0
          yro=0.d0
          zro=0.d0
        
        det1=(xk(j)-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(zk(j)-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(yk(j)-yk(i-1))*(zk(i)-zk(i-1))-
     & (zk(j)-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(xk(j)-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(yk(j)-yk(i-1))*(xk(i)-xk(i-1))
        
        det2=(xk(j+1)-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(zk(j+1)-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(yk(j+1)-yk(i-1))*(zk(i)-zk(i-1))-
     & (zk(j+1)-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(xk(j+1)-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(yk(j+1)-yk(i-1))*(xk(i)-xk(i-1))

        if (det1*det2.gt.0) ifun=0

        if (det1*det2.le.0) then 
         ! przeciecie sie plaszczyzny trojkata i-1,i,i+1
         ! z punktem wyzanaczonym przez prosta j,j+1
         ! dane wzorem na ro, F0_(i,1,i,i+1)(0)

       xji=xk(j+1)-xk(j)     ! wyznaczenie wektora od j do j+1
       yji=yk(j+1)-yk(j)
       zji=zk(j+1)-zk(j)
       
       det3=(xji-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(zji-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(yji-yk(i-1))*(zk(i)-zk(i-1))-
     & (zji-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(xji-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(yji-yk(i-1))*(xk(i)-xk(i-1))
       
       det0=(-xk(i-1))*(yk(i)-yk(i-1))*(zk(i+1)-zk(i-1))+
     & (xk(i)-xk(i-1))*(yk(i+1)-yk(i-1))*(-zk(i-1))+
     & (xk(i+1)-xk(i-1))*(-yk(i-1))*(zk(i)-zk(i-1))-
     & (-zk(i-1))*(yk(i)-yk(i-1))*(xk(i+1)-xk(i-1))-
     & (zk(i)-zk(i-1))*(yk(i+1)-yk(i-1))*(-xk(i-1))-
     & (zk(i+1)-zk(i-1))*(-yk(i-1))*(xk(i)-xk(i-1))
         
       ro=det1/(det3-det0)
       
       xro=xk(j)-ro*xji 
       yro=yk(j)-ro*yji
       zro=zk(j)-ro*zji
       
       f12xy=(xro-xk(i-1))*(yk(i)-yk(i-1))-(xk(i)-xk(i-1))*(yro-yk(i-1))
       f23xy=(xro-xk(i))*(yk(i+1)-yk(i))-(xk(i+1)-xk(i))*(yro-yk(i))
       f31xy=(xro-xk(i+1))*(yk(i-1)-yk(i+1))-(xk(i-1)-xk(i+1))*(yro-
     & yk(i+1))

        if ((f12xy*f23xy.le.0).or.(f23xy*f31xy.le.0)) then 
         ifun=0

        elseif ((f12xy*f23xy.gt.0).and.(f23xy*f31xy.gt.0)) then
          ifun=1

        endif
      endif 
      endif          
         enddo 


!########################################################################

       if (ifun.eq.0) then 

       do ii=i,n-1
        xk(ii)=xk(ii+1)
        yk(ii)=yk(ii+1)
        zk(ii)=zk(ii+1)
        na(ii)=na(ii+1)
       enddo
          
        xk(n)=0d0
        yk(n)=0d0
        zk(n)=0d0
        modify=modify+1
        n=n-1        
       elseif (ifun.eq.1) then
        i=i+1
       endif
!#######################################################################

      if (i.lt.n) then
       goto 566
      elseif ((i.eq.n).and.(modify.ge.1)) then
       i=2
       modify=0
       if(n.ge.3) goto 566
      endif

!#########################################################################

      if (i.eq.2) then
!-------------------------------------------------------------------------         
         if (knot.eq.0) knot1=1
         knot=1
!-------------------------------------------------------------------------       
         if ((kbinl.eq.0).and.(kbinr.eq.0)) then
            kbinl=0
            kbinr=0
         elseif (kendl.eq.0) then
            kendl=1     
            n=nres-kbinr      
            do i=1,n-1
             xk(i)=x0(i+kbegin-1)
             yk(i)=y0(i+kbegin-1)
             zk(i)=z0(i+kbegin-1)
            enddo
            n=n-1  
            kbinr=kbinr+1  
            i=2
            modify=0
            goto 566
         elseif (kendr.eq.0) then
            kendr=1      
            n=nres-koba  
            n=n-2     
            koba=koba+1 
            do i=1,n
             xk(i)=x0(i+kbegin-1+1)
             yk(i)=y0(i+kbegin-1+1)
             zk(i)=z0(i+kbegin-1+1)
            enddo
            i=2
            modify=0
            goto 566
         elseif (koba.gt.0) then
            koba=koba
         endif
       
!-------------------------------------------------------------------------       

      elseif (i.ne.2) then
!-------------------------------------------------------------------------             
         if (knot.eq.1) knot1=1
         knot=0
!-------------------------------------------------------------------------                
         if ((kbinl.lt.200).and.(kendl.eq.0)) then
            n=nres-kbinl
            n=n-1
            kbinl=kbinl+1
            do i=1,n
             xk(i)=x0(i+kbegin-1+kbinl)
             yk(i)=y0(i+kbegin-1+kbinl)
             zk(i)=z0(i+kbegin-1+kbinl)
            enddo
            i=2
            modify=0
            goto 566
         elseif ((kbinl.ge.200).and.(kendl.eq.0)) then
            kendl=1  
            n=nres-kbinr
            n=n-1
            kbinr=kbinr+1
            do i=1,n
             xk(i)=x0(i+kbegin-1)
             yk(i)=y0(i+kbegin-1)
             zk(i)=z0(i+kbegin-1)
            enddo
            i=2
            modify=0
            goto 566
         elseif ((kbinr.lt.200).and.(kendr.eq.0)) then
            n=nres-kbinr
            n=n-1        
            kbinr=kbinr+1  
            do i=1,n
             xk(i)=x0(i+kbegin-1)
             yk(i)=y0(i+kbegin-1)
             zk(i)=z0(i+kbegin-1)
            enddo
            i=2
            modify=0
            goto 566
         elseif ((kbinr.ge.200).and.(kendr.eq.0)) then
            kendr=1  
            n=nres-koba
            koba=koba+1
            n=n-2
            do i=1,n
             xk(i)=x0(i+kbegin-1+1)
             yk(i)=y0(i+kbegin-1+1)
             zk(i)=z0(i+kbegin-1+1)
            enddo
            i=2
            modify=0
            goto 566
         elseif (koba.le.1) then
            n=nres-(koba*2)
            koba=koba+1
            n=n-2
            do i=1,n-1
             xk(i)=x0(i+kbegin-1+koba+1)
             yk(i)=y0(i+kbegin-1+koba+1)
             zk(i)=z0(i+kbegin-1+koba+1)
            enddo
            i=2
            modify=0
            goto 566
         endif
!-------------------------------------------------------------------------       
      endif
!############################################################################      
      
       kprawy=nres-kbinr+kbegin-1
       klewy=kbinl+kbegin-1
       nres=nres1
       return
       end
