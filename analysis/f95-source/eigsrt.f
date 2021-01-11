      subroutine eigsrt(d,v,n,np)
      integer n,np
      real*8 d(np),v(np,np)
      integer i,j,k
      real*8 p
      do 13 i=1,n-1
       k=i
       p=d(i)
       do 11 j=i+1,n
       if(d(j).ge.p)then
          k=j
          p=d(j)
       endif
11     continue
       if(k.ne.i)then
        d(k)=d(i)
        d(i)=p
        do 12 j=1,n
          p=v(j,i)
          v(j,i)=v(j,k)
          v(j,k)=p
12      continue
       endif
13    continue
      return
      end
      
C  (C) Copr. 1986-92 Numerical Recipes Software W"..