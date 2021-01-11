      subroutine sort_set(n,q)
      implicit none
      integer n
      real*8 q(n)
      real*8 qi
      integer i,j
      
      do i=1,n
       do j=1,n
        if (q(j).gt.q(i)) then
         qi=q(i)
         q(i)=q(j)
         q(j)=qi
        endif
       enddo
      enddo
      return
      end 
