      subroutine norma(a)
      implicit none
      double precision a(3,3)
      integer i,j
      double precision sum
      do i=1,3
       sum = 0.0
       do j=1,3
        sum = sum + a(j,i)*a(j,i)
       enddo
       do j=1,3
          a(j,i) = a(j,i)/sqrt(sum)
       enddo
      enddo
      return
      end