      subroutine pvector(u,v,z,zmod)
      implicit none
      double precision u(3),v(3),z(3)
      double precision zmod
      z(1)=u(2)*v(3)-u(3)*v(2)
      z(2)=u(3)*v(1)-u(1)*v(3)
      z(3)=u(1)*v(2)-u(2)*v(1)
      zmod=sqrt(z(1)*z(1)+z(2)*z(2)+z(3)*z(3))
      return
      end