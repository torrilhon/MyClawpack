c
c
c =========================================================
      subroutine production(meqn,xi,yj,q,prod)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(meqn)
      dimension prod(meqn)
c     # source terms for 2d advaction equations
c
      tau = 0.1
c
      prod(1) = 5*xi*xi - 5*yj*yj
      prod(2) = -1.0/tau*q(2)
      prod(3) = -1.0/tau*q(3)

      return
      end
