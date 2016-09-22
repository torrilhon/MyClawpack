c
c
c =========================================================
      subroutine production(meqn,xi,yj,q,prod)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(meqn)
      dimension prod(meqn)

      common /cparam/ tau
c
      prod(1) = 0.d0
      prod(2) = 0.d0
      prod(3) = 0.d0
      prod(4) = -1.0/tau*q(4)
      prod(5) = -1.0/tau*q(5)
      prod(6) = -1.0/tau*q(6)
      prod(7) = -1.0/tau*q(7)

      return
      end
