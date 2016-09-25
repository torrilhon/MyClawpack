c
c
c =========================================================
      subroutine production(meqn,xi,yj,q,prod)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(meqn)
      dimension prod(meqn)

      common /cparam/ tau, p1, p2
c
c      tau = 0.1
c
      prod(1) = 0.d0
      prod(2) = -1.0/tau*q(2)
      prod(3) = -1.0/tau*q(3)

      return
      end
