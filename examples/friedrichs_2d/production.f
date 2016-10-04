c
c
c =========================================================
      subroutine production(meqn,xi,yj,t,q,prod)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(meqn)
      dimension prod(meqn)

      common /cparam/ tau, Re
c
      p  = q(1)
      ux = q(2)
      uy = q(3)
      sigxx = q(4)
      sigxy = q(5)
      sigyx = q(6)
      sigyy = q(7)

      alpha = 0.d0
      if (t>4.0) alpha = Re;
      
      prod(1) = 0.d0
      prod(2) = 0.d0
      prod(3) = 0.d0
      prod(4) = -1.0/tau*(sigxx - alpha*ux*ux)
      prod(5) = -1.0/tau*(sigxy - alpha*ux*uy)
      prod(6) = -1.0/tau*(sigyx - alpha*uy*ux)
      prod(7) = -1.0/tau*(sigyy - alpha*uy*uy)

      return
      end
