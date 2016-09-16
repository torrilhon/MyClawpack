c
c
c =========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      dimension  aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
c
c     # source terms for 2d advaction equations
c
      dimension qstar(meqn)
      common /cparam/ ubar,vbar
      tau = 0.1
c
c     # 2-stage Runge-Kutta method
c
      dt2    = dt/2.d0
      ndim   = 2

      do 10 i=1,mx
       do 10 j=1,my
         qstar(1) = q(1,i,j) - dt2/tau * q(1,i,j)
c
c        # second stage
c
         q(1,i,j) = q(1,i,j) - dt/tau * qstar(1)
   10    continue

      return
      end
