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
      tau = 100000.2
c
c     # 2-stage Runge-Kutta method
c
      dt2    = dt/2.d0
      ndim   = 2

      do 10 i=1,mx
       do 10 j=1,my
        do 10 m=1,meqn
         qstar(m) = q(m,i,j) - dt2/tau * q(m,i,j)
c
c        # second stage
c
         q(m,i,j) = q(m,i,j) - dt/tau * qstar(m)
   10    continue

      return
      end
