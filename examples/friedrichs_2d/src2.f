c
c
c =========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      dimension  aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      dimension qstar(meqn)
      dimension prod(meqn)
      common /cparam/ ubar,vbar
c
c     # 2-stage Runge-Kutta method
c
      dt2    = dt/2.0
      ndim   = 2

      do 10 i=1-mbc,mx+mbc
       xi = xlower + (i-0.5d0)*dx
       do 10 j=1-mbc,my+mbc
        yj = ylower + (j-0.5d0)*dy
        
        do m=1,meqn
         qstar(m) = q(m,i,j)
        end do
        call production(meqn,xi,yj,qstar,prod)
        do m=1,meqn
         qstar(m) = q(m,i,j) + dt2 * prod(m)
        end do

        call production(meqn,xi,yj,qstar,prod)
        do m=1,meqn
         q(m,i,j) = q(m,i,j) + dt * prod(m)
        end do

   10    continue

      return
      end
