c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux)
c     =====================================================
c
      implicit double precision (a-h,o-z)
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      dimension aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

c     # set all zero 
c     ---------------------------
      do 20 j=1,my
         yj = ylower + (j-0.5d0)*dy
         do 20 i=1,mx
            xi = xlower + (i-0.5d0)*dx
            do m=1,meqn
              q(m,i,j) = 0.d0
            end do
   20    continue

      return
      end
