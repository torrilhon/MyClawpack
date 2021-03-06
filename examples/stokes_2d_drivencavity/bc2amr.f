c
c  Modified to give inflow boundary conditions at the left edge
c
c ------------------------------------------------------------------
c
      subroutine bc2amr(val,aux,nrow,ncol,meqn,naux,
     1                  hx, hy, level, time,  
     2                  xlo_patch,  xhi_patch,  
     3                  ylo_patch, yhi_patch)
c
c
c :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     Take a grid patch with mesh widths hx,hy, of dimensions nrow by
c     ncol,  and set the values of any piece of
c     of the patch which extends outside the physical domain 
c     using the boundary conditions. 
c
c     ------------------------------------------------
c     # Standard boundary condition choices for amr2ez in clawpack
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (bottom), 4 (top):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary conditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     #            =  5  sphere bcs (left half maps to right half of same 
c     #                  side, and vice versa), as if domain folded in half
c     ------------------------------------------------
c
c     The corners of the grid patch are at 
c        (xlo_patch,ylo_patch)  --  lower left corner
c        (xhi_patch,yhi_patch) --  upper right corner
c
c     The physical domain itself is a rectangle bounded by
c        (xlower,ylower)  -- lower left corner
c        (xupper,yupper)  -- upper right corner
c     
c     the picture is the following: 
c
c               _____________________ (xupper,yupper)
c              |                     |  
c          ____|____ (xhi_patch,yhi_patch)   
c          |   |    |                |
c          |   |    |                |
c          |   |    |                |
c          |___|____|                |
c (xlo_patch,ylo_patch)              |
c              |                     |
c              |_____________________|
c   (xlower,ylower)
c        
c
c     Any cells that lie outside the physical domain are ghost cells whose
c     values should be set in this routine.  This is tested for by comparing
c     xlo_patch with xlower to see if values need to be set at the left, as in
c     the figure above, and similarly at the other boundaries.
c
c     Patches are guaranteed to have at least 1 row of cells filled
c     with interior values so it is possible to  extrapolate. 
c     Fix trimbd if you want more than 1 row pre-set.
c
c     Make sure the order the boundaries are specified is correct
c     so that diagonal corner cells are also properly taken care of.
c
c     Periodic boundaries are set before calling this routine, so if the
c     domain is periodic in one direction only you
c     can safely extrapolate in the other direction. 
c
c     Don't overwrite ghost cells in periodic directions!
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module, only: mthbc,xlower,ylower,xupper,yupper
      use amr_module, only: xperdom,yperdom,spheredom

      implicit double precision (a-h,o-z)

      dimension val(meqn,nrow,ncol), aux(naux,nrow,ncol)

      common /cparam/ tau, uW

      hxmarg = hx*.01
      hymarg = hy*.01

      if (xperdom .and. (yperdom .or. spheredom)) go to 499
c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      if (xlo_patch .ge. xlower-hxmarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 199
         endif
c
c     # number of grid cells from this patch lying outside physical domain:
      nxl = (xlower+hxmarg-xlo_patch)/hx
      if (nxl > 2) then
          write(6,*) '*** unexpected value nxl = ',nxl
          stop
          endif
c
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user definition:  (2 ghostcells required !!!)              ** LEFT **
      do j = 1,ncol
         do m=1,meqn
            val(m,nxl,j) = val(m,nxl+1,j)
         end do
         val(2,nxl,j) = -val(2,nxl,j)
         val(3,nxl,j) = -val(3,nxl,j)
         do m=1,meqn
            val(m,1,j) = val(m,nxl,j)
         end do
      enddo      
      go to 199
c
  110 continue
c     # zero-order extrapolation:      ** LEFT **
      do 115 j = 1,ncol
         do 115 i=1,nxl
            do 115 m=1,meqn
               val(m,i,j) = val(m,nxl+1,j)
  115       continue
      go to 199

  120 continue
c     # periodic:   handled elsewhere in amr
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):      ** LEFT **
      do 135 j = 1,ncol
         do 135 i=1,nxl
            do 135 m=1,meqn
               val(m,i,j) = val(m,2*nxl+1-i,j)
  135       continue
c     # negate the normal velocity:
      do 136 j = 1,ncol
         do 136 i=1,nxl
            val(2,i,j) = -val(2,i,j)
            val(3,i,j) = -val(3,i,j)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      if (xhi_patch .le. xupper+hxmarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 299
         endif
c
c     # number of grid cells lying outside physical domain:
      nxr = (xhi_patch - xupper + hxmarg)/hx
      ibeg = max0(nrow-nxr+1, 1)
c
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user definition:  (2 ghostcells required !!!)              ** RIGHT **
      do j = 1,ncol
         do m=1,meqn
            val(m,ibeg,j) = val(m,ibeg-1,j)
         end do
         val(2,ibeg,j) = -val(2,ibeg,j)
         val(3,ibeg,j) = -val(3,ibeg,j)
         do m=1,meqn
            val(m,nrow,j) = val(m,ibeg,j)
         end do
      enddo      
      go to 299

  210 continue
c     # zero-order extrapolation:      ** RIGHT **
      do 215 j = 1,ncol
         do 215 i=ibeg,nrow
            do 215 m=1,meqn
               val(m,i,j) = val(m,ibeg-1,j)
  215       continue
      go to 299
  
  220 continue
c     # periodic:   handled elsewhere in amr
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):      ** RIGHT **
      do 235 j = 1,ncol
         do 235 i=ibeg,nrow
            do 235 m=1,meqn
               val(m,i,j) = val(m,2*ibeg-1-i,j)
  235       continue
c     # negate the normal velocity:
      do 236 j = 1,ncol
         do 236 i=ibeg,nrow
            val(2,i,j) = -val(2,i,j)
            val(3,i,j) = -val(3,i,j)
  236    continue
      go to 299
   
  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      if (ylo_patch .ge. ylower-hymarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 399
         endif
c
c     # number of grid cells lying outside physical domain:
      nyb = (ylower+hymarg-ylo_patch)/hy
c
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user definition:   (2 ghostcells required !!!              ** BOTTOM **
      do i = 1,nrow
         do m=1,meqn
            val(m,i,nyb) = val(m,i,nyb+1)
         end do
         val(2,i,nyb) = -val(2,i,nyb)
         val(3,i,nyb) = -val(3,i,nyb)
         do m=1,meqn
            val(m,i,1) = val(m,i,nyb)
         end do
      enddo
      go to 399
c
  310 continue
c     # zero-order extrapolation:      ** BOTTOM **
      do 315 j=1,nyb
         do 315 i=1,nrow
            do 315 m=1,meqn
                val(m,i,j) = val(m,i,nyb+1)
  315       continue
      go to 399

  320 continue
c     # periodic:   handled elsewhere in amr
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):      ** BOTTOM **
      do 335 j=1,nyb
         do 335 i=1,nrow
            do 335 m=1,meqn
               val(m,i,j) =  val(m,i,2*nyb+1-j)
  335       continue
c     # negate the normal velocity:
      do 336 j=1,nyb
         do 336 i=1,nrow
            val(2,i,j) = -val(2,i,j)
            val(3,i,j) = -val(3,i,j)
  336    continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      if (yhi_patch .le. yupper+hymarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 499
         endif
c
c     # number of grid cells lying outside physical domain:
      nyt = (yhi_patch - yupper + hymarg)/hy
      jbeg = max0(ncol-nyt+1, 1)
c
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user definition:   (2 ghostcells required !!!              ** TOP **
      do i = 1,nrow
         do m=1,meqn
            val(m,i,jbeg) = val(m,i,jbeg-1)
         end do
         val(2,i,jbeg) = 2.d0-val(2,i,jbeg)
         val(3,i,jbeg) = -val(3,i,jbeg)
         do m=1,meqn
            val(m,i,ncol) = val(m,i,jbeg)
         end do
      enddo
      go to 499

  410 continue
c     # zero-order extrapolation:      ** TOP **
      do 415 j=jbeg,ncol
         do 415 i=1,nrow
            do 415 m=1,meqn
               val(m,i,j) =  val(m,i,jbeg-1)
  415       continue
      go to 499

  420 continue
c     # periodic:   handled elsewhere in amr
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):      ** TOP **
      do 435 j=jbeg,ncol
         do 435 i=1,nrow
            do 435 m=1,meqn
               val(m,i,j) =  val(m,i,2*jbeg-1-j)
  435       continue
c     # negate the normal velocity:
      do 436 j=jbeg,ncol
         do 436 i=1,nrow
            val(2,i,j) = 2.0*uW-val(2,i,j)
            val(3,i,j) = -val(3,i,j)
  436    continue
      go to 499

  499 continue

      return
      end



