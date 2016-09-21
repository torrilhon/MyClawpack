! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves,
!            s the speeds,
!
!            amdq = A^- Delta q,
!            apdq = A^+ Delta q,
!                   the decomposition of the flux difference
!                       f(qr(i-1)) - f(ql(i))
!                   into leftgoing and rightgoing parts respectively.
!

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr
! maux=0 and aux arrays are unused in this example.


    implicit double precision (a-h,o-z)
    common /cparam/ u,v

    dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension  apdq(meqn, 1-mbc:maxm+mbc)
    dimension  amdq(meqn, 1-mbc:maxm+mbc)

    do 30 i = 2-mbc, mx+mbc
        p  = ql(1,i) - qr(1,i-1)
        ux = ql(2,i) - qr(2,i-1)
        uy = ql(3,i) - qr(3,i-1)
        sigxx = ql(4,i) - qr(4,i-1)
        sigxy = ql(5,i) - qr(5,i-1)
        sigyx = ql(6,i) - qr(6,i-1)
        sigyy = ql(7,i) - qr(7,i-1)
                
        do k = 1, mwaves
          do m = 1, meqn
              wave(m,k,i) = 0.d0
          end do
        end do

        if (ixy == 1) then
          wm1 = 0.25*p - 0.5*ux/sqrt(2.d0) + 0.25*sigxx
          wp1 = 0.25*p + 0.5*ux/sqrt(2.d0) + 0.25*sigxx
          wm2 = 0.5*(uy-sigyx)
          wp2 = 0.5*(uy+sigyx)
!         # first wave
          wave(1,1,i) = wm1
          wave(2,1,i) = -sqrt(2.d0)*wm1
          wave(4,1,i) = wm1
!         # second wave
          wave(3,2,i) = wm2
          wave(6,2,i) = -wm2
!         # third wave
          wave(3,3,i) = wp2
          wave(6,3,i) = wp2
!         # fourth wave
          wave(1,4,i) = wp1
          wave(2,4,i) = sqrt(2.d0)*wp1
          wave(4,4,i) = wp1
        end if
        
        if (ixy == 2) then
          wm1 = 0.25*p - 0.5*uy/sqrt(2.d0) + 0.25*sigyy
          wp1 = 0.25*p + 0.5*uy/sqrt(2.d0) + 0.25*sigyy
          wm2 = 0.5*(ux-sigxy)
          wp2 = 0.5*(ux+sigxy)
!         # first wave
          wave(1,1,i) = wm1
          wave(3,1,i) = -sqrt(2.d0)*wm1
          wave(7,1,i) = wm1
!         # second wave
          wave(2,2,i) = wm2
          wave(5,2,i) = -wm2
!         # third wave
          wave(2,3,i) = wp2
          wave(5,3,i) = wp2
!         # fourth wave
          wave(1,4,i) = wp1
          wave(3,4,i) = sqrt(2.d0)*wp1
          wave(7,4,i) = wp1
        endif
    
!       # speeds
        s(1,i) = -sqrt(2.d0)
        s(2,i) = -1.d0
        s(3,i) = 1.d0
        s(4,i) = sqrt(2.d0)

!       # flux differences:
        do m = 1, meqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
        end do
        do k = 1, mwaves
          do m = 1, meqn
              amdq(m,i) = amdq(m,i) + dmin1(s(k,i),0.d0) * wave(m,k,i)
              apdq(m,i) = apdq(m,i) + dmax1(s(k,i),0.d0) * wave(m,k,i)
          end do
        end do
    
    30 END DO

    return
    end subroutine rpn2



