! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Riemann solver for the sample scalar equation
!  q_t + u*q_x + v*q_y = 0

! waves: 1
! equations: 1

! Conserved quantities:
!       1 q

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


!$$$      do 30 i = 2-mbc, mx+mbc-1
    do 30 i = 2-mbc, mx+mbc
!       # first wave
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        wave(2,1,i) = 0.d0
!       # second wave
        wave(1,2,i) = 0.d0
        wave(2,2,i) = ql(2,i) - qr(2,i-1)

        if (ixy == 1) then
            s(1,i) = 2*u
            s(2,i) = u/2
        else
            s(1,i) = v
            s(2,i) = v
        endif
    
!       # flux differences:
        do 31 m = 1, meqn
            amdq(m,i) = dmin1(s(1,i), 0.d0) * wave(m,1,i) + dmin1(s(2,i), 0.d0) * wave(m,2,i)
            apdq(m,i) = dmax1(s(1,i), 0.d0) * wave(m,1,i) + dmax1(s(2,i), 0.d0) * wave(m,2,i)
        31 END DO
    
    30 END DO

    return
    end subroutine rpn2
