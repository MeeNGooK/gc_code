! Assumption 1. u is treated as a constant when taking spatial derivatives
! Assumption 2. Time Independent
! Assumption 3. mu is an invariant. That is, u and mu are declared at initialization
! Assumption 4. The sign of sqrt{g} is not corrected
program main
    use calc_mod
    use fourier1
    use nc_reader
    use cubic_spline_mod
    use spline2_mod
    implicit none
    type(nc_data) :: nc

    real(8) :: m, q, eps, mu, v_perp,deltax   ! Constants
    real(8) :: x0(3)
    real(8) :: u,dudt, g, n_4   ! Updated scalars (u=velocity, mu=magnetic moment)
    real(8) :: b_para_n, b_n(3), bsup(3), bsub(3), b_abs, e_n(3), x(3), dxdt(3), grad_b(3)
    real(8) :: x_1(3), x_2(3), x_3(3), x_4(3), u_1, u_2, u_3, u_4
    real(8) :: n_1(3), n_2(3), n_3(3)  ! Note reference (reusable variables)
    real(8) :: r, z, zeta, dx
    integer :: sindex

!!!! Below are variables for storing spline coefficients
!!!! for spline3 coeffs

    real(8), allocatable :: b_abs_a(:,:), b_abs_b(:,:), b_abs_c(:,:), b_abs_d(:,:)
    real(8), allocatable :: b_abs_arr(:)
    real(8), allocatable :: g_a(:,:), g_b(:,:), g_c(:,:), g_d(:,:)
    real(8), allocatable :: g_arr(:)

    real(8), allocatable :: bsubs_a(:,:), bsubs_b(:,:), bsubs_c(:,:), bsubs_d(:,:)
    real(8), allocatable :: bsubs_arr(:)
    real(8), allocatable :: bsubu_a(:,:), bsubu_b(:,:), bsubu_c(:,:), bsubu_d(:,:)
    real(8), allocatable :: bsubu_arr(:)
    real(8), allocatable :: bsubv_a(:,:), bsubv_b(:,:), bsubv_c(:,:), bsubv_d(:,:)
    real(8), allocatable :: bsubv_arr(:)

    real(8), allocatable :: bsupu_a(:,:), bsupu_b(:,:), bsupu_c(:,:), bsupu_d(:,:)
    real(8), allocatable :: bsupu_arr(:)
    real(8), allocatable :: bsupv_a(:,:), bsupv_b(:,:), bsupv_c(:,:), bsupv_d(:,:)
    real(8), allocatable :: bsupv_arr(:)

    real(8), allocatable :: rmnc_a(:,:), rmnc_b(:,:), rmnc_c(:,:), rmnc_d(:,:)
    real(8), allocatable :: rmnc_arr(:)
    real(8), allocatable :: zmns_a(:,:), zmns_b(:,:), zmns_c(:,:), zmns_d(:,:)
    real(8), allocatable :: zmns_arr(:)


!!! for euler method setting

    real(8) :: dt
    INTEGER :: steps, istep, jstep, unit_out





    ! time measurement variables
    integer :: time_start, time_end, clock_rate
    real(8) :: elapsed

    INTEGER :: sig   ! s, theta, zeta order with eps(s, theta, zeta)=1 (for now)
    call system_clock(time_start, clock_rate)



    sig=1.0d0

    q = 1.602176634d-19  ! 전자 전하 electron charge
    m = 9.10938356d-31   ! 전자 질량 electron mass
    eps=m/q              ! epsilon

    ! Setting (arbitary values)@@@@@@@@@@@@@@@@@@@@@@@@@@@2
    x0=(/0.1, 0.1, 0.1/)  !(s, theta, zeta)
    u=1.0d5
    v_perp=1.0d5

    dt=1.0d-9
    steps=10000
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Header
    unit_out = 10
    open(unit=unit_out, file="dat.dat", status="replace", action="write", form="formatted")


    write(unit_out,*) "dt=", dt
    write(unit_out,*) "steps=", steps
    write(unit_out,*)
    write(unit_out,'(A)') "istep          p0              t0              rout            zout"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    call read_nc("wout_iota0_42.nc", nc)
    mu=0.5d0*m*(v_perp**2)/(compc(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bmnc))
    deltax=nc%phi(2)-nc%phi(1)

! Saving Spline3 Coeffs

    allocate(b_abs_a(nc%mnmax_nyq, nc%ns), b_abs_b(nc%mnmax_nyq, nc%ns-1), &
    b_abs_c(nc%mnmax_nyq, nc%ns), b_abs_d(nc%mnmax_nyq, nc%ns-1))
    allocate(g_a(nc%mnmax_nyq, nc%ns), g_b(nc%mnmax_nyq, nc%ns-1), &
    g_c(nc%mnmax_nyq, nc%ns), g_d(nc%mnmax_nyq, nc%ns-1))

    allocate(bsubs_a(nc%mnmax_nyq, nc%ns), bsubs_b(nc%mnmax_nyq, nc%ns-1), &
    bsubs_c(nc%mnmax_nyq, nc%ns), bsubs_d(nc%mnmax_nyq, nc%ns-1))
    allocate(bsubu_a(nc%mnmax_nyq, nc%ns), bsubu_b(nc%mnmax_nyq, nc%ns-1), &
    bsubu_c(nc%mnmax_nyq, nc%ns), bsubu_d(nc%mnmax_nyq, nc%ns-1))
    allocate(bsubv_a(nc%mnmax_nyq, nc%ns), bsubv_b(nc%mnmax_nyq, nc%ns-1), &
    bsubv_c(nc%mnmax_nyq, nc%ns), bsubv_d(nc%mnmax_nyq, nc%ns-1))

    allocate(bsupu_a(nc%mnmax_nyq, nc%ns), bsupu_b(nc%mnmax_nyq, nc%ns-1), &
    bsupu_c(nc%mnmax_nyq, nc%ns), bsupu_d(nc%mnmax_nyq, nc%ns-1))
    allocate(bsupv_a(nc%mnmax_nyq, nc%ns), bsupv_b(nc%mnmax_nyq, nc%ns-1), &
    bsupv_c(nc%mnmax_nyq, nc%ns), bsupv_d(nc%mnmax_nyq, nc%ns-1))


    allocate(rmnc_a(nc%mnmax, nc%ns), rmnc_b(nc%mnmax, nc%ns-1), &
    rmnc_c(nc%mnmax, nc%ns), rmnc_d(nc%mnmax, nc%ns-1))
    allocate(zmns_a(nc%mnmax, nc%ns), zmns_b(nc%mnmax, nc%ns-1), &
    zmns_c(nc%mnmax, nc%ns), zmns_d(nc%mnmax, nc%ns-1))



    allocate(b_abs_arr(nc%mnmax_nyq), g_arr(nc%mnmax_nyq), &
    bsubs_arr(nc%mnmax_nyq), bsubu_arr(nc%mnmax_nyq), bsubv_arr(nc%mnmax_nyq), &
    bsupu_arr(nc%mnmax_nyq), bsupv_arr(nc%mnmax_nyq), &
    rmnc_arr(nc%mnmax), zmns_arr(nc%mnmax))

    


    do jstep=1, nc%mnmax_nyq
        call spline_getter(deltax, nc%bmnc(jstep, :), b_abs_a(jstep, :), b_abs_b(jstep, :), &
        b_abs_c(jstep, :), b_abs_d(jstep, :))
        call spline_getter(deltax, nc%gmnc(jstep, :), g_a(jstep, :), g_b(jstep, :), &
        g_c(jstep, :), g_d(jstep, :))
        call spline_getter(deltax, nc%bsubsmns(jstep, :), bsubs_a(jstep, :), bsubs_b(jstep, :), &
        bsubs_c(jstep, :), bsubs_d(jstep, :))
        call spline_getter(deltax, nc%bsubumnc(jstep, :), bsubu_a(jstep, :), bsubu_b(jstep, :), &
        bsubu_c(jstep, :), bsubu_d(jstep, :))
        call spline_getter(deltax, nc%bsubvmnc(jstep, :), bsubv_a(jstep, :), bsubv_b(jstep, :), &
        bsubv_c(jstep, :), bsubv_d(jstep, :))
        call spline_getter(deltax, nc%bsupumnc(jstep, :), bsupu_a(jstep, :), bsupu_b(jstep, :), &
        bsupu_c(jstep, :), bsupu_d(jstep, :))
        call spline_getter(deltax, nc%bsupvmnc(jstep, :), bsupv_a(jstep, :), bsupv_b(jstep, :), &
        bsupv_c(jstep, :), bsupv_d(jstep, :))




    end do

    do jstep=1, nc%mnmax
        call spline_getter(deltax, nc%rmnc(jstep, :), rmnc_a(jstep, :), rmnc_b(jstep, :), &
        rmnc_c(jstep, :), rmnc_d(jstep, :))
        call spline_getter(deltax, nc%zmns(jstep, :), zmns_a(jstep, :), zmns_b(jstep, :), &
        zmns_c(jstep, :), zmns_d(jstep, :))
    end do





! Numerical integration using Euler method

    do istep=1,steps

        ! Getting range index of s
        sindex=int(x0(1)/deltax)+1
        if (sindex < 1) sindex = 1
        if (sindex > nc%ns-1) sindex = nc%ns-1
        dx=x0(1)-nc%phi(sindex)

        ! Evaluate spline3 at x0(1) for all required quantities
        
        !@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (mod(istep,1000)==0) then
            print*, "istep=", istep
            print*, "x=", x0
        end if
        !@@@@@@@@@@@@@@@@@@@@@@@@@@

        grad_b(1)=compc_ds(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bmnc)
        grad_b(2)=compc_dtheta(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bmnc)
        grad_b(3)=compc_dzeta(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bmnc)

        bsub(1)=comps(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bsubsmns)
        bsub(2)=compc(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bsubumnc)
        bsub(3)=compc(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bsubvmnc)
        bsup(1)=0.0d0
        bsup(2)=compc(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bsupumnc)
        bsup(3)=compc(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bsupvmnc)
        b_abs=compc(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%bmnc)
        g=compc(nc%phi, x0(1), x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, nc%gmnc)
        call crossp(grad_b, bsub, sig*g, n_1)   ! n_1 get
        n_2=-n_1*(mu/q)  ! 2행
        b_n=bsup-n_2*(mu*q)/(b_abs**2)  ! 3행
        b_para_n=b_abs   ! 4행
        e_n=grad_b*(mu/q) ! 5행
        call crossp(e_n, bsub, sig*g, n_3)   ! n_3 get, 6행
        n_4=b_n(1)*e_n(1)+b_n(2)*e_n(2)+b_n(3)*e_n(3)  ! 7행
        dxdt=((bsup-n_1*u*eps/(b_abs**2))*u/b_abs)+n_3/(b_abs**2)  ! 8행
        dudt=n_4*q/(m*b_abs)  ! 9행



        if (mod(istep,1000)==0) then
            print*, "grad_b=", grad_b
            print*, "bsub=", bsub
            print*, "bsup=", bsup
            print*, "b_abs=", b_abs
            print*, "g=", g
            print*, "n_1=", n_1
            print*, "n_2=", n_2
            print*, "b_n=", b_n
            print*, "b_para_n=", b_para_n
            print*, "e_n=", e_n
            print*, "n_3=", n_3
            print*, "n_4=", n_4
            print*, "dxdt=", dxdt
            print*, "dudt=", dudt
            print*, "---------------------"
        end if



        ! Update x0 and u
        u=u+dudt*dt
        x0=x0+dxdt*dt

        r=compc(nc%phi, x0(1), x0(2), x0(3), nc%xm, nc%xn, nc%rmnc)
        z=comps(nc%phi, x0(1), x0(2), x0(3), nc%xm, nc%xn, nc%zmns)
        write(unit_out,'(I8,3ES20.10)') istep,r,z,x0(3)

    end do

    print*, compc(nc%phi, 0.02d0, 0.1d0, 0.1d0, nc%xm_nyq, nc%xn_nyq, nc%gmnc)
    print*, nc%phi(2)-nc%phi(1)
    call system_clock(time_end)
    elapsed = real(time_end - time_start,8) / real(clock_rate,8)
    print *, 'Elapsed time (seconds): ', elapsed

end program main