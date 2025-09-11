! Assumption 1. u is treated as a constant when taking spatial derivatives
! Assumption 2. Time Independent
! Assumption 3. mu is an invariant. That is, u and mu are declared at initialization
! Assumption 4. The sign of sqrt{g} is not corrected
program main
    use calc_mod
    use fourier2
    use nc_reader
    use cubic_spline_mod
    use spline2_mod
    use device_coeffs
    use gpu_eval
    implicit none
    type(nc_data) :: nc

    real(8) :: m, q, eps, mu, v_perp,deltax   ! Constants
    real(8) :: x0(3)
    real(8) :: u,dudt, g, n_4   ! Updated scalars (u=velocity, mu=magnetic moment)
    real(8) :: b_para_n, b_n(3), bsup(3), bsub(3), b_abs, e_n(3), x(3), dxdt(3), grad_b(3)
    real(8) :: x_1(3), x_2(3), x_3(3), x_4(3),  u_1, u_2, u_3, u_4
    real(8) :: x_temp(3), u_temp
    real(8) :: n_1(3), n_2(3), n_3(3)  ! Note reference (reusable variables)
    real(8) :: r, z, zeta, dx
    integer :: sindex, zipper  ! zipper for writing(1 write for every 'zipper' steps)

!!!! Below are variables for storing spline coefficients
!!!! for spline3 coeffs

    real(8), allocatable :: b_abs_a(:,:), b_abs_b(:,:), b_abs_c(:,:), b_abs_d(:,:)
    real(8), allocatable :: b_abs_arr(:), b_abs_ds_arr(:)
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
    INTEGER :: steps, istep, jstep, kstep, lstep, unit_out





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
    x0=(/0.06, 0.1, 0.1/)  !(s, theta, zeta)
    u=1.0d5
    v_perp=0.5d4

    dt=1.0d-10
    steps=1000000
    zipper=100
    
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



    allocate(b_abs_arr(nc%mnmax_nyq), b_abs_ds_arr(nc%mnmax_nyq), g_arr(nc%mnmax_nyq), &
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

    ! end of saving spline3 coeffs.



! Numerical integration using rk4 method
    call init_device_coeffs(nc%ns, nc%mnmax, nc%mnmax_nyq, &
        b_abs_a,b_abs_b,b_abs_c,b_abs_d, g_a,g_b,g_c,g_d, &
        bsubs_a,bsubs_b,bsubs_c,bsubs_d, bsubu_a,bsubu_b,bsubu_c,bsubu_d, &
        bsubv_a,bsubv_b,bsubv_c,bsubv_d, bsupu_a,bsupu_b,bsupu_c,bsupu_d, &
        bsupv_a,bsupv_b,bsupv_c,bsupv_d, rmnc_a,rmnc_b,rmnc_c,rmnc_d, zmns_a,zmns_b,zmns_c,zmns_d)
    do istep=1,steps
        


        !!!!!!!! RK4 1st

        ! Getting range index of s
        sindex=int(x0(1)/deltax)+1
        if (sindex < 1) sindex = 1
        if (sindex > nc%ns-1) sindex = nc%ns-1
        dx=x0(1)-nc%phi(sindex)

        ! Evaluate spline3 at x0(1) for all required quantities
        call eval_at_sindex(sindex, dx, b_abs_arr, b_abs_ds_arr, g_arr, bsubs_arr, bsubu_arr, bsubv_arr, bsupu_arr, bsupv_arr, rmnc_arr, zmns_arr)
       
        !!!!!!!!!!!!! Writing r,z, zeta   !!!!!!!!!!!!!!!!
        r=compc_arg(x0(2), x0(3), nc%xm, nc%xn, rmnc_arr)
        z=comps_arg(x0(2), x0(3), nc%xm, nc%xn, zmns_arr)

        if (mod(istep,zipper)==0) then
            write(unit_out,'(I8,3ES20.10)') istep,r,z,x0(3)
        end if


        grad_b(1)=compc_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, b_abs_ds_arr)
        grad_b(2)=compc_dtheta_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        grad_b(3)=compc_dzeta_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)

        bsub(1)=comps_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, bsubs_arr)
        bsub(2)=compc_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, bsubu_arr)
        bsub(3)=compc_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, bsubv_arr)
        bsup(1)=0.0d0
        bsup(2)=compc_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, bsupu_arr)
        bsup(3)=compc_arg(x0(2),x0(3),nc%xm_nyq, nc%xn_nyq, bsupv_arr)
        b_abs=compc_arg(x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        g=compc_arg(x0(2), x0(3), nc%xm_nyq, nc%xn_nyq, g_arr)



        call crossp(grad_b, bsub, sig*g, n_1)   ! n_1 get
        n_2=-n_1*(mu/q)  ! 2행
        b_n=bsup-n_2*(mu*q)/(b_abs**2)  ! 3행
        b_para_n=b_abs   ! 4행
        e_n=grad_b*(mu/q) ! 5행
        call crossp(e_n, bsub, sig*g, n_3)   ! n_3 get, 6행
        n_4=b_n(1)*e_n(1)+b_n(2)*e_n(2)+b_n(3)*e_n(3)  ! 7행
        dxdt=((bsup-n_1*u_temp*eps/(b_abs**2))*u_temp/b_abs)+n_3/(b_abs**2)  ! 8행
        dudt=n_4*q/(m*b_abs)  ! 9행

        x_1=dxdt
        u_1=dudt

        !!!!!!!! RK4 2nd

        x_temp=x0+(dt/2.0d0)*x_1
        u_temp=u+(dt/2.0d0)*u_1
        
        ! Getting range index of s
        sindex=int(x_temp(1)/deltax)+1
        if (sindex < 1) sindex = 1
        if (sindex > nc%ns-1) sindex = nc%ns-1
        dx=x_temp(1)-nc%phi(sindex)

        ! Evaluate spline3 at x0(1) for all required quantities
        do kstep=1,nc%mnmax_nyq
            b_abs_arr(kstep)=b_abs_d(kstep,sindex)*(dx**3)+b_abs_c(kstep,sindex)*(dx**2)+b_abs_b(kstep,sindex)*dx+b_abs_a(kstep,sindex)
            b_abs_ds_arr(kstep)=3.0d0*b_abs_d(kstep,sindex)*(dx**2)+2.0d0*b_abs_c(kstep,sindex)*dx+b_abs_b(kstep,sindex)

            g_arr(kstep)=g_d(kstep,sindex)*(dx**3)+g_c(kstep,sindex)*(dx**2)+g_b(kstep,sindex)*dx+g_a(kstep,sindex)
            bsubs_arr(kstep)=bsubs_d(kstep,sindex)*(dx**3)+bsubs_c(kstep,sindex)*(dx**2)+bsubs_b(kstep,sindex)*dx+bsubs_a(kstep,sindex)
            bsubu_arr(kstep)=bsubu_d(kstep,sindex)*(dx**3)+bsubu_c(kstep,sindex)*(dx**2)+bsubu_b(kstep,sindex)*dx+bsubu_a(kstep,sindex)
            bsubv_arr(kstep)=bsubv_d(kstep,sindex)*(dx**3)+bsubv_c(kstep,sindex)*(dx**2)+bsubv_b(kstep,sindex)*dx+bsubv_a(kstep,sindex)

            bsupu_arr(kstep)=bsupu_d(kstep,sindex)*(dx**3)+bsupu_c(kstep,sindex)*(dx**2)+bsupu_b(kstep,sindex)*dx+bsupu_a(kstep,sindex) 
            bsupv_arr(kstep)=bsupv_d(kstep,sindex)*(dx**3)+bsupv_c(kstep,sindex)*(dx**2)+bsupv_b(kstep,sindex)*dx+bsupv_a(kstep,sindex) 

        end do

        do lstep=1,nc%mnmax
            rmnc_arr(lstep)=rmnc_d(lstep,sindex)*(dx**3)+rmnc_c(lstep,sindex)*(dx**2)+rmnc_b(lstep,sindex)*dx+rmnc_a(lstep,sindex)
            zmns_arr(lstep)=zmns_d(lstep,sindex)*(dx**3)+zmns_c(lstep,sindex)*(dx**2)+zmns_b(lstep,sindex)*dx+zmns_a(lstep,sindex)
        end do        



        grad_b(1)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_ds_arr)
        grad_b(2)=compc_dtheta_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        grad_b(3)=compc_dzeta_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)

        bsub(1)=comps_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubs_arr)
        bsub(2)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubu_arr)
        bsub(3)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubv_arr)
        bsup(1)=0.0d0
        bsup(2)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsupu_arr)
        bsup(3)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsupv_arr)
        b_abs=compc_arg(x_temp(2), x_temp(3), nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        g=compc_arg(x_temp(2), x_temp(3), nc%xm_nyq, nc%xn_nyq, g_arr)



        call crossp(grad_b, bsub, sig*g, n_1)   ! n_1 get
        n_2=-n_1*(mu/q)  ! 2행
        b_n=bsup-n_2*(mu*q)/(b_abs**2)  ! 3행
        b_para_n=b_abs   ! 4행
        e_n=grad_b*(mu/q) ! 5행
        call crossp(e_n, bsub, sig*g, n_3)   ! n_3 get, 6행
        n_4=b_n(1)*e_n(1)+b_n(2)*e_n(2)+b_n(3)*e_n(3)  ! 7행
        dxdt=((bsup-n_1*u_temp*eps/(b_abs**2))*u_temp/b_abs)+n_3/(b_abs**2)  ! 8행
        dudt=n_4*q/(m*b_abs)  ! 9행

        x_2=dxdt
        u_2=dudt


        !!!!!!!! RK4 3rd
        
        x_temp=x0+(dt/2.0d0)*x_2
        u_temp=u+(dt/2.0d0)*u_2
        
        ! Getting range index of s
        sindex=int(x_temp(1)/deltax)+1
        if (sindex < 1) sindex = 1
        if (sindex > nc%ns-1) sindex = nc%ns-1
        dx=x_temp(1)-nc%phi(sindex)

        ! Evaluate spline3 at x0(1) for all required quantities
        do kstep=1,nc%mnmax_nyq
            b_abs_arr(kstep)=b_abs_d(kstep,sindex)*(dx**3)+b_abs_c(kstep,sindex)*(dx**2)+b_abs_b(kstep,sindex)*dx+b_abs_a(kstep,sindex)
            b_abs_ds_arr(kstep)=3.0d0*b_abs_d(kstep,sindex)*(dx**2)+2.0d0*b_abs_c(kstep,sindex)*dx+b_abs_b(kstep,sindex)

            g_arr(kstep)=g_d(kstep,sindex)*(dx**3)+g_c(kstep,sindex)*(dx**2)+g_b(kstep,sindex)*dx+g_a(kstep,sindex)
            bsubs_arr(kstep)=bsubs_d(kstep,sindex)*(dx**3)+bsubs_c(kstep,sindex)*(dx**2)+bsubs_b(kstep,sindex)*dx+bsubs_a(kstep,sindex)
            bsubu_arr(kstep)=bsubu_d(kstep,sindex)*(dx**3)+bsubu_c(kstep,sindex)*(dx**2)+bsubu_b(kstep,sindex)*dx+bsubu_a(kstep,sindex)
            bsubv_arr(kstep)=bsubv_d(kstep,sindex)*(dx**3)+bsubv_c(kstep,sindex)*(dx**2)+bsubv_b(kstep,sindex)*dx+bsubv_a(kstep,sindex)

            bsupu_arr(kstep)=bsupu_d(kstep,sindex)*(dx**3)+bsupu_c(kstep,sindex)*(dx**2)+bsupu_b(kstep,sindex)*dx+bsupu_a(kstep,sindex) 
            bsupv_arr(kstep)=bsupv_d(kstep,sindex)*(dx**3)+bsupv_c(kstep,sindex)*(dx**2)+bsupv_b(kstep,sindex)*dx+bsupv_a(kstep,sindex) 

        end do

        do lstep=1,nc%mnmax
            rmnc_arr(lstep)=rmnc_d(lstep,sindex)*(dx**3)+rmnc_c(lstep,sindex)*(dx**2)+rmnc_b(lstep,sindex)*dx+rmnc_a(lstep,sindex)
            zmns_arr(lstep)=zmns_d(lstep,sindex)*(dx**3)+zmns_c(lstep,sindex)*(dx**2)+zmns_b(lstep,sindex)*dx+zmns_a(lstep,sindex)
        end do        



        grad_b(1)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_ds_arr)
        grad_b(2)=compc_dtheta_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        grad_b(3)=compc_dzeta_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)

        bsub(1)=comps_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubs_arr)
        bsub(2)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubu_arr)
        bsub(3)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubv_arr)
        bsup(1)=0.0d0
        bsup(2)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsupu_arr)
        bsup(3)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsupv_arr)
        b_abs=compc_arg(x_temp(2), x_temp(3), nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        g=compc_arg(x_temp(2), x_temp(3), nc%xm_nyq, nc%xn_nyq, g_arr)



        call crossp(grad_b, bsub, sig*g, n_1)   ! n_1 get
        n_2=-n_1*(mu/q)  ! 2행
        b_n=bsup-n_2*(mu*q)/(b_abs**2)  ! 3행
        b_para_n=b_abs   ! 4행
        e_n=grad_b*(mu/q) ! 5행
        call crossp(e_n, bsub, sig*g, n_3)   ! n_3 get, 6행
        n_4=b_n(1)*e_n(1)+b_n(2)*e_n(2)+b_n(3)*e_n(3)  ! 7행
        dxdt=((bsup-n_1*u_temp*eps/(b_abs**2))*u_temp/b_abs)+n_3/(b_abs**2)  ! 8행
        dudt=n_4*q/(m*b_abs)  ! 9행

        x_3=dxdt
        u_3=dudt


        !!!!!!!! RK4 4th
        
        
        x_temp=x0+dt*x_3
        u_temp=u+dt*u_3
        
        ! Getting range index of s
        sindex=int(x_temp(1)/deltax)+1
        if (sindex < 1) sindex = 1
        if (sindex > nc%ns-1) sindex = nc%ns-1
        dx=x_temp(1)-nc%phi(sindex)

        ! Evaluate spline3 at x0(1) for all required quantities
        do kstep=1,nc%mnmax_nyq
            b_abs_arr(kstep)=b_abs_d(kstep,sindex)*(dx**3)+b_abs_c(kstep,sindex)*(dx**2)+b_abs_b(kstep,sindex)*dx+b_abs_a(kstep,sindex)
            b_abs_ds_arr(kstep)=3.0d0*b_abs_d(kstep,sindex)*(dx**2)+2.0d0*b_abs_c(kstep,sindex)*dx+b_abs_b(kstep,sindex)

            g_arr(kstep)=g_d(kstep,sindex)*(dx**3)+g_c(kstep,sindex)*(dx**2)+g_b(kstep,sindex)*dx+g_a(kstep,sindex)
            bsubs_arr(kstep)=bsubs_d(kstep,sindex)*(dx**3)+bsubs_c(kstep,sindex)*(dx**2)+bsubs_b(kstep,sindex)*dx+bsubs_a(kstep,sindex)
            bsubu_arr(kstep)=bsubu_d(kstep,sindex)*(dx**3)+bsubu_c(kstep,sindex)*(dx**2)+bsubu_b(kstep,sindex)*dx+bsubu_a(kstep,sindex)
            bsubv_arr(kstep)=bsubv_d(kstep,sindex)*(dx**3)+bsubv_c(kstep,sindex)*(dx**2)+bsubv_b(kstep,sindex)*dx+bsubv_a(kstep,sindex)

            bsupu_arr(kstep)=bsupu_d(kstep,sindex)*(dx**3)+bsupu_c(kstep,sindex)*(dx**2)+bsupu_b(kstep,sindex)*dx+bsupu_a(kstep,sindex) 
            bsupv_arr(kstep)=bsupv_d(kstep,sindex)*(dx**3)+bsupv_c(kstep,sindex)*(dx**2)+bsupv_b(kstep,sindex)*dx+bsupv_a(kstep,sindex) 

        end do

        do lstep=1,nc%mnmax
            rmnc_arr(lstep)=rmnc_d(lstep,sindex)*(dx**3)+rmnc_c(lstep,sindex)*(dx**2)+rmnc_b(lstep,sindex)*dx+rmnc_a(lstep,sindex)
            zmns_arr(lstep)=zmns_d(lstep,sindex)*(dx**3)+zmns_c(lstep,sindex)*(dx**2)+zmns_b(lstep,sindex)*dx+zmns_a(lstep,sindex)
        end do        




        grad_b(1)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_ds_arr)
        grad_b(2)=compc_dtheta_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        grad_b(3)=compc_dzeta_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, b_abs_arr)

        bsub(1)=comps_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubs_arr)
        bsub(2)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubu_arr)
        bsub(3)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsubv_arr)
        bsup(1)=0.0d0
        bsup(2)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsupu_arr)
        bsup(3)=compc_arg(x_temp(2),x_temp(3),nc%xm_nyq, nc%xn_nyq, bsupv_arr)
        b_abs=compc_arg(x_temp(2), x_temp(3), nc%xm_nyq, nc%xn_nyq, b_abs_arr)
        g=compc_arg(x_temp(2), x_temp(3), nc%xm_nyq, nc%xn_nyq, g_arr)



        call crossp(grad_b, bsub, sig*g, n_1)   ! n_1 get
        n_2=-n_1*(mu/q)  ! 2행
        b_n=bsup-n_2*(mu*q)/(b_abs**2)  ! 3행
        b_para_n=b_abs   ! 4행
        e_n=grad_b*(mu/q) ! 5행
        call crossp(e_n, bsub, sig*g, n_3)   ! n_3 get, 6행
        n_4=b_n(1)*e_n(1)+b_n(2)*e_n(2)+b_n(3)*e_n(3)  ! 7행
        dxdt=((bsup-n_1*u_temp*eps/(b_abs**2))*u_temp/b_abs)+n_3/(b_abs**2)  ! 8행
        dudt=n_4*q/(m*b_abs)  ! 9행

        x_4=dxdt
        u_4=dudt

        ! Final update
        x0=x0+(dt/6.0d0)*(x_1+2.0d0*x_2+2.0d0*x_3+x_4)
        u=u+(dt/6.0d0)*(u_1+2.0d0*u_2+2.0d0*u_3+u_4)


        if (mod(istep,100000)==0) then
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

        

    end do

    call system_clock(time_end)
    elapsed = real(time_end - time_start,8) / real(clock_rate,8)
    print *, 'Elapsed time (seconds): ', elapsed

end program main