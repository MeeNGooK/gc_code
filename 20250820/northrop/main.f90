!가정 1. u는 공간 미분 시에는 상수 취급
! 가정 2. Time Independent
! 
program main
    use calc_mod
    use nc_reader
    use cubic_spline_mod
    implicit none
    type(nc_data) :: nc

    real(8) :: mu, m, q, eps   ! 상수
    real(8) :: x0(3), v0(3)  ! 초기 위치, 속도(실제 변위)
    real(8) :: b_para_n, b_n(3), b(3), b_abs, u, e_n(3), x(3), dxdt(3), grad_b(3), dudt
    real(8) :: x_1(3), x_2(3), x_3(3), x_4(3), u_1, u_2, u_3, u_4
    real(8) :: n_1(3), n_2(3), n_3(3), n_4(3)    ! 노트 참고(재활용 가능한 변수들)

    real(8) :: dt
    INTEGER :: steps, istep, unit_out

    real(8), ALLOCATABLE :: b_supu_coeff(:), b_supv_coeff(:), g_coeff(:)




    ! 시간 측정 변수
    integer :: time_start, time_end, clock_rate
    real(8) :: elapsed

    INTEGER :: sig   ! 일단 s, theta, zeta 순으로 eps(s, theta, zeta)=1인 걸로

    sig=1

    q = 1.602176634d-19  ! 전자 전하 electron charge
    m = 9.10938356d-31   ! 전자 질량 electron mass
    eps=m/q              ! epsilon
    
    ! 여기서부터는 설정값(arbitary values)@@@@@@@@@@@@@@@@@@@@@@@@@@@2
    x0=[0.4, 0.4, 0.4]  !(s, theta, zeta)
    v0=[0.1, 0.1, 0.1]  !(ds/dt, dtheta/dt, dzeta/dt)
    dt=0.00001
    steps=10000
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 쓰기 준비
    unit_out = 10
    open(unit=unit_out, file="dat.dat", status="replace", action="write", form="formatted")


    write(unit_out,*) "dt=", dt
    write(unit_out,*) "steps=", steps
    write(unit_out,*)
    write(unit_out,'(A)') "istep          p0              t0              rout            zout"

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    call read_nc("wout_iota0_42.nc", nc)
    allocate(b_supu_coeff(nc%ns), b_supv_coeff(nc%ns), g_coeff(nc%ns))




    


end program main