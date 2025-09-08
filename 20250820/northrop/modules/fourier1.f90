module fourier1
    use cubic_spline_mod
    implicit none
    ! 만약 필요하다면 sin계수의 미분도 만들자
    ! if needed, add some functions for sin coefficients
contains
    function compc(xs, s, theta, zeta,xm,xn,coeff) result(output)   ! cos 계수 받는 함수의 s, theta, zeta에서의 값
        integer, intent(in) :: xm(:), xn(:)   ! xs 자리에는 netcdf의 phi가 들어간다. (xs=nc.phi)
        real(8), intent(in) :: theta, zeta, s
        real(8), intent(in) :: xs(:)
        real(8), intent(in) :: coeff(:,:)        ! 2차원 계수 그대로 넣어줌(순서는 원래 계수, s 로 받음)
        real(8) :: output
        integer :: istep
        integer :: xm_l

        xm_l = size(xm)
        output = 0.0d0
        do istep = 1, xm_l
            output = output + spline_value(xs,coeff(istep, :), s)  * cos(xm(istep)*theta - xn(istep)*zeta*0.5)
        end do
    end function compc
    function compc_ds(xs, s, theta, zeta,xm,xn,coeff) result(output)
        integer, intent(in) :: xm(:), xn(:)
        real(8), intent(in) :: theta, zeta
        real(8), intent(in) :: s
        real(8), intent(in) :: xs(:)
        real(8), intent(in) :: coeff(:, :)
        real(8) :: output
        integer :: istep
        integer :: xm_l

        xm_l = size(xm)
        output = 0.0d0
        do istep = 1, xm_l
            output = output + spline_deriv(xs,coeff(istep, :), s) * cos(xm(istep)*theta - xn(istep)*zeta*0.5)
        end do
    end function compc_ds
    function compc_dtheta(xs,s,theta, zeta,xm,xn,coeff) result(output)
        integer, intent(in) :: xm(:), xn(:)
        real(8), intent(in) :: theta, zeta
        real(8), intent(in) :: s
        real(8), intent(in) :: xs(:)
        real(8), intent(in) :: coeff(:, :)
        real(8) :: output
        integer :: istep
        integer :: xm_l

        xm_l = size(xm)
        output = 0.0d0
        do istep = 1, xm_l
            output = output - spline_value(xs,coeff(istep, :), s) * xm(istep)*sin(xm(istep)*theta - xn(istep)*zeta*0.5)
        end do
    end function compc_dtheta
    function compc_dzeta(xs, s, theta, zeta,xm,xn,coeff) result(output)
        integer, intent(in) :: xm(:), xn(:)
        real(8), intent(in) :: theta, zeta
        real(8), intent(in) :: s
        real(8), intent(in) :: xs(:)
        real(8), intent(in) :: coeff(:, :)
        real(8) :: output
        integer :: istep
        integer :: xm_l

        xm_l = size(xm)
        output = 0.0d0
        do istep = 1, xm_l
            output = output - spline_value(xs,coeff(istep, :), s) * (xn(istep)*0.5)*sin(xm(istep)*theta - xn(istep)*zeta*0.5)
        end do
    end function compc_dzeta
    function comps(xs, s, theta2, zeta2,xm2,xn2,coeff2) result(output2)
        integer, intent(in) :: xm2(:), xn2(:)
        real(8), intent(in) :: theta2, zeta2
        real(8), intent(in) :: s
        real(8), intent(in) :: xs(:)
        real(8), intent(in) :: coeff2(:, :)
        real(8) :: output2
        integer :: istep2
        integer :: xm_l2

        xm_l2 = size(xm2)
        output2 = 0.0d0
        do istep2 = 1, xm_l2
            output2 = output2 + spline_value(xs,coeff2(istep2, :), s) * sin(xm2(istep2)*theta2 - xn2(istep2)*zeta2*0.5)
        end do
    end function comps
end module fourier1