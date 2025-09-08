module rk4
    use deriv
    implicit none
contains
    subroutine rk4_step(r, dt)
        real(8), dimension(3), intent(inout) :: r
        real(8), intent(in) :: dt
        real(8), dimension(3) :: k1, k2, k3, k4, r_temp, dot_r

        call compute_deriv(r, k1)
        r_temp = r + 0.5d0 * dt * k1
        call compute_deriv(r_temp, k2)
        r_temp = r + 0.5d0 * dt * k2
        call compute_deriv(r_temp, k3)
        r_temp = r + dt * k3
        call compute_deriv(r_temp, k4)

        dot_r = (k1 + 2.d0*k2 + 2.d0*k3 + k4) / 6.d0
        r = r + dt * dot_r
    end subroutine rk4_step
end module rk4