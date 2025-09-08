module calc_mod
    implicit none
contains

    subroutine crossp(a,b,j,c)
        real(8), intent(in) :: a(3), b(3)
        real(8), intent(in) :: j
        real(8), intent(out) :: c(3)
        c(1) = (a(2)*b(3) - a(3)*b(2))/j
        c(2) = (a(3)*b(1) - a(1)*b(3))/j
        c(3) = (a(1)*b(2) - a(2)*b(1))/j
    end subroutine crossp
end module calc_mod