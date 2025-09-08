module spline2_mod
  implicit none
contains

  subroutine spline_getter(dx, y, a_temp,b_temp,c_temp,d_temp)

    real(8), intent(in) :: dx, y(:)
    real(8) :: fval

    integer :: n, i, j
    real(8) :: h
    real(8), allocatable :: alpha(:), l(:), mu(:), z(:)
    real(8), intent(out) :: a_temp(:), b_temp(:), c_temp(:), d_temp(:)
    real(8), allocatable :: a(:), b(:), c(:), d(:)

    n = size(y)

    ! 등간격 간격 h
    h = dx
    allocate(alpha(n-1))
    allocate(l(n), mu(n), z(n))
    allocate(a(n), b(n-1), c(n), d(n-1))

    a = y

    ! alpha 계산
    do i = 2, n-1
       alpha(i) = (3.0d0/h) * (a(i+1)-2.0d0*a(i)+a(i-1))
    end do

    ! Thomas 알고리즘
    l(1) = 1.0d0
    mu(1) = 0.0d0
    z(1) = 0.0d0

    do i = 2, n-1
       l(i) = 4.0d0*h - h*mu(i-1)
       mu(i) = h/l(i)
       z(i) = (alpha(i) - h*z(i-1))/l(i)
    end do

    l(n) = 1.0d0
    z(n) = 0.0d0
    c(n) = 0.0d0

    do j = n-1, 1, -1
       c(j) = z(j) - mu(j)*c(j+1)
       b(j) = (a(j+1)-a(j))/h - h*(c(j+1)+2.0d0*c(j))/3.0d0
       d(j) = (c(j+1)-c(j))/(3.0d0*h)
    end do
    a_temp = a(1:n-1)
    b_temp = b
    c_temp = c(1:n-1)
    d_temp = d
    deallocate(a, b, c, d)


    deallocate(alpha, l, mu, z)

  end subroutine spline_getter
end module spline2_mod
