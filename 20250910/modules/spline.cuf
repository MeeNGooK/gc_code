module cubic_spline_mod
  implicit none
contains

  function spline_deriv(x, y, xval) result(dy)
    !-----------------------------------------------------------
    ! 입력:
    !   x(:)    - 등차수열 데이터 x값들
    !   y(:)    - 주어진 데이터 y값들
    !   xval    - 보간할 점 (x 범위 내에 있어야 함)
    !
    ! 출력:
    !   dy      - xval에서 cubic spline 보간의 도함수 값
    !-----------------------------------------------------------

    real(8), intent(in) :: x(:), y(:), xval
    real(8) :: dy

    integer :: n, i, j
    real(8) :: h
    real(8), allocatable :: alpha(:), l(:), mu(:), z(:)
    real(8), allocatable :: a(:), b(:), c(:), d(:)

    n = size(x)
    if (size(y) /= n) stop "Error: x and y size mismatch"

    ! 등간격 간격 h
    h = x(2) - x(1)

    allocate(alpha(n-1))
    allocate(l(n), mu(n), z(n))
    allocate(a(n), b(n-1), c(n), d(n-1))

    a = y

    ! alpha 계산 (등간격 h를 사용)
    do i = 2, n-1
       alpha(i) = (3.0d0/h) * (a(i+1)-2.0d0*a(i)+a(i-1))
    end do

    ! Thomas 알고리즘
    l(1) = 1.0d0
    mu(1) = 0.0d0
    z(1) = 0.0d0

    do i = 2, n-1
       l(i) = 2.0d0*(2.0d0*h) - h*mu(i-1)
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

    i = int( (xval - x(1))/h ) + 1
    if (i < 1) i = 1
    if (i > n-1) i = n-1

    ! 도함수
    dy = b(i) + 2.0d0*c(i)*(xval-x(i)) + 3.0d0*d(i)*(xval-x(i))**2

    deallocate(alpha, l, mu, z, a, b, c, d)

  end function spline_deriv
  function spline_value(x, y, xval) result(fval)
    !-----------------------------------------------------------
    ! 입력:
    !   x(:)    - 등차수열 데이터 x값들
    !   y(:)    - 주어진 데이터 y값들
    !   xval    - 보간할 점 (x 범위 내에 있어야 함)
    !
    ! 출력:
    !   fval    - xval에서 cubic spline 보간 함수값
    !-----------------------------------------------------------

    real(8), intent(in) :: x(:), y(:), xval
    real(8) :: fval

    integer :: n, i, j
    real(8) :: h, dx
    real(8), allocatable :: alpha(:), l(:), mu(:), z(:)
    real(8), allocatable :: a(:), b(:), c(:), d(:)

    n = size(x)
    if (size(y) /= n) stop "Error: x and y size mismatch"

    ! 등간격 간격 h
    h = x(2) - x(1)

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

    ! xval이 속하는 구간 index
    i = int( (xval - x(1))/h ) + 1
    if (i < 1) i = 1
    if (i > n-1) i = n-1

    dx = xval - x(i)

    ! 함수값 계산
    fval = a(i) + b(i)*dx + c(i)*dx**2 + d(i)*dx**3

    deallocate(alpha, l, mu, z, a, b, c, d)

  end function spline_value
end module cubic_spline_mod
