module vecmath

  use kinds

  implicit none

  private

  public :: vecmath_vector3d, &
    operator(+), operator(-), & 
    operator(*), operator(/), &
    operator(.dot.), operator(.cross.)

  type :: vecmath_vector3d
    real(kind=d) :: x, y, z
  contains
    procedure :: angle                   => vecmath_vector3d_angle
    procedure :: is_approximately_equal  => vecmath_vector3d_is_approximately_equal
    procedure :: is_approximately_zero   => vecmath_vector3d_is_approximately_zero
    procedure :: norm                    => vecmath_vector3d_norm
    procedure :: normalized              => vecmath_vector3d_normalized
  end type  vecmath_vector3d

  interface operator(+)
    module procedure vecmath_vector3d_plus_vecmath_vector3d
  end interface

  interface operator(-)
    module procedure vecmath_vector3d_minus_vecmath_vector3d
    module procedure vecmath_vector3d_negative
  end interface

  interface operator(*)
    module procedure vecmath_vector3d_times_real
  end interface

  interface operator(/)
    module procedure vecmath_vector3d_divided_by_real
  end interface

  interface operator(.dot.)
    module procedure vecmath_vector3d_dot_vecmath_vector3d
  end interface

  interface operator(.cross.)
    module procedure vecmath_vector3d_cross_vecmath_vector3d
  end interface

  type(vecmath_vector3d), parameter :: &
    vecmath_vector3d_zero = vecmath_vector3d(0.0_d, 0.0_d, 0.0_d)

contains

  ! Angle (in radians) between this vector and the other vector.
  !
  ! cos(phi) = (v_ji . v_jk) / (  )
  pure function vecmath_vector3d_angle(this, other) result(f)
    class(vecmath_vector3d), intent(in) :: this
    type(vecmath_vector3d), intent(in) :: other
    real(kind=d) :: f

    real(kind=d) :: cos_phi

    ! avoid division by zero for the zero vector 
    if (this%is_approximately_zero() .or. other%is_approximately_zero()) then
      cos_phi = 0.0_d
    else
      cos_phi = (this .dot. other) / (this%norm() * other%norm())
      if (cos_phi < -1.0_d) then
        cos_phi = -1.0_d
      else if (cos_phi > 1.0_d) then
        cos_phi = 1.0_d
      end if
    end if

    f = acos(cos_phi)
  end function vecmath_vector3d_angle

  pure function vecmath_vector3d_is_approximately_equal(this, other, tolerance) result(f)
    class(vecmath_vector3d), intent(in) :: this
    type(vecmath_vector3d), intent(in) :: other
    real(kind=d), intent(in), optional :: tolerance
    logical :: f

    type(vecmath_vector3d) :: diff

    real(kind=d) :: min_norm, tol

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = eps_d
    end if

    diff = this - other
    
    if (this%norm() <= other%norm()) then
      min_norm = this%norm()
    else
      min_norm = other%norm()
    end if

    f = diff%norm() <= tol * min_norm
  end function vecmath_vector3d_is_approximately_equal

  pure function vecmath_vector3d_is_approximately_zero(this, tolerance) result(f)
    class(vecmath_vector3d), intent(in) :: this
    real(kind=d), intent(in), optional :: tolerance
    logical :: f

    real(kind=d) :: tol

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = eps_d
    end if

    f = this%is_approximately_equal(vecmath_vector3d_zero, tol)
  end function vecmath_vector3d_is_approximately_zero

  pure function vecmath_vector3d_norm(this) result(f)
    class(vecmath_vector3d), intent(in) :: this
    real(kind=d) :: f

    f = sqrt(this .dot. this)
  end function vecmath_vector3d_norm

  pure function vecmath_vector3d_normalized(this) result(f)
    class(vecmath_vector3d), intent(in) :: this
    type(vecmath_vector3d) :: f

    f = this / this%norm()
  end function vecmath_vector3d_normalized

  pure function vecmath_vector3d_plus_vecmath_vector3d(v1, v2) result(f)
    type(vecmath_vector3d), intent(in) :: v1, v2
    type(vecmath_vector3d) :: f

    f%x = v1%x + v2%x
    f%y = v1%y + v2%y
    f%z = v1%z + v2%z
  end function vecmath_vector3d_plus_vecmath_vector3d

  pure function vecmath_vector3d_minus_vecmath_vector3d(v1, v2) result(f)
    type(vecmath_vector3d), intent(in) :: v1, v2
    type(vecmath_vector3d) :: f

    f%x = v1%x - v2%x
    f%y = v1%y - v2%y
    f%z = v1%z - v2%z
  end function vecmath_vector3d_minus_vecmath_vector3d

  pure function vecmath_vector3d_negative(v) result(f)
    type(vecmath_vector3d), intent(in) :: v
    type(vecmath_vector3d) :: f

    f%x = -v%x
    f%y = -v%y
    f%z = -v%z
  end function vecmath_vector3d_negative

  pure function vecmath_vector3d_times_real(v, r) result(f)
    type(vecmath_vector3d), intent(in) :: v
    real(kind=d), intent(in) :: r
    type(vecmath_vector3d) :: f

    f%x = v%x * r
    f%y = v%y * r
    f%z = v%z * r
  end function vecmath_vector3d_times_real

  pure function vecmath_vector3d_divided_by_real(v, r) result(f)
    type(vecmath_vector3d), intent(in) :: v
    real(kind=d), intent(in) :: r
    type(vecmath_vector3d) :: f

    f = v * (1.0_d / r)
  end function vecmath_vector3d_divided_by_real

  pure function vecmath_vector3d_dot_vecmath_vector3d(v1, v2) result(f)
    type(vecmath_vector3d), intent(in) :: v1, v2
    real(kind=d) :: f

    f = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
  end function vecmath_vector3d_dot_vecmath_vector3d

  pure function vecmath_vector3d_cross_vecmath_vector3d(v1, v2) result(f)
    type(vecmath_vector3d), intent(in) :: v1, v2
    type(vecmath_vector3d) :: f

    f%x = v1%y * v2%z - v1%z * v2%y
    f%y = v1%z * v2%x - v1%x * v2%z
    f%z = v1%x * v2%y - v1%y * v2%x
  end function vecmath_vector3d_cross_vecmath_vector3d

end module vecmath
