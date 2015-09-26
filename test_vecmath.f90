program test_vecmath

  use kinds
  use vecmath
  use test

  implicit none

  real(kind=d), parameter :: tolerance = 0.001_d

  type(vecmath_vector3d) :: zero_3d, &
    e_i_3d, e_j_3d, e_k_3d,                      &
    minus_e_i_3d, minus_e_j_3d, minus_e_k_3d,    &
    v

        zero_3d = vecmath_vector3d( 0.0_d,  0.0_d,  0.0_d)
         e_i_3d = vecmath_vector3d( 1.0_d,  0.0_d,  0.0_d)
         e_j_3d = vecmath_vector3d( 0.0_d,  1.0_d,  0.0_d)
         e_k_3d = vecmath_vector3d( 0.0_d,  0.0_d,  1.0_d)
   minus_e_i_3d = vecmath_vector3d(-1.0_d,  0.0_d,  0.0_d)
   minus_e_j_3d = vecmath_vector3d( 0.0_d, -1.0_d,  0.0_d)
   minus_e_k_3d = vecmath_vector3d( 0.0_d,  0.0_d, -1.0_d)
              v = vecmath_vector3d( 2.0_d,  3.0_d,  6.0_d)

  print *, "Testing..."

  print *, ""
  print *, "vector3d_angle()"
  call test_vector3d_angle()

  print *, ""
  print *, "vector3d_is_approximately_equal()"
  call test_vector3d_is_approximately_equal()

  print *, ""
  print *, "vector3d_norm()"
  call test_vector3d_norm()

  print *, ""
  print *, "vector3d_plus_vector3d()"
  call test_vector3d_plus_vector3d()
  
contains

  subroutine test_vector3d_angle()
    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      zero_3d%angle(e_i_3d), pi_d / 2.0_d, tolerance, "zero e_i"                                   &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_i_3d%angle(zero_3d), pi_d / 2.0_d, tolerance, "e_i zero"                                   &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      zero_3d%angle(zero_3d), pi_d / 2.0_d, tolerance, "zero zero"                                 &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_i_3d%angle(e_j_3d), pi_d / 2.0_d, tolerance, "e_i e_j"                                     &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_i_3d%angle(e_k_3d), pi_d / 2.0_d, tolerance, "e_i e_k"                                     &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_j_3d%angle(e_k_3d), pi_d / 2.0_d, tolerance, "e_j e_k"                                     &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_i_3d%angle(minus_e_i_3d), pi_d, tolerance, "e_i -e_i"                                      &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_j_3d%angle(minus_e_j_3d), pi_d, tolerance, "e_j -e_j"                                      &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_k_3d%angle(minus_e_k_3d), pi_d, tolerance, "e_k -e_k"                                      &
    )

  end subroutine test_vector3d_angle

  subroutine test_vector3d_is_approximately_equal()
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      zero_3d%is_approximately_equal(zero_3d, tolerance),                                          &
      .true.,                                                                                      &
      "zero"                                                                                       &
    )
    
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      e_i_3d%is_approximately_equal(e_i_3d, tolerance),                                            &
      .true.,                                                                                      &
      "x"                                                                                          &
    )
    
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      e_j_3d%is_approximately_equal(e_j_3d, tolerance),                                            &
      .true.,                                                                                      &
      "y"                                                                                          &
    )
    
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      e_k_3d%is_approximately_equal(e_k_3d, tolerance),                                            &
      .true.,                                                                                      &
      "z"                                                                                          &
    )

    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      v%is_approximately_equal(v, tolerance),                                                      &
      .true.,                                                                                      &
      "xyz"                                                                                        &
    )

    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      zero_3d%is_approximately_equal(e_i_3d, tolerance),                                           &
      .false.,                                                                                     &
      "false"                                                                                      &
    )
  end subroutine test_vector3d_is_approximately_equal


  subroutine test_vector3d_norm()
    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      zero_3d%norm(), 0.0_d, tolerance, "zero"                                                     &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_i_3d%norm(), 1.0_d, tolerance, "x"                                                         &
    )

    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_j_3d%norm(), 1.0_d, tolerance, "y"                                                         &
    )
    
    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      e_k_3d%norm(), 1.0_d, tolerance, "z"                                                         &
    )
    
    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      minus_e_i_3d%norm(), 1.0_d, tolerance, "-x"                                                  &
    )
    
    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      minus_e_j_3d%norm(), 1.0_d, tolerance, "-y"                                                  &
    )
    
    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      minus_e_k_3d%norm(), 1.0_d, tolerance, "-z"                                                  &
    )
    
    call test_mod_assert_equal_real                                                                &
    (                                                                                              &
      v%norm(), 7.0_d, tolerance, "xyz"                                                            &
    )
  end subroutine test_vector3d_norm


  subroutine test_vector3d_plus_vector3d()
    type(vecmath_vector3d) :: result

    result = zero_3d + zero_3d
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(zero_3d, tolerance),                                           &
      .true.,                                                                                      &
      "zero"                                                                                       &
    )

    result = e_i_3d + e_i_3d
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(vecmath_vector3d(2.0_d, 0.0_d, 0.0_d), tolerance),             &
      .true.,                                                                                      &
      "x"                                                                                          &
    )

    result = e_j_3d + e_j_3d
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(vecmath_vector3d(0.0_d, 2.0_d, 0.0_d), tolerance),             &
      .true.,                                                                                      &
      "y"                                                                                          &
    )

    result = e_k_3d + e_k_3d
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(vecmath_vector3d(0.0_d, 0.0_d, 2.0_d), tolerance),             &
      .true.,                                                                                      &
      "z"                                                                                          &
    )

    result = e_i_3d + e_j_3d
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(vecmath_vector3d(1.0_d, 1.0_d, 0.0_d), tolerance),             &
      .true.,                                                                                      &
      "xy"                                                                                         &
    )

    result = e_i_3d + e_k_3d
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(vecmath_vector3d(1.0_d, 0.0_d, 1.0_d), tolerance),             &
      .true.,                                                                                      &
      "xz"                                                                                         &
    )

    result = e_j_3d + e_k_3d
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(vecmath_vector3d(0.0_d, 1.0_d, 1.0_d), tolerance),             &
      .true.,                                                                                      &
      "yz"                                                                                         &
    )

    result = v + v
    call test_mod_assert_equal                                                                     &
    (                                                                                              &
      result%is_approximately_equal(vecmath_vector3d(4.0_d,  6.0_d, 12.0_d), tolerance),           &
      .true.,                                                                                      &
      "xyz"                                                                                        &
    )
  end subroutine test_vector3d_plus_vector3d

end program test_vecmath
