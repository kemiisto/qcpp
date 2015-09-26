module kinds

  ! Intel implementation of ieee_selected_real_kind() is buggy.
  ! https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/293623
  !
  ! use, intrinsic :: ieee_arithmetic

  implicit none

  ! single-precision floating-point format:
  !   "single" in ieee 754-1985, "binary32" in ieee 754-2008; 
  !   occupies 4 bytes (32 bits) in computer memory;
  !   6 to 9 significant decimal digits precision;
  !   38 approximate exponent range.
  ! integer, parameter :: s = ieee_selected_real_kind(6, 37)
  integer, parameter :: s = selected_real_kind(6, 37)
  real(kind=s), parameter :: eps_s = epsilon(1.0_s)
  real(kind=s), parameter :: pi_s = 4.0_s * atan(1.0_s)

  ! double-precision floating-point format:
  !   "double" in ieee 754-1985, "binary64" in ieee 754-2008;
  !   occupies 8 bytes (64 bits) in computer memory;
  !   15-17 significant decimal digits precision;
  !   308 approximate exponent range.
  ! integer, parameter :: d = ieee_selected_real_kind(15, 307)
  integer, parameter :: d = selected_real_kind(15, 307)
  real(kind=d), parameter :: eps_d = epsilon(1.0_d)
  real(kind=d), parameter :: pi_d = 4.0_d * atan(1.0_d)
  
end module kinds