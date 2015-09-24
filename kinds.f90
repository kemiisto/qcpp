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

  ! double-precision floating-point format:
  !   "double" in ieee 754-1985, "binary64" in ieee 754-2008;
  !   occupies 8 bytes (64 bits) in computer memory;
  !   15-17 significant decimal digits precision;
  !   308 approximate exponent range.
  ! integer, parameter :: d = ieee_selected_real_kind(15, 307)
  integer, parameter :: d = selected_real_kind(15, 307)
  
end module kinds