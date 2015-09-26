module test

  use kinds

  implicit none

  private

  public :: test_mod_assert_equal, test_mod_assert_equal_real

  interface test_mod_assert_equal
    module procedure test_mod_assert_equal_integer
    module procedure test_mod_assert_equal_logical
  end interface test_mod_assert_equal

  interface test_mod_assert_equal_real
    module procedure test_mod_assert_equal_real_s
    module procedure test_mod_assert_equal_real_d
  end interface test_mod_assert_equal_real

contains

! ------------------------------------------------------------------------------

#define SPECIFIC_OPERATION actual == expected

#define SPECIFIC_PROCEDURE test_mod_assert_equal_integer
#define SPECIFIC_TYPE integer
#include "test_mod_assert_equal_generic.inc"
#undef SPECIFIC_PROCEDURE
#undef SPECIFIC_TYPE

#undef SPECIFIC_OPERATION

! ------------------------------------------------------------------------------

#define SPECIFIC_OPERATION actual .eqv. expected

#define SPECIFIC_PROCEDURE test_mod_assert_equal_logical
#define SPECIFIC_TYPE logical
#include "test_mod_assert_equal_generic.inc"
#undef SPECIFIC_PROCEDURE
#undef SPECIFIC_TYPE

#undef SPECIFIC_OPERATION

! ------------------------------------------------------------------------------

#define SPECIFIC_PROCEDURE test_mod_assert_equal_real_s
#define REALKIND s
#include "test_mod_assert_equal_real_generic.inc"
#undef SPECIFIC_PROCEDURE
#undef REALKIND

#define SPECIFIC_PROCEDURE test_mod_assert_equal_real_d
#define REALKIND d
#include "test_mod_assert_equal_real_generic.inc"
#undef SPECIFIC_PROCEDURE
#undef REALKIND

  subroutine print_passed_message(description)
    character(len=*), intent(in) :: description

    print '(a, a10, a)', "|- test ", description, " passed."
  end subroutine print_passed_message

  subroutine print_failed_message(description)
    character(len=*), intent(in) :: description

    print '(a, a10, a)', "|- test ", description, " FAILED!"
  end subroutine print_failed_message

end module test
