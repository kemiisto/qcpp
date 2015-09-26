! From Knuth D.E. The art of computer programming (vol II).
! |a - b| <= e * |a| and |a - b| <= e * |b|
! or 
! |a - b| <= e * min(|a|, |b|)
subroutine SPECIFIC_PROCEDURE(actual, expected, tolerance, description)
  real(kind=REALKIND), intent(in) :: actual
  real(kind=REALKIND), intent(in) :: expected
  real(kind=REALKIND), intent(in) :: tolerance
  character(len=*), intent(in) :: description

  real(kind=REALKIND) :: diff_abs, min_abs

  diff_abs = abs(actual - expected)
  
  if (abs(actual) <= abs(expected)) then
    min_abs = abs(actual)
  else
    min_abs = abs(expected)
  end if

  if (diff_abs <= tolerance * min_abs) then
    call print_passed_message(description)
  else
    call print_failed_message(description)
    print *, "     actual:", actual
    print *, "   expected:", expected  
  end if
end subroutine SPECIFIC_PROCEDURE
