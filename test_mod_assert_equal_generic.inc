elemental impure subroutine SPECIFIC_PROCEDURE(actual, expected, description)
  SPECIFIC_TYPE, intent(in) :: actual
  SPECIFIC_TYPE, intent(in) :: expected
  character(len=*), intent(in) :: description

  if (SPECIFIC_OPERATION) then
    call print_passed_message(description)
  else
    call print_failed_message(description)
    print *, "     actual:", actual
    print *, "   expected:", expected  
  end if
end subroutine SPECIFIC_PROCEDURE