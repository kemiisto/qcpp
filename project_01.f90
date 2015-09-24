program project_01

  use kinds
  use vecmath
  use chem

  implicit none

  character(len=256) :: inp_file_name
  character(len=256) :: out_file_name
  type(chem_mod_molecule) :: molecule
  type(chem_mod_atom), pointer :: atom
  integer :: number_of_atoms
  integer :: i, j, k, l
  integer :: atomic_number
  type(vecmath_vector3d) :: position

  if (command_argument_count() /= 2) then
    print *, "Provide input and output file names."
    stop
  end if

  call get_command_argument(1, inp_file_name)
  call get_command_argument(2, out_file_name)

  print *, " Input file: ", trim(inp_file_name)
  print *, "Output file: ", trim(out_file_name)

  open(unit=1, file=inp_file_name, action="read")
  print *, "Reading input file..."
  read (1,*) number_of_atoms
  call molecule%set_number_of_atoms(7)
  
  do i = 1, number_of_atoms
    read (1,*) atomic_number, position
    atom => molecule%atom(i)
    atom%atomic_number = atomic_number
    atom%position = position
    !molecule%atom(i)%atomic_number = atomic_number
    !molecule%atom(i)%position = position
  end do
  close(unit=1)
  print *, "Done!"

  open(unit=2, file=out_file_name, action="write")

  print *, "Doing calculations..."

  write (2, '(a,i2)') "Number of atoms:", molecule%number_of_atoms()
  write (2, '(a)') "Input Cartesian coordinates:"
  do i = 1, molecule%number_of_atoms()
    atom => molecule%atom(i)
    write (2, '(i3,3f21.12)')          &
      atom%atomic_number, &
      atom%position%x,    &
      atom%position%y,    &
      atom%position%z
  end do

  write (2, '(a)') "Interatomic distances (bohr):"
  do i = 2, molecule%number_of_atoms()
    do j = 1, i - 1
      write (2, '(2i2,f9.5)') i - 1, j - 1, molecule%distance(i, j)
    end do
  end do
  write (2, '(bn)')
  
  write (2, '(a)') "Bond angles:"
  do i = 1, molecule%number_of_atoms()
    do j = 1, i - 1
      do k = 1, j - 1
        if (molecule%distance(i, j) < 4.0_d .and. molecule%distance(j, k) < 4.0_d) then
          write (2, '(i2,a,i2,a,i2,f11.6)') i - 1, "-", j - 1, "-", k - 1, &
            molecule%angle(i, j, k) * (180.0_d / acos(-1.0_d))
        end if
      end do
    end do
  end do
  write (2, '(bn)')

  write (2, '(a)') "Out-of-plane angles:"
  do i = 1, molecule%number_of_atoms()
    do k = 1, molecule%number_of_atoms()
      do j = 1, molecule%number_of_atoms()
        do l = 1, j - 1
          ! i /= j /= k /= l
          ! can be written as
          ! i /= j .and. i /= k .and. i /= l .and.
          ! j /= k .and. j /= l .and.
          ! k /= l
          ! besides j /= l is already guaranteed by the loop do l = 1, j - 1
          if (i /= j .and. i /= k .and. i /= l .and. j /= k .and. k /= l                           &
            .and. molecule%distance(i, k) < 4.0_d                                                  &
              .and. molecule%distance(j, k) < 4.0_d                                                &
                .and. molecule%distance(l, k) < 4.0_d) then
                  write (2, '(i2,a,i2,a,i2,a,i2,f11.6)')       &
                    i - 1, "-", j - 1, "-", k - 1, "-", l - 1, &
                    molecule%out_of_plane_angle(i, j, k, l) * (180.0_d / acos(-1.0_d))

          end if
        end do 
      end do 
    end do 
  end do

  write (2, '(a)') "Torsional angles:"
  do i = 1, molecule%number_of_atoms()
    do j = 1, i - 1
      do k = 1, j - 1
        do l = 1, k - 1
          if (molecule%distance(i, j) < 4.0_d &
                .and. molecule%distance(j, k) < 4.0_d &
                  .and. molecule%distance(k, l) < 4.0_d) then
            write (2, '(i2,a,i2,a,i2,a,i2,f11.6)')       &
              i - 1, "-", j - 1, "-", k - 1, "-", l - 1, &
            molecule%dihedral_angle(i, j, k, l) * (180.0_d / acos(-1.0_d))
          end if
        end do
      end do
    end do
  end do

  close(unit=2)

  print *, "Done!"

end program project_01
