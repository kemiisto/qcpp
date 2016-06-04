program project_01

  use fcl_constants
  use fcl_kinds
  use fcl_lapack
  use fcl_vecmath_vector_3d
  use chem

  implicit none

  character(len=79), parameter :: separator = repeat("-", 79)
  character(len=1), parameter :: tab = char(9)
  integer, parameter :: out_file_unit = 2

  real(kind=dp), parameter :: m_to_angstrom = 1.0d+10

  real(kind=dp), parameter :: bohr_to_m = bohr_radius
  real(kind=dp), parameter :: bohr_to_cm = bohr_to_m / centi 
  real(kind=dp), parameter :: bohr_to_angstrom = bohr_to_m * m_to_angstrom

  real(kind=dp), parameter :: amu_to_kg = unified_atomic_mass_unit
  real(kind=dp), parameter :: amu_to_g = amu_to_kg * kilo
  
  character(len=256) :: inp_file_name
  character(len=256) :: out_file_name
  type(chem_mod_molecule) :: molecule
  type(chem_mod_atom), pointer :: atom
  type(fcl_vecmath_mod_vector3d) :: position
  real(kind=dp), dimension(3, 3) :: moment_of_inertia_tensor
  real(kind=dp), dimension(3) :: principal_moments_of_inertia
  real(kind=dp), dimension(3) :: rotational_constants

  call process_arguments()

  call chem_mod_read_molecule_from_file(molecule, inp_file_name)

  open(unit=out_file_unit, file=out_file_name, action="write")

  call write_geometry()
  call write_distances()
  call write_angles()
  call write_out_of_plane_angles()
  call write_torsional_angles()

  call translate_molecule_to_center_of_mass()
  call calculate_moment_of_inertia_tensor()
  call write_moment_of_inertia_tensor()
  call calculate_principal_moments_of_inertia()
  call write_principal_moments_of_inertia()
  
  call classify_rotor()
  call calculate_rotational_constants()
  call write_rotational_constants()

  close(unit=out_file_unit)

contains

  subroutine write_geometry()
    integer :: i

    print *, "Writing geometry..."
    write (out_file_unit, "(a,i2)") "Number of atoms:", molecule%number_of_atoms()
    write (out_file_unit, "(a)") "Input Cartesian coordinates:"

    do i = 1, molecule%number_of_atoms()
      atom => molecule%atom_pointer(i)
      write (out_file_unit, "(i2,3f21.12)")          &
        atom%atomic_number,   &
        atom%position%x(),    &
        atom%position%y(),    &
        atom%position%z()
    end do
  end subroutine write_geometry

  subroutine write_distances()
    integer :: i, j

    print *, "Writing distances..."

    write (out_file_unit, "(a)") "Interatomic distances (bohr):"
    do i = 2, molecule%number_of_atoms()
      do j = 1, i - 1
        write (out_file_unit, "(2i2,f9.5)") i - 1, j - 1, molecule%distance(i, j)
      end do
    end do
    write (out_file_unit, "(bn)")
  end subroutine write_distances

  subroutine write_angles()
    integer :: i, j, k

    print *, "Writing angles..."

    write (out_file_unit, "(a)") "Bond angles:"
    do i = 1, molecule%number_of_atoms()
      do j = i + 1, molecule%number_of_atoms()
        do k = j + 1, molecule%number_of_atoms()
          if (molecule%distance(i, j) < 4.0_dp .and. molecule%distance(j, k) < 4.0_dp) then
            write (out_file_unit, "(i2,a,i2,a,i2,f11.6)") i - 1, "-", j - 1, "-", k - 1, &
              molecule%angle(i, j, k) * (180.0_dp / acos(-1.0_dp))
          end if
        end do
      end do
    end do
    write (out_file_unit, "(bn)")
  end subroutine write_angles

  subroutine write_out_of_plane_angles()
    integer :: i, j, k, l

    print *, "Writing out-of-plane angles..."

    write (out_file_unit, "(a)") "Out-of-plane angles:"
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
            if (i /= j .and. i /= k .and. i /= l .and. j /= k .and. k /= l   &
              .and. molecule%distance(i, k) < 4.0_dp                         &
                .and. molecule%distance(j, k) < 4.0_dp                       &
                  .and. molecule%distance(l, k) < 4.0_dp) then
                    write (out_file_unit, "(i2,a,i2,a,i2,a,i2,f11.6)")       &
                      i - 1, "-", j - 1, "-", k - 1, "-", l - 1, &
                      molecule%out_of_plane_angle(i, j, k, l) * (180.0_dp / acos(-1.0_dp))
            end if
          end do 
        end do 
      end do 
    end do
    write (out_file_unit, "(bn)")
  end subroutine write_out_of_plane_angles

  subroutine write_torsional_angles()
    integer :: i, j, k, l

    print *, "Writing torsional angles..."

    write (out_file_unit, "(a)") "Torsional angles:"
    do i = 1, molecule%number_of_atoms()
      do j = 1, i - 1
        do k = 1, j - 1
          do l = 1, k - 1
            if (molecule%distance(i, j) < 4.0_dp &
                  .and. molecule%distance(j, k) < 4.0_dp &
                    .and. molecule%distance(k, l) < 4.0_dp) then
              write (out_file_unit, "(i2,a,i2,a,i2,a,i2,f11.6)")  &
                i - 1, "-", j - 1, "-", k - 1, "-", l - 1, &
              molecule%dihedral_angle(i, j, k, l) * (180.0_dp / acos(-1.0_dp))
            end if
          end do
        end do
      end do
    end do
    write (out_file_unit, "(bn)")
  end subroutine write_torsional_angles

  subroutine translate_molecule_to_center_of_mass()
    print *, "Translating molecule to center of mass..."

    position = molecule%center_of_mass()
    write (out_file_unit, "(a,3f13.8)") "Molecular center of mass:",  &
      position%x(), position%y(), position%z()
    write (out_file_unit, "(bn)")
    call molecule%translate(-position)
  end subroutine translate_molecule_to_center_of_mass

  subroutine write_moment_of_inertia_tensor()
    integer :: i

    print *, "Writing moment of inertia tensor..."
    
    write (out_file_unit, "(a)") "Moment of inertia tensor:"
    write (out_file_unit, "(bn)")
    write (out_file_unit, "(3i12)") 1, 2, 3
    write (out_file_unit, "(bn)")
    write (out_file_unit, "(i5,3f13.7)") (i, moment_of_inertia_tensor(i, :), i = 1, 3)
    write (out_file_unit, "(bn)")
  end subroutine write_moment_of_inertia_tensor

  subroutine write_principal_moments_of_inertia()
    print *, "Writing principal moments of inertia..."

    write (out_file_unit, "(a)") "Principal moments of inertia (amu * bohr^2):"
    write (out_file_unit, "(a,3f13.6)") tab, principal_moments_of_inertia
    write (out_file_unit, "(bn)")

    write (out_file_unit, "(a)") "Principal moments of inertia (amu * AA^2):"
    write (out_file_unit, "(a,3f13.6)") tab, principal_moments_of_inertia * bohr_to_angstrom ** 2
    write (out_file_unit, "(bn)")

    write (out_file_unit, "(a)") "Principal moments of inertia (g * cm^2):"
    write (out_file_unit, "(a,3es13.6)") tab, principal_moments_of_inertia * amu_to_g * bohr_to_cm ** 2
    write (out_file_unit, "(bn)")
  end subroutine write_principal_moments_of_inertia

  subroutine classify_rotor()
    associate (m => principal_moments_of_inertia)
      if (molecule%number_of_atoms() == 2) then
        write (out_file_unit, "(a)") "Molecule is diatomic."
      else if (m(1) < 1d-4) then
        write (out_file_unit, "(a)") "Molecule is linear."
      else if (abs(m(1) - m(2)) < 1d-4 .and. abs(m(2) - m(3)) < 1d-4) then
        write (out_file_unit, "(a)") "Molecule is a spherical top."
      else if (abs(m(1) - m(2)) < 1d-4 .and. abs(m(2) - m(3)) > 1d-4) then
        write (out_file_unit, "(a)") "Molecule is an oblate symmetric top."
      else if (abs(m(1) - m(2)) > 1d-4 .and. abs(m(2) - m(3)) < 1d-4) then
        write (out_file_unit, "(a)") "Molecule is a prolate symmetric top."
      else 
        write (out_file_unit, "(a)") "Molecule is an asymmetric top."
      end if
      write (out_file_unit, "(bn)")
    end associate
  end subroutine

  subroutine write_rotational_constants()
    write (out_file_unit, "(a)") "Rotational constants (MHz):"
    write (out_file_unit, "(3(a,a4,f10.3))") &
      tab, "A = ", rotational_constants(1) * speed_of_light_in_vacuum / mega, &
      tab, "B = ", rotational_constants(2) * speed_of_light_in_vacuum / mega, &
      tab, "C = ", rotational_constants(3) * speed_of_light_in_vacuum / mega
    write (out_file_unit, "(bn)")
    write (out_file_unit, "(a)") "Rotational constants (cm-1):"
    write (out_file_unit, "(3(a,a4,f6.4))") &
      tab, "A = ", rotational_constants(1) * centi, &
      tab, "B = ", rotational_constants(2) * centi, &
      tab, "C = ", rotational_constants(3) * centi
  end subroutine write_rotational_constants

  subroutine process_arguments()
    if (command_argument_count() /= 2) then
      print *, "Provide input and output file names."
      stop
    end if

    call get_command_argument(1, inp_file_name)
    call get_command_argument(2, out_file_name)

    print *, separator
    print *, " Input file: ", trim(inp_file_name)
    print *, "Output file: ", trim(out_file_name)
    print *, separator
  end subroutine

  subroutine calculate_moment_of_inertia_tensor()
    print *, "Calculating moment of inertia tensor..."
    moment_of_inertia_tensor = molecule%moment_of_inertia_tensor()
  end subroutine

  subroutine calculate_principal_moments_of_inertia()
    print *, "Calculating principal moments of inertia..."
    call fcl_lapack_dsyev(moment_of_inertia_tensor, principal_moments_of_inertia, .false.)
  end subroutine

  subroutine calculate_rotational_constants()
    print *, "Calculating rotational constants..."
    ! Convert moments of inertia to SI units
    principal_moments_of_inertia = principal_moments_of_inertia * amu_to_kg * bohr_to_m ** 2
    ! Calculate rotational constants in SI units
    rotational_constants = planck_constant / (8 * pi ** 2 * speed_of_light_in_vacuum * principal_moments_of_inertia)
  end subroutine

end program project_01
