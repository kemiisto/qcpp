program project_02

  use fcl_constants
  use fcl_kinds
  use fcl_lapack
  use fcl_vecmath_vector_3d
  use chem

  implicit none

  character(len=79), parameter :: separator = repeat("-", 79)
  integer, parameter :: inp_file_unit = 1
  integer, parameter :: out_file_unit = 2

  character(len=256) :: geometry_inp_file_name
  character(len=256) :: hessian_inp_file_name
  character(len=256) :: out_file_name

  type(chem_mod_molecule) :: molecule
  
  type(fcl_vecmath_mod_vector3d) :: position
  real(kind=d), dimension(:,:), allocatable :: hessian
  real(kind=d), dimension(:), allocatable :: hessian_eigenvalues, frequencies

  if (command_argument_count() /= 3) then
    print *, "Provide geometry input, Hessian input and output file names."
    stop
  end if

  call get_command_argument(1, geometry_inp_file_name)
  call get_command_argument(2, hessian_inp_file_name)
  call get_command_argument(3, out_file_name)

  print *, separator
  print *, "Geometry input file: ", trim(geometry_inp_file_name)
  print *, " Hessian input file: ", trim(hessian_inp_file_name)
  print *, "        Output file: ", trim(out_file_name)
  print *, separator

  open(unit=inp_file_unit, file=geometry_inp_file_name, action="read")
  call read_geometry()
  close(unit=inp_file_unit)

  allocate (hessian(molecule%number_of_atoms() * 3, molecule%number_of_atoms() * 3))
  allocate (hessian_eigenvalues(molecule%number_of_atoms() * 3))
  allocate (frequencies(molecule%number_of_atoms() * 3))

  call read_hessian()
  call mass_weight_hessian()

  print *, "Finding Hessian eigenvalues..."
  call fcl_lapack_dsyev(hessian, hessian_eigenvalues, .false.)

  open(unit=out_file_unit, file=out_file_name, action="write")
  call write_hessian_eigenvalues()
  
  print *, "Calculating frequencies..."
  ! Converting eigenvalues to SI units
  hessian_eigenvalues = hessian_eigenvalues * hartree_energy / unified_atomic_mass_unit / bohr_radius ** 2
  ! Calculating frequencies
  frequencies = sqrt(abs(hessian_eigenvalues)) / (2 * pi * speed_of_light_in_vacuum) * centi
  where (hessian_eigenvalues < 0) frequencies = -frequencies
  call write_frequencies()

  close(unit=out_file_unit)
  print *, separator

  deallocate (hessian)
  deallocate (hessian_eigenvalues)
  deallocate (frequencies)

contains

  subroutine read_geometry()
    integer :: number_of_atoms, i
    real(kind=d) :: atomic_number
    type(chem_mod_atom), pointer :: atom

    print *, "Reading geometry..."
    read (inp_file_unit,*) number_of_atoms
    call molecule%set_number_of_atoms(number_of_atoms)
    
    do i = 1, number_of_atoms
      read (inp_file_unit,*) atomic_number, position
      atom => molecule%atom_pointer(i)
      atom%atomic_number = int(atomic_number)
      atom%position = position
    end do
  end subroutine read_geometry

  subroutine read_hessian()
    integer :: number_of_atoms, i, j

    open(unit=inp_file_unit, file=hessian_inp_file_name, action="read")
    print *, "Reading Hessian..."
    read (inp_file_unit,*) number_of_atoms

    if (number_of_atoms /= molecule%number_of_atoms()) then
      print *, "Different number of atoms in geometry & Hessian inputs!"
      stop
    end if

    do i = 1, number_of_atoms * 3
      do j = 1, number_of_atoms * 3, 3
        read (inp_file_unit, *) hessian(i, j:j+2)
      end do
    end do

    close(unit=1)
  end subroutine read_hessian

  subroutine mass_weight_hessian()
    type(chem_mod_atom), pointer :: atom_i, atom_j
    real(kind=d) :: mass_i, mass_j

    integer :: i, j, k, l

    print *, "Mass-weighting the Hessian..."
    do i = 1, molecule%number_of_atoms()
      atom_i => molecule%atom_pointer(i)
      mass_i = chem_mod_atomic_masses(atom_i%atomic_number)
      do j = 1, molecule%number_of_atoms()
        atom_j => molecule%atom_pointer(j)
        mass_j = chem_mod_atomic_masses(atom_j%atomic_number)
        k = 1 + (i - 1) * 3
        l = 1 + (j - 1) * 3
        hessian(k:k+2, l:l+2) = hessian(k:k+2, l:l+2) / sqrt(mass_i * mass_j)
      end do
    end do

    close(unit=1)
  end subroutine mass_weight_hessian

  subroutine write_hessian_eigenvalues()
    integer :: i

    print *, "Writing eigenvalues..."
    write (out_file_unit, "(a)") "Hessian eigenvalues (hartree/amu-bohr^2):"
    do i = molecule%number_of_atoms() * 3, 1, -1 
      write (out_file_unit, "(i4,f21.10)") i - 1, hessian_eigenvalues(i)
    end do
    write (out_file_unit, "(bn)")
  end subroutine write_hessian_eigenvalues

  subroutine write_frequencies()
    integer :: i

    print *, "Writing frequencies..."
    write (out_file_unit, "(a)") "Harmonic vibrational frequencies (cm^-1):"
    do i = molecule%number_of_atoms() * 3, 1, -1 
      write (out_file_unit, "(i4, f10.4)") i - 1, frequencies(i)
    end do
  end subroutine write_frequencies

end program project_02