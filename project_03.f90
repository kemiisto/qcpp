program project_03
  
  use fcl_kinds
  use fcl_lapack
  use chem
  use iso_fortran_env
  
  implicit none
  
  character(len=79), parameter :: separator = repeat("-", 79)
  
  integer, parameter :: inp_file_unit = 1
  integer, parameter :: out_file_unit = 2
  
  character(len=*), parameter :: molecule_file_name = "geom.dat"
  character(len=*), parameter :: nuclear_repulsion_energy_file_name = "enuc.dat"
  character(len=*), parameter :: overlap_integrals_file_name = "s.dat"
  character(len=*), parameter :: kinetic_energy_integrals_file_name = "t.dat"
  character(len=*), parameter :: nuclear_attraction_integrals_file_name = "v.dat"
  character(len=*), parameter :: two_electron_integrals_file_name = "eri.dat"
  
  character(len=256) :: inp_folder_name
  character(len=256) :: inp_file_name
  character(len=256) :: out_file_name

  type(chem_mod_molecule) :: molecule

  real(kind=d), dimension(:, :), allocatable :: overlap_integrals, &
    overlap_eigenvectors, &
    kinetic_energy_integrals, &
    nuclear_attraction_integrals, &
    core_hamiltonian, &
    symmetric_orthogonalization_matrix, &
    fock_matrix, &
    coefficients_matrix, &
    density_matrix
  real(kind=d), dimension(:), allocatable :: overlap_eigenvalues
  real(kind=d), dimension(:), allocatable :: orbital_energies
  real(kind=d), dimension(:), allocatable :: two_electron_integrals
  real(kind=d) :: nuclear_repulsion_energy
  integer :: basis_set_size, temp, i, j, m, number_of_occupied_orbitals
  
  if (command_argument_count() /= 2) then
    print *, "Provide input folder and output file name."
    stop
  end if
  
  call get_command_argument(1, inp_folder_name)
  call get_command_argument(2, out_file_name)

  print *, separator
  print *, "Input folder: ", trim(inp_folder_name)
  print *, " Output file: ", trim(out_file_name)
  print *, separator

  inp_file_name = trim(inp_folder_name)//"/"//molecule_file_name
  print *, "Reading molecule..."
  print *, "  file: "//trim(inp_file_name)
  call chem_mod_read_molecule_from_file(molecule, inp_file_name)

  call read_nuclear_repulsion_energy()  
  call determine_basis_set_size()
    
  allocate (overlap_integrals(basis_set_size, basis_set_size))
  allocate (overlap_eigenvectors(basis_set_size, basis_set_size))
  allocate (symmetric_orthogonalization_matrix(basis_set_size, basis_set_size))
  allocate (overlap_eigenvalues(basis_set_size))
  allocate (kinetic_energy_integrals(basis_set_size, basis_set_size))
  allocate (nuclear_attraction_integrals(basis_set_size, basis_set_size))
  allocate (core_hamiltonian(basis_set_size, basis_set_size))
  allocate (fock_matrix(basis_set_size, basis_set_size))
  allocate (coefficients_matrix(basis_set_size, basis_set_size))
  allocate (orbital_energies(basis_set_size))
  allocate (density_matrix(basis_set_size, basis_set_size))
  temp = (basis_set_size * (basis_set_size + 1)) / 2
  allocate (two_electron_integrals((temp * (temp + 1)) / 2))
  two_electron_integrals = 0.0_d

  call read_one_electron_integrals(overlap_integrals_file_name, overlap_integrals)
  call read_one_electron_integrals(kinetic_energy_integrals_file_name, kinetic_energy_integrals)
  call read_one_electron_integrals(nuclear_attraction_integrals_file_name, nuclear_attraction_integrals)

  print *, "Forming the core Hamiltonian..."
  core_hamiltonian = kinetic_energy_integrals + nuclear_attraction_integrals

  open(unit=out_file_unit, file=out_file_name, action="write")
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a,f20.10)") "Nuclear repulsion energy = ", nuclear_repulsion_energy
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Overlap Integrals:"
  call write_square_matrix(overlap_integrals)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Kinetic-Energy Integrals:"
  call write_square_matrix(kinetic_energy_integrals)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Nuclear Attraction Integrals:"
  call write_square_matrix(nuclear_attraction_integrals)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Core Hamiltonian:"
  call write_square_matrix(core_hamiltonian)
  
  call read_two_electron_integrals(two_electron_integrals)
  
  overlap_eigenvectors = overlap_integrals
  call fcl_lapack_dsyev(overlap_eigenvectors, overlap_eigenvalues, .true.)
  overlap_eigenvalues = overlap_eigenvalues ** (-1.0_d/2)
  overlap_integrals = 0.0_d
  do i = 1, basis_set_size
    overlap_integrals(i, i) = overlap_eigenvalues(i)
  end do
  symmetric_orthogonalization_matrix = matmul(overlap_eigenvectors, overlap_integrals)
  symmetric_orthogonalization_matrix = matmul(symmetric_orthogonalization_matrix, transpose(overlap_eigenvectors))
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	S^-1/2 Matrix:"
  call write_square_matrix(symmetric_orthogonalization_matrix)
  
  fock_matrix = matmul(transpose(symmetric_orthogonalization_matrix), core_hamiltonian)
  fock_matrix = matmul(fock_matrix, symmetric_orthogonalization_matrix)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Initial F' Matrix:"
  call write_square_matrix(fock_matrix)
  
  coefficients_matrix = fock_matrix
  call fcl_lapack_dsyev(coefficients_matrix, orbital_energies, .true.)
  coefficients_matrix = matmul(symmetric_orthogonalization_matrix, coefficients_matrix)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Initial C Matrix:"
  call write_square_matrix(coefficients_matrix)
  
  ! FIX: read geometry and calculate number of occupied MOs.
  number_of_occupied_orbitals = molecule%number_of_electrons() / 2
  density_matrix = 0.0_d
  do i = 1, basis_set_size
    do j = 1, basis_set_size
      do m = 1, number_of_occupied_orbitals
        density_matrix(i, j) = density_matrix(i, j) + coefficients_matrix(i, m) * coefficients_matrix(j, m)
      end do
    end do
  end do
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Initial Density Matrix:"
  call write_square_matrix(density_matrix)

  close(unit=out_file_unit)

  deallocate (overlap_integrals)
  deallocate (overlap_eigenvectors)
  deallocate (symmetric_orthogonalization_matrix)
  deallocate (kinetic_energy_integrals)
  deallocate (nuclear_attraction_integrals)
  deallocate (core_hamiltonian)
  deallocate (two_electron_integrals)
  deallocate (fock_matrix)
  deallocate (coefficients_matrix)
  deallocate (orbital_energies)
  deallocate (density_matrix)

contains
  
  subroutine read_nuclear_repulsion_energy()
    inp_file_name = trim(inp_folder_name)//"/"//nuclear_repulsion_energy_file_name
    print *, "Reading nuclear repulsion energy..."
    print *, "  file: "//trim(inp_file_name)
    open(unit=inp_file_unit, file=trim(inp_file_name), action="read")
    read (inp_file_unit,*) nuclear_repulsion_energy
    close(unit=inp_file_unit)
    print *, "  nuclear repulsion energy:", nuclear_repulsion_energy
  end subroutine read_nuclear_repulsion_energy
  
  subroutine determine_basis_set_size()
    integer :: stat, i
    
    inp_file_name = trim(inp_folder_name)//"/"//overlap_integrals_file_name
    print *, "Determining the basis set size... "
    print *, "  file: "//trim(inp_file_name)
    open(unit=inp_file_unit, file=trim(inp_file_name), action="read")
    do
      read (inp_file_unit,*,iostat=stat) i
      if (stat == iostat_end) exit
    end do
    basis_set_size = i
    close(unit=inp_file_unit)
    print *, "  basis set size:", basis_set_size
  end subroutine determine_basis_set_size
  
  subroutine read_one_electron_integrals(integrals_file_name, integrals)
    character(len=*), intent(in) :: integrals_file_name
    real(kind=d), dimension(:,:), intent(out) :: integrals
    
    integer :: n, line, i, j
    real(kind=d) :: integral
    
    inp_file_name = trim(inp_folder_name)//"/"//integrals_file_name
    print *, "Reading one-electron integrals..."
    print *, "  file: "//trim(inp_file_name)
    open(unit=inp_file_unit, file=trim(inp_file_name), action="read")
    
    n = size(integrals, 1)
    ! The number of elements in the upper (lower) triangle is n (n + 1) / 2.
    do line = 1, (n * (n + 1)) / 2
      read (inp_file_unit,*) i, j, integral
!       print *, i, j, integral
      integrals(i, j) = integral
      integrals(j, i) = integral
    end do
    
    close(unit=inp_file_unit)
  end subroutine read_one_electron_integrals
  
  subroutine read_two_electron_integrals(two_electron_integrals)
    real(kind=d), dimension(:), intent(out) :: two_electron_integrals
    
    integer :: stat, i, j, k, l, ij, kl, ijkl
    real(kind=d) :: integral
    
    inp_file_name = trim(inp_folder_name)//"/"//two_electron_integrals_file_name
    print *, "Reading two-electron integrals..."
    print *, "  file: "//trim(inp_file_name)
    open(unit=inp_file_unit, file=trim(inp_file_name), action="read")

    do
      read (inp_file_unit,*,iostat=stat) i, j, k, l, integral
      if (stat == iostat_end) exit
      ij = compound_index(i, j)
      kl = compound_index(k, l)
      ijkl = compound_index(ij, kl)
!       print *, i, j, k, l, ij, kl, ijkl
      two_electron_integrals(ijkl) = integral      
    end do
    
    close(unit=inp_file_unit)
  end subroutine read_two_electron_integrals
  
  subroutine write_square_matrix(matrix)
    real(kind=d), dimension(:,:), intent(in) :: matrix
    
    integer, parameter :: max_cols = 10
    
    integer :: cols, rows, m, i, j
    character(len=*), parameter :: header_format = "(10i12)"
    character(len=*), parameter :: line_format = "(i5,10f12.7)"
    
    cols = size(matrix, 1)
    rows = cols
    
    do j = 1, cols, max_cols
      if (cols > max_cols) then
        m = j + max_cols - 1
      else
        m = j + cols - 1
      end if
      write (out_file_unit, "(bn)")
      write (out_file_unit, header_format) (i, i = j, m)
      write (out_file_unit, "(bn)")
      do i = 1, rows
        write (out_file_unit, line_format) i, matrix(i, j:m)
      end do
      cols = cols - max_cols
    end do
  end subroutine write_square_matrix
  
  function compound_index(i, j) result(res)
    integer, intent(in) :: i, j
    integer :: res
    
    if (i > j) then
      res = i * (i - 1) / 2 + j
    else
      res = j * (j - 1) / 2 + i
    end if
  end function compound_index
  
end program project_03