program project_03
  
  use fcl_kinds
  use fcl_lapack
  use fcl_util
  use chem
  use iso_fortran_env
  
  implicit none
  
  character(len=79), parameter :: separator = repeat("-", 79)

  integer, parameter :: max_columns = 10
  character(len=*), parameter :: format = "f12.7"
  
  integer, parameter :: inp_file_unit = 1
  integer, parameter :: out_file_unit = 2

  real(kind=dp), parameter :: delta_1 = 1d-12;
  real(kind=dp), parameter :: delta_2 = 1d-10;
  
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

  real(kind=dp), dimension(:, :), allocatable :: overlap_integrals, &
    overlap_eigenvectors, &
    kinetic_energy_integrals, &
    nuclear_attraction_integrals, &
    core_hamiltonian, &
    orthogonalization_matrix, &
    fock_matrix, &
    orthogonalized_fock_matrix, &
    coefficients_matrix, &
    density_matrix, &
    previous_density_matrix
  real(kind=dp), dimension(:), allocatable :: overlap_eigenvalues
  real(kind=dp), dimension(:), allocatable :: orbital_energies
  real(kind=dp), dimension(:), allocatable :: two_electron_integrals
  real(kind=dp) :: nuclear_repulsion_energy
  integer :: basis_set_size, number_of_occupied_orbitals, iteration
  integer :: i, j, k, l, m
  integer :: ij, kl, ijkl, ik, jl, ikjl
  real(kind=dp) :: electronic_energy, previous_electronic_energy, delta_electronic_energy, rms_density
  
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
  call allocate_arrays()

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
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, overlap_integrals, format, max_columns, &
    headers=.true., decorate=.false.)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Kinetic-Energy Integrals:"
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, kinetic_energy_integrals, format, max_columns, &
    headers=.true., decorate=.false.)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Nuclear Attraction Integrals:"
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, nuclear_attraction_integrals, format, max_columns, &
    headers=.true., decorate=.false.)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Core Hamiltonian:"
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, core_hamiltonian, format, max_columns, &
    headers=.true., decorate=.false.)
  
  call read_two_electron_integrals(two_electron_integrals)
  
  overlap_eigenvectors = overlap_integrals
  call fcl_lapack_dsyev(overlap_eigenvectors, overlap_eigenvalues, .true.)
  overlap_eigenvalues = overlap_eigenvalues ** (-1.0_dp / 2)
  overlap_integrals = 0.0_dp
  do i = 1, basis_set_size
    overlap_integrals(i, i) = overlap_eigenvalues(i)
  end do
  orthogonalization_matrix = matmul(overlap_eigenvectors, overlap_integrals)
  orthogonalization_matrix = matmul(orthogonalization_matrix, transpose(overlap_eigenvectors))
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	S^-1/2 Matrix:"
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, orthogonalization_matrix, format, max_columns, &
    headers=.true., decorate=.false.)
  
  fock_matrix = core_hamiltonian
  call build_orthogonalized_fock_matrix()
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Initial F' Matrix:"
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, orthogonalized_fock_matrix, format, max_columns, &
    headers=.true., decorate=.false.)
  
  call diagonalize_orthogonalized_fock_matrix()
  call build_coefficients_matrix()
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Initial C Matrix:"
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, coefficients_matrix, format, max_columns, &
    headers=.true., decorate=.false.)
  
  number_of_occupied_orbitals = molecule%number_of_electrons() / 2
  call build_density_matrix()
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(a)") "	Initial Density Matrix:"
  write (out_file_unit, "(bn)")
  call fcl_util_pretty_print(out_file_unit, density_matrix, format, max_columns, &
    headers=.true., decorate=.false.)
  write (out_file_unit, "(bn)")
  write (out_file_unit, "(bn)")

  write (out_file_unit, "(a5,a15,a20,a20,a20)") "Iter", "E(elec)", "E(tot)", "Delta(E)", "RMS(D)"
  iteration = 0
  previous_electronic_energy = 0.0_dp
  call calculate_electronic_energy()
  write (out_file_unit, "(i3.2,2f21.12)") iteration, electronic_energy, &
    electronic_energy + nuclear_repulsion_energy
  write (out_file_unit, "(bn)")

  do
    iteration = iteration + 1
    call build_fock_matrix()

    if (iteration == 1) then
      write (out_file_unit, "(a)") "  Fock Matrix:"
      write (out_file_unit, "(bn)")
      call fcl_util_pretty_print(out_file_unit, fock_matrix, format, max_columns, &
        headers=.true., decorate=.false.)
    endif
    
    call build_orthogonalized_fock_matrix()
    call diagonalize_orthogonalized_fock_matrix()
    call build_coefficients_matrix()
    previous_density_matrix = density_matrix
    call build_density_matrix()
    previous_electronic_energy = electronic_energy
    call calculate_electronic_energy()

    delta_electronic_energy = electronic_energy - previous_electronic_energy
    rms_density = 0.0_dp
    do i = 1, basis_set_size
      do j = 1, basis_set_size
        rms_density = rms_density + (density_matrix(i,j) - previous_density_matrix(i,j)) ** 2
      end do
    end do
    rms_density = sqrt(rms_density)
    write (out_file_unit, "(i3.2,4f21.12)") iteration, electronic_energy, &
      electronic_energy + nuclear_repulsion_energy, delta_electronic_energy, rms_density

    if (delta_electronic_energy < delta_1 .and. rms_density < delta_2) exit 
  end do

  call deallocate_arrays()
  close(unit=out_file_unit)

contains

  subroutine allocate_arrays()
    integer :: n
    allocate (overlap_integrals(basis_set_size, basis_set_size))
    allocate (overlap_eigenvectors(basis_set_size, basis_set_size))
    allocate (orthogonalization_matrix(basis_set_size, basis_set_size))
    allocate (overlap_eigenvalues(basis_set_size))
    allocate (kinetic_energy_integrals(basis_set_size, basis_set_size))
    allocate (nuclear_attraction_integrals(basis_set_size, basis_set_size))
    allocate (core_hamiltonian(basis_set_size, basis_set_size))
    allocate (fock_matrix(basis_set_size, basis_set_size))
    allocate (orthogonalized_fock_matrix(basis_set_size, basis_set_size))
    allocate (coefficients_matrix(basis_set_size, basis_set_size))
    allocate (orbital_energies(basis_set_size))
    allocate (density_matrix(basis_set_size, basis_set_size))
    allocate (previous_density_matrix(basis_set_size, basis_set_size))
    n = (basis_set_size * (basis_set_size + 1)) / 2
    allocate (two_electron_integrals((n * (n + 1)) / 2))
  end subroutine allocate_arrays

  subroutine deallocate_arrays()
    deallocate (overlap_integrals)
    deallocate (overlap_eigenvectors)
    deallocate (orthogonalization_matrix)
    deallocate (kinetic_energy_integrals)
    deallocate (nuclear_attraction_integrals)
    deallocate (core_hamiltonian)
    deallocate (fock_matrix)
    deallocate (orthogonalized_fock_matrix)
    deallocate (coefficients_matrix)
    deallocate (orbital_energies)
    deallocate (density_matrix)
    deallocate (previous_density_matrix)
    deallocate (two_electron_integrals)
  end subroutine deallocate_arrays
  
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
    real(kind=dp), dimension(:,:), intent(out) :: integrals
    
    integer :: n, line, i, j
    real(kind=dp) :: integral
    
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
    real(kind=dp), dimension(:), intent(out) :: two_electron_integrals
    
    integer :: stat, i, j, k, l, ij, kl, ijkl
    real(kind=dp) :: integral
    
    inp_file_name = trim(inp_folder_name)//"/"//two_electron_integrals_file_name
    print *, "Reading two-electron integrals..."
    print *, "  file: "//trim(inp_file_name)
    open(unit=inp_file_unit, file=trim(inp_file_name), action="read")

    two_electron_integrals = 0.0_dp

    do
      read (inp_file_unit,*,iostat=stat) i, j, k, l, integral
      if (stat == iostat_end) exit
      ij = compound_index(i, j)
      kl = compound_index(k, l)
      ijkl = compound_index(ij, kl)
      two_electron_integrals(ijkl) = integral      
    end do
    
    close(unit=inp_file_unit)
  end subroutine read_two_electron_integrals
  
  function compound_index(i, j) result(res)
    integer, intent(in) :: i, j
    integer :: res
    
    if (i > j) then
      res = i * (i - 1) / 2 + j
    else
      res = j * (j - 1) / 2 + i
    end if
  end function compound_index

  subroutine build_fock_matrix()
    do i = 1, basis_set_size
      do j = 1, basis_set_size
        fock_matrix(i, j) = core_hamiltonian(i, j)
        do k = 1, basis_set_size
          do l = 1, basis_set_size
            ij = compound_index(i, j)
            kl = compound_index(k, l)
            ijkl = compound_index(ij, kl)
            ik = compound_index(i, k)
            jl = compound_index(j, l)
            ikjl = compound_index(ik, jl)
            fock_matrix(i, j) = fock_matrix(i, j) + density_matrix(k, l) * &
              (2 * two_electron_integrals(ijkl) - two_electron_integrals(ikjl))
          end do
        end do
      end do
    end do      
  end subroutine build_fock_matrix

  subroutine build_orthogonalized_fock_matrix()
    orthogonalized_fock_matrix = matmul(transpose(orthogonalization_matrix), fock_matrix)
    orthogonalized_fock_matrix = matmul(orthogonalized_fock_matrix, orthogonalization_matrix)
  end subroutine build_orthogonalized_fock_matrix

  subroutine diagonalize_orthogonalized_fock_matrix()
    call fcl_lapack_dsyev(orthogonalized_fock_matrix, orbital_energies, .true.)
  end subroutine diagonalize_orthogonalized_fock_matrix

  subroutine build_coefficients_matrix
    coefficients_matrix = matmul(orthogonalization_matrix, orthogonalized_fock_matrix)
  end subroutine build_coefficients_matrix

  subroutine build_density_matrix()
    density_matrix = 0.0_dp
    do i = 1, basis_set_size
      do j = 1, basis_set_size
        do m = 1, number_of_occupied_orbitals
          density_matrix(i, j) = density_matrix(i, j) + coefficients_matrix(i, m) * coefficients_matrix(j, m)
        end do
      end do
    end do
  end subroutine build_density_matrix

  subroutine calculate_electronic_energy()
    electronic_energy = 0.0_dp
    do i = 1, basis_set_size
      do j = 1, basis_set_size
        electronic_energy = electronic_energy + density_matrix(i, j) * &
          (core_hamiltonian(i, j) + fock_matrix(i, j))
      end do
    end do
  end subroutine calculate_electronic_energy

end program project_03