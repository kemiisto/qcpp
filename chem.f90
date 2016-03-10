module chem

  use fcl_kinds
  use fcl_vecmath_vector_3d

  implicit none

  private

  public :: chem_mod_atom, &
    chem_mod_molecule, &
    chem_mod_atomic_masses, &
    chem_mod_read_molecule_from_file

  type chem_mod_atom
    integer :: atomic_number
    type(fcl_vecmath_mod_vector3d) :: position
  contains
    procedure :: set_position => chem_mod_atom_set_position
  end type chem_mod_atom

  type chem_mod_atom_pointer
    type(chem_mod_atom), pointer :: atom
  end type chem_mod_atom_pointer

  type chem_mod_molecule
    private
    type(chem_mod_atom_pointer), dimension(:), allocatable :: atoms
  contains
    procedure :: number_of_atoms              => chem_mod_molecule_number_of_atoms
    procedure :: set_number_of_atoms          => chem_mod_molecule_set_number_of_atoms
    procedure :: atom                         => chem_mod_molecule_atom
    procedure :: atom_pointer                 => chem_mod_molecule_atom_pointer
    procedure :: distance                     => chem_mod_molecule_distance
    procedure :: angle                        => chem_mod_molecule_angle
    procedure :: out_of_plane_angle           => chem_mod_molecule_out_of_plane_angle
    procedure :: dihedral_angle               => chem_mod_molecule_dihedral_angle
    procedure :: center_of_mass               => chem_mod_molecule_center_of_mass
    procedure :: translate                    => chem_mod_molecule_translate
    procedure :: moment_of_inertia_tensor     => chem_mod_molecule_moment_of_inertia_tensor
  end type chem_mod_molecule

  ! Atomic masses of the most abundant isotopes of the first 18 elements.
  real(kind=d), dimension(18), parameter :: chem_mod_atomic_masses = &
  [                    &
     1.00782503223_d,  &
     4.00260325413_d,  &
     7.01600343660_d,  &
     9.01218306500_d,  &
    11.00930536000_d,  &
    12.00000000000_d,  &
    14.00307400443_d,  &
    15.99491461957_d,  &
    18.99840316273_d,  &
    19.99244017620_d,  &
    22.98976928200_d,  &
    23.98504169700_d,  &
    26.98153853000_d,  &
    27.97692653465_d,  &
    30.97376199842_d,  &
    31.97207117440_d,  &
    34.96885268200_d,  &
    39.96238312370_d   &
  ]   

contains

  subroutine chem_mod_atom_set_position(this, v)
    class(chem_mod_atom), intent(inout) :: this
    type(fcl_vecmath_mod_vector3d), intent(in) :: v

    this%position = v
  end subroutine chem_mod_atom_set_position

  pure function chem_mod_molecule_number_of_atoms(this) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer :: f

    f = size(this%atoms)
  end function chem_mod_molecule_number_of_atoms

  subroutine chem_mod_molecule_set_number_of_atoms(this, n)
    class(chem_mod_molecule), intent(inout) :: this
    integer, intent(in) :: n

    integer :: i

    allocate(this%atoms(n))
    do i = 1, n
      allocate(this%atoms(i)%atom)
    end do
    ! deallocate in destructor?
  end subroutine chem_mod_molecule_set_number_of_atoms

  pure function chem_mod_molecule_atom(this, i) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i
    type(chem_mod_atom) :: f

    f = this%atoms(i)%atom
  end function chem_mod_molecule_atom

  function chem_mod_molecule_atom_pointer(this, i) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i
    type(chem_mod_atom), pointer :: f

    f => this%atoms(i)%atom
  end function chem_mod_molecule_atom_pointer

  pure function chem_mod_molecule_distance(this, i, j) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i, j
    real(kind=d) :: f

    type(fcl_vecmath_mod_vector3d) :: v_i, v_j

    v_i = this%atoms(i)%atom%position
    v_j = this%atoms(j)%atom%position

    v_i = v_i - v_j
    f = v_i%norm()
  end function chem_mod_molecule_distance

  ! Angle (in radians) between atoms i-j-k, where j is the central atom.
  !
  ! cos(phi_ijk) = e_ji . e_jk = (v_ji . v_jk) / (||v_ji|| ||v_jk||)
  pure function chem_mod_molecule_angle(this, i, j, k) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(kind=d) :: f

    type(fcl_vecmath_mod_vector3d) :: v_i, v_j, v_k, v_ji, v_jk

    v_i = this%atoms(i)%atom%position
    v_j = this%atoms(j)%atom%position
    v_k = this%atoms(k)%atom%position

    v_ji = v_i - v_j
    v_jk = v_k - v_j

    f = v_ji%angle(v_jk)
  end function chem_mod_molecule_angle
  
  ! Angle (in radians) for atom i out of the plane containing atoms j-k-l 
  ! with k being the central atom, connected to i.
  ! In other words, angle between a bond i-k and a plane defined by two bonds j-k and k-l.
  ! 
  ! sin(theta_ijkl) = (e_kj x e_kl) / sin(phi_jkl) . e_ki
  pure function chem_mod_molecule_out_of_plane_angle(this, i, j, k, l) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i, j, k, l
    real(kind=d) :: f

    type(fcl_vecmath_mod_vector3d) :: v_i, v_j, v_k, v_l, e_kj, e_kl, e_ki
    real(kind=d) :: sin_theta

    v_i = this%atoms(i)%atom%position
    v_j = this%atoms(j)%atom%position
    v_k = this%atoms(k)%atom%position
    v_l = this%atoms(l)%atom%position
    
    e_kj = v_j - v_k
    e_kj = e_kj%normalized()
    e_kl = v_l - v_k
    e_kl = e_kl%normalized()
    e_ki = v_i - v_k
    e_ki = e_ki%normalized()

    sin_theta = ((e_kj .cross. e_kl) / sin(this%angle(j, k, l))) .dot. e_ki
    ! The arguments of asin function should be inbetween -1.0 and +1.0.
    ! However, numerical precision in the calculation of the cross- and dot-products 
    ! earlier in the calculation can yield results slightly outside this domain,
    ! in which case asin will return NaN.
    ! This is obviously undesirable so we test the argument before calling asin.
    if (sin_theta < -1.0_d) then
      sin_theta = -1.0_d
    else if (sin_theta > 1.0_d) then
      sin_theta = 1.0_d
    end if

    f = asin(sin_theta)
  end function chem_mod_molecule_out_of_plane_angle

  ! Dihedral (or torsion) angle of four atoms i-j-k-l, i.e. angle between planes i-j-k and j-k-l.
  ! Calculated as angle between the planes normal vectors.
  pure function chem_mod_molecule_dihedral_angle(this, i, j, k, l) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i, j, k, l
    real(kind=d) :: f

    type(fcl_vecmath_mod_vector3d) :: v_i, v_j, v_k, v_l, &
      v_ij, v_jk, v_kl, n_ijk, n_jkl, &
      cross

    v_i = this%atoms(i)%atom%position
    v_j = this%atoms(j)%atom%position
    v_k = this%atoms(k)%atom%position
    v_l = this%atoms(l)%atom%position

    v_ij = v_j - v_i
    v_jk = v_k - v_j
    v_kl = v_l - v_k

    n_ijk = v_ij .cross. v_jk
    n_jkl = v_jk .cross. v_kl

    f = n_ijk%angle(n_jkl)

    ! sign
    cross = n_ijk .cross. n_jkl
    if ((cross .dot. v_ij) > 0.0_d) f = -f

  end function chem_mod_molecule_dihedral_angle

  pure function chem_mod_molecule_center_of_mass(this) result(f)
    class(chem_mod_molecule), intent(in) :: this
    type(fcl_vecmath_mod_vector3d) :: f

    integer :: i
    type(chem_mod_atom) :: atom
    real(kind=d) :: total_m, total_mx, total_my, total_mz, m

    total_m  = 0.0_d
    total_mx = 0.0_d
    total_my = 0.0_d
    total_mz = 0.0_d
    do i = 1, this%number_of_atoms()
      atom = this%atom(i)
      m = chem_mod_atomic_masses(atom%atomic_number)
      total_m  = total_m  + m
      total_mx = total_mx + m * atom%position%x()
      total_my = total_my + m * atom%position%y()
      total_mz = total_mz + m * atom%position%z()
    end do

    f = fcl_vecmath_mod_vector3d([total_mx / total_m, total_my / total_m, total_mz / total_m])
  end function chem_mod_molecule_center_of_mass

  subroutine chem_mod_molecule_translate(this, v)
    class(chem_mod_molecule), intent(inout) :: this
    type(fcl_vecmath_mod_vector3d), intent(in) :: v

    integer :: i
    type(chem_mod_atom), pointer :: atom

    do i = 1, this%number_of_atoms()
      atom => this%atom_pointer(i)
      call atom%set_position(atom%position + v)
    end do
  end subroutine chem_mod_molecule_translate

  pure function chem_mod_molecule_moment_of_inertia_tensor(this) result(f)
    class(chem_mod_molecule), intent(in) :: this
    real(kind=d), dimension(3, 3) :: f

    integer :: i
    type(chem_mod_atom) :: atom
    real(kind=d) :: m

    f = 0.0_d
    do i = 1, this%number_of_atoms()
      atom = this%atom(i)
      m = chem_mod_atomic_masses(atom%atomic_number)
      f(1, 1) = f(1, 1) + m * (atom%position%y() ** 2 + atom%position%z() ** 2)
      f(2, 2) = f(2, 2) + m * (atom%position%x() ** 2 + atom%position%z() ** 2)
      f(3, 3) = f(3, 3) + m * (atom%position%x() ** 2 + atom%position%y() ** 2)
      f(1, 2) = f(1, 2) + m * atom%position%x() * atom%position%y()
      f(1, 3) = f(1, 3) + m * atom%position%x() * atom%position%z()
      f(2, 3) = f(2, 3) + m * atom%position%y() * atom%position%z()
      f(2, 1) = f(1, 2)
      f(3, 1) = f(1, 3)
      f(3, 2) = f(2, 3)
    end do
  end function chem_mod_molecule_moment_of_inertia_tensor
  
  subroutine chem_mod_read_molecule_from_file(molecule, inp_file_name)
    type(chem_mod_molecule), intent(out) :: molecule
    character(len=*), intent(in) :: inp_file_name

    integer, parameter :: inp_file_unit = 1
    integer :: i, number_of_atoms
    real :: atomic_number
    type(chem_mod_atom), pointer :: atom
    type(fcl_vecmath_mod_vector3d) :: position
    
    open(unit=inp_file_unit, file=inp_file_name, action="read")
    
    read (inp_file_unit,*) number_of_atoms
    call molecule%set_number_of_atoms(number_of_atoms)
    
    do i = 1, number_of_atoms
      read (inp_file_unit,*) atomic_number, position
      atom => molecule%atom_pointer(i)
      atom%atomic_number = int(atomic_number)
      atom%position = position
    end do
    
    close(unit=inp_file_unit)
  end subroutine chem_mod_read_molecule_from_file

end module chem
