module chem

  use kinds
  use vecmath

  implicit none

  private

  public :: chem_mod_atom, chem_mod_molecule

  type chem_mod_atom
    integer :: atomic_number
    type(vecmath_vector3d) :: position
  end type chem_mod_atom

  type chem_mod_atom_pointer
    type(chem_mod_atom), pointer :: atom
  end type chem_mod_atom_pointer

  type chem_mod_molecule
    private
    type(chem_mod_atom_pointer), dimension(:), allocatable :: atoms
  contains
    procedure :: number_of_atoms     => chem_mod_molecule_number_of_atoms
    procedure :: set_number_of_atoms => chem_mod_molecule_set_number_of_atoms
    procedure :: atom                => chem_mod_molecule_atom
    procedure :: distance            => chem_mod_molecule_distance
    procedure :: angle               => chem_mod_molecule_angle
    procedure :: out_of_plane_angle  => chem_mod_molecule_out_of_plane_angle
    procedure :: dihedral_angle      => chem_mod_molecule_dihedral_angle
  end type chem_mod_molecule

contains

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

  function chem_mod_molecule_atom(this, i) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i
    type(chem_mod_atom), pointer :: f

    f => this%atoms(i)%atom
  end function chem_mod_molecule_atom

  pure function chem_mod_molecule_distance(this, i, j) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i, j
    real(kind=d) :: f

    type(vecmath_vector3d) :: v_i, v_j

    v_i = this%atoms(i)%atom%position
    v_j = this%atoms(j)%atom%position

    f = sqrt( (v_i%x - v_j%x)**2 + (v_i%y - v_j%y)**2 + (v_i%z - v_j%z)**2 )
  end function chem_mod_molecule_distance

  ! Angle (in radians) between atoms i-j-k, where j is the central atom.
  !
  ! cos(phi_ijk) = e_ji . e_jk
  pure function chem_mod_molecule_angle(this, i, j, k) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(kind=d) :: f

    type(vecmath_vector3d) :: v_i, v_j, v_k, v_ji, v_jk

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

    type(vecmath_vector3d) :: v_i, v_j, v_k, v_l, e_kj, e_kl, e_ki
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
  ! Calculated as angle between the planes unit normal vectors.
  pure function chem_mod_molecule_dihedral_angle(this, i, j, k, l) result(f)
    class(chem_mod_molecule), intent(in) :: this
    integer, intent(in) :: i, j, k, l
    real(kind=d) :: f

    type(vecmath_vector3d) :: v_i, v_j, v_k, v_l, &
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

end module chem
