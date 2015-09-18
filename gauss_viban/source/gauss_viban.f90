program gauss_viban
use routines
implicit none

integer,parameter :: fcheckio = 1                   ! IO unit of the formatted checkpoint file
integer,parameter :: prec = 5                       ! number of digits to be printed out
real(dp),parameter :: pi = 3.141592653589793238_dp  ! pi
real(dp),parameter :: c0 = 299792458.0_dp           ! the speed of light
real(dp),parameter :: Eh = 4.3597438e-18_dp         ! the Hartree energy
real(dp),parameter :: a0 = 5.291772083e-11_dp       ! the bohr radius
real(dp),parameter :: hbar = 1.054571596e-34_dp     ! Planck's constant by 2 Pi
real(dp),parameter :: me = 9.10938188e-31_dp        ! the electron mass
real(dp),parameter :: amu2au = 1822.88848325_dp     ! to convert from amu to atomic units
integer :: Natoms                                   ! number of atoms (N)
integer :: Ncoords                                  ! number of coordinates (3N)
integer :: Nmodes                                   ! number of normal modes (3N-6)
integer :: i, j, info                               ! loop indices etc.
integer :: lwork                                    ! for LAPACK
real(dp) :: total_mass                              ! total mass of the molecule
real(dp),dimension(3) :: com                        ! position of the center of mass
real(dp),dimension(3) :: moments                    ! eigenvalues of inertia tensor
real(dp),dimension(3,3) :: inert                    ! the inertia tensor
real(dp),dimension(3,3) :: prinaxes                 ! eigenvectors of inertia tensor
real(dp),dimension(:),allocatable :: masses         ! the atomic masses (length N)
real(dp),dimension(:),allocatable :: mass_vector    ! the atomic masses (length 3N)
real(dp),dimension(:),allocatable :: x              ! the initial atomic positions
real(dp),dimension(:),allocatable :: q              ! the mass-weighted initial atomic positions
real(dp),dimension(:),allocatable :: s              ! the transformed atomic positions
real(dp),dimension(:),allocatable :: lambdas        ! the eigenvalues of the hessian
real(dp),dimension(:),allocatable :: freqs          ! the normal mode frequencies
real(dp),dimension(:),allocatable :: work           ! working array for the LAPACK routine
real(dp),dimension(:),allocatable :: red_masses     ! the reduced masses of the normal modes
real(dp),dimension(:),allocatable :: force_consts   ! the force constants of the normal modes
real(dp),dimension(:,:),allocatable :: m_mat        ! the inverse mass-matrix
real(dp),dimension(:,:),allocatable :: f_cart       ! the hessian in cartesian coordinates
real(dp),dimension(:,:),allocatable :: f_mwc        ! the hessian in mass weighted cart. coordinates
real(dp),dimension(:,:),allocatable :: f_diag       ! diagonalized hessian
real(dp),dimension(:,:),allocatable :: f_int        ! the hessian in internal coodinates
real(dp),dimension(:,:),allocatable :: f_sub        ! submatrix of the hessian without translation and rotation
real(dp),dimension(:,:),allocatable :: l_init       ! initial displacements (L matrix)
real(dp),dimension(:,:),allocatable :: l_mwc        ! mass-weighted cartesian displacements
real(dp),dimension(:,:),allocatable :: l_cart       ! non-mass-weighted cartesian displacements
real(dp),dimension(:,:),allocatable :: D            ! the D matrix
real(dp),dimension(:,:),allocatable :: metric       ! metric matrix to test orthogonality
real(dp),dimension(:,:),allocatable :: modes        ! the normal modes (eigenvectors of the hessian)
real(dp),dimension(:,:),allocatable :: gauss_modes  ! the normal modes calculated by gaussian
real(dp),dimension(:,:),allocatable :: test_hessian
character(len=50) :: fcheckfile                     ! name of the checkpoint file containing the data
logical :: lerr
logical :: cleanup = .true.


! open the checkpoint file:
call get_command_argument(1, fcheckfile)
if (len_trim(fcheckfile) == 0) then
    write(*,*) "ERROR opening fcheck file, ", fcheckfile
    stop
endif
open(unit=fcheckio,file=fcheckfile,status='old',action='read')

! read the data from the checkpoint file:
call searchstring(fcheckio, "Number of atoms", .true., lerr)
read(fcheckio,'(55x,i6)') Natoms
Ncoords = 3 * Natoms
Nmodes = Ncoords - 6
allocate(masses(Natoms))
allocate(mass_vector(Ncoords))
allocate(x(Ncoords))
allocate(q(Ncoords))
allocate(s(Ncoords))
allocate(lambdas(Ncoords))
allocate(freqs(Ncoords))
allocate(red_masses(Ncoords))
allocate(force_consts(Ncoords))
allocate(m_mat(Ncoords,Ncoords))
allocate(f_cart(Ncoords,Ncoords))
allocate(f_mwc(Ncoords,Ncoords))
allocate(f_diag(Ncoords,Ncoords))
allocate(f_int(Ncoords,Ncoords))
allocate(f_sub(Nmodes,Nmodes))
allocate(D(Ncoords,Ncoords))
allocate(metric(Ncoords,Ncoords))
allocate(modes(Ncoords,Ncoords))
allocate(gauss_modes(Ncoords,Nmodes))
allocate(test_hessian(Nmodes,Nmodes))
call searchstring(fcheckio, "Current cartesian coordinates", .true., lerr)
read(fcheckio,*)
i = 1
call readchkvec(fcheckio, Ncoords, x, i)
call searchstring(fcheckio, "Real atomic weights", .true., lerr)
read(fcheckio,*)
i = 1
call readchkvec(fcheckio, Natoms, masses, i)
call searchstring(fcheckio, "Cartesian Force Constants", .true., lerr)
read(fcheckio,*)
call readchkmat(fcheckio, Ncoords, f_cart)
call searchstring(fcheckio, "Vib-Modes", .true., lerr)
read(fcheckio,*)
call readgaussmodes(fcheckio, Ncoords, gauss_modes)
masses = masses * amu2au
total_mass = sum(masses)
do i = 1, Natoms
    mass_vector(3*(i-1)+1) = masses(i)
    mass_vector(3*(i-1)+2) = masses(i)
    mass_vector(3*(i-1)+3) = masses(i)
enddo
m_mat = 0.0_dp
do i = 1, Ncoords
    m_mat(i,i) = 1.0_dp / sqrt(mass_vector(i))
enddo

! write the initial data:
write(*,*) "Number of atoms: ", Natoms
write(*,*) "masses:"
call write_vector(masses, prec, cleanup)
write(*,*) "total mass: ", total_mass
write(*,*) "positions:"
call write_vector(x, prec, cleanup)
write(*,*) "mass-weighted positions:"
q = x * mass_vector
call write_vector(q, prec, cleanup)
write(*,*) "inverse-mass matrix:"
call write_matrix(m_mat, prec, cleanup)
write(*,*) "non-mass-weighted initial cartesian hessian:"
call write_matrix(f_cart, prec, cleanup)

! calculate the mass-weighted hessian:
call ortho_trans(f_cart, m_mat, f_mwc)
write(*,*) "mass-weighted cartesian hessian:"
call write_matrix(f_mwc, prec, cleanup)

! diagonalize the mass-weighted hessian:
f_diag = f_mwc
lwork = -1
allocate(work(1))
call dsyev('N', 'U', Ncoords, f_diag, Ncoords, lambdas, work, lwork, info)
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
f_diag = f_mwc
call dsyev('N', 'U', Ncoords, f_diag, Ncoords, lambdas, work, lwork, info)
write(*,*) "status of diagonalization: ", info

! conversion from atomic units to cm-1:
!  nu_tilde = 1 / (2 pi c 100 cm/m) * sqrt(lambda * Eh / (a0**2 u))
write(*,*) "eigenvalues:"
do i = 1, Ncoords
    write(*,'(1x, es15.7, f10.4)') lambdas(i), &
        sqrt(abs(lambdas(i)) * Eh / (a0**2 * me)) / (2 * pi * c0 * 100.0_dp)
enddo

! calculate the center of mass:
call calc_com(x, masses, total_mass, com)
write(*,*) "center of mass:"
call write_vector(com, prec, cleanup)

! calculate and diagonalize the inertia tensor:
call calc_inert(x, masses, total_mass, inert)
write(*,*) "inertia tensor:"
call write_matrix(inert, prec, cleanup)
prinaxes = inert
call dsyev('V', 'U', 3, prinaxes, 3, moments, work, lwork, info)
write(*,*) "status of diagonalization: ", info
write(*,*) "moments:"
call write_vector(moments, prec, cleanup)
write(*,*) "principal axes:"
call write_matrix(prinaxes, prec, cleanup)

! calculate the D matrix
D = 0.0_dp
do i = 1, Ncoords
    do j = 7, Ncoords
        call random_number(D(i,j))
    enddo
enddo
do i = 1, 3
    do j = 1, Natoms
        D(3*(j-1)+i,i) = sqrt(masses(j))
    enddo
enddo
do i = 1, Natoms    ! sum over atoms
    D(3*(i-1)+1,4) =  0.0_dp
    D(3*(i-1)+2,4) = -dot_product(x(3*(i-1)+1:3*(i-1)+3), prinaxes(3,:)) * sqrt(masses(i))
    D(3*(i-1)+3,4) =  dot_product(x(3*(i-1)+1:3*(i-1)+3), prinaxes(2,:)) * sqrt(masses(i))
    D(3*(i-1)+1,5) =  dot_product(x(3*(i-1)+1:3*(i-1)+3), prinaxes(3,:)) * sqrt(masses(i))
    D(3*(i-1)+2,5) =  0.0_dp
    D(3*(i-1)+3,5) = -dot_product(x(3*(i-1)+1:3*(i-1)+3), prinaxes(1,:)) * sqrt(masses(i))
    D(3*(i-1)+1,6) = -dot_product(x(3*(i-1)+1:3*(i-1)+3), prinaxes(2,:)) * sqrt(masses(i))
    D(3*(i-1)+2,6) =  dot_product(x(3*(i-1)+1:3*(i-1)+3), prinaxes(1,:)) * sqrt(masses(i))
    D(3*(i-1)+3,6) =  0.0_dp
enddo

write(*,*) "the initial D matrix:"
call write_matrix(D, prec, cleanup)
metric = matmul(transpose(D), D)
write(*,*) "the metric of the initial D matrix:"
call write_matrix(metric, prec, cleanup)

! normalize the D matrix:
do i = 1, Ncoords
    D(:,i) = D(:,i) / norm2(D(:,i))
enddo
write(*,*) "the normalized D matrix:"
call write_matrix(D, prec, cleanup)
metric = matmul(transpose(D), D)
write(*,*) "the metric of the normalized D matrix:"
call write_matrix(metric, prec, cleanup)

! orthogonalize the D-matrix:
call gram_schmidt(D)
write(*,*) "the orthogonalized D matrix:"
call write_matrix(D, prec, cleanup)
metric = matmul(transpose(D), D)
write(*,*) "the metric of the orthogonalized D matrix:"
call write_matrix(metric, prec, cleanup)

! transform the positions and the hessian to internal coordinates:
s = matmul(D, q)
write(*,*) "mass-weighted coordinates transformed by D:"
call write_vector(s, prec, cleanup)
call ortho_trans(f_mwc, D, f_int)
write(*,*) "the hessian in internal coordinates:"
call write_matrix(f_int, prec, cleanup)

! translate and rotate the molecule and write updated properties
do i = 1, Natoms
    x(3*(i-1)+1:3*i) = x(3*(i-1)+1:3*i) - com
enddo
do i = 1, Natoms
    x(3*(i-1)+1:3*i) = matmul(transpose(prinaxes), x(3*(i-1)+1:3*i))
enddo
q = x * mass_vector
call calc_com(x, masses, total_mass, com)
call calc_inert(x, masses, total_mass, inert)
write(*,*) "new positions:"
call write_vector(x, prec, cleanup)
write(*,*) "new mass-weighted positions:"
call write_vector(q, prec, cleanup)
write(*,*) "new center of mass:"
call write_vector(com, prec, cleanup)
write(*,*) "new inertia tensor:"
call write_matrix(inert, prec, cleanup)

! diagonalize the hessian submatrix:
f_sub = f_int(7:Ncoords,7:Ncoords)
lambdas = 0.0_dp
call dsyev('V', 'U', Nmodes, f_sub, Nmodes, lambdas, work, lwork, info)
write(*,*) "status of diagonalization: ", info
write(*,*) "eigenvalues of the hessian submatrix (real frequencies)"
do i = 1, Ncoords
    write(*,'(1x, es15.7, f10.4)') lambdas(i), &
        sqrt(abs(lambdas(i)) * Eh / (a0**2 * me)) / (2 * pi * c0 * 100.0)
enddo
freqs = sqrt(lambdas / (4.0_dp * pi**2 * (me * a0 * c0 / hbar)**2))
write(*,*) "frequencies in atomic units:"
call write_vector(freqs, prec, cleanup)

! diagonalize the full internal hessian:
l_init = f_int
call dsyev('V', 'U', Ncoords, l_init, Ncoords, freqs, work, lwork, info)
write(*,*) "status of diagonalization: ", info
write(*,*) "initial displacements:"
call write_matrix(l_init, prec, cleanup)
metric = matmul(transpose(l_init), l_init)
write(*,*) "metric of the initial displacements:"
call write_matrix(metric, prec, cleanup)

! transform the displacements:
l_mwc = matmul(D, l_init)
write(*,*) "transformed mass-weighted displacements:"
call write_matrix(l_mwc, prec, cleanup)
metric = matmul(transpose(l_mwc), l_mwc)
write(*,*) "metric of the mass-weighted displacements:"
call write_matrix(metric, prec, cleanup)
l_cart = matmul(m_mat, l_mwc)
write(*,*) "non-mass-weighted cartesian displacements:"
call write_matrix(l_cart, prec, cleanup)
metric = matmul(transpose(l_cart), l_cart)
write(*,*) "metric of the cartesian displacements:"
call write_matrix(metric, prec, cleanup)

! normalize l_cart:
do i = 1, Ncoords
    red_masses(i) = 1.0_dp / norm2(l_cart(:,i))**2
    l_cart(:,i) = l_cart(:,i) * sqrt(red_masses(i))
enddo

write(*,*) "the normalized non-mass-weighted cartesian displacements:"
call write_matrix(l_cart, prec, cleanup)
metric = matmul(transpose(l_cart), l_cart)
write(*,*) "metric of the normalized cartesian displacements:"
call write_matrix(metric, prec, cleanup)


write(*,*) "non-mass-weighted cartesian displacements in Gaussian HPmodes style:"
do i = 1, Ncoords
    do j = 1, Ncoords
        write(*,'(1x, f9.5)',advance='no') l_cart(i,j)
    enddo
    write(*,*)
enddo


! calculate force constants:
force_consts = 4 * pi**2 * freqs**2 * red_masses

write(*,*) "reduced masses of the normal modes:"
call write_vector(red_masses, prec, cleanup)
write(*,*) "force constants of the normal modes:"
call write_vector(force_consts, prec, cleanup)

write(*,*) "modes calculated by gaussian"
call write_matrix(gauss_modes, prec, cleanup)
write(*,*) "modes calculated by gaussian in Gaussian HPmodes style:"
do i = 1, Ncoords
    do j = 1, Nmodes
        write(*,'(1x, f9.5)',advance='no') gauss_modes(i,j)
    enddo
    write(*,*)
enddo




end program gauss_viban
