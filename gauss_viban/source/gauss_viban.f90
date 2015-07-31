program gauss_viban
use routines
implicit none

integer,parameter :: fcheckio = 1                   ! IO unit of the formatted checkpoint file
integer :: Natoms                                   ! number of atoms
integer :: i, j, info                               ! loop indices etc.
integer :: lwork                                    ! for LAPACK
real(dp) :: total_mass                              ! total mass of the molecule
real(dp),dimension(3) :: com                        ! position of the center of mass
real(dp),dimension(3) :: moments                    ! eigenvalues of inertia tensor
real(dp),dimension(3,3) :: inert                    ! the inertia tensor
real(dp),dimension(3,3) :: prinaxes                 ! eigenvectors of inertia tensor
real(dp),dimension(:),allocatable :: masses         ! the atomic masses (length N)
real(dp),dimension(:),allocatable :: mass_vector    ! the atomic masses (length 3N)
real(dp),dimension(:),allocatable :: positions      ! the atomic positions
real(dp),dimension(:),allocatable :: freqs          ! the normal mode frequencies
real(dp),dimension(:),allocatable :: work           ! working array for the LAPACK routine
real(dp),dimension(:,:),allocatable :: f_cart       ! the hessian in cartesian coordinates
real(dp),dimension(:,:),allocatable :: f_mwc        ! the hessian in mass weighted cart. coordinates
real(dp),dimension(:,:),allocatable :: f_diag       ! diagonalized hessian
real(dp),dimension(:,:),allocatable :: f_int        ! the hessian in internal coodinates
real(dp),dimension(:,:),allocatable :: D            ! the D matrix
real(dp),dimension(:,:),allocatable :: metric       ! metric matrix to test orthogonality
real(dp),dimension(:,:),allocatable :: modes        ! the normal modes (eigenvectors of the hessian)
real(dp),dimension(:,:),allocatable :: test_hessian
character(len=50) :: fcheckfile                     ! name of the checkpoint file containing the data
logical :: lerr


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
allocate(masses(Natoms))
allocate(mass_vector(3 * Natoms))
allocate(positions(3 * Natoms))
allocate(freqs(3 * Natoms))
allocate(f_cart(3 * Natoms, 3 * Natoms))
allocate(f_mwc(3 * Natoms, 3 * Natoms))
allocate(f_diag(3 * Natoms, 3 * Natoms))
allocate(f_int(3 * Natoms, 3 * Natoms))
allocate(D(3 * Natoms, 3 * Natoms))
allocate(metric(3 * Natoms, 3 * Natoms))
allocate(modes(3 * Natoms, 3 * Natoms))
allocate(test_hessian(3*Natoms-6, 3*Natoms-6))
call searchstring(fcheckio, "Current cartesian coordinates", .true., lerr)
read(fcheckio,*)
i = 1
call readchkvec(fcheckio, size(positions), positions, i)
call searchstring(fcheckio, "Real atomic weights", .true., lerr)
read(fcheckio,*)
i = 1
call readchkvec(fcheckio, size(masses), masses, i)
call searchstring(fcheckio, "Cartesian Force Constants", .true., lerr)
read(fcheckio,*)
call readchkmat(fcheckio, size(f_cart, 1), f_cart)
total_mass = sum(masses)
do i = 1, Natoms
    mass_vector(3*(i-1)+1) = masses(i)
    mass_vector(3*(i-1)+2) = masses(i)
    mass_vector(3*(i-1)+3) = masses(i)
enddo

! write the initial data:
write(*,*) "Number of atoms: ", Natoms
write(*,*) "masses:"
call write_vector(masses)
write(*,*) "total mass: ", total_mass
write(*,*) "positions:"
call write_vector(positions)
write(*,*) "non-mass-weighted initial cartesian hessian:"
call write_matrix(f_cart)

! calculate the mass-weighted hessian:
do i = 1, 3 * Natoms
    do j = 1, 3 * Natoms
        f_mwc(i,j) = f_cart(i,j) / sqrt(mass_vector(i) * mass_vector(j))
    enddo
enddo
write(*,*) "mass-weighted cartesian hessian:"
call write_matrix(f_mwc)

! diagonalize the mass-weighted hessian:
f_diag = f_mwc
lwork = -1
allocate(work(1))
call dsyev('N', 'U', 3 * Natoms, f_diag, 3 * Natoms, freqs, work, lwork, info)
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
f_diag = f_mwc
call dsyev('N', 'U', 3 * Natoms, f_diag, 3 * Natoms, freqs, work, lwork, info)
write(*,*) "status of diagonalization: ", info

! conversion from atomic units to cm-1:
!  nu_tilde = 1 / (2 pi c 100 cm/m) * sqrt(lambda * Eh / (a0**2 u))
write(*,*) "eigenvalues:"
do i = 1, 3 * Natoms
    write(*,'(1x, es15.7, f10.4)') freqs(i), sqrt(abs(freqs(i)) * 9.375829435e29_dp) / 1.883651567e11_dp
enddo

! calculate the center of mass and shift the molecule:
call calc_com(positions, masses, total_mass, com)
write(*,*) "center of mass:"
call write_vector(com)
do i = 1, Natoms
    positions(3*(i-1)+1:3*i) = positions(3*(i-1)+1:3*i) - com
enddo

! calculate and diagonalize the inertia tensor and rotate the molecule:
call calc_inert(positions, mass_vector, total_mass, inert)
write(*,*) "inertia tensor:"
call write_matrix(inert)
prinaxes = inert
call dsyev('V', 'U', 3, prinaxes, 3, moments, work, lwork, info)
write(*,*) "status of diagonalization: ", info
write(*,*) "moments:"
call write_vector(moments)
write(*,*) "principal axes:"
call write_matrix(prinaxes)
do i = 1, Natoms
    positions(3*(i-1)+1:3*i) = matmul(transpose(prinaxes), positions(3*(i-1)+1:3*i))
enddo

! write updated properties
call calc_com(positions, masses, total_mass, com)
call calc_inert(positions, masses, total_mass, inert)
write(*,*) "new positions:"
call write_vector(positions)
write(*,*) "new center of mass:"
call write_vector(com)
write(*,*) "new inertia tensor:"
call write_matrix(inert)

! calculate the D matrix
D = 0.0_dp
do i = 1, 3 * Natoms
    do j = 7, 3 * Natoms
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
    D(3*(i-1)+2,4) = -dot_product(positions(3*i-2:3*i), prinaxes(3,:)) * sqrt(masses(3*(i-1)+2))
    D(3*(i-1)+3,4) =  dot_product(positions(3*i-2:3*i), prinaxes(2,:)) * sqrt(masses(3*(i-1)+3))
    D(3*(i-1)+1,5) =  dot_product(positions(3*i-2:3*i), prinaxes(3,:)) * sqrt(masses(3*(i-1)+1))
    D(3*(i-1)+2,5) =  0.0_dp
    D(3*(i-1)+3,5) = -dot_product(positions(3*i-2:3*i), prinaxes(1,:)) * sqrt(masses(3*(i-1)+3))
    D(3*(i-1)+1,6) = -dot_product(positions(3*i-2:3*i), prinaxes(2,:)) * sqrt(masses(3*(i-1)+1))
    D(3*(i-1)+2,6) =  dot_product(positions(3*i-2:3*i), prinaxes(1,:)) * sqrt(masses(3*(i-1)+2))
    D(3*(i-1)+3,6) =  0.0_dp
enddo

write(*,*) "the D matrix:"
call write_matrix(D)

! normalize the D matrix:
do i = 1, 3 * Natoms
    D(:,i) = D(:,i) / norm2(D(:,i))
enddo

write(*,*) "the normalized D matrix:"
call write_matrix(D)
metric = matmul(transpose(D), D)
write(*,*) "the metric of the D matrix:"
call write_matrix(metric)

! orthogonalize the D-matrix:
call gram_schmidt(D)

write(*,*) "the orthogonalized D matrix:"
call write_matrix(D)
metric = matmul(transpose(D), D)
write(*,*) "the metric of the orthogonalized D matrix:"
call write_matrix(metric)

! transform the hessian to internal coordinates:
f_int = matmul(transpose(D), matmul(f_mwc, D))
write(*,*) "the hessian in internal coordinates:"
call write_matrix(f_int)


! diagonalize the internal hessian:
f_diag = f_int
call dsyev('V', 'U', 3 * Natoms, f_diag, 3 * Natoms, freqs, work, lwork, info)
write(*,*) "status of diagonalization: ", info

write(*,*) "eigenvalues:"
do i = 1, 3 * Natoms
    write(*,'(1x, es15.7, f10.4)') freqs(i), sqrt(abs(freqs(i)) * 9.375829435e29_dp) / 1.883651567e11_dp
enddo

write(*,*) "eigenvectors of the hessian:"
call write_matrix(f_diag)

test_hessian = f_int(7:9,7:9)
write(*,*) "test hessian:"
call write_matrix(test_hessian)
call dsyev('V', 'U', 3, test_hessian, 3, moments, work, 30, info)
write(*,*) "status of diagonalization: ", info

write(*,*) "eigenvalues:"
do i = 1, 3
    write(*,'(1x, es15.7, f10.4)') moments(i), sqrt(abs(moments(i)) * 9.375829435e29_dp) / 1.883651567e11_dp
enddo

write(*,*) "eigenvectors of the hessian:"
call write_matrix(test_hessian)




end program gauss_viban
