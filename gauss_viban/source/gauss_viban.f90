program readprogram
use routines
implicit none

integer,parameter :: id1 = 1            ! IO unit of the input file
integer,parameter :: id2 = 2            ! IO unit of the input file
integer,parameter :: dim1 = 9           ! length of the coordinate vector
integer :: i, j, info
real(dp) :: total_mass                  ! total mass of the molecule
real(dp),dimension(3) :: com            ! center of mass
real(dp),dimension(3) :: moments        ! eigenvalues of inertia tensor
real(dp),dimension(dim1) :: masses      ! the atomic masses
real(dp),dimension(dim1) :: positions   ! the atomic positions
real(dp),dimension(dim1) :: eigvals
real(dp),dimension(3,3) :: inert        ! the inertia tensor
real(dp),dimension(3,3) :: prinaxes     ! eigenvectors of inertia tensor
real(dp),dimension(3,3) :: test_hessian ! 
real(dp),dimension(dim1,dim1) :: f_cart ! the hessian in cartesian coordinates
real(dp),dimension(dim1,dim1) :: f_mwc  ! the hessian in mass weighted cart. coordinates
real(dp),dimension(dim1,dim1) :: f_diag ! diagonalized hessian
real(dp),dimension(dim1,dim1) :: f_int  ! the hessian in internal coodinates

real(dp),dimension(dim1,dim1) :: D      ! the D matrix
real(dp),dimension(dim1,dim1) :: metric ! the metric matrix

real(dp),dimension(30) :: work


! define the atomic masses:
masses(1:3) = 15.9949146_dp    ! oxygen
masses(4:9) = 1.00782504_dp    ! hydrogen
total_mass = masses(1) + masses(4) + masses(7)
write(*,*) "total mass: ", total_mass


! read the cartesian hessian from the file:
open(unit=id1,file='pos_data',status='old',action='read')
i = 1
call readchkvec(id1, dim1, positions, i)
close(unit=id1)
open(unit=id2,file='fc_data',status='old',action='read')
call readchkmat(id2, dim1, f_cart)
close(unit=id2)

write(*,*) "positions:"
do i = 1, dim1
    write(*,'(1x, f8.4)') positions(i)
enddo
write(*,*) "non-mass-weighted initial cartesian hessian:"
call write_matrix(f_cart)


! calculate the mass-weighted hessian:
do i = 1, dim1
    do j = 1, dim1
        f_mwc(i,j) = f_cart(i,j) / sqrt(masses(i) * masses(j))
    enddo
enddo
write(*,*) "mass-weighted cartesian hessian:"
call write_matrix(f_mwc)


! diagonalize the mass-weighted hessian:
f_diag = f_mwc
call dsyev('N', 'U', dim1, f_diag, dim1, eigvals, work, 30, info)
write(*,*) "status of diagonalization: ", info

! conversion from atomic units to cm-1:
!  nu_tilde = 1 / (2 pi c 100 cm/m) * sqrt(lambda * Eh / (a0**2 u))
write(*,*) "eigenvalues:"
do i = 1, dim1
    write(*,'(1x, es15.7, f10.4)') eigvals(i), sqrt(abs(eigvals(i)) * 9.375829435e29_dp) / 1.883651567e11_dp
enddo


! calculate the center of mass
call calc_com(positions, masses, total_mass, com)
write(*,*) "center of mass:"
do i = 1, 3
    write(*,'(1x, es15.7)') com(i)
enddo


! calculate the inertia tensor:
call calc_inert(positions, masses, total_mass, inert)
write(*,*) "inertia tensor:"
call write_matrix(inert)


! diagonalize the inertia tensor
prinaxes = inert
call dsyev('V', 'U', 3, prinaxes, 3, moments, work, 30, info)
write(*,*) "status of diagonalization: ", info
write(*,*) "moments:"
do i = 1, 3
    write(*,'(1x, es15.8)') moments(i)
enddo
write(*,*) "principal axes:"
call write_matrix(prinaxes)


! translate the molecule to the com:
do i = 1, dim1 / 3
    positions(3*(i-1)+1:3*i) = positions(3*(i-1)+1:3*i) - com
enddo


! rotate the molecule to the inertia frame:
do i = 1, dim1 / 3
    positions(3*(i-1)+1:3*i) = matmul(transpose(prinaxes), positions(3*(i-1)+1:3*i))
enddo


write(*,*) "new positions:"
do i = 1, dim1
    write(*,'(1x, f8.4)') positions(i)
enddo
call calc_com(positions, masses, total_mass, com)
write(*,*) "new center of mass:"
do i = 1, 3
    write(*,'(1x, es15.7)') com(i)
enddo
call calc_inert(positions, masses, total_mass, inert)
write(*,*) "new inertia tensor:"
call write_matrix(inert)


! calculate the D matrix
D = 0.0_dp
do i = 1, dim1
    D(i,i) = 1.0_dp
enddo
do i = 1, dim1
    do j = 7, dim1
        call random_number(D(i,j))
    enddo
enddo

D(1,1) = sqrt(masses(1))
D(4,1) = sqrt(masses(4))
D(7,1) = sqrt(masses(7))
D(2,2) = sqrt(masses(2))
D(5,2) = sqrt(masses(5))
D(8,2) = sqrt(masses(8))
D(3,3) = sqrt(masses(3))
D(6,3) = sqrt(masses(6))
D(9,3) = sqrt(masses(9))

do i = 1, dim1 / 3  ! sum over atom
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

!write(*,*) "the D matrix:"
!call write_matrix(D)

! normalize the D matrix
do i = 1, dim1
    D(:,i) = D(:,i) / norm2(D(:,i))
enddo

write(*,*) "the normalized D matrix:"
call write_matrix(D)
metric = matmul(transpose(D), D)
write(*,*) "the metric matrix:"
call write_matrix(metric)


! orthogonalize the D-matrix:
call gram_schmidt(D)

write(*,*) "the orthogonalized D matrix:"
call write_matrix(D)
metric = matmul(transpose(D), D)
write(*,*) "the metric matrix:"
call write_matrix(metric)


! transform the hessian to internal coordinates:
f_int = matmul(transpose(D), matmul(f_mwc, D))
!f_int = matmul(D, matmul(f_mwc, transpose(D)))
write(*,*) "the hessian in internal coordinates:"
call write_matrix(f_int)


! diagonalize the internal hessian:
f_diag = f_int
call dsyev('V', 'U', dim1, f_diag, dim1, eigvals, work, 30, info)
write(*,*) "status of diagonalization: ", info

write(*,*) "eigenvalues:"
do i = 1, dim1
    write(*,'(1x, es15.7, f10.4)') eigvals(i), sqrt(abs(eigvals(i)) * 9.375829435e29_dp) / 1.883651567e11_dp
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




end program readprogram
