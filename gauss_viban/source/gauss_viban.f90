program readprogram
use miquel_routines
implicit none

integer,parameter :: id1 = 1            ! IO unit of the input file
integer,parameter :: id2 = 2            ! IO unit of the input file
integer,parameter :: dim1 = 9           ! length of the coordinate vector
integer :: i, j, info
real(8) :: total_mass                   ! total mass of the molecule
real(8),dimension(3) :: com             ! center of mass
real(8),dimension(dim1) :: masses       ! the atomic masses
real(8),dimension(dim1) :: positions    ! the atomic positions
real(8),dimension(dim1) :: eigvals
real(8),dimension(3,3) :: inert         ! the inertia tensor
real(8),dimension(dim1,dim1) :: f_cart  ! the hessian in cartesian coordinates
real(8),dimension(dim1,dim1) :: f_mwc   ! the hessian in mass weighted cart. coordinates
real(8),dimension(dim1,dim1) :: f_diag  ! the hessian in mass weighted cart. coordinates
real(8),dimension(30) :: work


! define the atomic masses:
masses(1) = 15.9949146
masses(2) = 15.9949146
masses(3) = 15.9949146
masses(4) = 1.00782504
masses(5) = 1.00782504
masses(6) = 1.00782504
masses(7) = 1.00782504
masses(8) = 1.00782504
masses(9) = 1.00782504
total_mass = 15.9949146 + 2 * 1.00782504
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
    write(*,'(1x, es15.7)') positions(i)
enddo
write(*,*) "non-mass-weighted initial hessian:"
do i = 1, dim1
    do j = 1, dim1
        write(*,'(1x, es15.8)',advance='no') f_cart(i,j)
    enddo
    write(*,*)
enddo



! calculate the mass-weighted hessian:
do i = 1, dim1
    do j = 1, dim1
        f_mwc(i,j) = f_cart(i,j) / sqrt(masses(i) * masses(j))
    enddo
enddo


! diagonalize the mass-weighted hessian:
f_diag = f_mwc
call dsyev('N', 'U', dim1, f_diag, dim1, eigvals, work, 30, info)
write(*,*) "status of diagonalization: ", info

write(*,*) "non-mass-weighted initial hessian:"
do i = 1, dim1
    do j = 1, dim1
        write(*,'(1x, es15.8)',advance='no') f_cart(i,j)
    enddo
    write(*,*)
enddo

write(*,*) "mass-weighted initial hessian:"
do i = 1, dim1
    do j = 1, dim1
        write(*,'(1x, es15.8)',advance='no') f_mwc(i,j)
    enddo
    write(*,*)
enddo

! conversion from atomic units to cm-1:
!  nu_tilde = 1 / (2 pi c 100 cm/m) * sqrt(lambda * Eh / (a0**2 u))
write(*,*) "eigenvalues:"
do i = 1, dim1
    write(*,'(1x, 2es15.7)') eigvals(i), sqrt(eigvals(i) * 9.375829435e29) / 1.883651567e11
enddo


! calculate the center of mass
com(1) = (positions(1) * masses(1) + positions(4) * masses(4) + positions(7) * masses(7)) / total_mass
com(2) = (positions(2) * masses(2) + positions(5) * masses(5) + positions(8) * masses(8)) / total_mass
com(3) = (positions(3) * masses(3) + positions(6) * masses(6) + positions(9) * masses(9)) / total_mass
write(*,*) "center of mass:"
do i = 1, 3
    write(*,'(1x, es15.7)') com(i)
enddo


! calculate the inertia tensor:
inert(1,1) = masses(1) * (positions(2)**2 + positions(3)**2) &
           + masses(4) * (positions(5)**2 + positions(6)**2) &
           + masses(7) * (positions(8)**2 + positions(9)**2)
inert(2,2) = masses(1) * (positions(1)**2 + positions(3)**2) &
           + masses(4) * (positions(4)**2 + positions(6)**2) &
           + masses(7) * (positions(7)**2 + positions(9)**2)
inert(3,3) = masses(1) * (positions(1)**2 + positions(2)**2) &
           + masses(4) * (positions(4)**2 + positions(5)**2) &
           + masses(7) * (positions(7)**2 + positions(8)**2)
inert(1,2) = -(masses(1) * positions(1) * positions(2) &
             + masses(4) * positions(4) * positions(5) &
             + masses(7) * positions(7) * positions(8))
inert(1,3) = -(masses(1) * positions(1) * positions(3) &
             + masses(4) * positions(4) * positions(6) &
             + masses(7) * positions(7) * positions(9))
inert(2,3) = -(masses(1) * positions(2) * positions(3) &
             + masses(4) * positions(5) * positions(6) &
             + masses(7) * positions(8) * positions(9))
inert(2,1) = inert(1,2)
inert(3,1) = inert(1,3)
inert(3,2) = inert(2,3)
write(*,*) "inertia tensor:"
do i = 1, 3
    do j = 1, 3
        write(*,'(1x, es15.8)',advance='no') inert(i,j)
    enddo
    write(*,*)
enddo









end program readprogram
