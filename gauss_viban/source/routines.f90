module routines
implicit none

contains


subroutine readchkmat(id, dim1, mat)
implicit none

integer,intent(in) :: id                            ! IO unit to read from
integer,intent(in) :: dim1                          ! dimensionality of the matrix
real(8),intent(inout),dimension(dim1,dim1) :: mat   ! the matrix

integer :: i, k                                     ! loop indices

k = 1
do i = 1, dim1
    call readchkvec(id, i, mat(1:i,i), k)
    mat(i,1:i-1) = mat(1:i-1,i)
enddo

end subroutine readchkmat

subroutine readchkvec(id, dim1, vec, k)
implicit none

integer,intent(in) :: id                        ! IO unit to read from
integer,intent(in) :: dim1                      ! length of the vector
real(8),intent(inout),dimension(dim1) :: vec    ! the vector
integer,intent(inout) :: k

integer :: j

j = 1
do while (j <= dim1)
    !write(*,*) j, k
    read(id,'(1x, e15.8)',advance='no') vec(j)
    if (mod(k, 5) == 0) then
        read(id,*)
        k=1
    else
        k = k + 1
    endif
    j = j + 1
enddo

end subroutine readchkvec

subroutine calc_com(x, m, mtot, com)
!
implicit none

real(8),intent(in),dimension(:) :: x    ! positions of the atoms
real(8),intent(in),dimension(:) :: m    ! masses of the atoms
real(8),intent(in) :: mtot              ! total mass
real(8),intent(out),dimension(3) :: com ! center of mass

integer :: i

com = 0.0
do i = 1, size(x), 3
    com(1) = com(1) + (m(i) * x(i))
    com(2) = com(2) + (m(i) * x(i+1))
    com(3) = com(3) + (m(i) * x(i+2))
enddo
com = com / mtot

end subroutine calc_com

subroutine calc_inert(x, m, mtot, inert)
!
implicit none

real(8),intent(in),dimension(:) :: x        ! positions of the atoms
real(8),intent(in),dimension(:) :: m        ! masses of the atoms
real(8),intent(in) :: mtot                  ! total mass
real(8),intent(out),dimension(3,3) :: inert ! inertia tensor

integer :: i, j, k
real(8) :: factor
real(8),dimension(3) :: com

call calc_com(x, m, mtot, com)

do i = 1, 3
    do j = 1, 3
        inert(i,j) = 0.0
        do k = 1, size(x), 3
            factor = 0.0
            if (i == j) then
                factor = norm2(x(k:k+2) - com)
            endif
            factor = factor - (x(k+i-1) - com(i)) * (x(k+j-1) - com(j))
            inert(i,j) = inert(i,j) + (m(k) * factor)
        enddo
    enddo
enddo

end subroutine calc_inert

subroutine write_matrix(matrix)
!
implicit none

real(8),intent(in),dimension(:,:) :: matrix

integer :: i, j

do i = 1, size(matrix, 1)
    do j = 1, size(matrix, 2)
        write(*,'(1x, es15.8)',advance='no') matrix(i,j)
    enddo
    write(*,*)
enddo

end subroutine write_matrix

end module routines
