module routines
implicit none

integer,parameter :: dp = selected_real_kind(14)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine searchstring(iounit, keyword, back, lerr)
! Go to the first instance of "keyword" found in file "id"
implicit none

integer,intent(in) :: iounit            ! IO unit of the file to search
character(len=*),intent(in) :: keyword  ! the string to search for
logical,intent(in) :: back              ! wether to start from the beginning
logical,intent(inout) :: lerr           ! status

character(len=128) :: line
integer :: kwidx, ierr
logical :: stat

lerr = .false.
ierr = 0

if (back) rewind(iounit)

inquire(iounit,opened=stat,iostat=ierr)

if (.not. stat) return

do while (ierr == 0)
    read(iounit, '(A)', iostat=ierr) line
    kwidx = index(line, keyword)
    if (kwidx /= 0) then
        backspace(iounit)
        lerr = .true.
        return
    endif
enddo

lerr = .false.

end subroutine searchstring

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readchkmat(id, dim1, mat)
implicit none

integer,intent(in) :: id                            ! IO unit to read from
integer,intent(in) :: dim1                          ! dimensionality of the matrix
real(dp),intent(inout),dimension(dim1,dim1) :: mat  ! the matrix

integer :: i, k                                     ! loop indices

k = 1
do i = 1, dim1
    call readchkvec(id, i, mat(1:i,i), k)
    mat(i,1:i-1) = mat(1:i-1,i)
enddo

end subroutine readchkmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readchkvec(id, dim1, vec, k)
implicit none

integer,intent(in) :: id                        ! IO unit to read from
integer,intent(in) :: dim1                      ! length of the vector
real(dp),intent(inout),dimension(dim1) :: vec   ! the vector
integer,intent(inout) :: k

integer :: j

j = 1
do while (j <= dim1)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_com(x, m, mtot, com)
!
implicit none

real(dp),intent(in),dimension(:) :: x    ! positions of the atoms
real(dp),intent(in),dimension(:) :: m    ! masses of the atoms
real(dp),intent(in) :: mtot              ! total mass
real(dp),intent(out),dimension(3) :: com ! center of mass

integer :: i

com = 0.0_dp
do i = 1, size(m)
    com(1) = com(1) + (m(i) * x(3*(i-1)+1))
    com(2) = com(2) + (m(i) * x(3*(i-1)+2))
    com(3) = com(3) + (m(i) * x(3*(i-1)+3))
enddo
com = com / mtot

end subroutine calc_com

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_inert(x, m, mtot, inert)
!
implicit none

real(dp),intent(in),dimension(:) :: x        ! positions of the atoms
real(dp),intent(in),dimension(:) :: m        ! masses of the atoms
real(dp),intent(in) :: mtot                  ! total mass
real(dp),intent(out),dimension(3,3) :: inert ! inertia tensor

integer :: i, j, k
real(dp) :: factor
real(dp),dimension(3) :: com

call calc_com(x, m, mtot, com)

do i = 1, 3
    do j = 1, 3
        inert(i,j) = 0.0_dp
        do k = 1, size(m)
            factor = 0.0_dp
            if (i == j) then
                factor = norm2(x(3*(k-1)+1:3*(k-1)+3) - com)
            endif
            factor = factor - (x(3*(k-1)+i) - com(i)) * (x(3*(k-1)+j) - com(j))
            inert(i,j) = inert(i,j) + (m(k) * factor)
        enddo
    enddo
enddo

end subroutine calc_inert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_vector(vector, dig)
!
implicit none

real(dp),intent(in),dimension(:) :: vector  ! the vector to write
integer,intent(in),optional :: dig          ! number of decimal places

integer :: i
character(len=20) :: format_string

format_string = '(1x, '
do i = 1, size(vector)
    write(*,'(1x, es15.8)') vector(i)
enddo

end subroutine write_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_matrix(matrix, dig)
!
implicit none

real(dp),intent(in),dimension(:,:) :: matrix    ! the matrix to write
integer,intent(in),optional :: dig              ! number of decimal places

integer :: i, j

do i = 1, size(matrix, 1)
    do j = 1, size(matrix, 2)
        write(*,'(1x, es15.8)',advance='no') matrix(i,j)
    enddo
    write(*,*)
enddo

end subroutine write_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gram_schmidt(matrix)
!
! this subroutine will orthogonalize the column vectors of the matrix 'matrix'
! by the Gram-Schmidt-Method. The first column vector will remain the same.
! in this special case we start with the seventh column
!
implicit none

! arguments to the routine
real(dp),dimension(:,:),intent(inout) :: matrix     ! the matrix to be orthogonalized
! internal variables
integer :: i, j                                     ! loop indices
integer :: n                                        ! the dimension of the matrix
real(dp) :: factor                                  ! temp variable to scale the vectors
real(dp),dimension(:),allocatable :: temp_vec       ! temporary vector
real(dp),dimension(:,:),allocatable :: new_matrix   ! the new orthogonal matrix

n = size(matrix, 1)
allocate(temp_vec(n))
allocate(new_matrix(n,n))

new_matrix = 0.0
new_matrix(:,1:6) = matrix(:,1:6)

do i = 7, n
    temp_vec = matrix(:,i)
    do j = 1, i-1
        factor = dot_product(temp_vec, new_matrix(:,j))
        temp_vec = temp_vec - factor * new_matrix(:,j)
    enddo
    factor = norm2(temp_vec)
    temp_vec = temp_vec / factor
    new_matrix(:,i) = temp_vec
enddo

matrix = new_matrix

end subroutine gram_schmidt

end module routines
