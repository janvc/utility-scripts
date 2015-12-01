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

subroutine readgaussmodes(id, dim1, modes)
!
implicit none

integer,intent(in) :: id                                ! IO unit to read from
integer,intent(in) :: dim1                              ! number of coordinates (3N, not modes!!)
real(dp),intent(out),dimension(dim1,dim1-6) :: modes    ! the array to write the modes into

integer :: i, k

k = 1
do i = 1, dim1 - 6
    call readchkvec(id, dim1, modes(:,i), k)
enddo

end subroutine readgaussmodes

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

integer :: i
real(dp),dimension(3) :: com

call calc_com(x, m, mtot, com)

inert = 0.0_dp
do i = 1, size(m)
    inert(1,1) = inert(1,1) + m(i) * ((x(3*(i-1)+2)-com(2))**2 + (x(3*(i-1)+3)-com(3))**2)
    inert(2,2) = inert(2,2) + m(i) * ((x(3*(i-1)+1)-com(1))**2 + (x(3*(i-1)+3)-com(3))**2)
    inert(3,3) = inert(3,3) + m(i) * ((x(3*(i-1)+1)-com(1))**2 + (x(3*(i-1)+2)-com(2))**2)
    inert(1,2) = inert(1,2) - m(i) * (x(3*(i-1)+1)-com(1)) * (x(3*(i-1)+2)-com(2))
    inert(1,3) = inert(1,3) - m(i) * (x(3*(i-1)+1)-com(1)) * (x(3*(i-1)+3)-com(3))
    inert(2,3) = inert(2,3) - m(i) * (x(3*(i-1)+2)-com(2)) * (x(3*(i-1)+3)-com(3))
enddo

inert(2,1) = inert(1,2)
inert(3,1) = inert(1,3)
inert(3,2) = inert(2,3)

end subroutine calc_inert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_vector(vector, dig, clean)
!
implicit none

real(dp),intent(in),dimension(:) :: vector  ! the vector to write
integer,intent(in),optional :: dig          ! number of decimal places
logical,intent(in),optional :: clean        ! write small numbers as zero

integer :: i                        ! loop index
integer :: fw                       ! field width
real(dp) :: threshold               ! the threshold for setting vector elements to 0
character(len=12) :: format_string  ! the format string to write the numbers with
character(len=8) :: zero_string     ! the format string for the zero

threshold = norm2(vector) / 1.0e10_dp

if (present(dig)) then
    fw = dig + 7
    write(zero_string,'("(",i2,"x,i1)")') fw
    if (fw < 10) then
        write(format_string,'("(1x,  es",i1,".",i1,")")') fw, dig
        write(zero_string,'("( ",i1,"x,i1)")') fw
    else if (dig > 10) then
        write(format_string,'("(1x,es",i2,".",i2,")")') fw, dig
    else
        write(format_string,'("(1x, es",i2,".",i1,")")') fw, dig
    endif
else
    fw = 15
    format_string = '(1x, es15.8)'
endif

do i = 1, fw + 1
    write(*,'("-")',advance='no')
enddo
write(*,*)

do i = 1, size(vector)
    if (present(clean) .and. abs(vector(i)) < threshold) then
        if (clean) then
            write(*,zero_string) 0
        else
            write(*,format_string) vector(i)
        endif
    else
        write(*,format_string) vector(i)
    endif
enddo

do i = 1, fw + 1
    write(*,'("-")',advance='no')
enddo
write(*,*)

end subroutine write_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_matrix(matrix, dig, clean)
!
implicit none

real(dp),intent(in),dimension(:,:) :: matrix    ! the matrix to write
integer,intent(in),optional :: dig              ! number of decimal places
logical,intent(in),optional :: clean            ! write small numbers as zero

integer :: i, j                     ! loop indices
integer :: fw                       ! field width
real(dp) :: threshold               ! the threshold for setting matrix elements to 0
character(len=12) :: format_string  ! the format string to write the numbers with
character(len=8) :: zero_string     ! the format string for the zero

threshold = matrix_norm(matrix) / 1.0e10_dp

if (present(dig)) then
    fw = dig + 7
    write(zero_string,'("(",i2,"x,i1)")') fw
    if (fw < 10) then
        write(format_string,'("(1x,  es",i1,".",i1,")")') fw, dig
        write(zero_string,'("( ",i1,"x,i1)")') fw
    else if (dig > 10) then
        write(format_string,'("(1x,es",i2,".",i2,")")') fw, dig
    else
        write(format_string,'("(1x, es",i2,".",i1,")")') fw, dig
    endif
else
    fw = 15
    format_string = '(1x, es15.8)'
endif

do i = 1, size(matrix, 2)
    do j = 1, fw + 1
        write(*,'("-")',advance='no')
    enddo
enddo
write(*,*)

do i = 1, size(matrix, 1)
    do j = 1, size(matrix, 2)
        if (present(clean) .and. abs(matrix(i,j)) < threshold) then
            if (clean) then
                write(*,zero_string,advance='no') 0
            else
                write(*,format_string,advance='no') matrix(i,j)
            endif
        else
            write(*,format_string,advance='no') matrix(i,j)
        endif
    enddo
    write(*,*)
enddo

do i = 1, size(matrix, 2)
    do j = 1, fw + 1
        write(*,'("-")',advance='no')
    enddo
enddo
write(*,*)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp) function matrix_norm(matrix)
!
! this function calculates the Frobenius norm of a matrix
!
implicit none

! arguments to the routine:
real(dp),dimension(:,:),intent(in) :: matrix    ! the matrix to calculate the norm of

! internal variables:
integer :: i, j     ! loop indices

matrix_norm = 0.0_dp

do i = 1, size(matrix, 1)
    do j = 1, size(matrix, 2)
        matrix_norm = matrix_norm + matrix(i,j)**2
    enddo
enddo

matrix_norm = sqrt(matrix_norm)

end function matrix_norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ortho_trans(M, U, A)
!
! this subroutine performs the orthogonal transformation A = U^T * M * U
! all matrices have to be quadratic and equal in size
!
implicit none

! arguments to the routine:
real(dp),dimension(:,:),intent(in) :: M     ! the matrix to be transformed
real(dp),dimension(:,:),intent(in) :: U     ! the transformation matrix
real(dp),dimension(:,:),intent(out) :: A    ! the resulting matrix

! internal variables:
integer :: i, j, k, l   ! loop indices
integer :: dim1         ! dimension of the participating matrices

dim1 = size(M, 1)

do i = 1, dim1
    do j = 1, dim1
        do k = 1, dim1
            do l = 1, dim1
                A(i,j) = A(i,j) + U(k,i) * M(k,l) * U(l,j)
            enddo
        enddo
    enddo
enddo

end subroutine ortho_trans

end module routines
