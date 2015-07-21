module miquel_routines
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

end module miquel_routines
