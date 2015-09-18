! Copyright 2015 Jan von Cosel
!
! This file is part of utility-scripts.
!
! utility-scripts is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! utility-scripts is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have recieved a copy of the GNU General Public License
! along with utility-scripts. If not, see <http://www.gnu.org/licenses/>.

program mctdh_spec_fft
implicit none

integer,parameter :: dp = selected_real_kind(14)    ! double precision kind type
integer,parameter :: inputio = 1                    ! IO unit for the input file
integer,parameter :: specio = 2                     ! IO unit for the autocorrelation output
integer,parameter :: autoio = 3                     ! IO unit for the spectrum output
real(dp),parameter :: fs2au = 41.34137333656_dp     ! conversion between femtoseconds and atomic units of time
real(dp),parameter :: pi = 3.1415926535897932_dp    ! PI
integer :: nlines                                   ! number of lines (= time points) in the input
integer :: N                                        ! number of time points
integer :: stat                                     ! IO status
integer :: i, m, k                                  ! loop index
real(dp) :: tau                                     ! damping time
real(dp) :: auto_real                               ! dummy variable for the real part of auto
real(dp) :: auto_imag                               ! dummy variable for the imaginary part of auto
real(dp) :: auto_abs                                ! dummy variable for the absolute value of auto
real(dp) :: dt                                      ! the time interval
real(dp) :: df                                      ! the frequency interval
real(dp) :: nyquist                                 ! the Nyquist frequency
real(dp) :: tmax                                    ! the maximum time
real(dp),dimension(:),allocatable :: time           ! the time axis
real(dp),dimension(:),allocatable :: double_time    ! the double time axis
real(dp),dimension(:),allocatable :: double_freq    ! the double frequency axis
complex(dp),dimension(:),allocatable :: auto_inp    ! the autocorrelation function read from the file
complex(dp),dimension(:),allocatable :: auto_double ! the mirrored autocorrelation function
complex(dp),dimension(:),allocatable :: spec        ! the spectrum
complex(dp),dimension(:),allocatable :: spec_shift  ! the shifted spectrum
character(len=50) :: inputfile                      ! input file with the autocorrelation function
character(len=50) :: cmd_arg                        ! dummy string for the command line arguments


call get_command_argument(1, inputfile)
if (len_trim(inputfile) == 0) then
    write(*,*) "ERROR opening input file, ", inputfile
    stop
endif
open(unit=inputio,file=inputfile,status='old',action='read')

call get_command_argument(2, cmd_arg)
read(cmd_arg,*,iostat=stat) nlines
if (stat /= 0) then
    write(*,*) "ERROR reading nlines, status:", stat
    stop
endif

call get_command_argument(3, cmd_arg)
read(cmd_arg,*,iostat=stat) tau
if (stat /= 0) then
    write(*,*) "ERROR reading tau, status:", stat
    stop
endif
tau = tau * fs2au

N = 1
do while (N < nlines)
    N = N * 2
enddo

allocate(time(N))
allocate(double_time(2 * N))
allocate(double_freq((2 * N) + 1))
allocate(auto_inp(N))
allocate(auto_double(2 * N))
allocate(spec(2 * N))
allocate(spec_shift((2 * N) + 1))

! skip the first line of the input:
read(inputio,*)
auto_inp = (0.0_dp, 0.0_dp)
do i = 1, nlines
    read(inputio,*) time(i), auto_real, auto_imag, auto_abs
    auto_inp(i) = dcmplx(auto_real, auto_imag)
enddo
close(inputio)

! calculate the required times and frequencies:
dt = (time(2) - time(1)) * fs2au
tmax = dt * real(N - 1, dp)
nyquist = pi / dt
df = nyquist / real(N, dp)
do i = 1, (2 * N) - 1
    double_time(i) = -tmax + real(i - 1, dp) * dt
enddo
do i = 1, (2 * N) + 1
    double_freq(i) = -nyquist + real(i - 1, dp) * df
enddo

! construct the mirrored autocorrelation function with cosine weighting and exponential damping:
auto_double = (0.0_dp, 0.0_dp)
do i = 0, N - 1
    auto_double(i + 1) = dconjg(auto_inp(N - i))    ! complex conjugate the negative part!!!
    auto_double(2 * N - (i + 1)) = auto_inp(N - i)
enddo
auto_double = auto_double * cos(pi * double_time / (2.0_dp * tmax)) * exp(-abs(double_time) / tau)

! perform discrete fourier transformation:
! (Numerical Recipes p. 497, eq 12.1.7)
do m = 0, (2 * N) - 1
    spec(m + 1) = (0.0_dp, 0.0_dp)
    do k = 0, (2 * N) - 1
        spec(m + 1) = spec(m + 1) + auto_double(k + 1) * exp(2.0_dp * pi * (0.0_dp, 1.0_dp) * real(k*m, dp) / real((2*N)-1, dp))
    enddo
enddo

spec_shift(1:N) = spec(N+1:2*N)
spec_shift(N+1:(2*N)+1) = spec(1:N+1)

! write the results to the output files:
open(unit=autoio,file="auto_out.dat",status='replace',action='write')
open(unit=specio,file="spec_out.dat",status='replace',action='write')

write(autoio,'("# Time [au] Re(Auto) Im(Auto)")')
write(specio,'("# Frequency [au] Re(Spec) Im(Spec)")')

do i = 1, 2 * N
    write(autoio,'(f20.10, 2f16.10)') double_time(i), real(auto_double(i)), aimag(auto_double(i))
enddo
do i = 1, (2 * N) + 1
    write(specio,'(f15.10, 2f17.10)') double_freq(i), real(spec_shift(i)), aimag(spec_shift(i))
enddo
close(autoio)
close(specio)

end program mctdh_spec_fft
