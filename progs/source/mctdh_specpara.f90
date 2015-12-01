program mctdh_specpara
use routines
implicit none

integer,parameter :: gsfcio = 1                     ! IO unit of the file containing the GS frequencies
integer,parameter :: esfcio = 2                     ! IO unit of the file containing the ES frequencies
integer,parameter :: shiftio = 3                    ! IO unit of the file containing the shift vector (k)
integer,parameter :: duschio = 4                    ! IO unit of the file containing the Duschinsky matrix (J)
real(dp),parameter :: pi = 3.141592653589793238_dp  ! pi
real(dp),parameter :: c0 = 299792458.0_dp           ! the speed of light
real(dp),parameter :: Eh = 4.3597438e-18_dp         ! the Hartree energy
real(dp),parameter :: a0 = 5.291772083e-11_dp       ! the bohr radius
!real(dp),parameter :: u  = 1.66053873e-27_dp        ! the atomic mass unit
!real(dp),parameter :: hbar = 1.054571596e-34_dp     ! Planck's constant by 2 Pi
real(dp),parameter :: me = 9.10938188e-31_dp        ! the electron mass
!real(dp),parameter :: amu2au = 1822.88848325_dp     ! to convert from amu to atomic units
integer :: n, m, o                                  ! loop indices
integer :: Nmodes                                   ! number of normal modes
integer :: stat                                     ! IO status of the read operation
integer :: Nspfs                                    ! number of SPFs
integer :: spfs_written                             ! number of modes written so far
integer,dimension(2) :: max_indices                 ! location of the maximum value of the coupling array
logical :: line_started
character(len=50) :: argument                       ! the command argument
real(dp),dimension(:),allocatable :: v1             ! ground state frequencies
real(dp),dimension(:),allocatable :: v2             ! excited state frequencies
real(dp),dimension(:),allocatable :: f1             ! ground state force constants
real(dp),dimension(:),allocatable :: f2             ! excited state force constants
real(dp),dimension(:),allocatable :: k              ! shift vector
real(dp),dimension(:),allocatable :: fp             ! new force constants
real(dp),dimension(:),allocatable :: kappa          ! first-order coefficients
real(dp),dimension(:,:),allocatable :: J            ! Duschinsky matrix
real(dp),dimension(:,:),allocatable :: phi          ! coupling matrix
real(dp),dimension(:,:),allocatable :: phi_sort     ! sorted coupling matrix (for mode combination)
logical,dimension(:),allocatable :: written         ! has the mode been written to sbasis?


! the number of modes must be given as the argument:
call get_command_argument(1, argument)
if (len_trim(argument) == 0) then
    write(*,*) "ERROR: no argument given"
    stop
endif
read(argument,*,iostat=stat) Nmodes
if (stat /= 0) then
    write(*,*) "ERROR reading Nmodes, status:", stat
    stop
endif

! allocate the arrays:
allocate(v1(Nmodes))
allocate(v2(Nmodes))
allocate(f1(Nmodes))
allocate(f2(Nmodes))
allocate(k(Nmodes))
allocate(fp(Nmodes))
allocate(kappa(Nmodes))
allocate(written(Nmodes))
allocate(J(Nmodes,Nmodes))
allocate(phi(Nmodes,Nmodes))

! open the data files:
open(unit=gsfcio,file='gs_freqs',status='old',action='read')
open(unit=esfcio,file='es_freqs',status='old',action='read')
open(unit=shiftio,file='Displacement_Vector.dat',status='old',action='read')
open(unit=duschio,file='Duschinsky_Matrix.dat',status='old',action='read')

! read the data:
do n = 1, Nmodes
    read(gsfcio,*) v1(n)
enddo
do n = 1, Nmodes
    read(esfcio,*) v2(n)
enddo
do n = 1, Nmodes
    read(shiftio,*) k(n)
enddo
do n = 1, Nmodes
    do m = 1, Nmodes
        read(duschio,*) J(m,n)
    enddo
enddo

! calculate the force constants:
f1 = v1**2 * 4.0_dp * pi**2 * c0**2 * 10000.0_dp * a0**2 * me / Eh
f2 = v2**2 * 4.0_dp * pi**2 * c0**2 * 10000.0_dp * a0**2 * me / Eh

! write the data:
write(*,*) "ground state frequencies"
call write_vector(v1)
write(*,*) "excited state frequencies"
call write_vector(v2)
write(*,*) "ground state force constants"
call write_vector(f1)
write(*,*) "excited state force constants"
call write_vector(f2)
write(*,*) "shift vector"
call write_vector(k)
write(*,*) "Duschinsky matrix"
call write_matrix(J)

! calculate the new force constants:
do m = 1, Nmodes
    fp(m) = 0.0_dp
    do n = 1, Nmodes
        fp(m) = fp(m) + f2(n) * J(n,m)**2
    enddo
enddo

! calculate the first-order coefficients:
do m = 1, Nmodes
    kappa(m) = 0.0_dp
    do n = 1, Nmodes
        kappa(m) = kappa(m) + f2(n) * k(n) * J(n,m)
    enddo
enddo

! calculate the couplings:
phi = 0.0_dp
do m = 1, Nmodes
    do o = m + 1, Nmodes
        phi(m,o) = 0.0_dp
        do n = 1, Nmodes
            phi(m,o) = phi(m,o) + f2(n) * J(n,m) * J(n,o)
        enddo
    enddo
enddo

! write the results:
write(*,*) "new force constants"
call write_vector(fp)
write(*,*) "first-order coefficients"
call write_vector(kappa)
write(*,*) "couplings"
call write_matrix(phi)

!
! write MCTDH input/operator
!

! sbasis-section:
phi_sort = abs(phi)
written = .false.
write(*,'("sbasis-section")')
Nspfs = Nmodes / 2
do n = 1, Nspfs
    write(*,'("    ")',advance='no')
    spfs_written = 0
    write_loop: do
        max_indices = maxloc(phi_sort)
        phi_sort(max_indices(1), max_indices(2)) = 0.0_dp
        do m = 1, 2
            if (.not. written(max_indices(m))) then
                write(*,'("q_",i3.3)',advance='no') max_indices(m)
                spfs_written = spfs_written + 1
                written(max_indices(m)) = .true.
                if (spfs_written < 2) then
                    write(*,'(", ")',advance='no')
                endif
                if (spfs_written == 2) then
                    write(*,'("  =  2")')
                    exit write_loop
                endif
            endif
        enddo
    enddo write_loop
enddo
line_started = .false.
do n = 1, Nmodes
    if (.not. written(n)) then
        if (line_started) then
            write(*,'(", ")',advance='no')
        else
            write(*,'("    ")',advance='no')
        endif
        write(*,'("q_",i3.3)',advance='no') n
        line_started = .true.
    endif
enddo
if (line_started) then
    write(*,'("  =  2")')
endif
write(*,'("end-sbasis-section")')
write(*,*)

! pbasis-section:
write(*,'("pbasis-section")')
do n = 1, Nmodes
    write(*,'("    q_",i3.3,"   ho   10   xi-xf", 2f8.1)') &
        n, -(kappa(n) / fp(n)) - (3.3_dp / fp(n)**0.25_dp), -(kappa(n) / fp(n)) + (3.3_dp / fp(n)**0.25_dp)
enddo
write(*,'("end-pbasis-section")')
write(*,*)

! init_wf-section:
write(*,'("init_wf-section")')
write(*,'("    build")')
do n = 1, Nmodes
    write(*,'("        q_",i3.3,"  eigenf  Eq_",i3.3,"  pop = 1")') n, n
enddo
write(*,'("    end-build")')
write(*,'("end-init_wf-section")')

! parameter-section:
write(*,'("parameter-section")')
do n = 1, Nmodes
    write(*,'("    mass_q_",i3.3,"  =  1.0")') n
enddo
do n = 1, Nmodes
    write(*,'("    f_",i3.3,"       = ",d15.8)') n, f1(n)
enddo
do n = 1, Nmodes
    write(*,'("    fp_",i3.3,"      = ",d15.8)') n, fp(n)
enddo
do n = 1, Nmodes
    do m = n + 1, Nmodes
        write(*,'("    phi_",i3.3,"_",i3.3," = ",d15.8)') n, m, phi(n,m)
    enddo
enddo
do n = 1, Nmodes
    write(*,'("    kappa_",i3.3,"   = ",d15.8)') n, kappa(n)
enddo
write(*,'("end-parameter-section")')
write(*,*)

! hamiltonian-section:
write(*,'("hamiltonian-section")')
write(*,'("modes")',advance='no')
do n = 1, Nmodes
    write(*,'(" | q_",i3.3)',advance='no') n
enddo
write(*,*)
do n = 1, Nmodes
    write(*,'("1.0         |",i0," KE")') n
enddo
do n = 1, Nmodes
    write(*,'("0.5*fp_",i3.3,"  |",i0," q^2")') n, n
enddo
do n = 1, Nmodes
    do m = n + 1, Nmodes
        write(*,'("phi_",i3.3,"_",i3.3," |",i0," q |",i0," q")') n, m, n, m
    enddo
enddo
do n = 1, Nmodes
    write(*,'("kappa_",i3.3,"   |",i0," q")') n, n
enddo
write(*,'("end-hamiltonian-section")')
write(*,*)

! additional hamiltonians for the GS modes:
do n = 1, Nmodes
    write(*,'("hamiltonian-section_Eq_",i3.3)') n
    write(*,'("usediag")')
    write(*,'("modes     | q_",i3.3)') n
    write(*,'("1.0       |1 KE")')
    write(*,'("0.5*f_",i3.3," |1 q^2")') n
    write(*,'("end-hamiltonian-section")')
    write(*,*)
enddo

write(*,*) "###################################################################"

! mode couplings in descending order:
phi_sort = abs(phi)

do n = 1, (Nmodes * (Nmodes - 1)) / 2
    max_indices = maxloc(phi_sort)
    write(*,*) max_indices(1), max_indices(2), phi(max_indices(1), max_indices(2))
    phi_sort(max_indices(1), max_indices(2)) = 0.0_dp
enddo



end program mctdh_specpara
