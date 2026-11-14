!------------------------------------------------------------------------------
! shtools_grid_bench.f90: SHTOOLS (Fortran 95) driver of the grid synthesis benchmark.
! See benchmark/grid/code/README.md for the definitions, the build and the limits.
!
! Cases (one case per run):
!   dh1       Driscoll-Healy grid, n x n with n = 2*nmax+2, MakeGridDH(sampling=1)
!   dh2       equiangular Driscoll-Healy grid, n x 2n, MakeGridDH(sampling=2)
!   gl        NOT RUN: SHTOOLS MakeGridGLQ places 2*lmax+1 longitudes (spacing 360/(2*lmax+1)), so it does not
!             produce the contract Gauss-Legendre grid (2*nmax+2 longitudes, spacing 360/(2*nmax+2)). Exit 4.
!   regional  NOT RUN: SHTOOLS has no local grid synthesis. Exit 3.
!
! Grid rows run from 90 deg (row 0) southward, spacing 180/n, as documented in MakeGridDH.f95 (SHTOOLS 4.14.1):
!   lat_i = 90 - i*180/n,  lon_j = j*360/nlong,  i = 0..n-1, j = 0..nlong-1 (0-based here, geocentric latitude).
!
! Definitions (benchmark/shs/CONTRACT.md):
!   V_m2s2 = GM/R * MakeGridDH(cilm)   (cilm fully normalised, 4 pi, no Condon-Shortley phase; r = R, so the
!                                       radius scaling (R/r)^n is 1; degree 0 included, C00 as in the file)
!   T_m2s2 = V - U(R, lat), U the GRS80 normal potential, closed form (shs_common.normal_potential_closed).
!
! Timing: the wall time (omp_get_wtime) covers MakeGridDH only. The GM/R scaling and U are computed after the
! timed repetitions. Peak memory: VmHWM of /proc/self/status (Linux/WSL), in MB; -1 when not available.
! MakeGridDH.f95 has no OpenMP parallel region (only threadprivate declarations), so it runs on one thread.
!
! Output (never a whole grid):
!   <out-dir>/shtools_<case>_rows.csv  row,lat_deg,lon_deg,V_m2s2,T_m2s2 ; row is 0-based, 0 = northernmost
!   <out-dir>/timings.csv              appended, engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb
!
! Usage:
!   shtools_grid_bench --gfc F --nmax N --case dh1|dh2 --threads K --repeat R --out-dir DIR [--rows i,j,k]
!
! Exit codes: 0 ok, 2 bad arguments or model mismatch, 3 case regional, 4 case gl.
!
! References: Wieczorek and Meschede (2018), SHTools (MakeGridDH, Driscoll and Healy 1994 sampling);
! Heiskanen and Moritz (1967), Physical Geodesy (normal potential).
!------------------------------------------------------------------------------
program shtools_grid_bench
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32
    use omp_lib, only: omp_get_wtime
    use SHTOOLS, only: MakeGridDH
    implicit none

    real(dp), parameter :: PI = 3.141592653589793238_dp
    real(dp), parameter :: GRS80_A = 6378137.0_dp
    real(dp), parameter :: GRS80_E2 = 0.00669438002290_dp
    real(dp), parameter :: GRS80_GM = 3.986005e14_dp
    real(dp), parameter :: GRS80_OMEGA = 7.292115e-5_dp
    integer(int32), parameter :: MAX_SEL = 64

    character(len=1024) :: gfc_path, out_dir, grid_case, rows_arg
    integer(int32) :: nmax, threads, nrep, sampling, n, nlong, nret, ios_grid
    logical :: have_gfc, have_nmax, have_case, have_out, have_rows

    real(dp) :: gm, radius
    integer(int32) :: max_degree
    real(dp), allocatable :: cilm(:,:,:), grid(:,:)
    real(dp), allocatable :: lat_deg(:), lon_deg(:), phi_c(:), u_row(:)
    integer(int32), allocatable :: sel(:)
    integer(int32) :: nsel
    real(dp) :: t0, wall, peak
    integer(int32) :: rep, i, k
    character(len=32) :: recurrence

    call parse_arguments()
    if (nrep < 1 .or. threads < 1 .or. nmax < 0) then
        write(*,'(A)') 'shtools_grid_bench: repeat and threads must be >= 1 and nmax >= 0'
        stop 2
    end if

    if (trim(grid_case) == 'regional') then
        write(*,'(A)') 'shtools_grid_bench: case regional is not run by this driver: SHTOOLS has no local &
                       &(regional) grid synthesis.'
        stop 3
    end if
    if (trim(grid_case) == 'gl') then
        write(*,'(A)') 'shtools_grid_bench: case gl is not run by this driver: SHTOOLS MakeGridGLQ gives 2*lmax+1&
                       & longitudes (spacing 360/(2*lmax+1)), not the contract grid with 2*nmax+2 longitudes.'
        stop 4
    end if
    if (trim(grid_case) == 'dh1') then
        sampling = 1
    else if (trim(grid_case) == 'dh2') then
        sampling = 2
    else
        write(*,'(A)') 'shtools_grid_bench: --case must be dh1, dh2, gl or regional'
        stop 2
    end if

    ! ---- model (not timed) ----
    call read_gfc_header(trim(gfc_path), gm, radius, max_degree)
    if (nmax > max_degree) then
        write(*,'(A,I0,A,I0)') 'shtools_grid_bench: nmax exceeds max_degree of the model: ', nmax, ' > ', max_degree
        stop 2
    end if
    call read_gfc_coefficients(trim(gfc_path), nmax, cilm)

    ! ---- grid geometry (MakeGridDH.f95: n = 2*lmax+2, nlong = sampling*n) ----
    n = 2 * nmax + 2
    nlong = sampling * n
    allocate(grid(n, nlong))
    grid = 0.0_dp
    allocate(lat_deg(n), lon_deg(nlong), phi_c(n), u_row(n))
    do i = 1, n
        lat_deg(i) = 90.0_dp - real(i - 1, dp) * 180.0_dp / real(n, dp)
        phi_c(i) = lat_deg(i) * PI / 180.0_dp
    end do
    do k = 1, nlong
        lon_deg(k) = real(k - 1, dp) * 360.0_dp / real(nlong, dp)
    end do

    ! ---- timed repetitions: MakeGridDH only ----
    call omp_check()
    recurrence = 'native'
    do rep = 1, nrep
        t0 = omp_get_wtime()
        call MakeGridDH(grid, nret, cilm, nmax, norm=1, sampling=sampling, csphase=1, exitstatus=ios_grid)
        wall = omp_get_wtime() - t0
        if (ios_grid /= 0 .or. nret /= n) then
            write(*,'(A,I0,A,I0)') 'shtools_grid_bench: MakeGridDH failed, exitstatus=', ios_grid, ' n=', nret
            stop 1
        end if
        peak = peak_mem_mb()
        call append_timing(trim(out_dir), trim(grid_case), trim(recurrence), threads, nmax, n, nlong, rep, wall, peak)
        write(*,'(A,I0,A,F12.6,A)') 'repeat ', rep, ': wall=', wall, ' s'
    end do

    ! ---- V = GM/R * S (in place), U per latitude row, selected rows ----
    grid = grid * (gm / radius)
    do i = 1, n
        u_row(i) = normal_potential_closed(radius, phi_c(i))
    end do

    allocate(sel(MAX_SEL))
    nsel = 0
    call select_rows(n)
    call write_rows(trim(out_dir), trim(grid_case), sel, nsel, grid, lat_deg, lon_deg, u_row, n, nlong)
    write(*,'(A,A,A,I0,A,I0,A,I0,A,I0)') 'shtools_grid_bench: ', trim(grid_case), ' nmax=', nmax, ' nlat=', n, &
        ' nlon=', nlong, ' rows written=', nsel

contains

    ! ------------------------------------------------------------------
    ! Arguments
    ! ------------------------------------------------------------------
    subroutine parse_arguments()
        character(len=1024) :: arg, val
        integer :: nargs, kk, ios
        have_gfc = .false.; have_nmax = .false.; have_case = .false.; have_out = .false.; have_rows = .false.
        threads = 1; nrep = 1; rows_arg = ''
        nargs = command_argument_count()
        kk = 1
        do while (kk <= nargs)
            call get_command_argument(kk, arg)
            if (kk + 1 > nargs) then
                write(*,'(A,A)') 'missing value for ', trim(arg)
                stop 2
            end if
            call get_command_argument(kk + 1, val)
            select case (trim(arg))
            case ('--gfc')
                gfc_path = val
                have_gfc = .true.
            case ('--nmax')
                read(val, *, iostat=ios) nmax
                have_nmax = (ios == 0)
            case ('--case')
                grid_case = val
                have_case = .true.
            case ('--threads')
                read(val, *) threads
            case ('--repeat')
                read(val, *) nrep
            case ('--out-dir')
                out_dir = val
                have_out = .true.
            case ('--rows')
                rows_arg = val
                have_rows = .true.
            case default
                write(*,'(A,A)') 'unknown argument ', trim(arg)
                stop 2
            end select
            kk = kk + 2
        end do
        if (.not. (have_gfc .and. have_nmax .and. have_case .and. have_out)) then
            write(*,'(A)') 'usage: shtools_grid_bench --gfc F --nmax N --case dh1|dh2|gl|regional --threads K &
                           &--repeat R --out-dir DIR [--rows i,j,k]'
            stop 2
        end if
    end subroutine parse_arguments

    subroutine omp_check()
        ! Printed once so that the log shows the thread count the driver was started with.
        write(*,'(A,I0)') 'shtools_grid_bench: threads requested (informational, MakeGridDH is serial): ', threads
    end subroutine omp_check

    ! ------------------------------------------------------------------
    ! Model header (same keys as the SHS drivers) and coefficients
    ! ------------------------------------------------------------------
    subroutine read_gfc_header(path, gm_out, r_out, maxdeg)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: gm_out, r_out
        integer(int32), intent(out) :: maxdeg
        character(len=512) :: line, key, value
        integer :: u, ios, p
        logical :: have_gm, have_r, have_max, ended

        have_gm = .false.; have_r = .false.; have_max = .false.; ended = .false.
        open(newunit=u, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A)') 'shtools_grid_bench: cannot open the .gfc file'
            stop 2
        end if
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:11) == 'end_of_head') then
                ended = .true.
                exit
            end if
            line = adjustl(line)
            p = index(line, ' ')
            if (p == 0) cycle
            key = line(1:p-1)
            value = adjustl(line(p+1:))
            select case (trim(key))
            case ('earth_gravity_constant', 'gravity_constant')
                read(value, *) gm_out
                have_gm = .true.
            case ('radius')
                read(value, *) r_out
                have_r = .true.
            case ('max_degree')
                read(value, *) maxdeg
                have_max = .true.
            end select
        end do
        close(u)
        if (.not. ended) then
            write(*,'(A)') 'shtools_grid_bench: no end_of_head line in the model file'
            stop 2
        end if
        if (.not. (have_gm .and. have_r .and. have_max)) then
            write(*,'(A)') 'shtools_grid_bench: model header lacks earth_gravity_constant, radius or max_degree'
            stop 2
        end if
    end subroutine read_gfc_header

    subroutine read_gfc_coefficients(path, nmax_in, c)
        character(len=*), intent(in) :: path
        integer(int32), intent(in) :: nmax_in
        real(dp), allocatable, intent(out) :: c(:,:,:)
        character(len=512) :: line
        integer(int32) :: n, m, nread
        real(dp) :: cval, sval
        integer :: u, ios

        allocate(c(2, nmax_in + 1, nmax_in + 1))
        c = 0.0_dp
        nread = 0
        open(newunit=u, file=path, status='old', action='read')
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:3) /= 'gfc' .or. line(4:4) /= ' ') cycle
            read(line(4:), *) n, m, cval, sval
            if (n > nmax_in .or. m > n) cycle
            c(1, n + 1, m + 1) = cval
            c(2, n + 1, m + 1) = sval
            nread = nread + 1
        end do
        close(u)
        if (nread /= (nmax_in + 1) * (nmax_in + 2) / 2) then
            write(*,'(A,I0,A,I0)') 'coefficient count ', nread, ' does not match nmax ', nmax_in
            stop 2
        end if
    end subroutine read_gfc_coefficients

    ! ------------------------------------------------------------------
    ! GRS80 normal potential (shs_common.normal_potential_closed), at r and geocentric latitude phi_c
    ! ------------------------------------------------------------------
    real(dp) function q_func(x, e_lin)
        real(dp), intent(in) :: x, e_lin
        q_func = 0.5_dp * ((1.0_dp + 3.0_dp * x * x / (e_lin * e_lin)) * atan(e_lin / x) - 3.0_dp * x / e_lin)
    end function q_func

    real(dp) function normal_potential_closed(r, phi_c_in)
        real(dp), intent(in) :: r, phi_c_in
        real(dp) :: e_lin, b_semi, z, d, u2, u, sin2_beta, flattening
        e_lin = GRS80_A * sqrt(GRS80_E2)
        b_semi = GRS80_A * sqrt(1.0_dp - GRS80_E2)
        z = r * sin(phi_c_in)
        d = r * r - e_lin * e_lin
        u2 = 0.5_dp * (d + sqrt(d * d + 4.0_dp * e_lin * e_lin * z * z))
        u = sqrt(u2)
        sin2_beta = z * z / u2
        flattening = 0.5_dp * GRS80_OMEGA * GRS80_OMEGA * GRS80_A * GRS80_A * &
                     (q_func(u, e_lin) / q_func(b_semi, e_lin)) * (sin2_beta - 1.0_dp / 3.0_dp)
        normal_potential_closed = GRS80_GM / e_lin * atan(e_lin / u) + flattening
    end function normal_potential_closed

    ! ------------------------------------------------------------------
    ! Peak memory (VmHWM, MB) and timings
    ! ------------------------------------------------------------------
    real(dp) function peak_mem_mb()
        character(len=256) :: line
        real(dp) :: kb
        integer :: u, ios
        peak_mem_mb = -1.0_dp
        open(newunit=u, file='/proc/self/status', status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:6) == 'VmHWM:') then
                read(line(7:), *, iostat=ios) kb
                if (ios == 0) peak_mem_mb = kb / 1024.0_dp
                exit
            end if
        end do
        close(u)
    end function peak_mem_mb

    subroutine append_timing(dir, gcase, rec, nthreads, nm, nl_, nlo, rep_i, wall_s, peak_mb)
        character(len=*), intent(in) :: dir, gcase, rec
        integer(int32), intent(in) :: nthreads, nm, nl_, nlo, rep_i
        real(dp), intent(in) :: wall_s, peak_mb
        character(len=1100) :: path
        logical :: exists
        integer :: u, ios

        path = trim(dir) // '/timings.csv'
        inquire(file=trim(path), exist=exists)
        open(newunit=u, file=trim(path), status='unknown', position='append', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(A)') 'shtools_grid_bench: cannot append to timings.csv (does --out-dir exist?)'
            stop 1
        end if
        if (.not. exists) write(u, '(A)') 'engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb'
        write(u, '(A,A,A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES20.10E3,A,F12.3)') &
            'shtools,', trim(gcase), ',', trim(rec), ',', nthreads, ',', nm, ',', nl_, ',', nlo, ',', rep_i, ',', &
            wall_s, ',', peak_mb
        close(u)
    end subroutine append_timing

    ! ------------------------------------------------------------------
    ! Row selection (1-based internally): global rule of the benchmark, or --rows (0-based on the command line)
    ! ------------------------------------------------------------------
    integer(int32) function nearest_row(target_deg)
        real(dp), intent(in) :: target_deg
        integer(int32) :: ii
        real(dp) :: dmin, dd
        nearest_row = 1
        dmin = huge(1.0_dp)
        do ii = 1, n
            dd = abs(lat_deg(ii) - target_deg)
            if (dd < dmin) then
                dmin = dd
                nearest_row = ii
            end if
        end do
    end function nearest_row

    subroutine add_row(v)
        integer(int32), intent(in) :: v
        integer(int32) :: kk
        do kk = 1, nsel
            if (sel(kk) == v) return
        end do
        if (nsel >= MAX_SEL) then
            write(*,'(A)') 'shtools_grid_bench: too many selected rows'
            stop 2
        end if
        nsel = nsel + 1
        sel(nsel) = v
    end subroutine add_row

    subroutine select_rows(nlat)
        integer(int32), intent(in) :: nlat
        real(dp) :: targets(9)
        integer(int32) :: kk, jj, tmp, v, ios, p
        character(len=1024) :: rest, token

        if (have_rows) then
            rest = rows_arg
            do
                p = index(rest, ',')
                if (p == 0) then
                    token = trim(rest)
                else
                    token = rest(1:p-1)
                end if
                read(token, *, iostat=ios) v
                if (ios /= 0 .or. v < 0 .or. v >= nlat) then
                    write(*,'(A)') 'shtools_grid_bench: --rows must be a comma-separated list of indices in [0, nlat-1]'
                    stop 2
                end if
                call add_row(v + 1)
                if (p == 0) exit
                rest = rest(p+1:)
            end do
        else
            targets = [0.0_dp, 30.0_dp, 60.0_dp, 85.0_dp, 89.5_dp, -30.0_dp, -60.0_dp, -85.0_dp, -89.5_dp]
            do kk = 1, size(targets)
                call add_row(nearest_row(targets(kk)))
            end do
            call add_row(1)
            call add_row(nlat)
        end if
        ! insertion sort, ascending
        do kk = 2, nsel
            jj = kk
            do while (jj > 1)
                if (sel(jj - 1) <= sel(jj)) exit
                tmp = sel(jj - 1)
                sel(jj - 1) = sel(jj)
                sel(jj) = tmp
                jj = jj - 1
            end do
        end do
    end subroutine select_rows

    subroutine write_rows(dir, gcase, s, ns, g, lat, lon, ur, nlat, nlo)
        character(len=*), intent(in) :: dir, gcase
        integer(int32), intent(in) :: s(:), ns, nlat, nlo
        real(dp), intent(in) :: g(:,:), lat(:), lon(:), ur(:)
        character(len=1100) :: path
        integer :: u, ios, kk, jj, ii
        real(dp) :: vv, tt

        path = trim(dir) // '/shtools_' // trim(gcase) // '_rows.csv'
        open(newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,'(A)') 'shtools_grid_bench: cannot write the rows CSV (does --out-dir exist?)'
            stop 1
        end if
        write(u, '(A)') 'row,lat_deg,lon_deg,V_m2s2,T_m2s2'
        do kk = 1, ns
            ii = s(kk)
            do jj = 1, nlo
                vv = g(ii, jj)
                tt = vv - ur(ii)
                write(u, '(I0,A,ES26.16E3,A,ES26.16E3,A,ES26.16E3,A,ES26.16E3)') &
                    ii - 1, ',', lat(ii), ',', lon(jj), ',', vv, ',', tt
            end do
        end do
        close(u)
    end subroutine write_rows

end program shtools_grid_bench
