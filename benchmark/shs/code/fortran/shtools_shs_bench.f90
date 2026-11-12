!------------------------------------------------------------------------------
! SHTOOLS (Fortran 95) driver for the spherical-harmonic synthesis benchmark.
!
! Engine 'shtools' (definitions: benchmark/shs/results/csv.md). It mirrors code/python/run_pyshtools.py:
! the same model (tide-free .gfc, fully normalised 4pi), the same radius scaling, the same
! gradient, the same geodetic -> geocentric conversion, U and gamma_SH. The formulas are
! ported line by line from code/python/shs_common.py; read that file for the reasoning.
!
! One run evaluates, at each point, V, T and dg in one pass and writes one result CSV:
!   id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal   (15 significant digits)
! plus a companion .md, and appends one row per repetition to the timings CSV.
!
! Library vs driver parallelism (checked in SHTOOLS 4.14.1 sources):
!   * MakeGridPoint and MakeGravGridPoint contain no OpenMP directives: one point is serial.
!   * Their Legendre helpers (PlmBar, PlmBar_d1) keep saved work arrays declared
!     THREADPRIVATE under OpenMP, so the OpenMP build (libSHTOOLS-mp) can be called from
!     several threads at once. The driver relies on this and checks it (see README).
!   * The driver parallelises the point loops with OpenMP: V over groups of points that share
!     a radius (each thread owns a radius-scaled copy of cilm), dg over points.
!
! Usage:
!   shtools_shs_bench --gfc <tide_free.gfc> --points <csv> --nmax N --out <result.csv>
!       --timings <timings.csv> --set colorado|latsweep --threads K --repeat R
!------------------------------------------------------------------------------
program shtools_shs_bench
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32, compiler_version
    use omp_lib
    use SHTOOLS, only: MakeGridPoint, MakeGravGridPoint
    implicit none

    ! Constants: same values as shs_common.py (GRS80) and pi as numpy.
    real(dp), parameter :: PI = 3.141592653589793238_dp
    real(dp), parameter :: GRS80_A = 6378137.0_dp
    real(dp), parameter :: GRS80_E2 = 0.00669438002290_dp
    real(dp), parameter :: GRS80_GM = 3.986005e14_dp
    real(dp), parameter :: GRS80_J2 = 1.08263e-3_dp
    real(dp), parameter :: GRS80_OMEGA = 7.292115e-5_dp
    real(dp), parameter :: MGAL_PER_MS2 = 1.0e5_dp
    character(len=*), parameter :: FUNCTIONAL_TAG = 'V_T_dg'
    character(len=*), parameter :: POINT_HEADER = 'id,lat_deg,lon_deg,h_m'
    character(len=*), parameter :: RESULT_HEADER = &
        'id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal'
    ! Library build flags. build_shtools_driver.sh replaces this default in its copy of the source.
    character(len=*), parameter :: LIB_FLAGS = "not recorded (compiled without build_shtools_driver.sh)"
    character(len=*), parameter :: TIMING_HEADER = &
        'engine,recurrence,set,nmax,n_points,threads,repeat,wall_s,cpu_s,notes'

    ! Arguments
    character(len=1024) :: gfc_path, points_path, out_path, timings_path, set_name
    integer(int32) :: nmax, threads, nrep

    ! Model
    character(len=256) :: modelname, tide_system, norm_key
    real(dp) :: gm, r0
    integer(int32) :: max_degree
    real(dp), allocatable :: cilm(:,:,:)

    ! Points
    integer(int32) :: npts
    character(len=64), allocatable :: ids(:)
    real(dp), allocatable :: lat_deg(:), lon_deg(:), h_m(:)
    real(dp), allocatable :: r(:), phi_c(:), lat_c_deg(:)

    ! Radius groups (points with identical r share one radius-scaled cilm)
    integer(int32), allocatable :: order(:), gstart(:)
    integer(int32) :: ngroups

    ! Functionals
    real(dp), allocatable :: v(:), vr(:), vlat(:), vlon(:)
    real(dp), allocatable :: t_vals(:), dg_vals(:)

    real(dp) :: t0, load_s, w0, wall, c0, c1, cpu
    integer(int32) :: rep, i
    character(len=256) :: compiler_str, build_flags
    character(len=1024) :: args_text

    call parse_arguments()
    if (nrep < 1 .or. threads < 1 .or. nmax < 0) then
        write(*,'(A)') 'repeat, threads must be >= 1 and nmax >= 0'
        stop 2
    end if
    call omp_set_num_threads(threads)

    ! ---- model load (timed separately, reported in notes) ----
    t0 = omp_get_wtime()
    call read_gfc_header(gfc_path, modelname, gm, r0, max_degree, norm_key, tide_system)
    if (nmax > max_degree) then
        write(*,'(A,I0,A,I0)') 'nmax exceeds max_degree of the model: ', nmax, ' > ', max_degree
        stop 2
    end if
    call read_gfc_coefficients(gfc_path, nmax, cilm)
    load_s = omp_get_wtime() - t0

    ! ---- points ----
    call read_points(points_path, npts, ids, lat_deg, lon_deg, h_m)
    allocate(r(npts), phi_c(npts), lat_c_deg(npts))
    do i = 1, npts
        call geodetic_to_geocentric(lat_deg(i), h_m(i), r(i), phi_c(i))
        lat_c_deg(i) = phi_c(i) * 180.0_dp / PI
    end do

    ! Group points by identical radius (same grouping as np.unique in the python runner).
    allocate(order(npts), gstart(npts + 1))
    call sort_by_radius(r, order)
    ngroups = 0
    do i = 1, npts
        if (i == 1) then
            ngroups = 1
            gstart(1) = 1
        else if (r(order(i)) /= r(order(i-1))) then
            ngroups = ngroups + 1
            gstart(ngroups) = i
        end if
    end do
    gstart(ngroups + 1) = npts + 1

    allocate(v(npts), vr(npts), vlat(npts), vlon(npts))

    ! ---- timed repetitions: synthesis only ----
    do rep = 1, nrep
        w0 = omp_get_wtime()
        call cpu_time(c0)
        call synth_potential(nmax, cilm, gm, r0, r, lat_c_deg, lon_deg, order, gstart, &
                             ngroups, v)
        call synth_gradient(nmax, cilm, gm, r0, r, lat_c_deg, lon_deg, vr, vlat, vlon, &
                            phi_c)
        call cpu_time(c1)
        cpu = c1 - c0
        wall = omp_get_wtime() - w0

        call get_build_info(compiler_str, build_flags)
        call append_timing(timings_path, rep, npts, cpu, wall, load_s, ngroups, &
                           compiler_str, build_flags)
        write(*,'(A,I0,A,F10.3,A,F10.3)') 'repeat ', rep, ': wall=', wall, ' s cpu=', cpu
    end do

    ! ---- functionals of the last repetition (same as the python runner) ----
    allocate(t_vals(npts), dg_vals(npts))
    do i = 1, npts
        t_vals(i) = v(i) - normal_potential_closed(r(i), phi_c(i))
    end do
    dg_vals = gravity_disturbance_mgal(vr, vlat, vlon, r, phi_c)

    call write_result_csv(out_path, npts, ids, lat_deg, lon_deg, h_m, v, t_vals, dg_vals)
    call write_result_md(out_path, npts, nmax, set_name, modelname, gm, r0, tide_system, &
                         load_s, ngroups, compiler_str, build_flags)

    write(*,'(A)') 'first rows:'
    do i = 1, min(3, npts)
        write(*,'(A)') trim(ids(i)) // ',' // trim(fmt15(lat_deg(i))) // ',' // &
            trim(fmt15(lon_deg(i))) // ',' // trim(fmt15(h_m(i))) // ',' // &
            trim(fmt15(v(i))) // ',' // trim(fmt15(t_vals(i))) // ',' // trim(fmt15(dg_vals(i)))
    end do

contains

    ! ------------------------------------------------------------------
    ! Arguments
    ! ------------------------------------------------------------------
    subroutine parse_arguments()
        character(len=1024) :: arg, val
        integer(int32) :: k, nargs
        logical :: have_gfc, have_points, have_out, have_timings, have_nmax, have_set

        have_gfc = .false.; have_points = .false.; have_out = .false.
        have_timings = .false.; have_nmax = .false.; have_set = .false.
        threads = 1; nrep = 1; nmax = -1; set_name = ''

        nargs = command_argument_count()
        args_text = ''
        do k = 1, nargs
            call get_command_argument(k, arg)
            args_text = trim(args_text) // ' ' // trim(arg)
        end do
        k = 1
        do while (k <= nargs)
            call get_command_argument(k, arg)
            if (k + 1 > nargs) then
                write(*,'(A,A)') 'missing value for ', trim(arg)
                stop 2
            end if
            call get_command_argument(k + 1, val)
            select case (trim(arg))
            case ('--gfc')
                gfc_path = val; have_gfc = .true.
            case ('--points')
                points_path = val; have_points = .true.
            case ('--out')
                out_path = val; have_out = .true.
            case ('--timings')
                timings_path = val; have_timings = .true.
            case ('--nmax')
                read(val, *) nmax; have_nmax = .true.
            case ('--threads')
                read(val, *) threads
            case ('--repeat')
                read(val, *) nrep
            case ('--set')
                set_name = val; have_set = .true.
            case default
                write(*,'(A,A)') 'unknown argument ', trim(arg)
                stop 2
            end select
            k = k + 2
        end do

        if (.not. (have_gfc .and. have_points .and. have_out .and. have_timings &
                   .and. have_nmax .and. have_set)) then
            write(*,'(A)') 'usage: shtools_shs_bench --gfc F --points CSV --nmax N --out CSV ' // &
                           '--timings CSV --set colorado|latsweep [--threads K] [--repeat R]'
            stop 2
        end if
        if (trim(set_name) /= 'colorado' .and. trim(set_name) /= 'latsweep') then
            write(*,'(A)') '--set must be colorado or latsweep'
            stop 2
        end if
    end subroutine parse_arguments

    ! ------------------------------------------------------------------
    ! Model: header (GM, R, max degree, norm, tide system) and coefficients
    ! ------------------------------------------------------------------
    subroutine read_gfc_header(path, name, gm_out, r_out, maxdeg, norm, tide)
        character(len=*), intent(in) :: path
        character(len=*), intent(out) :: name, norm, tide
        real(dp), intent(out) :: gm_out, r_out
        integer(int32), intent(out) :: maxdeg
        character(len=512) :: line, key, value
        integer :: u, ios, p
        logical :: have_gm, have_r, have_max, have_norm, have_tide

        have_gm = .false.; have_r = .false.; have_max = .false.
        have_norm = .false.; have_tide = .false.
        name = 'unknown'
        open(newunit=u, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,A)') 'cannot open model file ', trim(path)
            stop 2
        end if
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) then
                write(*,'(A)') 'no end_of_head line found in model file'
                stop 2
            end if
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            if (line(1:11) == 'end_of_head') exit
            p = scan(trim(line), ' '//char(9))
            if (p == 0) cycle
            key = line(1:p-1)
            value = adjustl(line(p+1:))
            select case (trim(key))
            case ('modelname')
                name = trim(value)
            case ('earth_gravity_constant', 'gravity_constant')
                gm_out = parse_real(value); have_gm = .true.
            case ('radius')
                r_out = parse_real(value); have_r = .true.
            case ('max_degree')
                read(value, *) maxdeg; have_max = .true.
            case ('norm')
                norm = trim(value); have_norm = .true.
            case ('tide_system')
                tide = trim(value); have_tide = .true.
            end select
        end do
        close(u)

        if (.not. (have_gm .and. have_r .and. have_max .and. have_norm .and. have_tide)) then
            write(*,'(A)') 'model header lacks a required key'
            stop 2
        end if
        if (trim(norm) /= 'fully_normalized') then
            write(*,'(A,A)') 'norm is not fully_normalized: ', trim(norm)
            stop 2
        end if
        if (trim(tide) /= 'tide_free') then
            write(*,'(A,A)') 'tide_system is not tide_free: ', trim(tide)
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
    ! Points CSV: id,lat_deg,lon_deg,h_m (same checks as shs_common.read_points)
    ! ------------------------------------------------------------------
    subroutine read_points(path, n, idv, lat, lon, h)
        character(len=*), intent(in) :: path
        integer(int32), intent(out) :: n
        character(len=64), allocatable, intent(out) :: idv(:)
        real(dp), allocatable, intent(out) :: lat(:), lon(:), h(:)
        character(len=512) :: line
        integer :: u, ios, nlines, p1, p2, p3

        open(newunit=u, file=path, status='old', action='read')
        read(u, '(A)') line
        if (trim(adjustl(line)) /= POINT_HEADER) then
            write(*,'(A,A)') 'points header must be ' // POINT_HEADER // ', got ', trim(line)
            stop 2
        end if
        nlines = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            nlines = nlines + 1
        end do
        rewind(u)
        read(u, '(A)') line
        allocate(idv(nlines), lat(nlines), lon(nlines), h(nlines))
        n = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            n = n + 1
            p1 = index(line, ',')
            p2 = index(line(p1+1:), ',') + p1
            p3 = index(line(p2+1:), ',') + p2
            idv(n) = adjustl(line(1:p1-1))
            lat(n) = parse_real(line(p1+1:p2-1))
            lon(n) = parse_real(line(p2+1:p3-1))
            h(n) = parse_real(line(p3+1:))
            if (abs(lat(n)) > 90.0_dp) then
                write(*,'(A)') 'latitude outside [-90, 90]'
                stop 2
            end if
        end do
        close(u)
        if (n < 1) stop 2
    end subroutine read_points

    ! ------------------------------------------------------------------
    ! Geometry (shs_common.geodetic_to_geocentric)
    ! ------------------------------------------------------------------
    subroutine geodetic_to_geocentric(lat_d, h_in, r_out, phic_out)
        real(dp), intent(in) :: lat_d, h_in
        real(dp), intent(out) :: r_out, phic_out
        real(dp) :: phi, sinphi, nrad, p, z

        phi = lat_d * PI / 180.0_dp
        sinphi = sin(phi)
        nrad = GRS80_A / sqrt(1.0_dp - GRS80_E2 * sinphi * sinphi)
        p = (nrad + h_in) * cos(phi)
        z = (nrad * (1.0_dp - GRS80_E2) + h_in) * sinphi
        r_out = sqrt(p * p + z * z)
        phic_out = atan2(z, p)
    end subroutine geodetic_to_geocentric

    ! Indices sorted by r (insertion sort; point counts are a few hundred to thousands).
    subroutine sort_by_radius(rv, idx)
        real(dp), intent(in) :: rv(:)
        integer(int32), intent(out) :: idx(:)
        integer(int32) :: i, j, key

        do i = 1, size(rv)
            idx(i) = i
        end do
        do i = 2, size(rv)
            key = idx(i)
            j = i - 1
            do while (j >= 1)
                if (rv(idx(j)) <= rv(key)) exit
                idx(j + 1) = idx(j)
                j = j - 1
            end do
            idx(j + 1) = key
        end do
    end subroutine sort_by_radius

    ! ------------------------------------------------------------------
    ! Normal field (shs_common.py): zonal series, gamma_SH, closed-form U
    ! ------------------------------------------------------------------
    real(dp) function ellipsoid_zonal(k)
        integer(int32), intent(in) :: k
        real(dp) :: c20, e2, sgn
        real(dp), parameter :: SQRT5 = 2.236067977499789696_dp

        if (k == 0) then
            ellipsoid_zonal = 1.0_dp
            return
        end if
        c20 = -GRS80_J2 / sqrt(5.0_dp)
        if (k == 1) then
            ellipsoid_zonal = c20
            return
        end if
        e2 = GRS80_E2
        if (mod(k, 2) == 0) then
            sgn = 1.0_dp
        else
            sgn = -1.0_dp
        end if
        ellipsoid_zonal = (sgn * 3.0_dp * e2 ** k / &
                          (real(2*k+1, dp) * real(2*k+3, dp) * sqrt(real(4*k+1, dp)))) &
                          * (1.0_dp - real(k, dp) - 5.0_dp * SQRT5 * real(k, dp) * c20 / e2)
    end function ellipsoid_zonal

    ! gamma_SH(r, phi') [m/s^2]: port of normal_gravity_gamma (degree 0..20 Legendre, zonal sum).
    real(dp) function normal_gravity_gamma(rr, phic)
        real(dp), intent(in) :: rr, phic
        real(dp) :: t, u, q, sr, sp, qn, dvr, dvp, dvcr, dvcp, gx, gz, nn, a1, b1
        real(dp) :: p0(0:20), p1(0:20)
        integer(int32) :: n, k
        real(dp) :: c

        t = sin(phic)
        u = cos(phic)
        if (abs(phic) == 0.5_dp * PI) then
            t = sign(1.0_dp, phic)
            u = 0.0_dp
        end if
        q = GRS80_A / rr

        p0(0) = 1.0_dp
        p0(1) = sqrt(3.0_dp) * t
        p1(0) = 0.0_dp
        p1(1) = sqrt(3.0_dp) * u
        do n = 2, 20
            nn = real(n, dp)
            p0(n) = sqrt(4.0_dp * nn * nn - 1.0_dp) / nn * t * p0(n-1) &
                    - (nn - 1.0_dp) / nn * sqrt((2.0_dp * nn + 1.0_dp) / (2.0_dp * nn - 3.0_dp)) * p0(n-2)
            a1 = sqrt((2.0_dp * nn + 1.0_dp) * (2.0_dp * nn - 1.0_dp) / ((nn - 1.0_dp) * (nn + 1.0_dp)))
            if (n >= 3) then
                b1 = sqrt((2.0_dp * nn + 1.0_dp) * nn * (nn - 2.0_dp) / &
                          ((nn - 1.0_dp) * (nn + 1.0_dp) * (2.0_dp * nn - 3.0_dp)))
            else
                b1 = 0.0_dp
            end if
            p1(n) = a1 * t * p1(n-1) - b1 * p1(n-2)
        end do

        sr = 0.0_dp
        sp = 0.0_dp
        qn = 1.0_dp
        do k = 0, 10
            n = 2 * k
            c = ellipsoid_zonal(k)
            nn = real(n, dp)
            sr = sr + qn * (nn + 1.0_dp) * c * p0(n)
            sp = sp + qn * c * sqrt(nn * (nn + 1.0_dp) / 2.0_dp) * p1(n)
            qn = qn * q * q
        end do
        dvr = -GRS80_GM / (rr * rr) * sr
        dvp = GRS80_GM / rr * sp
        dvcr = GRS80_OMEGA ** 2 * rr * u * u
        dvcp = -GRS80_OMEGA ** 2 * rr * rr * u * t
        gx = (dvp + dvcp) / rr
        gz = dvr + dvcr
        normal_gravity_gamma = sqrt(gx * gx + gz * gz)
    end function normal_gravity_gamma

    real(dp) function q_func(x, e_lin)
        real(dp), intent(in) :: x, e_lin
        q_func = 0.5_dp * ((1.0_dp + 3.0_dp * x * x / (e_lin * e_lin)) * atan(e_lin / x) &
                           - 3.0_dp * x / e_lin)
    end function q_func

    ! Gravitational part of the GRS80 normal potential, closed form (shs_common.normal_potential_closed).
    real(dp) function normal_potential_closed(rr, phic)
        real(dp), intent(in) :: rr, phic
        real(dp) :: e_lin, b_semi, z, d, u2, uu, sin2b, q0, flat

        e_lin = GRS80_A * sqrt(GRS80_E2)
        b_semi = GRS80_A * sqrt(1.0_dp - GRS80_E2)
        z = rr * sin(phic)
        d = rr * rr - e_lin * e_lin
        u2 = 0.5_dp * (d + sqrt(d * d + 4.0_dp * e_lin * e_lin * z * z))
        uu = sqrt(u2)
        sin2b = z * z / u2
        q0 = q_func(b_semi, e_lin)
        flat = 0.5_dp * GRS80_OMEGA ** 2 * GRS80_A ** 2 * (q_func(uu, e_lin) / q0) * &
               (sin2b - 1.0_dp / 3.0_dp)
        normal_potential_closed = GRS80_GM / e_lin * atan(e_lin / uu) + flat
    end function normal_potential_closed

    ! dg = |g| - gamma_SH in mGal (shs_common.gravity_disturbance_mgal), vectorised over points.
    function gravity_disturbance_mgal(vr_in, vlat_in, vlon_in, rr, phic) result(dg)
        real(dp), intent(in) :: vr_in(:), vlat_in(:), vlon_in(:), rr(:), phic(:)
        real(dp) :: dg(size(rr))
        real(dp) :: cosp, sinp, w2, vcr, vclat, gx, gy, gz, g
        integer(int32) :: i

        w2 = GRS80_OMEGA ** 2
        do i = 1, size(rr)
            cosp = cos(phic(i))
            sinp = sin(phic(i))
            vcr = w2 * rr(i) * cosp * cosp
            vclat = -w2 * rr(i) * rr(i) * cosp * sinp
            gx = (vlat_in(i) + vclat) / rr(i)
            gy = -vlon_in(i) / (rr(i) * cosp)
            gz = vr_in(i) + vcr
            g = sqrt(gx * gx + gy * gy + gz * gz)
            dg(i) = (g - normal_gravity_gamma(rr(i), phic(i))) * MGAL_PER_MS2
        end do
    end function gravity_disturbance_mgal

    ! ------------------------------------------------------------------
    ! Synthesis
    ! ------------------------------------------------------------------

    ! V = GM/r * MakeGridPoint(cilm scaled by (R/r)^n), one radius-scaled copy per radius group.
    ! Groups are distributed over threads; each thread owns its scaled copy.
    subroutine synth_potential(lmax, c, gmv, r0v, rv, latc, lon, idx, gs, ng, vout)
        integer(int32), intent(in) :: lmax, ng
        real(dp), intent(in) :: c(:,:,:), gmv, r0v, rv(:), latc(:), lon(:)
        integer(int32), intent(in) :: idx(:), gs(:)
        real(dp), intent(out) :: vout(:)
        real(dp), allocatable :: scaled(:,:,:), sclv(:)
        real(dp) :: rk
        integer(int32) :: g, k, l, m, j, p

        !$omp parallel private(scaled, sclv, g, rk, l, m, j, p)
        allocate(scaled(2, lmax + 1, lmax + 1), sclv(0:lmax))
        !$omp do schedule(dynamic, 1)
        do g = 1, ng
            rk = rv(idx(gs(g)))
            do l = 0, lmax
                sclv(l) = (r0v / rk) ** l
            end do
            do m = 1, lmax + 1
                do l = 1, lmax + 1
                    scaled(1, l, m) = c(1, l, m) * sclv(l - 1)
                    scaled(2, l, m) = c(2, l, m) * sclv(l - 1)
                end do
            end do
            do k = gs(g), gs(g + 1) - 1
                p = idx(k)
                vout(p) = gmv / rk * MakeGridPoint(scaled, lmax, latc(p), lon(p), 1, 1)
            end do
        end do
        !$omp end do
        deallocate(scaled, sclv)
        !$omp end parallel
    end subroutine synth_potential

    ! Gradient of V in the geocentric spherical frame, converted to dV/dr, dV/dphi', dV/dlambda
    ! exactly as run_pyshtools.py does (omega = 0: no centrifugal term here).
    subroutine synth_gradient(lmax, c, gmv, r0v, rv, latc, lon, vr_out, vlat_out, vlon_out, phic)
        integer(int32), intent(in) :: lmax
        real(dp), intent(in) :: c(:,:,:), gmv, r0v, rv(:), latc(:), lon(:), phic(:)
        real(dp), intent(out) :: vr_out(:), vlat_out(:), vlon_out(:)
        real(dp) :: b(3)
        integer(int32) :: i

        !$omp parallel do schedule(dynamic, 4) private(b)
        do i = 1, size(rv)
            b = MakeGravGridPoint(c, lmax, gmv, r0v, rv(i), latc(i), lon(i), 0.0_dp)
            vr_out(i) = b(1)
            vlat_out(i) = -rv(i) * b(2)
            vlon_out(i) = rv(i) * cos(phic(i)) * b(3)
        end do
        !$omp end parallel do
    end subroutine synth_gradient

    ! ------------------------------------------------------------------
    ! Output
    ! ------------------------------------------------------------------
    character(len=32) function fmt15(x)
        real(dp), intent(in) :: x
        character(len=32) :: buf

        write(buf, '(ES23.14E2)') x
        fmt15 = adjustl(buf)
    end function fmt15

    real(dp) function parse_real(text)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: tmp
        integer :: i

        tmp = text
        do i = 1, len(tmp)
            if (tmp(i:i) == 'D' .or. tmp(i:i) == 'd') tmp(i:i) = 'E'
        end do
        read(tmp, *) parse_real
    end function parse_real

    subroutine write_result_csv(path, n, idv, lat, lon, h, vv, tt, dg)
        character(len=*), intent(in) :: path
        integer(int32), intent(in) :: n
        character(len=64), intent(in) :: idv(:)
        real(dp), intent(in) :: lat(:), lon(:), h(:), vv(:), tt(:), dg(:)
        integer :: u, i

        open(newunit=u, file=path, status='replace', action='write')
        write(u, '(A)') RESULT_HEADER
        do i = 1, n
            write(u, '(A,A,A,A,A,A,A,A,A,A,A,A,A)') trim(idv(i)), ',', trim(fmt15(lat(i))), ',', &
                trim(fmt15(lon(i))), ',', trim(fmt15(h(i))), ',', trim(fmt15(vv(i))), ',', &
                trim(fmt15(tt(i))), ',', trim(fmt15(dg(i)))
        end do
        close(u)
    end subroutine write_result_csv

    subroutine get_build_info(compiler, flags)
        character(len=*), intent(out) :: compiler, flags
        compiler = compiler_version()
        flags = LIB_FLAGS
        call sanitize(compiler)
        call sanitize(flags)
    end subroutine get_build_info

    ! Commas would break the CSV notes column; replace them.
    subroutine sanitize(s)
        character(len=*), intent(inout) :: s
        integer :: i
        do i = 1, len(s)
            if (s(i:i) == ',') s(i:i) = ';'
        end do
    end subroutine sanitize

    subroutine append_timing(path, rep, n, cpu, wall, load, ng, compiler, flags)
        character(len=*), intent(in) :: path, compiler, flags
        integer(int32), intent(in) :: rep, n, ng
        real(dp), intent(in) :: cpu, wall, load
        logical :: exists
        integer(8) :: fsize
        integer :: u
        character(len=2048) :: notes, line
        character(len=32) :: s_nmax, s_n, s_thr, s_rep, s_wall, s_cpu, s_load, s_ng

        inquire(file=path, exist=exists, size=fsize)
        write(s_nmax, '(I0)') nmax
        write(s_n, '(I0)') n
        write(s_thr, '(I0)') threads
        write(s_rep, '(I0)') rep
        write(s_wall, '(F20.6)') wall
        write(s_cpu, '(F20.6)') cpu
        write(s_load, '(F20.3)') load
        write(s_ng, '(I0)') ng
        notes = 'load_s=' // trim(adjustl(s_load)) // '; distinct_radii=' // trim(adjustl(s_ng)) // &
            '; parallel=driver_openmp(V:radius_groups,dg:points); library_point_routines=serial' // &
            '; threads=' // trim(adjustl(s_thr)) // '; platform=Windows; compiler=' // trim(compiler) // &
            '; library_flags=' // trim(flags) // '; functional=' // FUNCTIONAL_TAG
        call sanitize(notes)
        line = 'shtools,native,' // trim(set_name) // ',' // trim(adjustl(s_nmax)) // ',' // &
            trim(adjustl(s_n)) // ',' // trim(adjustl(s_thr)) // ',' // trim(adjustl(s_rep)) // ',' // &
            trim(adjustl(s_wall)) // ',' // trim(adjustl(s_cpu)) // ',' // trim(notes)
        open(newunit=u, file=path, status='unknown', position='append', action='write')
        if (.not. exists .or. fsize == 0) write(u, '(A)') TIMING_HEADER
        write(u, '(A)') trim(line)
        close(u)
    end subroutine append_timing

    subroutine write_result_md(csv_path, n, lmax, set, name, gmv, r0v, tide, load, ng, compiler, flags)
        character(len=*), intent(in) :: csv_path, set, name, tide, compiler, flags
        integer(int32), intent(in) :: n, lmax, ng
        real(dp), intent(in) :: gmv, r0v, load
        character(len=1024) :: md_path, csv_base
        integer :: u, p
        character(len=32) :: buf
        character(len=40) :: stamp
        integer :: vals(8)

        md_path = csv_path
        p = len_trim(md_path)
        if (p > 4 .and. md_path(p-3:p) == '.csv') then
            md_path(p-2:p) = 'md'
        else
            md_path = trim(md_path) // '.md'
        end if
        p = index(csv_path, '/', back=.true.)
        if (p == 0) p = index(csv_path, '\', back=.true.)
        csv_base = csv_path(p+1:)

        call date_and_time(values=vals)
        write(stamp, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') vals(1), vals(2), &
            vals(3), vals(5), vals(6), vals(7)

        open(newunit=u, file=md_path, status='replace', action='write')
        write(u, '(A)') '# ' // trim(csv_base)
        write(u, '(A)') ''
        write(u, '(A)') 'Engine: **SHTOOLS 4.14.1 (Fortran, native, single algorithm)**. Generated by ' // &
            'shtools_shs_bench on ' // trim(stamp) // '. V, T and dg are computed in the same run.'
        write(u, '(A)') ''
        write(u, '(A)') '## Columns'
        write(u, '(A)') ''
        write(u, '(A)') '| column | unit | meaning |'
        write(u, '(A)') '|---|---|---|'
        write(u, '(A)') '| id | - | point identifier from the input CSV |'
        write(u, '(A)') '| lat_deg | deg | geodetic latitude on GRS80 (input, echoed) |'
        write(u, '(A)') '| lon_deg | deg | longitude (input, echoed) |'
        write(u, '(A)') '| h_m | m | ellipsoidal height on GRS80 (input, echoed) |'
        write(u, '(A)') '| V_m2s2 | m^2/s^2 | gravitational potential V of the model at the point (no normal potential ' // &
            'subtracted; degree 0 included, C00 as in the file) |'
        write(u, '(A)') '| T_m2s2 | m^2/s^2 | disturbing potential T = V - U at the point (V: model, ' // &
            'U: GRS80 normal potential, closed form) |'
        write(u, '(A)') '| dg_mgal | mGal | gravity disturbance dg = \|g\| - gamma_SH at the point ' // &
            '(g includes the centrifugal term; gamma_SH: project series) |'
        write(u, '(A)') ''
        write(u, '(A)') '## Command and parameters'
        write(u, '(A)') ''
        write(u, '(A)') '```'
        write(u, '(A)') 'code/fortran/shtools_shs_bench ' // trim(args_text)
        write(u, '(A)') '```'
        write(u, '(A)') ''
        write(u, '(A,I0,A)') '- nmax: `', lmax, '`'
        write(u, '(A,A,A)') '- set: `', trim(set), '`'
        write(u, '(A,I0,A)') '- threads: `', threads, '`'
        write(u, '(A,I0,A)') '- repeat: `', nrep, '`'
        write(u, '(A)') '- normalization: `4pi (fully normalised)`'
        write(u, '(A)') '- csphase: `1 (no Condon-Shortley phase)`'
        write(u, '(A)') '- functionals: `V, T, dg from the same synthesis pass`'
        write(u, '(A)') '- gravity disturbance: `|g| - gamma_SH with centrifugal term, omega=7.292115e-5`'
        write(u, '(A)') '- library: `SHTOOLS 4.14.1 (OpenMP build, libSHTOOLS-mp)`'
        write(u, '(A)') ''
        write(u, '(A)') '## Model'
        write(u, '(A)') ''
        write(buf, '(ES23.10E2)') gmv
        write(u, '(A,A,A)') '- modelname: `', trim(name), '`'
        write(u, '(A,A,A)') '- GM [m^3/s^2]: `', trim(adjustl(buf)), '`'
        write(buf, '(F20.6)') r0v
        write(u, '(A,A,A)') '- R [m]: `', trim(adjustl(buf)), '`'
        write(u, '(A,A,A)') '- tide_system: `', trim(tide), '`'
        write(buf, '(F20.3)') load
        write(u, '(A,A,A)') '- load_s: `', trim(adjustl(buf)), '`'
        write(u, '(A)') ''
        write(u, '(A,I0,A)') 'Points: ', n, '.'
        write(u, '(A)') ''
        write(u, '(A)') '## Notes'
        write(u, '(A)') ''
        write(u, '(A)') '- V = GM/r * sum_n (R/r)^n S_n from MakeGridPoint on the radius-scaled coefficients: ' // &
            'the model potential itself, no normal potential subtracted, degree 0 included (C00 as in the file).'
        write(u, '(A)') '- T = V - U, U = GRS80 level-ellipsoid gravitational potential (closed form), ' // &
            'with the V of this run.'
        write(u, '(A)') '- Functionals evaluated at the actual point (geodetic h), not on the ellipsoid.'
        write(buf, '(I0)') ng
        write(u, '(A,A,A)') '- pyshtools-style radius scaling: one scaled copy of cilm per distinct radius (', &
            trim(adjustl(buf)), ' radii); the copy is made inside the timed block. Each OpenMP thread owns its copy.'
        write(u, '(A)') '- Gradient from MakeGravGridPoint (omega=0); centrifugal term added in code.'
        write(u, '(A)') '- Parallelism: MakeGridPoint and MakeGravGridPoint are serial inside SHTOOLS (no OpenMP ' // &
            'directives). The driver runs the radius groups (V) and the points (gradient) in parallel with ' // &
            'OpenMP; the Legendre work arrays of the library are THREADPRIVATE in the OpenMP build.'
        write(u, '(A)') '- Compiler of the driver: ' // trim(compiler) // '. Library flags: ' // trim(flags) // '.'
        write(u, '(A)') ''
        close(u)
        write(*,'(A,A)') 'wrote ', trim(md_path)
    end subroutine write_result_md

end program shtools_shs_bench
