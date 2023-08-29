!!$ include "mkl_omp_offload.f90"
!!$ include "mkl.fi"
!!$ include "blas.f90"

program solve_complex_matrix 

! Decide whether to use 32- or 64-bit integer type
!#if defined(MKL_ILP64)
!!$  use onemkl_lapack_omp_offload_ilp64   ! 64-bit
!#else
!!!$  use onemkl_lapack_omp_offload_lp64    ! 32-bit
!#endif

    implicit none

    integer                    :: n = 1424, batch_size = 1, nrhs = 1, cycles = 2
    integer                    :: lda, stride_a, stride_ipiv
    integer                    :: ldb, stride_b
    integer,       allocatable :: ipiv(:,:), info_rf(:), info_rs(:)

    complex*16, allocatable :: a(:,:), b(:,:), a_orig(:,:), b_orig(:,:), x(:)
    complex*16, allocatable :: a2d(:,:), a2d_orig(:,:)
    real (kind=8)        :: residual, residual2, threshold = 1.0e-9
    real (kind=8)        :: a1, a2, b1, b2
    complex*16              :: aa, alpha, beta

    integer (kind=8) :: start_time, end_time, clock_precision
    real    (kind=8) :: cycle_time, total_time = 0.0d0

    integer               :: i, j, k, c, allocstat, stat
    character (len = 132) :: allocmsg
    character (len =  32) :: arg1, arg2

    alpha = (1.0, 0.0)
    beta = (0.0, 0.0)
    ! Simple command-line parser with no error checks
    do i = 1, command_argument_count(), 2
        call get_command_argument(i, arg1)
        select case (arg1)
            case ('-n')
                call get_command_argument(i+1, arg2)
                read(arg2, *, iostat=stat) n
            case ('-b')
                call get_command_argument(i+1, arg2)
                read(arg2, *, iostat=stat) batch_size
            case ('-r')
                call get_command_argument(i+1, arg2)
                read(arg2, *, iostat=stat) nrhs
            case ('-c')
                call get_command_argument(i+1, arg2)
                read(arg2, *, iostat=stat) cycles
            case default
                print *, 'Unrecognized command-line option:', arg1
                stop
        end select
    enddo
    print *, 'Matrix dimensions:', n
    print *, 'Batch size:', batch_size
    print *, 'Number of RHS:', nrhs
    print *, 'Number of test cycles:', cycles

    lda         = n
    stride_a    = n * lda
    stride_ipiv = n
    ldb         = n
    stride_b    = n * nrhs



    !! batch version
    ! Allocate memory for linear algebra computations
    allocate (a(stride_a, batch_size), b(n, batch_size*nrhs), &
              stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)

#if !defined(_OPENMP)
    allocate(ipiv(stride_ipiv, batch_size), &
             stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)
#endif

    allocate(info_rf(batch_size), info_rs(batch_size), &
             stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)

    ! Allocate memory for error checking
    allocate (a_orig(stride_a, batch_size), b_orig(n, batch_size*nrhs), x(n), &
              stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)

    !! single matrix version
    allocate (a2d(n,n), a2d_orig(n,n), &
              stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)


    call system_clock(count_rate = clock_precision)
    call random_seed()

    call init_from_file(a2d,n) 
    !call init_from_file_batch(a,n,batch_size) 
    !call init_random(a,n,batch_size) 

    do c = 1, cycles
        call system_clock(start_time)   ! Start timer

        !call matrix_inverse_omp_onemkl_batch(a,n,batch_size)
        !call matrix_inverse_wrapper_c_batch(a,n,batch_size) !! not implemented yet

        call matrix_inverse_omp_onemkl(a2d,n)
        !call matrix_inverse_wrapper_c(a2d,n)

        call system_clock(end_time)   ! Stop timer

        call check_inverse(a2d,a2d_orig,n)
        !call check_inverse_batch(a,a_orig,n,batch_size)

        cycle_time = dble(end_time - start_time) / dble(clock_precision)
        print *, 'Computation completed successfully', cycle_time, 'seconds'
        ! Do not include the JIT compilation overhead of the first cycle.
            if (c > 1) then
                total_time = total_time + cycle_time
            endif
    enddo

    print *, 'Total time:', total_time, 'seconds'

    ! Clean up
#if defined(_OPENMP)
    deallocate (a, b, a_orig, b_orig, x, info_rf, info_rs)
#else
    deallocate (a, b, a_orig, b_orig, x, ipiv, info_rf, info_rs)
#endif
end program solve_complex_matrix

