!$ include "mkl_omp_offload.f90"
!!$ include "mkl_lapack.h"
!!$ include "mkl.fi"
!!$ include "blas.f90"
program solve_batched_linear_systems
! Decide whether to use 32- or 64-bit integer type
#if defined(MKL_ILP64)
!$  use onemkl_lapack_omp_offload_ilp64   ! 64-bit
#else
!$  use onemkl_lapack_omp_offload_lp64    ! 32-bit
#endif
    implicit none
    integer                    :: n = 4096, batch_size = 8, nrhs = 1, cycles = 5
    integer                    :: lda, stride_a, stride_ipiv
    integer                    :: ldb, stride_b
    integer,       allocatable :: ipiv(:,:), info_rf(:), info_rs(:)
#if defined(DC)
    complex*16, allocatable :: a(:,:), b(:,:), a_orig(:,:), b_orig(:,:), x(:), aread(:,:)
    real (kind=8)        :: residual, residual2, threshold = 1.0e-9
    complex*16              :: aa, alpha, beta
#else
    real (kind=8), allocatable :: a(:,:), b(:,:), a_orig(:,:), b_orig(:,:), x(:)
    real (kind=8)              :: residual, threshold = 1.0d-9
#endif
    real (kind=8)        :: a1, a2, b1, b2
    integer (kind=8) :: start_time, end_time, clock_precision
    real    (kind=8) :: cycle_time, total_time = 0.0d0
    integer               :: i, j, k, c, allocstat, stat, u
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
 
    ! Allocate memory for linear algebra computations
    allocate (a(stride_a, batch_size), b(n, batch_size*nrhs), &
#if !defined(_OPENMP)
              ipiv(stride_ipiv, batch_size),                  &
#endif
              info_rf(batch_size), info_rs(batch_size),       &
              aread(12450,12450),                             &
              stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)
 
    ! Allocate memory for error checking
    allocate (a_orig(stride_a, batch_size), b_orig(n, batch_size*nrhs), x(n), &
              stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)
 
    call system_clock(count_rate = clock_precision)
    call random_seed()

    print *,"Reading data!"
    open(u,file='KKR_matrix_for_tests/12450x12450.dat',form="unformatted")
    read(u) aread
    close(u)
    print *,"shape of a is ",shape(aread)
    do c = 1, batch_size
    do i = 1, 12450
    do j = 1, 12450
        a(i+(j-1)*12450,c) = aread(i,j)
    enddo
    enddo
    enddo

 

    print*, 'Got the data!'

    do c = 1, cycles

#if defined(DC)

        do i = 1, n
        do j = 1, batch_size*nrhs
          call random_number(b1)
          call random_number(b2)
          b(i,j) = (2.5,2.5) - cmplx(5.0*b1, 5.0*b2)
        enddo
        enddo

        a_orig = a
        b_orig = b

#else
        ! Initialize the matrices with a random number in the interval (-0.5, 0.5)
        call random_number(a)
        a = 0.5 - a
        ! Make diagonal band values larger to ensure well-conditioned matrices
        do i = 1, n
            a(i+(i-1)*lda,:) = a(i+(i-1)*lda,:) + 50.0
            if (i .ne. 1) a(i+(i-1)*lda-1,:) = a(i+(i-1)*lda-1,:) + 20.0
            if (i .ne. n) a(i+(i-1)*lda+1,:) = a(i+(i-1)*lda+1,:) + 20.0
        enddo

        ! Initialize the RHS with a random number in the interval (-2.5, 2.5)
        call random_number(b)
        b = 2.5 - (5.0 * b)
        a_orig = a
        b_orig = b
#endif
!        anorm = zlarge(norm,n,n,a,lda,rwork)
!        CALL zgecon(norm,n,a,lda,anorm,rcond,work,rwork,info)
!        print *, 'Condition number: ', rcond
        call system_clock(start_time)   ! Start timer
        print *, "Offloading!"
        ! Compute the LU factorizations and solve the linear systems using OpenMP offload.
        ! On entry, "a" contains the input matrices. On exit, it contains the factored matrices.
        !$omp target data map(to:a) map(tofrom: b) map(from:info_rf, info_rs)    &
        !$                          map(alloc:ipiv(1:stride_ipiv, 1:batch_size))
            !$omp dispatch
#if defined(DC)
            call zgetrf_batch_strided(n, n, a, lda, stride_a, ipiv, stride_ipiv, batch_size, info_rf)
#else
            call dgetrf_batch_strided(n, n, a, lda, stride_a, ipiv, stride_ipiv, batch_size, info_rf)
#endif
            !$omp dispatch
#if defined(DC)
            call zgetrs_batch_strided('N', n, nrhs, a, lda, stride_a, ipiv, stride_ipiv, &
                                       b, ldb, stride_b, batch_size, info_rs)
#else
            call dgetrs_batch_strided('N', n, nrhs, a, lda, stride_a, ipiv, stride_ipiv, &
                                       b, ldb, stride_b, batch_size, info_rs)
#endif
        !$omp end target data
        call system_clock(end_time)   ! Stop timer
        if (any(info_rf .ne. 0)) then
            print *, 'Error: getrf_batch_strided returned with errors.'
            stop
        elseif (any(info_rs .ne. 0)) then
            print *, 'Error: getrs_batch_strided returned with errors.'
            stop
        else
            ! Compute a_orig*b and compare result to saved RHS
            do i = 1, batch_size
                do j = 1, nrhs
                    x = (0.0, 0.0)
#if defined(DC)
                    call zgemv('N', n, n, alpha, a_orig(:,i), lda, b(:,(i-1)*nrhs+j), 1, beta, x, 1)
#else
                    call dgemv('N', n, n, 1.0d0, a_orig(:,i), lda, b(:,(i-1)*nrhs+j), 1, 0.0d0, x, 1)
#endif
                    ! Check relative residual
#if defined(DC)
                    !residual = dznrm2(n, b_orig(:,(i-1)*nrhs+j)-x(:), 1) / dznrm2(n, b_orig(:,(i-1)*nrhs+j), 1)
                    residual=0.0d0
                    residual2=0.0d0
                    !print *,'b_orig=',b_orig
                    !print *,'x=',x
                    do k=1,n
                      aa=b_orig(k,(i-1)*nrhs+j)-x(k)
                      residual=residual+real(aa)**2+aimag(aa)**2
                      aa=b_orig(k,(i-1)*nrhs+j)
                      residual2=residual2+real(aa)**2+aimag(aa)**2
                    enddo
                    residual=residual/residual2
                    print *,'residual=',residual,'residual2=',residual2
#else
                    residual = norm2(b_orig(:,(i-1)*nrhs+j) - x(:)) / norm2(b_orig(:,(i-1)*nrhs+j))
#endif
                    if (residual > threshold) then
                       print *, 'Warning: relative residual of ', residual
                    endif
                enddo
            enddo
            cycle_time = dble(end_time - start_time) / dble(clock_precision)
            print *, 'Computation completed successfully', cycle_time, 'seconds'
            ! Do not include the JIT compilation overhead of the first cycle.
            if (c > 1) then
                total_time = total_time + cycle_time
            endif
        endif
    enddo
    print *, 'Total time:', total_time, 'seconds'
    ! Clean up
#if defined(_OPENMP)
    deallocate (a, b, a_orig, b_orig, x, info_rf, info_rs, aread)
#else
    deallocate (a, b, a_orig, b_orig, x, ipiv, info_rf, info_rs, aread)
#endif
end program solve_batched_linear_systems
