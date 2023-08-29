!$ include "mkl_omp_offload.f90"

subroutine matrix_inverse_omp_onemkl_batch(A_i,N,batch_size)

! Decide whether to use 32- or 64-bit integer type
#if defined(MKL_ILP64)
!$  use onemkl_lapack_omp_offload_ilp64   ! 64-bit
#else
!$  use onemkl_lapack_omp_offload_lp64    ! 32-bit
#endif
      
  implicit none
  integer :: N,batch_size
  integer :: i,j,allocstat,stat
  complex*16, intent(inout) :: A_i(N*N,batch_size)
  complex*16 :: adi(N*N,batch_size)
  complex*16 :: Ainv(N*N,batch_size)
  !A_i is a flatterned N*N 2D matrix in the first dimension and the batch_size
  !number of batches in the second dimension
  integer :: lda, stride_a, stride_ipiv
  integer, allocatable :: ipiv(:,:), info_rf(:), info_rs(:)
  character (len = 132) :: allocmsg
  character (len =  32) :: arg1, arg2
  lda=N
  adi=A_i
  stride_a=N*lda
  stride_ipiv=N
  print *,"shape of a is ",shape(adi)
  ! Allocate memory for linear algebra computations
#if !defined(_OPENMP)
    allocate(ipiv(stride_ipiv, batch_size), &
             stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)
#endif
    allocate(info_rf(batch_size), info_rs(batch_size), &
             stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)

  !how to set the remaining ipiv stride_ipiv info_rf?
  !how to set zgetri openmp strided version?

  ! Compute the LU factorizations and inverse the matrix using OpenMP offload.
  ! On entry, "adi" contains the input matrices. On exit, adi contains the inversed matrices.
  !$omp target data map(tofrom:adi) map(tofrom:Ainv) map(from:info_rf, info_rs)    &
  !$                          map(alloc:ipiv(1:stride_ipiv, 1:batch_size))
    !$omp dispatch
    call zgetrf_batch_strided(lda, lda, adi, lda, stride_a, ipiv, stride_ipiv, batch_size, info_rf)
    !$omp dispatch
    call zgetri_oop_batch_strided(lda, adi, lda, stride_a, ipiv, stride_ipiv, Ainv, lda, stride_a, batch_size, info_rs)
  !$omp end target data
  A_i=Ainv
end subroutine
