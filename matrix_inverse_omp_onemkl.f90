!$ include "mkl_omp_offload.f90"

subroutine matrix_inverse_omp_onemkl(A_i,N)

! Decide whether to use 32- or 64-bit integer type
#if defined(MKL_ILP64)
!$  use onemkl_lapack_omp_offload_ilp64   ! 64-bit
#else
!$  use onemkl_lapack_omp_offload_lp64    ! 32-bit
#endif
      
  implicit none
  integer :: N
  integer :: i,j,allocstat,stat
  complex*16, intent(inout) :: A_i(N,N)
  complex*16 :: adi(N,N)
  complex*16 :: work(N,N) ! work is a workspace array
  !A_i is a (N,N) 2D matrix 
  integer :: lda, lwork
  integer :: info_rf, info_rs
  integer, allocatable :: ipiv(:)
  character (len = 132) :: allocmsg
  character (len =  32) :: arg1, arg2
  lda=N
  adi=A_i
  print *,"shape of adi is ",shape(adi)
  lwork=N
  ! Allocate memory for linear algebra computations
#if !defined(_OPENMP)
    allocate(ipiv(N), &
             stat = allocstat, errmsg = allocmsg)
    if (allocstat > 0) stop trim(allocmsg)
#endif

  ! Compute the LU factorizations and inverse the matrix using OpenMP offload.
  ! On entry, "adi" contains the input matrices. On exit, adi contains the inversed matrices.
  !$omp target data map(tofrom:adi) map(tofrom:work) map(from:lwork, info_rf, info_rs)    &
  !$                          map(alloc:ipiv(1:N)) 
    !$omp dispatch
    call zgetrf(N, N, adi, lda, ipiv, info_rf)
    !$omp dispatch
    call zgetri(N, adi, lda, ipiv, work, lwork, info_rs)
  !$omp end target data
  A_i=adi
end subroutine
