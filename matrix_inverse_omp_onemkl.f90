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
  integer :: i,j,stat,allocstat=0
  complex*16, intent(inout) :: A_i(N,N)
  !A_i is a (N,N) 2D matrix 
  integer :: lda, lwork
  integer :: info_rf=1, info_rs=1
  character (len = 132) :: allocmsg
  character (len =  32) :: arg1, arg2

  integer, allocatable :: ipiv(:) ! all mkl used matrix need allocatable attribute to work properly with openmp dispatch
  complex*16, allocatable :: adi(:,:)
  complex*16, allocatable :: work(:) ! work is a workspace array

  lda=N
  adi=A_i
  lwork=N*640 ! offload disabled and set lwork=-1 to get the optimized lwork from work(1)
  !lwork=N*4
  ! Allocate memory for linear algebra computations

  allocate(ipiv(N), stat = allocstat, errmsg = allocmsg)
  allocate(adi(N,N), stat = allocstat, errmsg = allocmsg)
  allocate(work(N*N), stat = allocstat, errmsg = allocmsg)
  print *,"shape of adi: ",shape(adi)
  print *,"shape of work, ipiv: ",shape(work),shape(ipiv)

  ! Compute the LU factorizations and inverse the matrix using OpenMP offload.
  ! On entry, "adi" contains the input matrices. On exit, adi contains the inversed matrices.
  !$omp target data map(tofrom:adi) map(alloc:work) map(tofrom:lwork, info_rf, info_rs)    &
  !$                          map(alloc:ipiv(1:N)) 
    !$omp dispatch
    call zgetrf(lda, lda, adi, lda, ipiv, info_rf)
    !$omp dispatch
    call zgetri(lda, adi, lda, ipiv, work, lwork, info_rs)
  !$omp end target data
  !print *,"optimized lwork value from work(1)=",work(1)
  A_i=adi
end subroutine
