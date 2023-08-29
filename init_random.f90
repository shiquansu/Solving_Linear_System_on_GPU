
subroutine init_random(A_i,N,batch_size)
  implicit none
  integer :: N,lda,batch_size
  integer :: i,j,k
  complex*16, allocatable, intent(inout) :: A_i(:,:)
  real(8) :: a1,a2

  lda=N
  ! Initialize the matrices with a random number in the interval (-0.5, 0.5)
  call random_seed()
  do i = 1, N*N
  do j = 1, batch_size
    call random_number(a1)
    call random_number(a2)
    A_i(i,j) = (0.5,0.5) - cmplx(a1, a2)
  enddo
  enddo

  ! Make diagonal band values larger to ensure well-conditioned matrices
  do i = 1, N 
    A_i(i+(i-1)*lda,:) = A_i(i+(i-1)*lda,:) + (200000.0,200000.0)
    !if (i .ne. 1) A_i(i+(i-1)*lda-1,:) = A_i(i+(i-1)*lda-1,:) + (2000.0,10.0)
    !if (i .ne. n) A_i(i+(i-1)*lda+1,:) = A_I(i+(i-1)*lda+1,:) + (2000.0,10.0)
  enddo

end subroutine
