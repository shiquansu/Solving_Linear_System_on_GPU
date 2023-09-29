subroutine check_inverse_batch(A_inv,A,N,batch_size)
  implicit none
  integer :: N,lda,batch_size,stride_a
  integer :: i,j
  complex*16 :: A_inv(N*N,batch_size),A(N*N,batch_size),I_unit(N*N,batch_size)
  complex*16 :: alpha=(1.0d0,0.0d0), beta=(-1.0d0,0.0d0), e2=(0.0d0,0.0d0)

  real*8  :: residual, threshold = 1.0d-9
  lda=N
  stride_a=lda*N
  I_unit=(0.0d0,0.0d0)

  ! set I_unit as batch unit matrix
  do i = 1, N
    I_unit(i+(i-1)*lda,:) = (1.0d0,0.0d0)
  enddo

  call zgemm_batch_strided('N','N', N, N, N, alpha, A, lda, stride_a, A_inv, lda, stride_a, beta, I_unit, lda, stride_a, batch_size)

  do j=1,batch_size
  do i=1,stride_a
    e2=e2+I_unit(i,j)*I_unit(i,j)
  enddo
  enddo

  residual=abs(e2)/real(N*batch_size)
  if (residual > threshold) then
    print *, 'Warning: relative residual of ', residual 
  else
    print *, 'A*Inv(A)=I, PASSED'
  endif
end subroutine

