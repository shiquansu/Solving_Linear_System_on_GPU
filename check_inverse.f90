subroutine check_inverse(A_inv,A,N)
  implicit none
  integer :: N,lda
  integer :: i,j
  complex*16 :: A_inv(N,N),A(N,N),I_unit(N,N)
  complex*16 :: alpha=(1.0d0,0.0d0), beta=(-1.0d0,0.0d0), e2=(0.0d0,0.0d0)
  real*8  :: residual, threshold = 1.0d-9
  lda=N
  I_unit=(0.0d0,0.0d0)
  do i=1,N
    I_unit(i,i)=(1.0d0,0.0d0)
  enddo
  call zgemm('N','N',N,N,N,alpha,A,lda,A_inv,lda,beta,I_unit,lda)

  do j=1,N
  do i=1,N
    e2=e2+I_unit(i,j)*I_unit(i,j)
  enddo
  enddo
  residual=abs(e2)/real(N)
  if (residual > threshold) then
    print *, 'Warning: relative residual of ', residual 
  else
    print *, 'A*Inv(A)=I, PASSED'
  endif
end subroutine

