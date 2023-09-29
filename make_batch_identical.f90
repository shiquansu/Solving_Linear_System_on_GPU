subroutine make_batch_identical(a2d, a2d_batch, N, batch_size)
        !   init_from_file_batch(A_i,N,batch_size)
  implicit none
  integer :: N,batch_size
  integer :: i,j,k,u
  !complex*16, allocatable, intent(inout) :: A_i(:,:)
  complex*16 :: a2d(N,N)
  complex*16 :: a2d_batch(N*N,batch_size)
  !a2d_batch is a flatterned N*N 2D matrix in the first dimension and the batch_size
  !number of batches in the second dimension

  print *,"shape of a2d in make_batch_identical is ",shape(a2d)
  print *,"shape of a2d_batch in make_batch_identical is ",shape(a2d_batch)

  do k = 1, batch_size
  do i = 1, N
  do j = 1, N
    a2d_batch(i+(j-1)*N,k) = a2d(i,j)
  enddo
  enddo
  enddo

end subroutine
