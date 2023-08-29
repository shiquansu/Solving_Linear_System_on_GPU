subroutine init_from_file_batch(A_i,N,batch_size)
  implicit none
  integer :: N,batch_size
  integer :: i,j,k,u
  !complex*16, allocatable, intent(inout) :: A_i(:,:)
  complex*16 :: A_i(N*N,batch_size)
  complex*16, allocatable :: aread(:,:)
  !A_i is a flatterned N*N 2D matrix in the first dimension and the batch_size
  !number of batches in the second dimension
  allocate( aread(N,N) )

  print *, 'inside subroutine init_from_file batch_size=', batch_size
  open(u,file='/nfs/site/home/shiquans/testLU/data/1424x1424.dat',form="unformatted")
  !open(u,file='/nfs/site/home/shiquans/testLU/data/12450x12450.dat',form="unformatted")
  read(u) aread
  close(u)
  print *,"shape of aread is ",shape(aread)

  do k = 1, batch_size
  do i = 1, N
  do j = 1, N
    A_i(i+(j-1)*N,k) = aread(i,j)  
  enddo
  enddo
  enddo

  print *, 'A_i(1,1)=',A_i(1,1)
  print *, 'A_i(1+n,1)=',A_i(1+n,1)
  print *, 'A_i(n*n-2,1)=',A_i(n*n-2,1)
  print *, 'A_i(n*n,1)=',A_I(n*n,1)

  deallocate( aread )
end subroutine
