subroutine init_from_file(A_i,N)
  implicit none
  integer :: N
  integer :: i,j,k,u
  !complex*16, allocatable, intent(inout) :: A_i(:,:)
  complex*16 :: A_i(N,N)
  complex*16, allocatable :: aread(:,:)
  !A_i is a flatterned N*N 2D matrix in the first dimension and the batch_size
  !number of batches in the second dimension
  allocate( aread(N,N) )

  print *, 'inside subroutine init_from_file'
  open(u,file='/nfs/site/home/shiquans/testLU/data/1424x1424.dat',form="unformatted")
  !open(u,file='/nfs/site/home/shiquans/testLU/data/12450x12450.dat',form="unformatted")
  read(u) aread
  close(u)
  print *,"shape of aread is ",shape(aread)

  do i = 1, N
  do j = 1, N
    A_i(i,j) = aread(i,j)  
  enddo
  enddo

  print *, 'A_i(1,1)=',A_i(1,1)
  print *, 'A_i(2,1)=',A_i(2,1)
  print *, 'A_i(n,n-2)=',A_i(n,n-2)
  print *, 'A_i(n,n)=',A_i(n,n)

  deallocate( aread )
end subroutine
