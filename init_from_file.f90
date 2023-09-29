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

  print *,"shape of aread is ",shape(aread)
  !open(u,file='/nfs/site/home/shiquans/testLU/data/1424x1424.dat',form="unformatted")
  open(u,file='/nfs/site/home/shiquans/testLU/data/12450x12450.dat',form="unformatted")
  read(u) aread
  close(u)

  do i = 1, N
  do j = 1, N
    A_i(i,j) = aread(i,j)  
  enddo
  enddo

  deallocate( aread )
end subroutine
