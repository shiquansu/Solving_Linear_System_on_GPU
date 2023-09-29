subroutine init_trivial(A_i,N)
  implicit none
  integer :: N
  integer :: i,j,k,u
  !complex*16, allocatable, intent(inout) :: A_i(:,:)
  complex*16 :: A_i(N,N)
  !A_i is a N*N 2D matrix, trivial to solve by 2x2 along diagonal 

  do i = 1, N
  do j = 1, N
    A_i(i,j) = (0.0d0,0.0d0)  
  enddo
  enddo
  do i = 1, N
    A_i(i,i) = (2.0d0,0.0d0)
  enddo
  !A_i(1,2)=(0.0d0,0.0d0)
  !A_i(2,1)=(0.0d0,0.0d0)
  A_i(1,1)=(3.0d0,0.0d0)
  A_i(2,2)=(2.0d0,0.0d0)
  A_i(1,2)=(1.0d0,0.0d0)
  A_i(2,1)=(5.0d0,0.0d0)

  print *,"shape of A_i of init_trivial is ",shape(A_i)

end subroutine
