subroutine init_trivial(A_i,N)
  implicit none
  integer :: N
  integer :: i,j,k,u
  !complex*16, allocatable, intent(inout) :: A_i(:,:)
  complex*16 :: A_i(N,N)
  real(8)::a1,a2
  !A_i is a N*N 2D matrix, initial by random number, making diagonal dominated 

  call random_seed()
  do i = 1, N
  do j = 1, N
    call random_number(a1)
    call random_number(a2)
    A_i(i,j) = (0.5,0.5) - cmplx(a1, a2)
  enddo
  enddo
  do i = 1, N
    A_i(i,i) = A_i(i,i)+(2000000.0d0,10.0d0)
  enddo

  print *,"shape of A_i of init_random is ",shape(A_i)

end subroutine
