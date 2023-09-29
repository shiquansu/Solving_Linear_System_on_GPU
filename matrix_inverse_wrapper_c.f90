module getri_wrapper
  use iso_c_binding
  interface
     !
     ! single matrix version 
     integer(c_int) function zgetr_i(a, lda, n) bind(C, name="zgetr_i")
       use iso_c_binding
       implicit none
       type(c_ptr),value :: a
       integer(c_int), value :: lda
       integer(c_int),value :: n
     end function zgetr_i
     !
     ! TODO batch version
  end interface
end module getri_wrapper

subroutine matrix_inverse_wrapper_c(a_i,n)
use iso_c_binding
use getri_wrapper
implicit none
complex*16 ::a_i(n,n)

type(c_ptr) :: a_ptr
complex*16, allocatable, target :: a(:,:)
integer :: n, lda, info, u
real*8 :: timer1, timer2

lda=n
allocate(a(n,n))
a=a_i
a_ptr = c_loc(a)
info = 0
print *,"Got the data, direct inversion subroutine starts, shape of a is ",shape(a)
info = zgetr_i(a_ptr,lda,n)
a_i=a
deallocate(a)


end subroutine
