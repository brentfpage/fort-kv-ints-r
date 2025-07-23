function fact(n)
!   use omp_lib
    implicit none
    integer :: n, i
    real :: fact
    fact = 1.0
    do i=2,n
      fact = fact * i
    enddo
end function fact
    
