!******************************************************************
!*  Inversion of a symmetric matrix by Cholesky decomposition.    *
!*  The matrix must be positive definite.                         * 
!* -------------------------------------------------------------- *
!* REFERENCE:                                                     *
!*             From a Java Library Created by Vadim Kutsyy,       *
!*             "http://www.kutsyy.com".                           *
!* -------------------------------------------------------------- * 
!* SAMPLE RUN:                                                    *
!*                                                                *
!* Inversion of a square real symetric matrix by Cholevsky method *
!* (The matrix must positive definite).                           *
!*                                                                *
!* Size = 4                                                       *
!*                                                                *
!* Matrix A:                                                      *
!*   5.000000  -1.000000  -1.000000  -1.000000                    *
!*  -1.000000   5.000000  -1.000000  -1.000000                    *
!*  -1.000000  -1.000000   5.000000  -1.000000                    *
!*  -1.000000  -1.000000  -1.000000   5.000000                    *
!*                                                                *
!* Determinant = 432.000000                                       *
!*                                                                *
!* Matrix Inv(A):                                                 *
!*   0.250000   0.083333   0.083333   0.083333                    *
!*   0.083333   0.250000   0.083333   0.083333                    *
!*   0.083333   0.083333   0.250000   0.083333                    *
!*   0.083333   0.083333   0.083333   0.250000                    *
!*                                                                *
!*                      F90 Release By Jean-Pierre Moreau, Paris. *
!* -------------------------------------------------------------- *
!* Release 1.1 : added verification Inv(A) * A = I.               *
!******************************************************************
module solve_eg

  implicit none

contains

! ------------------------------------------------
!     Cholesky decomposition.

!     input    n  size of matrix
!     input    A  Symmetric positive def. matrix
!     output  aa  lower deomposed matrix
!     uses        choldc1(int,MAT,VEC)
! ------------------------------------------------
Subroutine choldc(n,A,aa)
  integer n
  real*8 A(0:n-1,0:n-1), aa(0:n-1,0:n-1)
  integer i,j, ialloc 
  real*8, pointer :: p(:)
  allocate(p(0:n-1),stat=ialloc)

  aa = A
 
  call choldc1(n, aa, p)

  do i = 0, n-1
    aa(i,i) = p(i)
    do j = i + 1, n-1
	  aa(i,j) = 0.d0
    end do
  end do
  deallocate(p)
  return
End

! -----------------------------------------------------
!     Inverse of Cholesky decomposition.

!     input    n  size of matrix
!     input    A  Symmetric positive def. matrix
!     output  aa  inverse of lower decomposed matrix
!     uses        choldc1(int,MAT,VEC)         
! -----------------------------------------------------
Subroutine choldcsl(n,A,aa)
  integer n
  real*8 A(0:n-1,0:n-1), aa(0:n-1,0:n-1)
  integer i,j,k, ialloc
  real*8 sum
  real*8, pointer :: p(:)
  allocate(p(0:n-1),stat=ialloc)

  aa = A

  call choldc1(n, aa, p)

  do i = 0, n-1
    aa(i,i) = 1.d0 / p(i)
    do j = i + 1, n-1
      sum = 0.d0
      do k = i, j-1
	    sum = sum - aa(j,k) * aa(k,i)
      end do
      aa(j,i) = sum / p(j)
    end do
  end do
  deallocate(p)
  return
End
 
! ----------------------------------------------------------------------
! Computation of Determinant of the matrix using Cholesky decomposition

! input    n  size of matrix
! input    a  Symmetric positive def. matrix
! return      det(a)
! uses        choldc(int,MAT,MAT)
! ----------------------------------------------------------------------
real*8 Function choldet(n,a)
  integer n
  real*8 a(0:n-1,0:n-1)
  real*8, pointer :: c(:,:)
  real*8 d
  integer i, ialloc
  allocate(c(0:n-1,0:n-1),stat=ialloc)
  d=1.d0
  call choldc(n,a,c)
  do i = 0, n-1
    d = d * c(i,i)
  end do
  choldet = d * d
  deallocate(c)
  return
End

 
! ---------------------------------------------------
!   Matrix inverse using Cholesky decomposition

!   input    n  size of matrix
!   input	 A  Symmetric positive def. matrix
!   output  aa  inverse of A
!   uses        choldc1(int,MAT,VEC)
! ---------------------------------------------------
Subroutine cholsl(n,A,aa)
  integer n
  real*8 A(0:n-1,0:n-1), aa(0:n-1,0:n-1)
  integer i,j,k

  call choldcsl(n,A,aa)

  do i = 0, n-1
    do j = i + 1, n-1
	  aa(i,j) = 0.d0
    end do
  end do

  do i = 0, n-1
    aa(i,i) = aa(i,i) * aa(i,i)
    do k = i + 1, n-1
      aa(i,i) = aa(i,i) + aa(k,i) * aa(k,i)
    end do

    do j = i + 1, n-1
      do k = j, n-1
        aa(i,j) = aa(i,j) + aa(k,i) * aa(k,j)
      end do
    end do
  end do
  do i = 0,  n-1
    do j = 0, i-1
      aa(i,j) = aa(j,i)
    end do
  end do
  return
End

! -------------------------------------------------
! main method for Cholesky decomposition.
!
! input         n  size of matrix
! input/output  a  matrix
! output        p  vector of resulting diag of a
! author:       <Vadum Kutsyy, kutsyy@hotmail.com>
! -------------------------------------------------
Subroutine choldc1(n,a,p)
  integer n
  real*8 a(0:n-1,0:n-1), p(0:n-1)
  integer i,j,k
  real*8 sum
  do i = 0, n-1
    do j = i, n-1
      sum = a(i,j)
      do k = i - 1, 0, -1
        sum = sum - a(i,k) * a(j,k)
      end do
      if (i.eq.j) then
        if (sum <= 0.d0)  &
          print *,' the matrix is not positive definite!'
        p(i) = dsqrt(sum)
      else
        a(j,i) = sum / p(i)
      end if
	end do
  end do
  return
End

! print a square real matrix A of size n with caption s
! (n items per line).
Subroutine MatPrint(s,n,A)
  character*(*) s
  integer n
  real*8 A(0:n-1,0:n-1)
  integer i,j
  print *,' '
  print *,' ', s
  do i=0, n-1
    write(*,10) (A(i,j),j=0,n-1)
  end do
  return
10 format(10F10.6)
End

Integer Function Check_Matrix(n,A)
  integer n
  real*8 A(0:n-1,0:n-1)
  integer i,j,k
  real*8 sum
  Check_Matrix=1
  do i = 0, n-1
    do j = i, n-1
      sum = A(i,j)
      do k = i - 1, 0, -1
        sum = sum - a(i,k) * a(j,k)
      end do 
      if (i.eq.j) then
        if (sum <= 0.d0) Check_Matrix=0
      end if
	end do
  end do
  return
End

!*******************************************
!*    MULTIPLICATION OF TWO SQUARE REAL    *                                     
!*    MATRICES                             *
!* --------------------------------------- *                                     
!* INPUTS:    A  MATRIX N*N                *                                     
!*            B  MATRIX N*N                *                                     
!*            N  INTEGER                   *                                     
!* --------------------------------------- *                                     
!* OUTPUTS:   C  MATRIX N*N PRODUCT A*B    *                                     
!*                                         *
!*******************************************
Subroutine MATMULT(n,A,B,C)
  integer n
  real*8 A(0:n-1,0:n-1), B(0:n-1,0:n-1), C(0:n-1,0:n-1)
  real*8 SUM
  integer I,J,K
  do I=0, n-1                                                                  
    do J=0, n-1
      SUM= 0.d0                                                                
      do K=0, n-1
        SUM=SUM+A(I,K)*B(K,J)
	  end do  	                                               
      C(I,J)=SUM                                                            
    end do
  end do
  return	                                                             
End

!end of file choles.f90


end module solve_eg
