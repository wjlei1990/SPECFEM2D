program test_main

	double precision, dimension(2,200) :: A
	double precision, dimension(200,3) :: B
	double precision, dimension(2,3) :: C,D

	A=1.0
	B=1.0
	!A=reshape((/1,0,0,0,0,1/),(/2,3/))
	!B=reshape((/1,2,3,4,5,6/),(/3,2/))
  C=matmul(A,B)
  !do i=1,3
	!	do j=1,2
	!		print *, B(i,j)
!			print *, A(i,j)
!		enddo
!	enddo
	print *, C
	print *,A

	D(1,:)=C(2,:)
	D(2,:)=C(1,:)
!	print *,D

end program test_main
