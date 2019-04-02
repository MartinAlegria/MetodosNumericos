PROGRAM SECOND_PARTIAL

	IMPLICIT NONE
	CALL gauss_elimination()
	!CALL lu_decomp()
	CALL gauss_seidel()


END PROGRAM SECOND_PARTIAL

!********** System of Linear Equations **********!

SUBROUTINE gauss_elimination()
	IMPLICIT NONE
	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l 
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:,:), allocatable :: calcs
	DOUBLE PRECISION, dimension (:), allocatable :: results
	DOUBLE PRECISION :: temps

	open(unit = 10, file = "test.txt")
	read(10,*)n
	write(*,*)n
	n_1 = n+1
	allocate ( matrix(n,n_1) )
	allocate ( calcs(n,n_1) )
	allocate ( results(n) )    
	do i=1,n
		read(10,*)matrix(i,:)
		write(*,*)matrix(i,:)
	enddo

	calcs = matrix
	close(10) 
	!********** READ FILE **********!


	!********** TURN MATRIX INTO LOWER TRIANGLE MATRIX **********!
	do i=1,n-1
		matrix(i,:) = matrix(i,:)/calcs(i,i)
		calcs = matrix
		do j=i+1,n
			calcs = matrix
			do k=1,n+1
				matrix(j,k) = ((-1)*calcs(j,i) * matrix(i,k)) + matrix(j,k)
			enddo
		enddo
		calcs = matrix
		write(*,*) "------------------- New matrix -------------------"
		do l=1,n
			write(*,*)matrix(l,:)
		enddo
	enddo

	write(*,*) "-------------------Final matrix -------------------"
	do i=1,n
			write(*,*)matrix(i,:)
	enddo
	!********** TURN MATRIX INTO LOWER TRIANGLE MATRIX **********!

	!********** BACKWARDS SUBSTITUTION **********!
	results(n) = (matrix(n,n_1)/matrix(n,n))!--- Get Xsub1 first

	do i=n-1, 1, -1
		temps = 0
		do j=i+1, n
			temps = temps + matrix(i,j)*results(j)
		enddo
		results(i) = (matrix(i,n_1) - temps)/matrix(i,i)
	enddo

	!********** EXPORT TO CSV **********!
	open (unit = 2, file = "gauss_elimination.csv")
	write(2,*)"Xsub",(",",j,j=1,n)
	write(2,*)"-", (",",results(i), i=1,n)
	close(2)
	!********** EXPORT TO CSV **********!

	!********** BACKWARDS SUBSTITUTION **********! !*********** DONE ***********!
END SUBROUTINE gauss_elimination

SUBROUTINE lu_decomp()
	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l 
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:,:), allocatable :: lm
	DOUBLE PRECISION, dimension (:,:), allocatable :: um
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: temps

	open(unit = 10, file = "test.txt")
	read(10,*)n
	n_1 = n+1
	allocate ( matrix(n,n) )
	allocate ( lm(n,n) ) 
	allocate ( um(n,n) ) 
	allocate ( y(n) )
	allocate ( x(n) )
	allocate ( b(n) )

	write(*,*) "************* Read Matrix *************"   
	do i=1,n
		read(10,*)matrix(i,:),b(i)
		write(*,*)matrix(i,:)
	enddo

 
	close(10)
	!********** READ FILE **********!

	lm(:,1) = matrix(:,1)
	um(1,:) = matrix(1,:)/lm(1,1)
	um(1,1) = 1;

	do k=2,n
		do j=2,n
			
			do i=j,n
					lm(i,j) = matrix(i,j) - dot_product(lm(i,1:j-1), um(1:j-1,j))
			enddo
			um(k,j) = (matrix(k,j) - dot_product(lm(k,1:k-1), um(1:k-1,j)) )/lm(k,k)
		enddo
	enddo

	do i=n-1,1,-1
		do j=n,i+1,-1
			lm(i,j) = 0
		enddo
	enddo

	do i=2,n
		do j=1,i-1
			um(i,j)=0
		enddo
	enddo

	write(*,*) "************* L matrix *************"
	do i=1,n
		write(*,*)lm(i,:)
	enddo
	write(*,*) "************* U matrix *************"
	do i=1,n
		write(*,*)um(i,:)
	enddo
	write(*,*) "************* B matrix *************"
	write(*,*) b(:)

	! ******* FOWARD SUBSTITUTION
	do i=1,n
		temps = b(i)
		do j=1,i-1
			temps = temps - lm(i,j)*y(j)
		enddo
		y(i) = temps/lm(i,i)
	enddo

	write(*,*) "Y MATRIX *********************"
	write(*,*) y(:)

	x(n) = (y(n)/um(n,n))!--- Get Xsub1 first
	write(*,*) x(n)
	do i=n-1, 1, -1
		temps = 0
		do j=i+1, n
			temps = temps + um(i,j)*x(j)
		enddo
		x(i) = (y(i) - temps)/um(i,i)
	enddo


	!********** EXPORT TO CSV **********!
	open (unit = 2, file = "lu_decomp.csv")
	write(2,*)"Xsub",(",",j,j=1,n)
	write(2,*)"-", (",",x(i), i=1,n)
	close(2)
	!********** EXPORT TO CSV **********! !*********** DONE ***********!
END SUBROUTINE lu_decomp

SUBROUTINE gauss_seidel()

	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l 
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: past_x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: temps, sum, error

	open(unit = 10, file = "test.txt")
	read(10,*)n
	n_1 = n+1
	allocate ( matrix(n,n) )
	allocate ( x(n) )
	allocate ( b(n) )

	write(*,*) "************* Read Matrix *************"   
	do i=1,n
		read(10,*)matrix(i,:),b(i)
		write(*,*)matrix(i,:)
	end do

	x(:) = 0

	do

	!********** EVALUATE EQUATIONS **********!
	do i=1,n
		write(*,*) "************* for 1 *************" 
		temps = 0
		past_x = x
		do j=1,n
			write(*,*) "************* for 2 *************" 
			if (j /= i) then
				temps = temps + matrix(i,j)*x(j)
			endif
		end do
		x(i) = (b(i) - temps)/matrix(i,i)
	end do

	!********** EVALUATE EQUATIONS **********!

	write(*,*) x(:)

	!********** CALCULATE AND COMPARE ERRORS **********!
	sum = 0
	do i=1,n
		error = x(i) - past_x(i)
		sum = sum + error
	enddo

	write(*,*) sum

	if (abs(sum)<= 0.0001) then
		EXIT
	endif

	!********** CALCULATE AND COMPARE ERRORS **********!

 	end do

 	write(*,*) "@@@@ FINAL -------", x(:)

	close(10)
	!********** READ FILE **********!
END SUBROUTINE gauss_seidel

!********** Interpolation **********!

SUBROUTINE power_series()
END SUBROUTINE power_series

SUBROUTINE lagrange()
END SUBROUTINE lagrange

SUBROUTINE newton()
END SUBROUTINE newton

!********* FILE ************!



