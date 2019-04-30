PROGRAM SECOND_PARTIAL

	IMPLICIT NONE
	INTEGER :: selection
	CHARACTER::loop
	!CALL gauss_elimination()
	!CALL lu_decomp()
	!CALL gauss_seidel()
	!CALL power_series()
	!CALL lagrange()
	!CALL newton()
	DO
      write (*,*) "WELCOME!"
      write (*,*) "########### BE SURE TO READ THE README.txt FILE BEFORE USING THE PROGRAM ###############"
      write (*,*) "WHAT METHOD DO YOU WANT TO USE [TYPE THE NUMBER OF THE DESIRED METHOD]"
      write (*,*) "########### SYSTEM OF EQUATIONS ###############"
      write (*,*) "1) GAUSS ELIMINATION"
      write (*,*) "2) LU DECOMPOSITION"
      write (*,*) "3) GAUSS-SEIDEL"
      write (*,*) "########### INTERPOLATION ###############"
      write (*,*) "4) POWER SERIES"
      write (*,*) "5) LAGRANGE"
      write (*,*) "6) NEWTON"
			write (*,*) "7) REGRESSION"
			write (*,*) "8) TRAPEZOIDAL"


      read(*,*)selection

      SELECT CASE(selection)
         CASE(1)
            CALL gauss_elimination()
         CASE(2)
            CALL lu_decomp()
         CASE(3)
            CALL gauss_seidel()
         CASE(4)
            CALL power_series()
         CASE(5)
          	CALL lagrange()
         CASE(6)
          	CALL newton()
					CASE(7)
						CALL regression()
					CASE(8)
						CALL trapezoidal()
      END SELECT

      write (*,*) "DO YOU WANT TO TRY ANOTHER METHOD ? [Y/N]"
      read(*,*) loop
      if(loop == "N" .OR. loop == "n") then
         exit
      end if

  END DO

END PROGRAM SECOND_PARTIAL
SUBROUTINE trapezoidal()
	INTEGER :: n,i,j,k,l,n_1
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	!DOUBLE PRECISION, dimension (:), allocatable :: poly
	DOUBLE PRECISION :: res=0.0

	open(unit = 10, file = "Regression.txt")
	read(10,*)n,r
	allocate ( y(n) )
	allocate ( x(n) )
	!allocate ( poly(n) )

	write(*,*) "************* Read Data *************"
	do i=1,n
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo

	!write(*,*) "************* Number of regression: ", r
	close(10)
	n_1=n-1
	do i=1,n_1
		res=((x(i+1)-x(i))/2)*(y(i)+y(i+1))+res
	end do
	write(*,*) res

END SUBROUTINE trapezoidal
SUBROUTINE regression()
	INTEGER :: n,i,j,k,l
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: poly
	DOUBLE PRECISION :: sumx=0.0, sumy=0.0,r, sumxy=0.0,sumx2=0.0,a0,a1

	open(unit = 10, file = "Regression.txt")
	read(10,*)n,r
	allocate ( y(n) )
	allocate ( x(n) )
	allocate ( poly(n) )

	write(*,*) "************* Read Data *************"
	do i=1,n
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo

	write(*,*) "************* Number of regression: ", r
	close(10)
IF (r<=1) then
	write(*,*) "Linear regression"
	do i=1,n
		sumx=sumx+x(i)
		sumy=sumy+y(i)
		sumxy=sumxy+x(i)*y(i)
		sumx2=sumx2+x(i)**2
	end do
a1=(n*sumxy-sumx*sumy)/(n*sumx2-(sumx**2))
a0=(sumy/n)-a1*(sumx/n)
write(*,*) a0
write(*,*) a1

ELSE
	write(*,*) "Cuadratic or cubic regression"
END IF

END SUBROUTINE regression

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
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
	open (unit = 2, file = "gauss_elimination.csv")
	write(2,*)"Xsub",(",",j,j=1,n)
	write(2,*)"-", (",",results(i), i=1,n)
	close(2)
	!********** EXPORT TO CSV **********!

	!********** BACKWARDS SUBSTITUTION **********! !*********** DONE ***********!  !####### DONE -------------->
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
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
	write(2,*)"Xsub",(",",j,j=1,n)
	write(2,*)"-", (",",x(i), i=1,n)
	close(2)
	!********** EXPORT TO CSV **********! !*********** DONE ***********!			!####### DONE -------------->
END SUBROUTINE lu_decomp

SUBROUTINE gauss_seidel()

	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: past_x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: temps, sum, error,tol

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

	close(10)

	write(*,*) "Enter your tolerance"
	read(*,*)tol

	do

	!********** EVALUATE EQUATIONS **********!
	do i=1,n
		temps = 0
		past_x = x
		do j=1,n
			if (j /= i) then
				temps = temps + matrix(i,j)*x(j)
			endif
		end do
		x(i) = (b(i) - temps)/matrix(i,i)
	end do

	!********** EVALUATE EQUATIONS **********!


	!********** CALCULATE AND COMPARE ERRORS **********!
	sum = 0
	do i=1,n
		error = x(i) - past_x(i)
		sum = sum + error
	enddo

	if (abs(sum)<= tol) then
		EXIT
	endif

	!********** CALCULATE AND COMPARE ERRORS **********!

 	end do

 	!********** EXPORT TO CSV **********!
	open (unit = 2, file = "gauss_seidel.csv")
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
	write(2,*)"Xsub",(",",j,j=1,n)
	write(2,*)"-", (",",x(i), i=1,n)
	close(2)
	!********** EXPORT TO CSV **********!
	!********** READ FILE **********!		!####### DONE -------------->
END SUBROUTINE gauss_seidel

!********** Interpolation **********!

SUBROUTINE power_series()
	!********** READ FILE **********!
	INTEGER :: n,i,j,k,l
	DOUBLE PRECISION, dimension (:,:), allocatable :: lm
	DOUBLE PRECISION, dimension (:,:), allocatable :: um
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: up, down,r, res

	open(unit = 10, file = "interpolation.txt")
	read(10,*)n,r
	allocate ( y(n) )
	allocate ( x(n) )
	allocate ( b(n) )
	allocate ( matrix(n,n) )
	allocate ( lm(n,n) )
	allocate ( um(n,n) )

	write(*,*) "************* Read Data *************"
	do i=1,n
		read(10,*)x(i),b(i)
		write(*,*)x(i),b(i)
	enddo

	write(*,*) "************* Number to interpolate: ", r
	write(*,*) "************************************************"

	do i=1,n
		do j=1,n
			if (j/=1) then
				matrix(i,j) = x(i)**(j-1)
			endif
		enddo
	end do

	matrix(:,1) = 1

	write(*,*) "************* MATRIX TO SOLVE *************"
	do i=1,n
		write(*,*) matrix(i,:)
	enddo

	close(10)

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

	res = 0
	do i=1,n
		if (i /=1) then
			res = res + (x(i) *(r**(i-1)))
		else
			res = res + x(i)
		endif
	enddo

	write(*,*)res

	!********** EXPORT TO CSV **********!
	open (unit = 2, file = "power_series.csv")
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
	write(2,*) "POLYNOMIAL:"
	write(2,*) x(1), "+",(x(j), "x^", (j-1),j=2,n)
	write(2,*)"RESULT:"
	write(2,*) res
	close(2)
	!********** EXPORT TO CSV **********!


END SUBROUTINE power_series

SUBROUTINE lagrange()
	!********** READ FILE **********!
	INTEGER :: n,i,j,k,l
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: poly
	DOUBLE PRECISION :: up, down,r, res

	open(unit = 10, file = "interpolation.txt")
	read(10,*)n,r
	allocate ( y(n) )
	allocate ( x(n) )
	allocate ( poly(n) )

	write(*,*) "************* Read Data *************"
	do i=1,n
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo

	write(*,*) "************* Number to interpolate: ", r

	close(10)

	do i=1, n
		up = 1
		down = 1
		do j=1, n
			if(j /= i) then
				up = up*(r-x(j))
				down = down*(x(i)-x(j))
			end if
		end do!DO DEL J
		poly(i) = up/down
	end do



	res = dot_product(poly(:),y(:))
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "

	!********** EXPORT TO CSV **********!
	open (unit = 2, file = "lagrange.csv")
	write(2,*) "RESPUESTA:"
	write(2,*) res
	close(2)
	!********** EXPORT TO CSV **********!

	!********** READ FILE **********!				!####### DONE -------------->
END SUBROUTINE lagrange

SUBROUTINE newton()

	!********** READ FILE **********!
	INTEGER :: n,n_1,i,j,k,l
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:,:), allocatable :: diff
	DOUBLE PRECISION :: temp, fin,r

	open(unit = 10, file = "interpolation.txt")
	read(10,*)n,r
	n_1 = n-1
	allocate ( y(n) )
	allocate ( x(n) )
	allocate ( diff(n,n) )

	write(*,*) "************* Read Data *************"
	do i=1,n
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo
	write(*,*) "************* Number to interpolate: ", r

	close(10)

	diff(:,1) = y(:)

	do j = 2,n
		do i=1, (n-j+1)
			diff(i,j) = (diff(i + 1, j - 1) - diff(i, j - 1)) / (x(i + j - 1) - x(i));
		enddo
	enddo


	temp = 1
	fin = 0

	do i=2,n
		do j=1,i-1
			temp = temp*(r-x(j))
		enddo
		fin = fin + temp*diff(1,i)
	enddo

	fin = fin + y(1)
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "

	!********** EXPORT TO CSV **********!
	open (unit = 2, file = "newton_div_diff.csv")
	write(2,*)"DIVIDED DIFFERENCES"
	do i=1, n
		write(2,*) ( diff(i,j),"," ,j=2,n)
	enddo
	write(2,*) "RESPUESTA:"
	write(2,*) fin
	close(2)
	!********** EXPORT TO CSV **********!


END SUBROUTINE newton

!********* FILE ************!
