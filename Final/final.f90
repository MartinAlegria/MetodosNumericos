PROGRAM SECOND_PARTIAL

	IMPLICIT NONE
	INTEGER :: selection
	CHARACTER::loop

	write(*,*) " ############### NUMERICAL METHODS ###############"
	write(*,*) " PLEASE BE SURE TO READ THE README.txt FILE BEFORE USING THE PRORGAM "
	write(*,*) " WHAT DO YOU WANT TO SOLVE ? "
	write(*,*) " 1.) NON LINEAR EQUATIONS"
	write(*,*) " 2.) SYSTEMS OF LINEAR EQUATIONS"
	write(*,*) " 3.) INTERPOLATION"
	write(*,*) " 4.) REGRESSION"
	write(*,*) " 5.) NUMERICAL INTEGRATION"
	write(*,*) " 6.) ORDINARY DIFFERENTIAL EQUATIONS"
  read(*,*) selection
	SELECT CASE(selection)
		CASE(1) ! *********** NON LINEAR EQUATIONS!
			write(*,*) "WHAT METHOD DO YOU WANT TO USE? "
			write(*,*) " 1.) BISECTION"
			write(*,*) " 2.) FALSE POSITION"
			write(*,*) " 3.) NEWTHON RAPHSON"
			write(*,*) " 4.) SECANT"
			read(*,*) selection
			SELECT CASE(selection)
			CASE(1)
				CALL bisection()
			CASE(2)
				CALL false_position()
			CASE(3)
				CALL newton_raph()
			CASE(4)
				CALL secant()

			END SELECT
		CASE(2)! ***********  SYSTEMS OF LINEAR EQUATIONS!


		CASE(3)! *********** INTERPOLATION!


		CASE(4)! *********** REGRESSION!


		CASE(5)! *********** NUMERICAL INTEGRATION!


		CASE(6)! *********** ORDINARY DIFFERENTIAL EQUATIONS!


	END SELECT


	!CALL simpson1_3()
	!CALL simpson3_8()
	!CALL trapezoid()
	!CALL euler()
	!CALL mod_euler()
	!CALL rk3()
	!CALL rk4()

	write(*,*) "DO YOU WANT TO USE DATA OR "

END PROGRAM SECOND_PARTIAL

!********** Non-Linear Equations **********!

SUBROUTINE bisection()
  IMPLICIT NONE
  DOUBLE PRECISION :: a,b,c,f
  REAL :: tol, iters, start
  LOGICAL:: interval
  interval = .FALSE.

  print *," #################### BISECTION ####################"
  CALL intervalos(a,b)

  DO WHILE (interval .EQV. .FALSE.)
     write(*,*) f(a)
     write(*,*) f(b)
     if((f(a)*f(b)) .LT. 0) then
        print *,"Perfect!"
        interval = .TRUE.
     else
        print *,"No interval, try again"
        print *, "------------------------------------"
        CALL intervalos(a,b)
     end if
  END DO

  print *,"Input the tolerance:"
  read (*,*) tol
  print *, "Input the number of iterations"
  read (*,*) iters

  DO
     c = (a+b)/2.0
     if(f(c) .GT. 0) then
        b = c
     else
        a = c
     end if

     if (abs(b-a) <= tol ) then
        write (*,*) "The result is: ", c
        write (*,*) "Iters:", start
        EXIT
     end if

     if (start == iters) then
        write (*,*) "The result after all iterations is: ", c
        write (*,*) "Iters:", start
        EXIT
     end if
     start = start +1
 END DO  !---------------> BISECTION!
END SUBROUTINE bisection

SUBROUTINE false_position()
   IMPLICIT NONE
   DOUBLE PRECISION :: a,b,c,f
   REAL:: tol,iters,start
   LOGICAL interval
   interval = .FALSE.

   print *," #################### FALSE POSITION METHOD ####################"
   CALL intervalos(a,b)

   DO WHILE (interval .EQV. .FALSE.)
   write(*,*) f(a)
   write(*,*) f(b)
   if((f(a)*f(b)) .LT. 0) then
      print *,"Perfect!"
      interval = .TRUE.
   else
      print *,"No interval, try again"
      print *, "---------------------------------------"
      CALL intervalos(a,b)
   end if
   END DO

   print *,"Input the tolerance"
   read (*,*) tol
   print *, "Input the number of iterations"
   read (*,*) iters

   DO
   c=b-(f(b)*(b-a))/(f(b)-f(a))
   if(f(c) .GT. 0)then
      b = c
   else
      a = c
   end if

   if( abs(f(c)) <=  tol )then
      write (*,*) "The result is: ", c
      write (*,*) "Iters:", start
      EXIT
   end if

   if (start == iters) then
      write (*,*) "The result After all the iterations is: ", c
      write (*,*) "Iters:", start
      EXIT
   end if

   start = start + 1
   END DO !---------------> FALSE POSITION!
END SUBROUTINE false_position

SUBROUTINE newton_raph()
   IMPLICIT NONE
   DOUBLE PRECISION:: guess,p,f,f_prime
   REAL:: tol,iters,start

   print *," #################### NEWTON-RAPHSON METHOD ####################"
   write (*,*) "Enter your guess:"
   read (*,*) guess
   print *,"Enter the tolerance:"
   read (*,*) tol
   print *, "Enter the number of iterations desired:"
   read (*,*) iters
   p = guess
   start = 0

   DO
      write (*,*) "P: ", p
      write (*,*) "F(P):", f(p)
      write (*,*) "F'(P):", f_prime (p)
      if ( abs(f(p)) <= tol) then
         write (*,*) "The result is:", p
         write (*,*) "Iters: ", start
         EXIT
      end if

      if (start == iters) then
         write (*,*) "The result after all iterations is:", p
         write (*,*) "Iters: ", start
         EXIT
      end if
      p = p - (f(p)/ f_prime(p))
      write (*,*) "New P: ", p
      write (*,*) "Iters: ", start
      write (*,*) "*************************"
      start = start + 1

   END DO !---------------> NEWTON-RAPHSON!
END SUBROUTINE newton_raph

SUBROUTINE secant()
   IMPLICIT NONE
   DOUBLE PRECISION:: guess, p, p_1,f,temp
   REAL:: tol,iters,start

   print *," #################### SECANT METHOD ####################"

   write (*,*) "Enter your guess (CANT BE ZERO [0]):"
   read (*,*) guess
   print *,"Enter the tolerance:"
   read (*,*) tol
   print *, "Enter the number of iterations desired:"
   read (*,*) iters

   p = guess
   start  = 0

   write (*,*) "P: ", p
   write (*,*) "F(P):", f(p)
   if (abs(f(p)) <= tol) then
      write (*,*) "The result is:", p
      write (*,*) "Iters: ", start
   end if
   p_1 = p
   p = p * 0.99

   DO
      write (*,*) "P: ", p
      write (*,*) "F(P):", f(p)
      if (abs(f(p)) <= tol) then
         write (*,*) "The result is:", p
         write (*,*) "Iters: ", start
         EXIT
      end if

      if (start == iters) then
         write (*,*) "The result is:", p
         write (*,*) "Iters: ", start
         EXIT
      end if
      temp = p
      p = p - (f(p)*(p-p_1))/(f(p)-f(p_1))
      p_1 = temp
      write (*,*) "New P: ", p
      write (*,*) "Iters: ", start
      write (*,*) "*************************"
      start = start + 1

   END DO !---------------> SECANT!
END SUBROUTINE secant

!********** System of Linear Equations **********!

SUBROUTINE gauss_elimination()
	IMPLICIT NONE
	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:,:), allocatable :: calcs
	DOUBLE PRECISION, dimension (:), allocatable :: results
	DOUBLE PRECISION :: temps

	open(unit = 10, file = "system_eq.txt")
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

	!********** BACKWARDS SUBSTITUTION **********! !*********** DONE ***********!  !####### DONE --------------> !---------------> GAUSS ELIMINATION!
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

	open(unit = 10, file = "system_eq.txt")
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
	!********** EXPORT TO CSV **********! !*********** DONE ***********!			!####### DONE --------------> !---------------> LU-DECOMPOSITION!
END SUBROUTINE lu_decomp

SUBROUTINE gauss_seidel()

	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: past_x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: temps, sum, error,tol

	open(unit = 10, file = "system_eq.txt")
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
	!********** READ FILE **********!		!####### DONE --------------> !---------------> GAUSS-SEIDEL!
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
	!********** EXPORT TO CSV **********! !---------------> POWER SERIES!
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
	write(2,*) "RESULT:"
	write(2,*) res
	close(2)
	!********** EXPORT TO CSV **********!

	!********** READ FILE **********!				!####### DONE --------------> !---------------> LAGRANGE!
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

	do j=2,n
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
	!********** EXPORT TO CSV **********! !---------------> NEWTON DIVIDED DIFERENCES!
END SUBROUTINE newton

!********** REGRESSION **********!
SUBROUTINE ExponentialR()

	!********** READ FILE **********!
	INTEGER :: row,i,j,k,l,contpoint
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION :: columns, ymed, st, sr, sumx, sumy, sumxy, sumx2

	open(unit = 10, file = "file.txt")
	read(10,*)row,columns
	allocate ( y(row) )
	allocate ( x(row) )

	write(*,*) "************* Read Data *************"
	do i=1,row
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo

	close(10)
			ymed = 0
			St = 0
			Sr = 0

			sumx = 0
			sumy = 0
			sumxy = 0
			sumx2 = 0
			do i = 1, row
				sumx = sumx + x(i)
				sumy = sumy + log(y(i))
				sumx2 = sumx2 + x(i)**2
				sumxy = sumxy + x(i)*log((y(i)))

			end do

			xmed = sumx/row
			ymed = sumy/row
			a1 = (row*sumxy-sumx*sumy)/(row*sumx2-sumx**2)
			a0 = ymed - a1*xmed
			print*, '        sumx                    sumy                    sumx2                    sumxy '
			print*, sumx,sumy,sumx2,sumxy

			print *, 'original values for '
			print*, 'a0= ',a0,' a1 = ', a1
			write (15,*) ' original values for a0 = ',a0,' a1 = ',a1
			do i = 1, row
				St = St + (log(y(i) - ymed)**2)
				Sr = Sr + (log(y(i)) - a1*x(i) - a0)**2
			end do
			R2 = (St - Sr)/St
			R = sqrt (R2)
			Print *, '        Sr                    St                    R2                    R '
			print*, Sr,St,R2,R
			a0 = exp(a0)
			Write (*,*) ' Values for the logarithmic equation: a0 = ',a0,' a1= ',a1
			Write (*,*) 'So the equations is: ',a0,'*e',a1
	!
	! Evaluation of points
	!
		contpoint=1
			Do while (contpoint == 1)
				Print *, 'press 1 if you want to calculate the regression value for a point, 0 if not'
				read *, contpoint
				if (contpoint == 1) then
					print *, 'what value do you want to calculate'
					read *, xint
					yint = a0*exp(a1*xint)
					write (*,*) 'the value for x = ',xint, ' is ',yint
					write (15,*) xint,yint
				endif
			end do
END SUBROUTINE ExponentialR


SUBROUTINE LogarithmicR()

	!********** READ FILE **********!
	INTEGER :: row,i,j,k,l,contpoint
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION :: columns, ymed, st, sr, sumx, sumy, sumxy, sumx2

	open(unit = 10, file = "file.txt")
	read(10,*)row,columns
	allocate ( y(row) )
	allocate ( x(row) )

	write(*,*) "************* Read Data *************"
	do i=1,row
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo

	close(10)
			ymed = 0
			St = 0
			Sr = 0

			sumx = 0
			sumy = 0
			sumxy = 0
			sumx2 = 0
			do i = 1, row
				sumx = sumx + log10(x(i))
				sumy = sumy + log10(y(i))
				sumx2 = sumx2 + log10(x(i))**2
				sumxy = sumxy + log10(x(i))*log10((y(i)))

			end do

			xmed = sumx/row
			ymed = sumy/row
			a1 = (row*sumxy-sumx*sumy)/(row*sumx2-sumx**2)
			a0 = ymed - a1*xmed

			print*, '        log(sumx)                   log(sumy)                    log(sumx2)                  log(sumxy) '
			print*, sumx,sumy,sumx2,sumxy

			print *, 'original values for '
			print*, 'log(a0)= ',10**a0,' a1 = ', a1
			write (15,*) ' original values for a0 = ',10**a0,' a1 = ',a1
			do i = 1, row
				St = St + (log10(y(i) - ymed)**2)
				Sr = Sr + (log10(y(i)-(a0+a1*log10(x(i)))**2))
			end do

			R2 = (St - Sr)/St
			R = sqrt (R2)
			Print *, '        Sr                    St                    R2                    R '
			print*, Sr,St,R2,R
			a0 = 10**a0
			Write (*,*) ' Values for the logarithmic equation: a0 = ',a0,' a1= ',a1
			Write (*,*) 'So the equations is: ',a0,'*e',a1
			Write (*,*) 'Sr = ',Sr, ' St =  ',St, ' R2 = ',R2, ' R = ',R
			Write (*,*) '           x              y'
	!
	! Evaluation of points
	!
		contpoint=1
			Do while (contpoint == 1)
				Print *, 'press 1 if you want to calculate the regression value for a point, 0 if not'
				read *, contpoint
				if (contpoint == 1) then
					print *, 'what value do you want to calculate'
					read *, xint
					yint = a0*exp(a1*xint)
					write (*,*) 'the value for x = ',xint, ' is ',yint
					write (15,*) xint,yint
				endif
			end do
END SUBROUTINE LogarithmicR

SUBROUTINE PolynomialRegression()

END SUBROUTINE PolynomialRegression
!********** NUMERICAL INTEGRATION **********!

SUBROUTINE trapezoid()

	INTEGER:: op,lol
	DOUBLE PRECISION:: upper_int, lower_int, diff,h,f,res
	DOUBLE PRECISION, dimension (:), allocatable :: f_x
	DOUBLE PRECISION, dimension (:), allocatable :: x

	write (*,*) " ################# TRAPEZOID #################"

	write (*,*) "WILL YOU BE USING SCATTERED DATA OR A FUNCTION ?"
	write (*,*) "1) DATA"
	write (*,*) "2) FUNCTION"

	read(*,*)op

	if (op == 1) then

		open(unit = 10, file = "file.txt")
		read(10,*)inters

		allocate(x(inters))
		allocate(f_x(inters))

		do i=1,inters
			read(10,*)x(i),f_x(i)
			write(*,*)x(i),f_x(i)
		enddo

		close(10)

		upper_int = x(inters)
		lower_int = x(1)
		diff = upper_int-lower_int
		h = (diff)/inters

		res = f_x(1) + f_x(inters)
		lol = inters - 1
		do i=2,lol
				res = res + 2 * f_x(i)
		end do

		res = res * (h/2)
		write(*,*) "RESULT = ", real(res)

	else

		CALL intervalos(lower_int,upper_int)
		diff = upper_int-lower_int
		write (*,*) "HOW MANY TRAPEZOIDS DO YOU WANT ?"
		write (*,*) "TIP: SOMETIMES USING LARGER INTERVALS GIVE A BETTER RESULT"
		read(*,*)inters
		h = diff/inters
		allocate( x(inters) )
		allocate( f_x(inters) )
		i = 0

		res = f(lower_int) + f(upper_int)
		lol = inters - 1
		do i=1,lol
				res = res + 2 * f(lower_int +i*h)
		end do

		res = res * (h/2)
		write(*,*) "RESULT = ", real(res)

	endif

END SUBROUTINE trapezoid

SUBROUTINE simpson1_3()

	INTEGER:: op, inters, tol, i
	DOUBLE PRECISION:: upper_int, lower_int, diff,h,f,res
	DOUBLE PRECISION, dimension (:), allocatable :: f_x
	DOUBLE PRECISION, dimension (:), allocatable :: x

	write (*,*) " ################# SIMPSON 1/3 #################"

	write (*,*) "WILL YOU BE USING SCATTERED DATA OR A FUNCTION ?"
	write (*,*) "1) DATA"
	write (*,*) "2) FUNCTION"

	read(*,*)op

	if (op == 1) then

		open(unit = 10, file = "file.txt")
		read(10,*)inters

		allocate(x(inters))
		allocate(f_x(inters))

		do i=1,inters
			read(10,*)x(i),f_x(i)
			write(*,*)x(i),f_x(i)
		enddo

		close(10)

		upper_int = x(inters)
		lower_int = x(1)

		diff = upper_int-lower_int
		h = (diff)/inters

		res = 0
		DO i=1,inters
			if (i==0 .OR. i==inters) then
				res = res + f_x(i)
			else if (MOD(i,2) /= 0) then
				res = res + 4*f_x(i)
			else
				res = res + 2*f_x(i)
			endif
		END DO

		res = res * (h/3)
		write(*,*) "THE ANSWER IS = ", real(res)
	else
		CALL intervalos(lower_int,upper_int)
		diff = upper_int-lower_int
		write (*,*) "HOW MANY INTERVALS DO YOU WANT ?"
		write (*,*) "TIP: SOMETIMES USING LARGER INTERVALS GIVE A BETTER RESULT"
		read(*,*)inters
		h = diff/inters
		allocate( x(inters) )
		allocate( f_x(inters) )
		i = 0
		! --  CREATE AN ARRAY WITH ALL THE VALUES OF X WE WILL EVALUATE -- !
		DO i=1,inters
			x(i) = lower_int + i*h
		END DO
		! --  CREATE AN ARRAY WITH THE EVALUATIONS OF THE VALUES OF X -- !
		DO i=1,inters
			f_x(i) = f(x(i))
		END DO
		res = 0
		DO i=1,inters
			if (i==0 .OR. i==inters) then ! -- SUM THE FIRST AND LAST TERMS --!
				res = res + f_x(i)
			else if (MOD(i,2) /= 0) then ! -- IF i IS NOT PAIR THE EVALUATION IS MULTIPLIED BY 4, ELSE IT IS MULT BY 2 --!
				res = res + 4*f_x(i)
			else
				res = res + 2*f_x(i)
			endif
		END DO

		res = res * (h/3)
		write(*,*) "THE ANSWER IS = ", real(res)
	endif

END SUBROUTINE simpson1_3

SUBROUTINE simpson3_8()
	INTEGER:: op, inters, tol, i
	DOUBLE PRECISION:: upper_int, lower_int, diff,h,f,res
	DOUBLE PRECISION, dimension (:), allocatable :: f_x
	DOUBLE PRECISION, dimension (:), allocatable :: x

	write (*,*) " ################# SIMPSON 3/8 #################"

	write (*,*) "WILL YOU BE USING SCATTERED DATA OR A FUNCTION ?"
	write (*,*) "1) DATA"
	write (*,*) "2) FUNCTION"

	read(*,*)op

	if (op == 1) then
		open(unit = 10, file = "file.txt")
		read(10,*)inters

		allocate(x(inters))
		allocate(f_x(inters))

		do i=1,inters
			read(10,*)x(i),f_x(i)
			write(*,*)x(i),f_x(i)
		enddo

		close(10)

		upper_int = x(inters)
		lower_int = x(1)

		diff = upper_int-lower_int
		h = (diff)/inters

		res = 0
		DO i=1,inters
			if (i==0 .OR. i==inters) then
				res = res + f_x(i)
			else if (MOD(i,3) == 0) then
				res = res + 2*f_x(i)
			else
				res = res + 3*f_x(i)
			endif
		END DO

		res = res * (3*h/8)
		write(*,*) "THE ANSWER IS = ", real(res)

	else
		CALL intervalos(lower_int,upper_int)
		diff = upper_int-lower_int
		write (*,*) "HOW MANY INTERVALS DO YOU WANT ?"
		write (*,*) "TIP: SOMETIMES USING LARGER INTERVALS GIVE A BETTER RESULT"
		read(*,*)inters
		h = diff/inters
		allocate( x(inters) )
		allocate( f_x(inters) )
		i = 0
		! --  CREATE AN ARRAY WITH ALL THE VALUES OF X WE WILL EVALUATE -- !
		DO i=1,inters
			x(i) = lower_int + (i*h)
		END DO
		! --  CREATE AN ARRAY WITH THE EVALUATIONS OF THE VALUES OF X -- !
		DO i=1,inters
			f_x(i) = f(x(i))
		END DO

		res = 0
		DO i=1,inters
			if (i==0 .OR. i==inters) then
				res = res + f_x(i)
			else if (MOD(i,3) == 0) then
				res = res + 2*f_x(i)
			else
				res = res + 3*f_x(i)
			endif
		END DO

		res = res * (3*h/8)
		write(*,*) "THE ANSWER IS = ", real(res)
	endif

END SUBROUTINE simpson3_8

!********** DIFFERENTIAL EQUATIONS **********!

SUBROUTINE euler()

	DOUBLE PRECISION:: y, initial_x, h, aprox_x,df

	y = 0.00000
	initial_x = 0.000000
	h = 120
	aprox_x = 120

	write (*,*) " ################# EULER METHOD #################"

	DO WHILE( initial_x < aprox_x)
		write(*,*)  df(initial_x, y)
		y = y + h * df(initial_x, y)
		initial_x = initial_x + h
	END DO

	write(*,*) "X: ",initial_x, "Y:", y

END SUBROUTINE euler

SUBROUTINE mod_euler()

	DOUBLE PRECISION:: y, initial_x, h, aprox_x,df,k1,k2

	y = 0.000000
	initial_x = 0.000000
	h = 120
	aprox_x = 120

	write (*,*) " ################# MODIFIED EULER METHOD #################"

	DO WHILE( initial_x < aprox_x)
		!write(*,*) initial_x, y
		k1 = df(initial_x,y)
		k2 = df(initial_x+1, y+ (h*k1))
		y = y + (h/2)*(k1+k2)
		initial_x = initial_x + h
	END DO

	write(*,*) "X: ",initial_x, "Y:", y

END SUBROUTINE mod_euler

SUBROUTINE rk3()

	DOUBLE PRECISION:: y, initial_x,initial_y, h, aprox_x,df,k1,k2,k3,n,n1
	INTEGER::numinterval

	initial_y = 0.000000
	initial_x = 0.000000
	numinterval = 120
	aprox_x = 120
	h = (aprox_x-initial_x)/numinterval
	n1 = n+1

	write(*,*)h

	y = initial_y
	write (*,*) " ################# RUNGE KUTTA 3rd ORDER #################"

	do i = 1, numinterval
		!write(*,*) initial_x, y
		k1 = df(initial_x,y)
		k2 = df( initial_x+(h/2), y+(h/2)*k1 )
		k3 = df( (intial_x + h), y-k1*h+2*h*k2 )
		!,yold-h*k1+2*h*k2
		y = y + (h/6)*(k1+4*k2+k3)
		!write(*,*)initial_x,y
		initial_x = initial_x + h

	END DO

	write(*,*) "RES:", initial_x, y


END SUBROUTINE rk3

SUBROUTINE rk4()

	DOUBLE PRECISION:: y, initial_x,initial_y, h, aprox_x,df,k1,k2,k3,k4,n,n1
	INTEGER::numinterval

	initial_y = 0.000000
	initial_x = 0.000000
	numinterval = 120
	aprox_x = 120
	h = (aprox_x-initial_x)/numinterval
	n1 = n+1

	y = initial_y
	write (*,*) " ################# RUNGE KUTTA 4th ORDER #################"

	do i = 1, numinterval
		!write(*,*) initial_x, y
		k1 = df(initial_x,y)
		k2 = df( initial_x+(h/2), y+(h/2)*k1 )
		k3 = df( initial_x+(h/2), y+(h/2)*k2 )
		k4 = df( initial_x +h, y +k3*h)
		!,yold-h*k1+2*h*k2
		y = y + (h/6)*(k1+2*k2+2*k3+k4)
		!write(*,*)initial_x,y
		initial_x = initial_x + h

	END DO

	write(*,*) "RES:", initial_x, y

END SUBROUTINE rk4

!********** UTILS **********!

SUBROUTINE intervalos(a,b)
  IMPLICIT NONE
  DOUBLE PRECISION :: a,b
  print *,"Lower bound of the interval:"
  read (*,*) a
  print *,"Upper bound of the interval:"
  read (*,*) b
END SUBROUTINE intervalos

DOUBLE PRECISION FUNCTION f(x)
  IMPLICIT NONE
  DOUBLE PRECISION :: x
  f = (1/(1+ x*x))
END FUNCTION f

DOUBLE PRECISION FUNCTION f_prime(x)
  IMPLICIT NONE
  DOUBLE PRECISION :: x
  f_prime = 2 - 6*x + 12*(x**2)
END FUNCTION f_prime

DOUBLE PRECISION FUNCTION df(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION :: x,y
  df = (10)/( (30000)*(  -(x**2)/(3600) + (x/30) + 1   ))
END FUNCTION df
