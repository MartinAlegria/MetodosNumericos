PROGRAM SECOND_PARTIAL

	IMPLICIT NONE
	INTEGER :: selection,nl,se,inter,reg,numi,ordi
	INTEGER::loop
	loop = 1

	do while (loop == 1)
		write(*,*) " ############### NUMERICAL METHODS ###############"
	write(*,*) " PLEASE BE SURE TO READ THE README.txt FILE BEFORE USING THE PRORGAM "
	write(*,*) " WHAT DO YOU WANT TO SOLVE ? "
	write(*,*) " 1.) NON LINEAR EQUATIONS"
	write(*,*) " 2.) SYSTEMS OF LINEAR EQUATIONS"
	write(*,*) " 3.) INTERPOLATION"
	write(*,*) " 4.) REGRESSION"
	write(*,*) " 5.) NUMERICAL INTEGRATION"
	write(*,*) " 6.) ORDINARY DIFFERENTIAL EQUATIONS"

	read(*,*)selection

	SELECT CASE(selection)
		CASE(1) ! *********** NON LINEAR EQUATIONS!
			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) BISECTION"
			write(*,*) " 2.) FALSE POSITION"
			write(*,*) " 3.) NEWTON-RALPHSON"
			write(*,*) " 4.) SECANT"
			read(*,*)nl

			SELECT CASE(nl)
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
			
			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) GAUSSIAN ELIMINATION"
			write(*,*) " 2.) LU DECOMPOSITION"
			write(*,*) " 3.) GAUSS-SEIDEL"
			read(*,*)se
			SELECT CASE(se)
			CASE(1)
				CALL gauss_elimination()
			CASE(2)
				CALL lu_decomp()
			CASE(3)
				CALL gauss_seidel()
			END SELECT

		CASE(3)! *********** INTERPOLATION!
		
			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) POWER SERIES"
			write(*,*) " 2.) LAGRANGE"
			write(*,*) " 3.) NEWTON DIVIDED DIFFERENCES"
			read(*,*)inter
			SELECT CASE(inter)
			CASE(1)
				CALL power_series()
			CASE(2)
				CALL lagrange()
			CASE(3)
				CALL newton()
			END SELECT

		CASE(4)! *********** REGRESSION!

			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) POLYNOMIAL"
			write(*,*) " 2.) EXPONENTIAL"
			write(*,*) " 3.) LOGARITHMIC"
			read(*,*)reg
			SELECT CASE(reg)
			CASE(1)
				CALL PolynomialRegression()
			CASE(2)
				CALL ExponentialR()
			CASE(3)
				CALL LogarithmicR()
			END SELECT

		CASE(5)! *********** NUMERICAL INTEGRATION!

			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) TRAPEZOID"
			write(*,*) " 2.) SIMPSON 1/3"
			write(*,*) " 3.) SIMPSON 3/8"
			read(*,*)numi
		SELECT CASE(numi)
		CASE(1)
			CALL trapezoid()
		CASE(2)
			CALL simpson1_3()
		CASE(3)
			CALL simpson3_8()
		END SELECT

		CASE(6)! *********** ORDINARY DIFFERENTIAL EQUATIONS!

			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) EULER"
			write(*,*) " 2.) MODIFIED EULER"
			write(*,*) " 3.) RANGE KUTTEN 3rd ORDER"
			write(*,*) " 4.) RANGE KUTTEN 4th ORDER"
		SELECT CASE(ordi)
			CASE(1)
				CALL euler()
			CASE(2)
				CALL mod_euler()
			CASE(3)
				CALL rk3()
			CASE(4)
				CALL rk4()
		END SELECT
	END SELECT	

	write(*,*) "DO YOU WANT TO TRY ANOTHER METHOD ?"
	write(*,*)"1 == yes, 0 == no"
	READ(*,*)loop
	end do

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
  
  ! CHECA SI EL INTERVALO ES VALIDO
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
  !AQUI SE HACE EL METODO E ITERA
  DO
     c = (a+b)/2.0
     if(f(c) .GT. 0) then
        b = c
     else
        a = c
     end if

     !SI CONVERGE SE SALE
     if (abs(b-a) <= tol ) then
        write (*,*) "The result is: ", c
				write (*,*) "Iters:", start	
				write(*,*) "ERROR -- ", abs(b-a)
        EXIT
     end if

     if (start == iters) then
        write (*,*) "The result after all iterations is: ", c
				write (*,*) "Iters:", start
				write(*,*) "ERROR -- ", abs(b-a)
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

	! CHECA SI EL INTERVALO ES VALIDO
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

   ! EL METODO
   DO
   c=b-(f(b)*(b-a))/(f(b)-f(a))
   if(f(c) .GT. 0)then
      b = c
   else
      a = c
   end if

   !SI CONVERGE SE SALE
   if( abs(f(c)) <=  tol )then
      write (*,*) "The result is: ", c
			write (*,*) "Iters:", start
			write(*,*) "ERROR -- ", abs(f(c))
      EXIT
   end if

   if (start == iters) then
      write (*,*) "The result After all the iterations is: ", c
			write (*,*) "Iters:", start
			write(*,*) "ERROR -- ", abs(f(c))
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
   !EL METODO
   DO
      write (*,*) "P: ", p
      write (*,*) "F(P):", f(p)
      write (*,*) "F'(P):", f_prime (p)
      if ( abs(f(p)) <= tol) then
         write (*,*) "The result is:", p
				 write (*,*) "Iters: ", start
				 write(*,*) "ERROR -- ", abs(f(p))
         EXIT
      end if

      if (start == iters) then
         write (*,*) "The result after all iterations is:", p
				 write (*,*) "Iters: ", start
				 write(*,*) "ERROR -- ", abs(f(p))
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
				 write(*,*) "ERROR -- ", abs(f(p))
         EXIT
      end if

      if (start == iters) then
         write (*,*) "The result is:", p
				 write (*,*) "Iters: ", start
				 write(*,*) "ERROR -- ", abs(f(p))
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
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix,clone
	DOUBLE PRECISION, dimension (:,:), allocatable :: calcs
	DOUBLE PRECISION, dimension (:), allocatable :: results
	DOUBLE PRECISION :: temps
	character (len=20) :: file_read, file_write

	write (*,*) " ################# GAUSSIAN ELIMINATION #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write
	open(unit = 10, file = file_read)
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

	clone = matrix

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
	open (unit = 2, file = file_write)
	write(2,*)"MATRIX TO BE SOLVED:"
	do i=1,n
		write(2,*)clone(i,:)
	enddo
	write(2,*)"Xsub",(",",j,j=1,n)
	write(2,*)"-", (",",results(i), i=1,n)
	close(2)
	!********** EXPORT TO CSV **********!

	!********** BACKWARDS SUBSTITUTION **********! !*********** DONE ***********!  !####### DONE --------------> !---------------> GAUSS ELIMINATION!
END SUBROUTINE gauss_elimination

SUBROUTINE lu_decomp()
	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l, choice 
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix,clone
	DOUBLE PRECISION, dimension (:,:), allocatable :: lm
	DOUBLE PRECISION, dimension (:,:), allocatable :: um
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: temps, val
	character (len=20) :: file_read, file_write
	
	write (*,*) " ################# LU - DECOMPOSITION #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	open(unit = 10, file = file_read)
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

	clone = matrix

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

	choice = 1
	do while (choice == 1)
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
		open (unit = 2, file = file_write)
		write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
		do i=1,n
			write(2,*)clone(i,:)
		enddo
		write(2,*)"Xsub",(",",j,j=1,n)
		write(2,*)"-", (",",x(i), i=1,n)
		close(2)
		!********** EXPORT TO CSV **********!

		write(*,*) "DO YOU WANT TO TRY ANOTHER RIGH HAND SIDE ?"
		write(*,*) "[1 = YES, 0 = NO]"
		read(*,*)choice

		if(choice == 1) then
			do i = 1, n
				write(*,*) "Enter the ", i, "value"
				read(*,*)val
				b(i) = val
			end do
		
			write(*,*)"THIS IS YOUR VECTOR"
			write(*,*)b(:)

		end if 
		
	end do

 !*********** DONE ***********!			!####### DONE --------------> !---------------> LU-DECOMPOSITION!
END SUBROUTINE lu_decomp

SUBROUTINE gauss_seidel()

	!********** READ FILE **********!
	INTEGER :: n, n_1,i,j,k,l 
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix,clone
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: past_x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: temps, sum, error,tol
	character (len=20) :: file_read, file_write

	write (*,*) " ################# GAUSS_SEIDEL #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	open(unit = 10, file = file_read)
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

	clone = matrix

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
	open (unit = 2, file = file_write)
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
	do i=1,n
		write(2,*)clone(i,:)
	enddo
	write(2,*)"Xsub",(",",j,j=1,n)
	write(2,*)"-", (",",x(i), i=1,n)
	close(2)
	!********** EXPORT TO CSV **********!
	!********** READ FILE **********!		!####### DONE --------------> !---------------> GAUSS-SEIDEL!
END SUBROUTINE gauss_seidel

!********** Interpolation **********!

SUBROUTINE power_series()
	!********** READ FILE **********!
	INTEGER :: n,i,j,k,l ,choice
	DOUBLE PRECISION, dimension (:,:), allocatable :: lm
	DOUBLE PRECISION, dimension (:,:), allocatable :: um
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: b
	DOUBLE PRECISION :: up, down,r, res

	character (len=20) :: file_read, file_write

	write (*,*) " ################# POWER SERIES #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write


	open(unit = 10, file = file_read)
	read(10,*)n
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

	open (unit = 2, file = file_write)
	write(2,*) "POLYNOMIAL:"
	write(2,*) x(1), "+",(x(j), "x^", (j-1),j=2,n)

	choice = 1

	do while(choice == 1)
	
		write(*,*) "AT WHAT POINT DO YOU WANT TO EVALUATE"
		read(*,*) r

		!EVALUA LA ECUACION
		res = 0
		do i=1,n
			if (i /=1) then
				res = res + (x(i) *(r**(i-1)))
			else
				res = res + x(i)
			endif
		enddo

		write(2,*)	"X: ", r , "Y: ", res
		write(*,*) "Point ", r, " RES = ", res

		write(*,*) "DO YOU WANT TO EVALUATE ANOTHER POINT ?"
		write(*,*) "[ 1 == YES, 0 == NO]"
		read(*,*) choice


	end do

	!********** EXPORT TO CSV **********!
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
	close(2)
	!********** EXPORT TO CSV **********! !---------------> POWER SERIES!
END SUBROUTINE power_series

SUBROUTINE lagrange()
	!********** READ FILE **********!
	INTEGER :: n,i,j,k,l,choie
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: poly
	DOUBLE PRECISION :: up, down,r, res

	character (len=20) :: file_read, file_write

	write (*,*) " ################# LAGRANGE #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	open(unit = 10, file = file_read)
	read(10,*)n
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
	open (unit = 2, file = file_write)

	choice = 1
	do while(choice == 1)
	write(*,*) "AT WHAT POINT DO YOU WANT TO EVALUATE"
	read(*,*) r

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
		write(2,*)	"X: ", r , "Y: ", res
		write(*,*) "Point ", r, " RES = ", res

		write(*,*) "DO YOU WANT TO EVALUATE ANOTHER POINT ?"
		write(*,*) "[ 1 == YES, 0 == NO]"
		read(*,*) choice
	end do
	
	write(*,*) "************* RESULTS EXPORTED TO CSV ************* "
	close(2)
	!********** EXPORT TO CSV **********!

	!********** READ FILE **********!				!####### DONE --------------> !---------------> LAGRANGE!
END SUBROUTINE lagrange

SUBROUTINE newton()

	!********** READ FILE **********!
	INTEGER :: n,n_1,i,j,k,l,choice
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:,:), allocatable :: diff
	DOUBLE PRECISION :: temp, fin,r
	character (len=20) :: file_read, file_write

	write (*,*) " ################# NEWTON DIVIDED DIFFERENCES #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	open(unit = 10, file = file_read)
	read(10,*)n
	n_1 = n-1
	allocate ( y(n) )
	allocate ( x(n) )
	allocate ( diff(n,n) )

	write(*,*) "************* Read Data *************"   
	do i=1,n
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo

	close(10)

	diff(:,1) = y(:)

	do j=2,n
		do i=1, (n-j+1)
			diff(i,j) = (diff(i + 1, j - 1) - diff(i, j - 1)) / (x(i + j - 1) - x(i));
		enddo
	enddo

	open (unit = 2, file = file_write)


	choice = 1
	do while (choice == 1)
		
		write(*,*) "WHICH NUMBER DO YOU WANT TO INTERPOLATE:"
		read(*,*)r 
		temp = 1
		fin = 0
	
		do i=2,n
			do j=1,i-1
				temp = temp*(r-x(j))
			enddo
			fin = fin + temp*diff(1,i)
		enddo
	
		fin = fin + y(1)
		write(2,*) "X = ", r, " Y = ", fin
		write(*,*) "DO YOU WANT TO TRY ANOTHER NUMBER ?:"
		write(*,*) "1 == YES, 0 == NO"
		read(*,*)choice

	end do

	!********** EXPORT TO CSV **********!
	write(2,*)"DIVIDED DIFFERENCES"
	do i=1, n
		write(2,*) ( diff(i,j),"," ,j=2,n)
	enddo
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
	character (len=20) :: file_read, file_write

	write (*,*) " ################# EXPONENTIAL REGRESSION #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	open(unit = 10, file = file_read)
	read(10,*)row
	allocate ( y(row) )
	allocate ( x(row) )

	write(*,*) "************* Read Data *************"
	do i=1,row
		write(*,*)i
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
			Write (*,*) 'So the equations are: ',a0,'*e',a1
	!
	! Evaluation of points

		open(unit = 2, file = file_write)
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
				endif
				write (2,*) "X: ",xint,"Y: ",yint
			end do
			close(2)
END SUBROUTINE ExponentialR

SUBROUTINE LogarithmicR()

	!********** READ FILE **********!
	INTEGER :: row,i,j,k,l,contpoint
	DOUBLE PRECISION, dimension (:), allocatable :: y
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION :: columns, ymed, st, sr, sumx, sumy, sumxy, sumx2

		character (len=20) :: file_read, file_write

	write (*,*) " ################# LOGARITHMIC REGRESSION #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	open(unit = 10, file = file_read)
	read(10,*)row
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

		open(unit = 2, file = file_write)

			Do while (contpoint == 1)
				Print *, 'press 1 if you want to calculate the regression value for a point, 0 if not'
				read *, contpoint
				if (contpoint == 1) then
					print *, 'what value do you want to calculate'
					read *, xint
					yint = a0*exp(a1*xint)
					write (*,*) 'the value for x = ',xint, ' is ',yint
					write (15,*) xint,yint
					write (2,*) "X: ",xint,"Y: ",yint
				endif
			end do
			close(2)
END SUBROUTINE LogarithmicR

SUBROUTINE PolynomialRegression()

	INTEGER :: row,i,j,k,l,point, degree,size,lol,n,choice,r
	DOUBLE PRECISION, dimension (:), allocatable :: y, y_otra
	DOUBLE PRECISION, dimension (:), allocatable :: x
	DOUBLE PRECISION, dimension (:), allocatable :: a,util,fill
	DOUBLE PRECISION, dimension (:,:), allocatable :: matrix,clone
	DOUBLE PRECISION, dimension (:,:), allocatable :: lm
	DOUBLE PRECISION, dimension (:,:), allocatable :: um
	DOUBLE PRECISION :: columns, ymed, st, sr, sumx, sumy, sumxy, sumx2,sum,temps,res
	character (len=20) :: file_read, file_write

	write (*,*) " ################# POLYNOMIAL REGRESSION #################"

	write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
	write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
	read(*,*)file_read
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	open(unit = 10, file = file_read)
	read(10,*)row
	allocate ( y(row) )
	allocate ( x(row) )

	write(*,*) "************* Read Data *************"
	do i=1,row
		read(10,*)x(i),y(i)
		write(*,*)x(i),y(i)
	enddo
	close(10)

	n = row
	write(*,*) "INPUT THE DEGREE"
	read(*,*)degree
	size = degree+1
	lol = 2*degree

	allocate( matrix(size,size) )
	allocate( lm(size,size) )
	allocate( um(size,size) )
	allocate(a(size))
	allocate(y_otra(size))
	allocate( util(size))
	allocate( fill(lol))


	!MATRIX
	do i = 1,lol
		sum = 0
		do j = 1, row
			sum = sum + x(j)**i
		end do
		fill(i) = sum
	end do

	write(*,*)fill(:)

	matrix(1,1) = row
	matrix(1,2:) = fill(1:)
	do i = 2, size
		matrix(i,:) = fill((i-1):)
	end do

	do i = 1, size
		write(*,*)matrix(i,:)
	end do

	!right hand
	sum = 0
	do i = 1, size
		do j = 1, row
			sum = sum + (  (x(i)**(j-1) )*y(j) )
		end do
		util(i) = sum
	end do

	write(*,*)util(:)

	lm(:,1) = matrix(:,1)
	um(1,:) = matrix(1,:)/lm(1,1)
	um(1,1) = 1;

	clone = matrix

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

	do i=1,size
		temps = util(i)
		do j=1,i-1
			temps = temps - lm(i,j)*y(j)
		enddo
		y_otra(i) = temps/lm(i,i)
	enddo

	write(*,*) "Y MATRIX *********************"
	write(*,*) y_otra(:)

	a(size) = (y_otra(size)/um(size,size))!--- Get Xsub1 first
		write(*,*) a(size)
		do i=size-1, 1, -1
			temps = 0
			do j=i+1, n
				temps = temps + um(i,j)*a(j)
			enddo
			a(i) = (y_otra(i) - temps)/um(i,i)
		enddo

		write(*,*) "RESULTS --- ", a(:)

		choice = 1
		open(unit = 2, file = file_write)

		do while(choice == 1)
			write(*,*) "WHAT POINT DO YOU WANT TO EVALUATE:"
			read(*,*)r
			res = 0
			do i = 1, size
				res = res + a(i)*(r**i)
			end do
			write(*,*) "X: ", r, "Y: ", res
			write(2,*) "X: ", r, "Y: ", res
			write(*,*) "DO YOU WANT TO EVALUATE ANOTHER POINT ? "
			write(*,*) "1 == YES, 0 == NO"
			read(*,*)choice
		end do
		close(2	)

END SUBROUTINE PolynomialRegression

!********** NUMERICAL INTEGRATION **********!

SUBROUTINE trapezoid()

	INTEGER:: op,lol
	DOUBLE PRECISION:: upper_int, lower_int, diff,h,f,res
	DOUBLE PRECISION, dimension (:), allocatable :: f_x
	DOUBLE PRECISION, dimension (:), allocatable :: x
	character (len=20) :: file_read, file_write


	write (*,*) " ################# TRAPEZOID #################"

	write (*,*) "WILL YOU BE USING SCATTERED DATA OR A FUNCTION ?"
	write (*,*) "1) DATA"
	write (*,*) "2) FUNCTION"

	read(*,*)op

	if (op == 1) then
		
		write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
		write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
		read(*,*)file_read
		write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
		read(*,*)file_write
		open(unit = 10, file = file_read)
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
		open(unit = 10, file = file_write)
		res = res * (h/2)	
		write(*,*) "RESULT = ", real(res)
		write(2,*) "RESULT = ", real(res)
		close(2)

	else
		write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
		read(*,*)file_write
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

		open(unit = 10, file = file_write)
		write(*,*) "RESULT = ", real(res)
		write(2,*) "RESULT = ", real(res)
		close(2)

	endif

END SUBROUTINE trapezoid

SUBROUTINE simpson1_3()

	INTEGER:: op, inters, i
	DOUBLE PRECISION:: tol,orig,iters,upper_int, lower_int, diff,h,f,res, old, new, maxiters,c,c1,c2,err
	LOGICAL:: converged
	DOUBLE PRECISION, dimension (:), allocatable :: f_x
	DOUBLE PRECISION, dimension (:), allocatable :: x
	character (len=20) :: file_read, file_write

	write (*,*) " ################# SIMPSON 1/3 #################"

	write (*,*) "WILL YOU BE USING SCATTERED DATA OR A FUNCTION ?"
	write (*,*) "1) DATA"
	write (*,*) "2) FUNCTION"

	read(*,*)op

	if (op == 1) then

		write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
		write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
		read(*,*)file_read
		write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
		read(*,*)file_write
		open(unit = 10, file = file_read)
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
		open(unit = 10, file = file_write)
		write(*,*) "THE ANSWER IS = ", real(res)
		write(2,*) "THE ANSWER IS = ", real(res)
		close(2)
		
	else
		CALL intervalos(lower_int,upper_int)
		diff = upper_int-lower_int
		converged = .FALSE.
		inters = 2
		h = diff/inters
		write (*,*) "HOW MANY ITERATIONS:"
		read(*,*) maxiters
		write (*,*) "TOLERANCE:"
		read(*,*) tol
		iters = 1
		new = 0
		orig = (h/3)*(f(lower_int)+4*f((upper_int+lower_int)/2)+f(upper_int))
		old = orig
		write(*,*) "INITIAL: ", orig

		DO WHILE (converged .eqv. .FALSE. .and. iters < maxiters)
			inters = inters +2
			h = diff/inters
			c = lower_int
			do i = 1, inters/2
				c1 = c+h
				c2 = c+2*h
				new = new + f(c)+4*f(c1) + f(c2)
				c = c2
			end do
			new = (h/3)*new
			err = abs((new-old)/new)
			if (err <= tol) then
				write(*,*) "CONVERGED --"
				converged = .TRUE.
			endif
			iters = iters + 1
			old = new
			new = 0
		END DO
		open(unit = 10, file = file_write)
		write(*,*) "THE ANSWER IS = ", real(old)
		write(2,*) "THE ANSWER IS = ", real(old)
		close(2)
	endif

END SUBROUTINE simpson1_3

SUBROUTINE simpson3_8()
	INTEGER:: op, inters, i
	DOUBLE PRECISION:: tol,orig,iters,upper_int, lower_int, diff,h,f,res, old, new, maxiters,c,c1,c2,c3,err
	LOGICAL:: converged
	DOUBLE PRECISION, dimension (:), allocatable :: f_x
	DOUBLE PRECISION, dimension (:), allocatable :: x
	character (len=20) :: file_read, file_write

	write (*,*) " ################# SIMPSON 3/8 #################"

	write (*,*) "WILL YOU BE USING SCATTERED DATA OR A FUNCTION ?"
	write (*,*) "1) DATA"
	write (*,*) "2) FUNCTION"

	read(*,*)op

	if (op == 1) then
		write(*,*) "INPUT THE NAME OF THE FILE TO BE USED"
		write(*,*) "*IMPORTANT!!!* REMEMBER THAT YOUR FILE HAS TO HAVE THE FORMAT GIVEN IN THE MANUAL"
		read(*,*)file_read
		write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
		read(*,*)file_write
		open(unit = 10, file = file_read)
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
		open(unit = 10, file = file_write)
		write(*,*) "THE ANSWER IS = ", real(res)
		write(2,*) "THE ANSWER IS = ", real(res)
		close(2)

	else
		CALL intervalos(lower_int,upper_int)
		diff = upper_int-lower_int
		converged = .FALSE.
		inters = 3
		h = diff/inters
		write (*,*) "HOW MANY ITERATIONS:"
		read(*,*) maxiters
		write (*,*) "TOLERANCE:"
		read(*,*) tol
		iters = 1
		new = 0
		orig = (3*h/8)*(f(lower_int)+3*f((2*upper_int+lower_int)/3)+3*f((upper_int+2*lower_int)/3)+f(upper_int))
		write(*,*) "INITIAL: ", orig
		old = orig
		write(*,*)iters, " ", maxiters
		DO WHILE ((converged .eqv. .FALSE.) .and. (iters .lt. maxiters))
			inters = inters +3
			h = diff/inters
			c = lower_int
			do i = 1, inters/3
				c1 = c+h
				c2 = c+2*h
				c3 = c+3*h
				new = new +f(c)+3*f((2*c3+c)/3)+3*f((c3+2*c)/3)+f(c3)
				c = c3
			end do
			new = (3*h/8)*new
			err = abs((new-old)/new)
			write(*,*) err
			if (err <= tol) then
				converged = .TRUE.
				write(*,*) "CONVERGED"
			endif
			iters = iters + 1
			old = new
			new = 0
		END DO
		open(unit = 10, file = file_write)
		write(*,*) "THE ANSWER IS = ", real(old)
		write(2,*) "THE ANSWER IS = ", real(old)
		close(2)
	endif

END SUBROUTINE simpson3_8

!********** DIFFERENTIAL EQUATIONS **********!

SUBROUTINE euler()

	DOUBLE PRECISION:: y, initial_x, h, aprox_x,df, old, new,errs,tol
	LOGICAL:: converged
	character (len=20) ::file_write

	write (*,*) " ################# EULER METHOD #################"
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	write (*,*) " GIVE ME THE INITIAL Y "
	read(*,*)y
	write (*,*) " GIVE ME THE INITIAL X "
	read(*,*)initial_x
	write (*,*) " GIVE ME THE H "
	read(*,*)h
	write (*,*) " GIVE ME THE APROX X "
	read(*,*)aprox_x
	write (*,*) " GIVE ME THE TOLERANCE "
	read(*,*)tol

	old = 0	
	converged = .FALSE.
	DO WHILE ( (initial_x < aprox_x ).and.( converged .eqv. .FALSE.))
		write(*,*)  df(initial_x, y)
		y = y + h * df(initial_x, y)
		initial_x = initial_x + h

		new = y
		err = abs( (new-old)/new )
		write(*,*) "ERR = ", err, " ---",tol
		WRITE(*,*) "X: ", initial_x, " Y:", y

		if (err <= tol) then
			write(*,*) "CONVERGED"
			converged = .TRUE.
			EXIT
		endif

		old = new
	end do

	write(*,*) "SOLUTION *****************"
	write(*,*) "X: ",initial_x, "Y: ", y

	open(unit = 2, file = file_write)
	write(2,*) "SOLUTIONS"
	write(2,*) "X: ",initial_x, "Y: ", y
	write(*,*) "***************** EXPORTED TO FILE *****************"
	close(2)

END SUBROUTINE euler

SUBROUTINE mod_euler()

	DOUBLE PRECISION:: y, initial_x, h, aprox_x,df,k1,k2,tol,old,new,errs
	LOGICAL:: converged
	character (len=20) ::file_write

	write (*,*) " ################# MOD EULER METHOD #################"
	write(*,*) "INPUT THE NAME OF THE FILE TO BE EXPORTED"
	read(*,*)file_write

	write (*,*) " GIVE ME THE INITIAL Y "
	read(*,*)y
	write (*,*) " GIVE ME THE INITIAL X "
	read(*,*)initial_x
	write (*,*) " GIVE ME THE H "
	read(*,*)h
	write (*,*) " GIVE ME THE APROX X "
	read(*,*)aprox_x
	write (*,*) " GIVE ME THE TOLERANCE "
	read(*,*)tol


	old = 0	
	converged = .FALSE.
	DO WHILE ( (initial_x < aprox_x ).and.( converged .eqv. .FALSE.))
		!write(*,*) initial_x, y
		k1 = df(initial_x,y)
		k2 = df(initial_x+1, y+ (h*k1))
		y = y + (h/2)*(k1+k2)
		initial_x = initial_x + h

		new = y
		err = abs( (new-old)/new )
		write(*,*) "ERR = ", err, " ---",tol
		WRITE(*,*) "X: ", initial_x, " Y:", y

		if (err <= tol) then
			write(*,*) "CONVERGED"
			converged = .TRUE.
			EXIT
		endif

		old = new

	END DO
	write(*,*) "SOLUTION *****************"
	write(*,*) "X: ",initial_x, "Y: ", y

	open(unit = 2, file = file_write)
	write(2,*) "SOLUTIONS"
	write(2,*) "X: ",initial_x, "Y: ", y
	write(*,*) "***************** EXPORTED TO FILE *****************"
	close(2)

END SUBROUTINE mod_euler

SUBROUTINE rk3()

	DOUBLE PRECISION:: y, initial_x,initial_y, h, aprox_x,df,k1,k2,k3,n,n1
	LOGICAL:: coverged
	INTEGER::numinterval
	write (*,*) " ################# RUNGE KUTTA 3rd ORDER #################"
	
	write (*,*) " GIVE ME THE INITIAL Y "
	read(*,*)initial_y
	write (*,*) " GIVE ME THE INITIAL X "
	read(*,*)initial_x
	write (*,*) " GIVE ME THE H "
	read(*,*)h
	write (*,*) " GIVE ME THE APROX X "
	read(*,*)aprox_x

	y = initial_y!NO QUITAR ESTO O SE
	numinterval = (aprox_x - initial_x)/h

	write(*,*)h

	y = initial_y

	do i = 1, numinterval
		!write(*,*) initial_x, y
		k1 = df(initial_x,y)
		k2 = df( initial_x+(h/2), y+(h/2)*k1 )
		k3 = df( (initial_x + h), y-k1*h+2*h*k2 )
		!,yold-h*k1+2*h*k2
		y = y + (h/6)*(k1+4*k2+k3)
		write(*,*)initial_x,y
		initial_x = initial_x + h

	END DO

	write(*,*) "RES:", initial_x, y
	
END SUBROUTINE rk3

SUBROUTINE rk4()

	DOUBLE PRECISION:: y, initial_x,initial_y, h, aprox_x,df,k1,k2,k3,k4,n,n1
	INTEGER::numinterval


	write (*,*) " ################# RUNGE KUTTA 4th ORDER #################"
	write (*,*) " GIVE ME THE INITIAL Y "
	read(*,*)initial_y
	write (*,*) " GIVE ME THE INITIAL X "
	read(*,*)initial_x
	write (*,*) " GIVE ME THE H "
	read(*,*)h
	write (*,*) " GIVE ME THE APROX X "
	read(*,*)aprox_x

	y = initial_y!NO QUITAR ESTO O SE

	numinterval = (aprox_x-initial_x)/h
	write(*,*)numinterval

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
	f = 10*x + 25
END FUNCTION f

DOUBLE PRECISION FUNCTION f_prime(x)
  IMPLICIT NONE
  DOUBLE PRECISION :: x
  f_prime = 10
END FUNCTION f_prime

DOUBLE PRECISION FUNCTION df(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION :: x,y
  df = (-3200)/(2*(4.D0*DATAN(1.D0))*(x))
END FUNCTION df




