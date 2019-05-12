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

	read(*,*)selection

	SELECT CASE(selection)
		CASE(1) ! *********** NON LINEAR EQUATIONS!
			Write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) BISECTION"
			write(*,*) " 2.) FALSE POSITION"
			write(*,*) " 3.) NEWTON-RALPHSON"
			Write(*,*) " 4.) SECANT"
			

		CASE(2)! ***********  SYSTEMS OF LINEAR EQUATIONS!
			
			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) GAUSSIAN ELIMINATION"
			write(*,*) " 2.) LU DECOMPOSITION"
			write(*,*) " 3.) GAUSS-SEIDEL"

		CASE(3)! *********** INTERPOLATION!
		
			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) POWER SERIES"
			write(*,*) " 2.) LAGRANGE"
			write(*,*) " 3.) NEWTON DIVIDED DIFFERENCES"

		CASE(4)! *********** REGRESSION!

		CASE(5)! *********** NUMERICAL INTEGRATION!

			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) TRAPEZOID"
			write(*,*) " 2.) SIMPSON 1/3"
			write(*,*) " 3.) SIMPSON 3/8"

		CASE(6)! *********** ORDINARY DIFFERENTIAL EQUATIONS!

			write(*,*) " WHAT METHOD DO YOU WANT TO USE ? "
			write(*,*) " 1.) EULER"
			write(*,*) " 2.) MODIFIED EULER"
			write(*,*) " 3.) RANGE KUTTEN 3rd ORDER"
			write(*,*) " 4.) RANGE KUTTEN 4th ORDER"
	END SELECT

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
   !EL METODO
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

	open(unit = 10, file = "third.csv")
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

	!EVALUA LA ECUACION
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

	INTEGER:: op, inters, i
	DOUBLE PRECISION:: tol,orig,iters,upper_int, lower_int, diff,h,f,res, old, new, maxiters,c,c1,c2,err
	LOGICAL:: converged
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
		write(*,*) "THE ANSWER IS = ", real(old)
	endif

END SUBROUTINE simpson1_3

SUBROUTINE simpson3_8()
	INTEGER:: op, inters, i
	DOUBLE PRECISION:: tol,orig,iters,upper_int, lower_int, diff,h,f,res, old, new, maxiters,c,c1,c2,c3,err
	LOGICAL:: converged
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
		write(*,*) "THE ANSWER IS = ", real(old)
	endif

END SUBROUTINE simpson3_8

!********** DIFFERENTIAL EQUATIONS **********!

SUBROUTINE euler()

	DOUBLE PRECISION:: y, initial_x, h, aprox_x,df, old, new,errs,tol
	LOGICAL:: converged

	y = 0.00000
	initial_x = 0.000000
	h = 1
	aprox_x = 120
	tol = 0.0

	old = 0
	write (*,*) " ################# EULER METHOD #################"
	converged = .FALSE.

	DO WHILE( initial_x < aprox_x)
		write(*,*)  df(initial_x, y)
		y = y + h * df(initial_x, y)
		initial_x = initial_x + h

		new = y
		err = abs( (new-old)/new )
		write(*,*) "ERR = ", err, " ---",tol

		if (err <= tol) then
			write(*,*) "CONVERGED"
			converged = .TRUE.
		endif

		old = new

	END DO

	write(*,*) "X: ",initial_x, "Y:", y

END SUBROUTINE euler

SUBROUTINE mod_euler()

	DOUBLE PRECISION:: y, initial_x, h, aprox_x,df,k1,k2,tol

	y = 0.000000
	initial_x = 0.000000
	h = 120
	aprox_x = 120

	write (*,*) " ################# MODIFIED EULER METHOD #################"
	write (*,*) " GIVE ME THE TOLERANCE:"
	read(*,*) tol

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
  f = (  1/SQRT(2*(4.D0*DATAN(1.D0))) * exp(-0.5*(x**2)))
END FUNCTION f

DOUBLE PRECISION FUNCTION f_prime(x)
  IMPLICIT NONE
  DOUBLE PRECISION :: x
  f_prime = 30194.50584*(x**2)+1.5
END FUNCTION f_prime

DOUBLE PRECISION FUNCTION df(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION :: x,y
  df = (-3200)/(2*(4.D0*DATAN(1.D0))*(x))
END FUNCTION df




