PROGRAM SECOND_PARTIAL

	IMPLICIT NONE
	INTEGER :: selection
	CHARACTER::loop

	CALL simpson1_3()
	CALL simpson3_8()

	

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

!********** NUMERICAL INTEGRATION **********!

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

		write(*,*) "Intervals" , x(:)
		write(*,*) "FX" , f_x(:)
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






