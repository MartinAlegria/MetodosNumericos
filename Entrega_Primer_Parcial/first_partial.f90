PROGRAM first_partial
  IMPLICIT NONE
  DOUBLE PRECISION :: a,b,c,f,guess,p,f_prime,temp, p_1
  REAL :: tol, iters, start
  INTEGER:: selection
  LOGICAL:: interval
  CHARACTER::loop
      
   DO
      write (*,*) "WELCOME!"
      write (*,*) "WHAT METHOD DO YOU WANT TO USE [TYPE THE NUMBER OF THE DESIRED METHOD]"
      write (*,*) "1) BISECTION"
      write (*,*) "2) FALSE POSITION"
      write (*,*) "3) NEWTON RAPHSON"
      write (*,*) "4) SECANT"
      read (*,*) selection

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
      
      write (*,*) "DO YOU WANT TO TRY ANOTHER METHOD ? [Y/N]"
      read(*,*) loop
      if(loop == "N" .OR. loop == "n") then
         exit
      end if

   END DO
  
END PROGRAM first_partial

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
 END DO 
  
  
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
   END DO

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
   
   END DO

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
      
   END DO

END SUBROUTINE secant

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
  f = 5 + 2*x - 3*(x**2) + 4*(x**3)
END FUNCTION f

DOUBLE PRECISION FUNCTION f_prime(x)
  IMPLICIT NONE
  DOUBLE PRECISION :: x
  f_prime = 2 - 6*x + 12*(x**2)
END FUNCTION f_prime
