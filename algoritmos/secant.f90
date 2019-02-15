! SECANT METHOD
! GUESS -> INITIAL GUESS
PROGRAM secant
  IMPLICIT NONE
  REAL :: guess, p,p_1, f, tol, iters, start, temp

  write (*,*) "Enter your guess:"
  read (*,*) guess
  print *,"Enter the tolerance:"
  read (*,*) tol
  print *, "Enter the number of iterations desired:"
  read (*,*) iters

  p = guess
  start  = 0

  write (*,*) "P: ", p
  write (*,*) "F(P):", f(p)
  if (f(p) <= tol) then
     write (*,*) "The result is:", p
     write (*,*) "Iters: ", start
  end if
  p_1 = p
  p = p * 0.99
  
  DO
     write (*,*) "P: ", p
     write (*,*) "F(P):", f(p)
     if (f(p) <= tol .AND. f(p) >=0) then
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
  
END PROGRAM secant

REAL FUNCTION f(x)
  IMPLICIT NONE
  REAL :: x
  f = 5*x**5 - 4*x**4 +2*x**3 -3*x**2 + 8*x - 15
END FUNCTION f
