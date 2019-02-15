PROGRAM newton_raph
  IMPLICIT NONE
  REAL :: guess, p ,f,f_prime,tol,iters,start

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
     write (*,*) "F'(P):", f_prime(p)
     if (f(p) <= tol .AND. f(p) >= 0) then
        write (*,*) "The result is:", p
        write (*,*) "Iters: ", start
        EXIT
     end if

     if (start == iters) then
        write (*,*) "The result is:", p
        write (*,*) "Iters: ", start
        EXIT
     end if
     p = p - (f(p)/ f_prime(p))
     write (*,*) "New P: ", p
     write (*,*) "Iters: ", start
     write (*,*) "*************************"
     start = start + 1
     
  END DO
  
END PROGRAM newton_raph

REAL FUNCTION f(x)
  IMPLICIT NONE
  REAL :: x
  f = 5 + 2*x - 3*(x**2) + 4*(x**3)
END FUNCTION f

REAL FUNCTION f_prime(x)
  IMPLICIT NONE
  REAL :: x
  f_prime = 2 - 6*x + 12*(x**2)
END FUNCTION f_prime
