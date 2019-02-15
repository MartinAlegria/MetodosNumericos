PROGRAM bisection
  IMPLICIT NONE
  REAL :: a,b,c,f
  REAL :: tol,iters,start
  LOGICAL :: interval
  interval = .FALSE. 
  CALL intervalos(a,b)

  DO WHILE (interval .EQV. .FALSE.)
     write(*,*) f(a)
     write(*,*) f(b)
     if((f(a)*f(b)) .LT. 0) then
        print *,"Perfecto!"
        interval = .TRUE.
     else
        print *,"No hay intervalo, intentalo otra vez"
        print *, "------------------------------------"
        CALL intervalos(a,b)
     end if
  END DO

  print *,"Ingresa la tolerancia:"
  read (*,*) tol
  print *, "Ingresa el numero de iteraciones"
  read (*,*) iters

  CALL binary_search(a,b,tol,iters,f,start)
  
END PROGRAM bisection


SUBROUTINE intervalos(a,b)
  IMPLICIT NONE
  REAL :: a,b
  print *,"Extremo menor del intervalo:"
  read (*,*) a
  print *,"Extremo mayor del intervalo:"
  read (*,*) b
END SUBROUTINE intervalos

SUBROUTINE binary_search(a,b,tol,iters,f,start)
  IMPLICIT NONE
  REAL :: a,b,mid,iters,f
  REAL :: tol, start
  DO
     mid = (a+b)/2.0
     if(f(mid) .GT. 0) then
        a = mid
     else
        b = mid
     end if

     if ((b-a) <= tol .AND. (b-a) >= 0) then
        write (*,*) "The result is: ", mid
        write (*,*) "Iters:", start
        EXIT
     end if

     if (start == iters) then
        write (*,*) "The result is: ", mid
        write (*,*) "Iters:", start
        EXIT
     end if
     start = start +1
 END DO 
END SUBROUTINE binary_search

REAL FUNCTION f(x)
  IMPLICIT NONE
  REAL :: x
 
  f = (-0.6*x**2)+(2.4*x) + 5.5

END FUNCTION f
