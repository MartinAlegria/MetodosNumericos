PROGRAM false_position

  IMPLICIT NONE
  REAL :: a,b,c,f,tol,iters,start
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

  CALL fp(a,b,c,f,tol,iters,start)
  
END PROGRAM false_position


SUBROUTINE intervalos(a,b)
  IMPLICIT NONE
  REAL :: a,b
  print *,"Extremo menor del intervalo:"
  read (*,*) a
  print *,"Extremo mayor del intervalo:"
  read (*,*) b
END SUBROUTINE intervalos

SUBROUTINE fp(a,b,c,f,tol,iters,start)
  IMPLICIT NONE
  REAL:: a,b,c,f,tol,iters,start

  DO
     c= a + ( (f(a)*(b-a))/(f(a)-f(b)) )
     if(f(c) .GT. 0)then
        b = c
     else
        c = c
     end if

     if( (b-a) <= tol .AND. (b-a) >= 0 )then
        write (*,*) "The result is: ", c
        write (*,*) "Iters:", start
        EXIT
     end if

     if (start == iters) then
        write (*,*) "The result is: ", c
        write (*,*) "Iters:", start
        EXIT
     end if

     start = start + 1
  END DO
  
END SUBROUTINE fp

REAL FUNCTION f(x)
  IMPLICIT NONE
  REAL :: x
 
  f = 0.654*x*(1-exp(-(135/x)))-35

END FUNCTION f
