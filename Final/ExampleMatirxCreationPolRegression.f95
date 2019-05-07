 !
 !
 ! 1 THIS CAN BE USED TO GENERATE THE SYSTEM OF EQUATIONS FOR POLYNOMIAL REGRESSION
 ! 2 "poldegree" is the degree of the polynomial
 ! 3 "poldegree2" is 2 times poldegree + 1
 ! 4 row is the number of points, in this example is not related to the size of tha matrix to be solved
 ! that size would be related to "poldegree"
 ! 5 this code works for any polynomial degree, including linear
 ! 6 rememember that it doesn't make a lot of sense doing regressions higher than 5th degree.
 ! 7 This code could be translated to any programming language like Python.
 ! 8 This code doesn't include the "reading" of the data, that would be similar to the one you developed for
 ! Interpolation methods.
 !********** READ FILE **********!
PROGRAM POLREGRESSION

 INTEGER :: row,i,j, poldegree,deg1,poldegree2,row_1
 DOUBLE PRECISION, dimension (:), allocatable :: y
 DOUBLE PRECISION, dimension (:,:), allocatable :: c
 DOUBLE PRECISION, dimension (:), allocatable :: b
 DOUBLE PRECISION, dimension (:,:), allocatable :: A


 open(unit = 10, file = "file.txt")
 read(10,*)row,r
 row_1=row+1
 allocate ( y(deg1) )
 allocate ( c(row,row) )
 allocate(  b(deg1))

 write(*,*) "PASAME EL DIGRI"
 read(*,*)poldegree
 deg1 = poldegree+1
 allocate(A(deg1,deg1))
row_1=row-1

 write(*,*) "************* Read Data *************"
 do i=1,row_1
   read(10,*)c(0,i),c(1,i)
   write(*,*)c(0,i),c(1,i)
 enddo

 close(10)

    poldegree2 = poldegrre*2 + 1

		do i = 1,poldegree + 1
			b(i) = 0
		end do

		do i = 1,poldegree2
			y(i) = 0
		end do

		do i = 1, poldegree2!Aqui hace los valores de la matrix
			do j=1,row
				y(i) = y(i)+ c(j,1)**(i-1)
			end do
		end do

		do i = 1, poldegree + 1
			do j=1,row
				b(i) = b(i)+ c(j,1)**(i-1)*c(j,2)
			end do
		end do

		Do i = 1, poldegree+1	! here we calculate the coefficients for the matrix
			Do j= 1, poldegree+1
				A(i,j) = y(i+j-1)
			End do
		End do

    	print *, 'this is the system of equations'
    	write (15,*) 'This is the system of equations to be solved'
   		do i = 1,poldegree+1
      		print *, (A(i,j),j=1,poldegree+1), B(i)! This is the system of equations
      		write (15,*) (A(i,j),j=1,poldegree+1), B(i)
    	end do

    END PROGRAM POLREGRESSION
