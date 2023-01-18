program test
    use util
    implicit none
 
    integer, parameter :: matrix_size = 3
    double precision, dimension(matrix_size,matrix_size) :: A, B
    double precision, dimension(matrix_size) :: alphar, alphai, beta, alphar_d, alphai_d, beta_d
    double precision, parameter :: M = 1.1d5
    double precision, parameter :: eps = 1.1d-8
    double precision, parameter :: d = 1.d-2
    real, dimension(matrix_size) :: alphar_s, alphai_s, beta_s
    character(len=40) :: arg, dir
 
    integer :: i, niter, hdeflation, tdeflation
 
    ! Get the name of the directory of the executable
    call get_command_argument(0,arg)
    dir = arg( 1: index( arg, "/", .TRUE. ) )
 
    ! Test diagonizable matrices
    open (unit = 1, file = trim(dir) // "data_special.csv")

    call dlaset( 'G', matrix_size, matrix_size, 0.0d0, 1.0d0, A, matrix_size )
    call dlaset( 'G', matrix_size, matrix_size, 0.0d0, 1.0d0, B, matrix_size )

    A(2,1) = eps
    A(3,2) = eps
    A(1,1) = 1.0d0
    A(2,2) = (1 + d)
    A(3,3) = (1 + 2*d)
    A(1,2) = M
    A(2,3) = M

    call test_eigenvalues_double(A,B, 3, 1, alphar_d, alphai_d, beta_d, niter)
    call sort_eigenvalues(alphar_d, alphai_d, beta_d)

    A(1,2) = A(1,2)/M
    A(2,2) = A(2,2)/M
    A(3,2) = A(3,2)/M
    B(2,2) = B(2,2)/M


    A(1,3) = A(1,3)/M
    A(2,3) = A(2,3)/M
    A(3,3) = A(3,3)/M
    B(3,3) = B(3,3)/M

    write(*,*) "eigenvalues in double precision: "
    do i=1,matrix_size
        write(*,*) complex( alphar_d(i)/beta_d(i), alphai_d(i)/beta_d(i) )
    end do

    do hdeflation = 1,3
        do tdeflation = 1,1

            call test_eigenvalues_real(A,B, hdeflation, tdeflation, alphar_s, alphai_s, beta_s, niter)
            call sort_eigenvalues_real(alphar_s, alphai_s, beta_s)

            write(1,'(I8,",",ES16.6)',advance="no") niter, maxval( abs(alphar_s/beta_s - alphar_d/beta_d)/abs(alphar_d/beta_d) )
            if(hdeflation .lt. 3 .or. tdeflation .lt. 2) write(1,'(A1)',advance="no") ","

            write(*,*) "eigenvalues in single precision: ", hdeflation
            do i=1,matrix_size
                write(*,*) complex( alphar_s(i)/beta_s(i), alphai_s(i)/beta_s(i) )
            end do
        end do
    end do
    write(1,*) '\n'
 
    close(1)
 
 contains
 
 
 
 end program
 