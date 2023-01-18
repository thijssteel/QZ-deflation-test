program test
   use util
   implicit none

   integer, parameter :: matrix_size = 50
   double precision, dimension(matrix_size,matrix_size) :: A, B
   double precision, dimension(matrix_size) :: alphar, alphai, beta, alphar_d, alphai_d, beta_d
   real, dimension(matrix_size) :: alphar_s, alphai_s, beta_s
   character(len=40) :: arg, dir

   integer :: i, niter, seed_size, hdeflation, tdeflation, iexp
   integer, allocatable :: seed( : )

   ! Get the name of the directory of the executable
   call get_command_argument(0,arg)
   dir = arg( 1: index( arg, "/", .TRUE. ) )


   call random_seed( size = seed_size)
   allocate(seed(seed_size))
   do i=1,seed_size
      seed(i) = 0
   end do
   seed(1) = 1302
   call random_seed( put = seed )

   ! Test well conditioned matrices
   open (unit = 1, file = trim(dir) // "data_well_conditioned.csv")
   do iexp = 1,10000

      call random_number(alphar)
      do i = 1,matrix_size
         beta(i) = 1.0d0
      end do
      beta(matrix_size) = 1.0
      call sort_eigenvalues(alphar, alphai, beta)

      call generate_pencil(A,B, alphar, alphai, beta, 1.0d0)
      ! call generate_graded_pencil(A,B)
      call test_eigenvalues_double(A,B, 1, 1, alphar_d, alphai_d, beta_d, niter)
      call sort_eigenvalues(alphar_d, alphai_d, beta_d)

      do hdeflation = 1,3
         do tdeflation = 1,1

            call test_eigenvalues_real(A,B, hdeflation, tdeflation, alphar_s, alphai_s, beta_s, niter)
            call sort_eigenvalues_real(alphar_s, alphai_s, beta_s)

            write(1,'(I8,",",ES16.6)',advance="no") niter, maxval( abs(alphar_s/beta_s - alphar_d/beta_d)/abs(alphar_d/beta_d) )
            if(hdeflation .lt. 3 .or. tdeflation .lt. 2) write(1,'(A1)',advance="no") ","
         end do
      end do
      write(1,*) '\n'

   end do
   close(1)

   ! Test ill conditioned matrices matrices
   open (unit = 1, file = trim(dir) // "data_ill_conditioned.csv")
   do iexp = 1,10000

      call random_number(alphar)
      do i = 1,matrix_size
         beta(i) = 1.0d0
      end do
      call sort_eigenvalues(alphar, alphai, beta)

      call generate_pencil(A,B, alphar, alphai, beta, 1.0d3)
      ! call generate_graded_pencil(A,B)
      call test_eigenvalues_double(A,B, 1, 1, alphar_d, alphai_d, beta_d, niter)
      call sort_eigenvalues(alphar_d, alphai_d, beta_d)

      do hdeflation = 1,3
         do tdeflation = 1,1

            call test_eigenvalues_real(A,B, hdeflation, tdeflation, alphar_s, alphai_s, beta_s, niter)
            call sort_eigenvalues_real(alphar_s, alphai_s, beta_s)

            write(1,'(I8,",",ES16.6)',advance="no") niter, maxval( abs(alphar_s/beta_s - alphar_d/beta_d)/abs(alphar_d/beta_d) )
            if(hdeflation .lt. 3 .or. tdeflation .lt. 2) write(1,'(A1)',advance="no") ","
         end do
      end do
      write(1,*) '\n'

   end do
   close(1)

   ! Test with beta not equal to one to try to induce a difference between elementwise and strict
   open (unit = 1, file = trim(dir) // "data_beta_ne_one.csv")
   do iexp = 1,10000

      call random_number(alphar)
      call random_number(beta)
      do i = 1,matrix_size
         beta(i) = 1.0d2**beta(i)
      end do
      ! call sort_eigenvalues(alphar, alphai, beta)

      call generate_pencil(A,B, alphar, alphai, beta, 1.0d2)
      ! call generate_graded_pencil(A,B)
      call test_eigenvalues_double(A,B, 1, 1, alphar_d, alphai_d, beta_d, niter)
      call sort_eigenvalues(alphar_d, alphai_d, beta_d)

      do hdeflation = 1,3
         do tdeflation = 1,1

            call test_eigenvalues_real(A,B, hdeflation, tdeflation, alphar_s, alphai_s, beta_s, niter)
            call sort_eigenvalues_real(alphar_s, alphai_s, beta_s)

            write(1,'(I8,",",ES16.6)',advance="no") niter, maxval( abs(alphar_s/beta_s - alphar_d/beta_d)/abs(alphar_d/beta_d) )
            if(hdeflation .lt. 3 .or. tdeflation .lt. 2) write(1,'(A1)',advance="no") ","
         end do
      end do
      write(1,*) '\n'

   end do
   close(1)

   ! Test graded matrices
   open (unit = 1, file = trim(dir) // "data_graded.csv")
   do iexp = 1,10000

      call generate_graded_pencil(A,B,.false.)
      call test_eigenvalues_double(A,B, 1, 1, alphar_d, alphai_d, beta_d, niter)
      call sort_eigenvalues(alphar_d, alphai_d, beta_d)

      do hdeflation = 1,3
         do tdeflation = 1,1

            call test_eigenvalues_real(A,B, hdeflation, tdeflation, alphar_s, alphai_s, beta_s, niter)
            call sort_eigenvalues_real(alphar_s, alphai_s, beta_s)

            write(1,'(I8,",",ES16.6)',advance="no") niter, maxval( abs(alphar_s/beta_s - alphar_d/beta_d)/abs(alphar_d/beta_d) )
            if(hdeflation .lt. 3 .or. tdeflation .lt. 2) write(1,'(A1)',advance="no") ","
         end do
      end do
      write(1,*) '\n'

   end do
   close(1)



   deallocate(seed)

contains



end program
