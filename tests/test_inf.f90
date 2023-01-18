program test
   use util
   implicit none

   integer, parameter :: matrix_size = 50
   double precision, dimension(matrix_size,matrix_size) :: A, B
   double precision, dimension(matrix_size) :: alphar, alphai, beta, alphar_d, alphai_d, beta_d
   real, dimension(matrix_size) :: alphar_s, alphai_s, beta_s
   character(len=40) :: arg, dir

   integer :: i, niter, seed_size, hdeflation, tdeflation, iexp, n_inf, n_inf_s, n_inf_s2, nexp
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

   write(*,*) "testing matrix with B constructed to have 6 infinite eigenvalues"

   ! Test pencils with set number of infinite eigenvalues
   do tdeflation = 1,3
      call random_seed( put = seed )
      nexp = 1000
      n_inf_s = 0
      n_inf_s2 = 0
      do iexp = 1,nexp

         n_inf = 6
         call generate_inf_pencil(A,B, n_inf)
         call test_eigenvalues_double(A,B, 1, tdeflation, alphar_d, alphai_d, beta_d, niter)
         call sort_eigenvalues(alphar_d, alphai_d, beta_d)

         do i = 1, matrix_size
            if( beta_d(i) == 0.0d0) then
               n_inf_s = n_inf_s + 1
            end if
         end do

      end do
      write(*,*) "Deflation type: ", tdeflation, "number of infinite eigenvalues detected: ",real(n_inf_s)/nexp
   end do

   write(*,*)
   write(*,*) "testing matrix with very large alpha (no infinite eigenvalues)"

   ! Test pencils with large, but finite eigenvalues
   do tdeflation = 1,3
      call random_seed( put = seed )
      nexp = 1000
      n_inf_s = 0
      n_inf_s2 = 0
      do iexp = 1,nexp

         do i = 1,matrix_size
            alphar(i) = 1.0d1**((20.0d0*i)/matrix_size)
            beta(i) = 1.0d0
         end do
         call generate_pencil(A,B, alphar, alphai, beta, 10.0d0)
         call sort_eigenvalues(alphar, alphai, beta)

         call test_eigenvalues_double(A,B, 1, tdeflation, alphar_d, alphai_d, beta_d, niter)
         call sort_eigenvalues(alphar_d, alphai_d, beta_d)

         do i = 1, matrix_size
            if( beta_d(i) == 0.0d0) then
               n_inf_s = n_inf_s + 1
            end if
         end do

      end do
      write(*,*) "Deflation type: ", tdeflation, "number of infinite eigenvalues detected: ",real(n_inf_s)/nexp
   end do

   write(*,*)
   write(*,*) "testing matrix with very small beta (2 eigenvalues numerically infinite)"

   ! Test pencils with large, but finite eigenvalues
   do tdeflation = 1,3
      call random_seed( put = seed )
      nexp = 1000
      n_inf_s = 0
      n_inf_s2 = 0
      do iexp = 1,nexp

         do i = 1,matrix_size
            alphar(i) = 1.0d0
            beta(i) = 1.0d1**((-16.0d0*(i-1))/(matrix_size-1))
         end do
         call generate_pencil(A,B, alphar, alphai, beta, 1.0d0)
         call sort_eigenvalues(alphar, alphai, beta)

         call test_eigenvalues_double(A,B, 1, tdeflation, alphar_d, alphai_d, beta_d, niter)
         call sort_eigenvalues(alphar_d, alphai_d, beta_d)

         do i = 1, matrix_size
            if( beta_d(i) == 0.0d0) then
               n_inf_s = n_inf_s + 1
            end if
         end do

      end do
      write(*,*) "Deflation type: ", tdeflation, "number of infinite eigenvalues detected: ",real(n_inf_s)/nexp
   end do

   write(*,*)
   write(*,*) "testing matrix with very small alpha and beta (2 eigenvalues numerically infinite)"

   ! Test pencils with large, but finite eigenvalues
   do tdeflation = 1,3
      call random_seed( put = seed )
      nexp = 1000
      n_inf_s = 0
      n_inf_s2 = 0
      do iexp = 1,nexp

         do i = 1,matrix_size
            alphar(i) = 1.0d1**((-8.0d0*i)/matrix_size)
            beta(i) = 1.0d1**((-16.0d0*i)/matrix_size)
         end do
         call generate_pencil(A,B, alphar, alphai, beta, 10.0d0)
         call sort_eigenvalues(alphar, alphai, beta)

         call test_eigenvalues_double(A,B, 1, tdeflation, alphar_d, alphai_d, beta_d, niter)
         call sort_eigenvalues(alphar_d, alphai_d, beta_d)

         do i = 1, matrix_size
            if( beta_d(i) == 0.0d0) then
               n_inf_s = n_inf_s + 1
            end if
         end do

      end do
      write(*,*) "Deflation type: ", tdeflation, "number of infinite eigenvalues detected: ",real(n_inf_s)/nexp
   end do
   write(*,*)

   deallocate(seed)

end program
