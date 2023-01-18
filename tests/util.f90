module util

contains

subroutine test_eigenvalues_double(H,T, hdeflation, tdeflation, ar, ai, be, niter)
   double precision, dimension(:,:), intent(in) :: H, T
   double precision, dimension(:), intent(out) :: ar, ai, be
   integer, intent(in) :: hdeflation, tdeflation
   integer, intent(out) :: niter

   double precision, allocatable :: work(:), T_copy(:,:), H_copy(:,:)
   double precision :: dum(1)
   integer :: lwork, n, info

   n = size(T, 1)
   lwork = n
   allocate( work(lwork) )
   allocate( H_copy(n,n) )
   allocate( T_copy(n,n) )
   H_copy = H
   T_copy = T
   deflation = 3
   call DHGEQZV1( 'E', 'N', 'N', n, 1, n, H_copy, n, T_copy, n, ar, ai, be, dum, 1, dum, 1, work, lwork,&
                  info, hdeflation, tdeflation, niter )
   deallocate(work,H_copy,T_copy)

end subroutine

subroutine test_eigenvalues_real(H,T, hdeflation, tdeflation, ar, ai, be, niter)
   double precision, dimension(:,:), intent(in) :: H, T
   real, dimension(:), intent(out) :: ar, ai, be
   integer, intent(in) :: hdeflation, tdeflation
   integer, intent(out) :: niter
   
   real, allocatable :: work(:), T_copy(:,:), H_copy(:,:)
   real :: dum(1)
   integer :: lwork, n, info

   n = size(T, 1)
   lwork = n
   allocate( work(lwork) )
   allocate( H_copy(n,n) )
   allocate( T_copy(n,n) )
   H_copy = real(H)
   T_copy = real(T)
   call SHGEQZV1( 'E', 'N', 'N', n, 1, n, H_copy, n, T_copy, n, ar, ai, be, dum, 1, dum, 1, work, lwork,&
                  info, hdeflation, tdeflation, niter )
   deallocate(work,H_copy,T_copy)

end subroutine

subroutine generate_pencil(H,T, ar, ai, be, condition)
   double precision, dimension(:,:), intent(out) :: H, T
   double precision, dimension(:), intent(in) :: ar, ai, be
   double precision, intent(in) :: condition

   double precision, allocatable :: Q(:,:)
   double precision, allocatable :: work(:)
   double precision, allocatable :: S(:)
   double precision :: dum(1)
   integer :: lwork, n, info, jcol, jrow

   n = size(T, 1)

   ! Workspace query
   call dgeqrf( n, n, T, n, dum, dum, -1,info)
   lwork = int(dum(1))
   call dormqr( 'L', 'T',n,n,n,T,n,dum,H,n,dum,-1,info)
   lwork = max(lwork, int(dum(1))) + 1
   allocate( work(lwork) )

   ! Initialize
   call dlaset( 'all', n,n, 0.0d0, 1.0d0, H, n )
   call dlaset( 'all', n,n, 0.0d0, 1.0d0, T, n )

   ! Set diagonals of H and T to desired eigenvalues (TODO, complex)
   do jcol = 1,n
      H(jcol,jcol) = ar(jcol)
      T(jcol,jcol) = be(jcol)
   end do

   ! Multiply H,T with random matrices from left and right
   allocate(Q(n,n))
   call random_number(Q)
   call dgeqrf( n, n, Q, n, work, work(n+1), lwork-n,info)
   call dormqr( 'L', 'N',n,n,n,Q,n,work,H,n,work(n+1),lwork-n,info)
   call dormqr( 'L', 'N',n,n,n,Q,n,work,T,n,work(n+1),lwork-n,info)

   call random_number(Q)
   call dgeqrf( n, n, Q, n, work, work(n+1), lwork-n,info)
   call dormqr( 'R', 'N',n,n,n,Q,n,work,H,n,work(n+1),lwork-n,info)
   call dormqr( 'R', 'N',n,n,n,Q,n,work,T,n,work(n+1),lwork-n,info)
   deallocate(Q)

   ! Multiply H,T with diagonal scaling to set eigenvector condition
   allocate(S(n))
   do jcol = 1,n
      S(jcol) = condition**((jcol-1)/real(n-1))
   end do

   do jcol = 1,n
      call dscal( n, S(jcol), H(1,jcol), 1 )
      call dscal( n, S(jcol), T(1,jcol), 1 )

      call dscal( n, S(jcol), H(jcol,1), n )
      call dscal( n, S(jcol), T(jcol,1), n )
   end do
   deallocate(S)

   ! Multiply H,T with random matrices from left and right
   allocate(Q(n,n))
   call random_number(Q)
   call dgeqrf( n, n, Q, n, work, work(n+1), lwork-n,info)
   call dormqr( 'L', 'N',n,n,n,Q,n,work,H,n,work(n+1),lwork-n,info)
   call dormqr( 'L', 'N',n,n,n,Q,n,work,T,n,work(n+1),lwork-n,info)

   call random_number(Q)
   call dgeqrf( n, n, Q, n, work, work(n+1), lwork-n,info)
   call dormqr( 'R', 'N',n,n,n,Q,n,work,H,n,work(n+1),lwork-n,info)
   call dormqr( 'R', 'N',n,n,n,Q,n,work,T,n,work(n+1),lwork-n,info)
   deallocate(Q)

   ! Hessenberg-triangular reduction
   call dgeqrf( n, n, T, n, work, work(n+1), lwork-n,info)
   call dormqr( 'L', 'T',n,n,n,T,n,work,H,n,work(n+1),lwork-n,info)
   call dgghrd( 'N','N', N, 1, N, H, n, T, n, dum, 1, dum, 1, info )

   deallocate(work)

end subroutine

subroutine generate_graded_pencil(H,T, includeT)
   double precision, dimension(:,:), intent(out) :: H, T
   logical :: includeT

   double precision, allocatable :: work(:)
   double precision, allocatable :: S(:)
   double precision :: dum(1)
   integer :: lwork, n, info, jcol

   n = size(T, 1)

   ! Workspace query
   call dgeqrf( n, n, T, n, dum, dum, -1,info)
   lwork = int(dum(1))
   call dormqr( 'L', 'T',n,n,n,T,n,dum,H,n,dum,-1,info)
   lwork = max(lwork, int(dum(1))) + 1
   allocate( work(lwork) )

   ! Initialize
   call random_number(H)
   call random_number(T)

   ! Scale rows to get a graded matrix
   allocate(S(n))
   S(1) = 1.0d0
   do jcol = 2,n
      S(jcol) = 1.0d1**(-3.0d0*(jcol)/n)
   end do

   do jcol = 1,n
      call dscal( n, S(jcol), H(1,jcol), 1 )
      call dscal( n, 1/S(jcol), T(1,jcol), 1 )

      call dscal( n, S(jcol), H(jcol,1), n )
      call dscal( n, 1/S(jcol), T(jcol,1), n )
   end do
   deallocate(S)


   if(.not. includeT) call dlaset( 'all', n,n, 0.0d0, 1.0d0, T, n )

   ! Hessenberg-triangular reduction
   call dgeqrf( n, n, T, n, work, work(n+1), lwork-n,info)
   call dormqr( 'L', 'T',n,n,n,T,n,work,H,n,work(n+1),lwork-n,info)
   call dgghrd( 'N','N', N, 1, N, H, n, T, n, dum, 1, dum, 1, info )

   deallocate(work)

end subroutine

subroutine generate_graded_inf_pencil(H,T)
   double precision, dimension(:,:), intent(out) :: H, T

   double precision, allocatable :: work(:)
   double precision, allocatable :: S(:)
   double precision :: dum(1)
   integer :: lwork, n, info, jcol

   n = size(T, 1)

   ! Workspace query
   call dgeqrf( n, n, T, n, dum, dum, -1,info)
   lwork = int(dum(1))
   call dormqr( 'L', 'T',n,n,n,T,n,dum,H,n,dum,-1,info)
   lwork = max(lwork, int(dum(1))) + 1
   allocate( work(lwork) )

   ! Initialize
   call random_number(H)
   call random_number(T)

   ! Scale rows to get a graded matrix
   allocate(S(n))
   S(1) = 1.0d0
   do jcol = 2,n
      S(jcol) = 1.0d1**(-8.0d0*(jcol)/n)
   end do

   do jcol = 1,n
      ! call dscal( n, S(jcol), H(1,jcol), 1 )
      call dscal( n, S(jcol), T(1,jcol), 1 )

      ! call dscal( n, S(jcol), H(jcol,1), n )
      call dscal( n, S(jcol), T(jcol,1), n )
   end do
   deallocate(S)

   ! Hessenberg-triangular reduction
   call dgeqrf( n, n, T, n, work, work(n+1), lwork-n,info)
   call dormqr( 'L', 'T',n,n,n,T,n,work,H,n,work(n+1),lwork-n,info)
   call dgghrd( 'N','N', N, 1, N, H, n, T, n, dum, 1, dum, 1, info )

   deallocate(work)

end subroutine

subroutine generate_inf_pencil(H,T, n_inf)
   double precision, dimension(:,:), intent(out) :: H, T
   integer, intent(inout) :: n_inf

   double precision, allocatable :: work(:)
   double precision, allocatable :: S(:)
   double precision :: dum(1)
   integer :: lwork, n, info, m1, m2

   n = size(T, 1)

   ! Workspace query
   call dgeqrf( n, n, T, n, dum, dum, -1,info)
   lwork = int(dum(1))
   call dormqr( 'L', 'T',n,n,n,T,n,dum,H,n,dum,-1,info)
   lwork = max(lwork, int(dum(1))) + 1
   allocate( work(lwork) )

   ! Initialize
   m1 = (n-n_inf+1)/2
   m2 = n-m1
   call normal_random_number(H)
   T = 0.0d0
   call normal_random_number(T(1:m1,1:m2))
   call normal_random_number(T(m1+1:n,m2+1:n))

   ! Hessenberg-triangular reduction
   call dgeqrf( n, n, T, n, work, work(n+1), lwork-n,info)
   call dormqr( 'L', 'T',n,n,n,T,n,work,H,n,work(n+1),lwork-n,info)
   call dgghrd( 'N','N', N, 1, N, H, n, T, n, dum, 1, dum, 1, info )

   deallocate(work)

end subroutine

subroutine sort_eigenvalues( ar, ai, be )
   double precision, dimension(:), intent(inout) :: ar, ai, be

   double precision :: swap
   integer :: k, j, n

   n = size( ar )

   ! Bubble sort
   do k = n,1,-1
      do j = 1,k-1
         if( ar(j)/be(j) .GT. ar(j+1)/be(j+1) ) then
            swap = ar(j)
            ar(j) = ar(j+1)
            ar(j+1) = swap

            swap = ai(j)
            ai(j) = ai(j+1)
            ai(j+1) = swap

            swap = be(j)
            be(j) = be(j+1)
            be(j+1) = swap

         end if
      end do
   end do

end subroutine

subroutine sort_eigenvalues_real( ar, ai, be )
   real, dimension(:), intent(inout) :: ar, ai, be

   real :: swap
   integer :: k, j, n

   n = size( ar )

   ! Bubble sort
   do k = n,1,-1
      do j = 1,k-1
         if( ar(j)/be(j) .GT. ar(j+1)/be(j+1) ) then
            swap = ar(j)
            ar(j) = ar(j+1)
            ar(j+1) = swap

            swap = ai(j)
            ai(j) = ai(j+1)
            ai(j+1) = swap

            swap = be(j)
            be(j) = be(j+1)
            be(j+1) = swap

         end if
      end do
   end do

end subroutine

subroutine normal_random_number( X )
   implicit none
   double precision, intent(out) :: X(:,:)
   integer :: n, m, i, j
   double precision, parameter :: pi = 4.d0*atan(1.d0)
   double precision :: u(2)

   m = size(X, 1)
   n = size(X, 2)

   do i = 1,m
      do j = 1,n
         call random_number(u)
         X(i,j) = sqrt(-2.d0*log(u(1)))*cos(2.d0*pi*u(2))
      end do
   end do

end subroutine

end module