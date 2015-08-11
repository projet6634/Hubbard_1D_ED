module m_combination
    public :: combinations
    public :: choose
    private


contains

    !> Find all the combinations of k numbers in (1,2,3,...,N).
    !! @param N
    !! @param k
    !! @return results list of combinations
    !! @return N_c the number of combinations
    subroutine combinations(n, k, co, N_c)
      integer, intent(in) :: n, k
      integer, dimension(:,:), allocatable, intent(out) :: co
      integer, intent(out) :: N_c

      integer :: i, j, s, ix, kx, hm, t
      integer :: err

      hm = choose(n, k, err)
      N_c = hm
      if ( err /= 0 ) then
         return
      end if

      allocate(co(hm, k))
      do i = 0, hm-1
         ix = i; kx = k
         do s = 0, n-1
            if ( kx == 0 ) exit
            t = choose(n-(s+1), kx-1)
            if ( ix < t ) then
               co(i+1,kx) = s+1
               kx = kx - 1
            else
               ix = ix - t
            end if
         end do
      end do

    end subroutine combinations

    !> Calculate the number of combinations : nCk
    !! @param N
    !! @param k
    !! @return (N k)
    function choose(n, k, err)
      integer :: choose
      integer, intent(in) :: n, k
      integer, optional, intent(out) :: err

      integer :: imax, i, imin, ie

      ie = 0
      if ( (n < 0 ) .or. (k < 0 ) ) then
         choose = 0
         ie = 1
      else
         if ( n < k ) then
            choose = 0
         else if ( n == k ) then
            choose = 1
         else
            imax = max(k, n-k)
            imin = min(k, n-k)
            choose = 1
            do i = imax+1, n
               choose = choose * i
            end do
            do i = 2, imin
               choose = choose / i
            end do
         end if
      end if
      if ( present(err) ) err = ie
    end function choose
end module m_combination
