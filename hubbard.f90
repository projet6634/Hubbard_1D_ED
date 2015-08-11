module m_hubbard
    use m_combination
    public :: construct_basis
    private

contains 

    !> Find occupation number basis
    subroutine construct_basis(N, V, N_d)
        integer(kind=4), intent(in) :: N
        integer(kind=4), intent(out) :: N_d
        integer, intent(out), dimension(:,:), allocatable :: V 

        integer, dimension(:,:), allocatable :: occ_ups, occ_downs

        integer i,j,k,a,b,n_up,N_c

        N_d = dimension_hilbert(N)
        allocate(V(N_d, 2*N))
        V = 0.0

        i=1
        do n_up=0,N
            call combinations(N, n_up, occ_ups, N_c)
            call combinations(N, N-n_up, occ_downs, N_c)

            ! For each configurations of ups and downs,
            do j=1,N_c
                do k=1,N_c
                    ! Allocate the occupation to V(i,:)

                    do a=1,n_up
                        V(i, occ_ups(j,a)) = 1
                    end do
                    do a=1,N-n_up
                        V(i, N+occ_downs(k,a)) = 1
                    end do

                    i = i+1
                end do
            end do
        end do

    end subroutine construct_basis

    ! Dimension of the Hilbert space
    ! sum_i=1...N  (N i)^2
    integer(kind=4) function dimension_hilbert(N) result(N_d)
        integer(kind=4), intent(in) :: N
        integer :: i, c
        N_d = 0
        do i=0,N
            c = choose(N, i) 
            N_d = N_d + c*c
        end do
    end function dimension_hilbert
end module m_hubbard
