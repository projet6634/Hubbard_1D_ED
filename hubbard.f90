module m_hubbard
    use m_combination
    use m_datetime
    implicit none
    integer, parameter:: dp=kind(0.d0)

    public :: init_hubbard, solve
    
    integer, public :: N ! Number of sites
    integer, public :: N_d ! Dimension of the Hilbert space (number of basis)
    real(dp), public :: U ! Hubbard U (in unit of t)
    real(dp) :: t=-1.0_dp ! sign of the hopping parameter

    integer, dimension(:,:), allocatable, public :: V ! Set of occupation number basis vectors
    real(dp), dimension(:,:), allocatable, public :: H ! H matrix elements (in unit of t)
    real(dp), dimension(:), allocatable, public :: E ! Eigenvalues
                                                   
    real(dp), dimension(:,:,:), allocatable, public :: A ! Hopping op 

    private

contains 
    subroutine init_hubbard(Nin, Uin)
        integer, intent(in) :: Nin
        real(dp), intent(in) :: Uin

        N = Nin
        U = Uin

        write(*,*) "start init_hubbard"
        call timestamp

        call calc_dimension_hilbert

        call construct_basis

        call construct_hopping_op

        call construct_hamiltonian

        write(*,*) "end init_hubbard"
        call timestamp
    end subroutine init_hubbard

    !> Dimension of the Hilbert space
    subroutine calc_dimension_hilbert
        integer :: i,c

        N_d = 0
        do i=0,N
            c = choose(N, i) 
            N_d =N_d + c*c
        end do
    end subroutine calc_dimension_hilbert

    !> Find occupation number basis
    subroutine construct_basis
        integer, dimension(:,:), allocatable :: occ_ups, occ_downs

        integer i,j,k,a,b,n_up,N_c

        allocate(V(N_d, 2*N))
        V = 0

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

    subroutine construct_hopping_op
        integer i,i2,j,k,l,sum_n
        logical orth

        allocate(A(2*N, N_d, N_d))


        do i=1,2*N
            do j=1,N_d
                do k=1,N_d
                   
                    ! Periodic boundary condition
                    if (mod(i,N)==0) then
                        i2=i-N+1
                    else
                        i2=i+1
                    end if
                    

                    if (V(k,i2)*(1-V(k,i))==0) then
                        A(i,j,k)=0
                    else
                        ! Orthogonality test
                        orth = .false.
                        do l=1,2*N
                            if (l==i) then
                                orth = V(j,l) .ne. (V(k,l)+1)
                            else if (l==i2) then
                                orth = V(j,l) .ne. (V(k,l)-1)
                            else
                                orth = V(j,l) .ne. V(k,l)
                            end if

                            if (orth) then
                                exit
                            end if
                        end do

                        if (orth) then
                            A(i,j,k)=0
                        else 
                            ! Calculate the phase factor
                            sum_n = 0
                            if (i2>1) then
                                do l=1,i2-1
                                    sum_n = sum_n + V(k,l)
                                end do
                            end if

                            if (i>1) then
                                do l=1,i-1
                                    sum_n = sum_n + V(k,l)
                                end do
                                ! Correction due to PBC
                                if (i2<i) then
                                    sum_n = sum_n - 1
                                end if
                            end if

                            if (mod(sum_n,2)==0) then
                                A(i,j,k) = 1
                            else
                                A(i,j,k) = -1
                            end if

                        end if

                    end if
                    
                end do
            end do
        end do
    end subroutine construct_hopping_op

    subroutine construct_hamiltonian

        integer :: i,j,k

        allocate(H(N_d, N_d))

        H = 0.0_dp

        ! Hopping
        do i=1,2*N
            do j=1,N_d
                do k=1,N_d
                    if (N==2) then
                        H(j,k) = H(j,k) + t * A(i,j,k) 
                    else
                        H(j,k) = H(j,k) + t * (A(i,j,k) + A(i,k,j))
                    end if
                end do
            end do
        end do

        ! Coulomb 
        do j=1,N_d
            do i=1,N
                H(j,j) = H(j,j) + U*V(j,i)*V(j,i+N)
            end do
        end do

    end subroutine construct_hamiltonian

    subroutine solve
        real(dp) :: DUMMY(1,1), WORK(3*N_d-1)
        integer :: info
        allocate(E(N_d))

        ! Diagonalize
        call DSYEV('V', 'U', N_d, H, N_d, E, WORK, 3*N_d-1, info)
    end subroutine
end module m_hubbard
