
! 1D single band hubbard model
program main
    use m_combination
    use m_hubbard
    implicit none

    integer(kind=4) :: N, &   ! Number of sites
                       N_d    ! Hilbert space dimension
    real(kind=8), dimension(:,:), allocatable :: H ! H matrix elements (in unit of t)
    real :: U = 1.0 ! U/t

    integer, dimension(:,:), allocatable :: V

    integer :: i

    call construct_basis(3, V, N_d)

    print *, "N_d= ", N_d
    do i=1,N_d
        print *,V(i,:)
    end do
    ! Construct the Hubbard Hamiltonian
    
    ! Diagonalize

end program main



