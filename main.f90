
! 1D single band hubbard model
program main
    use m_hubbard
    use m_datetime
    implicit none

    integer :: N_sites = 2, i,j,k, &
                basisfile=11, &
                eigfile=12, &
                eigvfile=13, &
                n_gs
                
    real(kind=8) :: U_over_t = 3, &
                    E_gs, &
                    tol=1e-5, &
                    gs_coeff_cut=0.1

    ! Read input parameters
    read(*,*) N_sites, U_over_t

    write(*,*) "[Input parameters]"
    write(*,*) "Number of sites = ", N_sites
    write(*,*) "U/t = ", U_over_t

    write(*,*) "[Initialization]"
    write(*,*) "start init_hubbard"
    call timestamp

    call init_hubbard(N_sites, U_over_t)

    write(*,*) "end init_hubbard"
    call timestamp


    write(*,*) "Hilbert space dimension ", N_d

    open(basisfile, file="basis", status="replace")
    open(eigfile, file="eigenvalues", status="replace")
    open(eigvfile, file="eigenvectors", status="replace")

    do i=1,N_d
        write(basisfile,"(I10,4x)", advance='no') i
        do j=1,2*N
            write(basisfile, "(I2)", advance='no') V(i,j)
        end do
        write(basisfile,*) ''
    end do

    write(*,*) "start diagonalization"
    call timestamp

    call solve

    write(*,*) "end of diagonalization"
    call timestamp

    write(*,*) "write eigenvalues, eigenvectors"
    do i=1,N_d
        write(eigfile,*) E(i)
    end do

    write(*,*) "writing eigenvectors"
    do i=1,N_d
        do j=1,N_d
            write(eigvfile,"(F5.2)",advance='no') H(i,j)
        end do
        write(eigvfile, *) ''
    end do

    write(*,*) "analyze_result"
    call timestamp

    n_gs = 1
    ! Assuming the eigenvalues are in ascending order
    E_gs = E(1)
    do i=2,N_d
        if (abs(E(i)-E_gs)<tol) then
            n_gs = n_gs + 1
        else
            exit
        end if
    end do
    write(*,*) NEW_LINE('a')
    write(*,*) "[Ground state]"
    write(*,*) "Number of ground states : ", n_gs

    do i=1,n_gs
        write(*,*) "### Ground state ", i
        do j=1,N_d
            if (H(j,i)>gs_coeff_cut) then
                write(*,"(6x,A6,F5.3,A6,I3)") "coeff:", H(j,i), ", idx:", j
                write(*,"(6x,A6)",advance='no') "n_up " 
                do k=1,N
                    write(*,"(I2)", advance='no') V(j,k)
                end do
                write(*,*) ''
                write(*,"(6x,A6)",advance='no') "n_down " 
                do k=N+1,2*N
                    write(*,"(I2)", advance='no') V(j,k)
                end do
                write(*,*) NEW_LINE('a')
            end if
        end do
    end do

    close(basisfile)
    close(eigfile)
    close(eigvfile)

    write(*,*) "End of run"
    call timestamp
end program main

