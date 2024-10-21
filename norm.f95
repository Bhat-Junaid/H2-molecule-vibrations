! This program calculates the rotational constants for the v = 0 and v = 1 states
PROGRAM rotational_cons
    IMPLICIT NONE
    INTEGER, PARAMETER :: N = 2001
    REAL(KIND=8) :: B, x_0(N), y_0(N), x_1(N), y_1(N), r_0(N), r_1(N)
    REAL(KIND=8) :: coeff_0, coeff_1, f_0(N), f_1(N), i_0, i_1, B_0, B_1
    INTEGER :: i

    ! Load data from files for v = 0 and v = 1
    OPEN(UNIT=1, FILE="Wavefunction_0_int.txt", STATUS="OLD")
    DO i = 1, N
        READ(1, *) x_0(i), y_0(i)    
        r_0(i) = 1.0D0 / (x_0(i)**2)
    ENDDO
    CLOSE(1)

    OPEN(UNIT=2, FILE="Wavefunction_1_int.txt", STATUS="OLD")
    DO i = 1, N
        READ(2, *) x_1(i), y_1(i)
        r_1(i) = 1.0D0 / (x_1(i)**2)
    ENDDO
    CLOSE(2)

    ! Calculation of the Normalization coefficients
    coeff_0 = 0.0D0
    DO i = 1, N
        coeff_0 = coeff_0 + y_0(i)**2
    ENDDO

    coeff_1 = 0.0D0
    DO i = 1, N
        coeff_1 = coeff_1 + y_1(i)**2
    ENDDO

    ! Calculation of the averages in order to get B_0 and B_1
    i_0 = 0.0D0
    DO i = 1, N
        f_0(i) = y_0(i)**2 * r_0(i) / coeff_0
        i_0 = i_0 + f_0(i)
    ENDDO

    i_1 = 0.0D0
    DO i = 1, N
        f_1(i) = y_1(i)**2 * r_1(i) / coeff_1
        i_1 = i_1 + f_1(i)
    ENDDO
    ! calculation of actual values of B_0 and B_1
    B = 60.8D0               
    B_0 = B * i_0
    B_1 = B * i_1

    WRITE(*, '(A,F12.6,A,F12.6)') "B0 is ", B_0, " and B1 is ", B_1

END PROGRAM rotational_cons
! B0 is    59.814183 and B1 is    57.705493