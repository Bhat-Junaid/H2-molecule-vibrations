PROGRAM SCHROD
    IMPLICIT NONE
    INTEGER (KIND=4), PARAMETER :: N=601
    INTEGER (KIND=4):: NCMAX, IFLAG, M, JLEV, NCOUNT, NODE, I
    REAL (KIND=8), DIMENSION(0:N):: V, F, RR
    REAL (KIND=8):: RN, A, B, EPS, TOL, FINIT, ENERG, H, HC, E, DE, R0
    REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: R0_VALUES, ENERG_VALUES, NCOUNT_VALUES
    
    DATA R0/0.1D0/, RN/6.6D0/, A/1.3924D0/, B/1.4840D-03/, EPS/1.0D-10/, TOL/1.0D-4/, &
         FINIT/0.D0/, ENERG/-0.1205D0/, NCMAX/30/, IFLAG/1/, M/11/, JLEV/0/
    
    ALLOCATE(R0_VALUES(100), ENERG_VALUES(100), NCOUNT_VALUES(100)) ! Allocate arrays for values
    
    DO I = 1, 100
        R0 = R0 + 0.01   ! Vary R0 values from 0.1 to 1.0
        H = (RN - R0) / DFLOAT(N)  ! calculation of the step
        HC = H * H / B  ! expression h2/b
        
        ! calculation of the vector V(0:N) containing the potential
        CALL POT(N, R0, H, A, B, JLEV, RR, V)
        
        ! putting the trial energy value
        E = ENERG
        NCOUNT = 0
        
        DO WHILE (NCOUNT <= NCMAX)
            ! propagation; NCOUNT giving the number of iterations
            CALL SHOOT(N, A, H, HC, FINIT, EPS, E, V, F, NODE, DE, M, IFLAG)
            NCOUNT = NCOUNT + 1
            IF (ABS(DE/E).LE.TOL) EXIT ! Exit loop if convergence criterion met
            E = E + DE
        END DO
        
        WRITE(*, *) "For R0 =", R0, " - ENERG =", E, " - NCOUNT =", NCOUNT
        
        R0_VALUES(I) = R0
        ENERG_VALUES(I) = E
        NCOUNT_VALUES(I) = NCOUNT
    END DO
    
    ! Save values of R0, ENERG, and NCOUNT to a file
    OPEN(UNIT=1, FILE="R0_Energ_Ncount_12.txt", STATUS="UNKNOWN")
    DO I = 1, 100
        WRITE(1, 100) R0_VALUES(I), ENERG_VALUES(I), NCOUNT_VALUES(I)
    END DO
    CLOSE(1)
    STOP
    
    100 FORMAT(3E15.7)
    
    END

    ! Rest of the subroutine definitions remain unchanged
    SUBROUTINE POT(N,R0,H,A,B,JLEV,RR,V)
        IMPLICIT NONE
        INTEGER (KIND=4):: N,JLEV,I
        REAL (KIND=8), DIMENSION(0:N):: V,RR
        REAL (KIND=8):: R0,H,A,B
        DO I=0,N
        RR(I)=R0+I*H
        V(I)=DEXP(-2.D0*A*(RR(I)-1.D0))-2.D0*DEXP(-A*(RR(I)-1.D0))
        ENDDO
        RETURN
        END
        SUBROUTINE SHOOT(N,A,H,HC,FINIT,EPS,E,V,F,NODE,DE,M,IFLAG)
        IMPLICIT NONE
        INTEGER (KIND=4):: N,NODE,M,IFLAG,ISTEP,IINIT,IFIN,I
        REAL (KIND=8), DIMENSION(0:N):: V,F
        REAL (KIND=8):: A,B,H,HC,FINIT,EPS,E,DE,F0,F1,COEFF,AN,FSM,FM
        ISTEP=1 ! propagation do the right
        IINIT=1
        IFIN=N-1
        IF(IFLAG.EQ.0) IFIN=M
        F0=FINIT
        F1=EPS
        CALL PROPAG(N,HC,ISTEP,IINIT,IFIN,IFLAG,F0,F1,E,V,F,M)
        FM=F(M)
        ISTEP=-1 ! propagation to the left
        IINIT=N-1
        IFIN=M+1
        F0=FINIT
        F1=EPS
        CALL PROPAG(N,HC,ISTEP,IINIT,IFIN,IFLAG,F0,F1,E,V,F,M)
        COEFF=F(M)/FM
        DO I=0,M-1
        F(I)=F(I)*COEFF
        ENDDO
        NODE=0 ! calculation of nodes:
        AN=0.D0 ! calculation of the norm for the wavefunction
        DO I=1,N ! F0=0
        AN=AN+F(I)*F(I)
        IF(F(I)*F(I-1).LT.0.D0) NODE=NODE+1
        ENDDO
        FSM=F(M+1)+F(M-1)-2.D0*F(M) ! correction of the energy:
        DE=F(M)*((V(M)-E)*F(M)-FSM/HC)/AN
        RETURN
        END
        SUBROUTINE PROPAG(N,HC,ISTEP,IINIT,IFIN,IFLAG,F0,F1,E,V,F,M)
        IMPLICIT NONE
        INTEGER (KIND=4):: N,ISTEP,IINIT,IFIN,IFLAG,M,I,I1
        REAL (KIND=8), DIMENSION(0:N):: V,F
        REAL (KIND=8):: HC,F0,F1,E
        F(IINIT-ISTEP)=F0
        F(IINIT)=F1
        DO I=IINIT,IFIN,ISTEP
        I1=I+ISTEP
        F(I1)=(2.0D0+HC*(V(I)-E))*F(I)-F(I-ISTEP)
        IF(IFLAG.EQ.ISTEP) THEN
        IF(F(I1).LT.F(I)) THEN
        M=I1
        RETURN
        ENDIF
        ENDIF
        ENDDO
        RETURN
        END