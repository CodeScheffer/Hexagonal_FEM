MODULE Corpotri_module
    USE precision_module, ONLY : DP
    USE Material_Module, ONLY : MontaCEPT
    USE Elementotri_Module, ONLY : MontaXYtri, MontaJtri, MontadNtri, CorrigedNtri, MontaBtri, Montaglstri

    IMPLICIT NONE

    REAL(KIND=DP), PARAMETER :: PI = ACOS(-1.0_DP)

CONTAINS

! Monta as funções de interpolação do elemento Q4

FUNCTION MontaNtri(r, s) RESULT(N)

    REAL(KIND=DP), INTENT(IN) :: r, s

    REAL(KIND=DP) :: N(2,8)
    REAL(KIND=DP) :: N1, N2, N3, N4

    N1 = 0.25_DP * (1.0_DP - r) * (1.0_DP - s)
    N2 = 0.25_DP * (1.0_DP + r) * (1.0_DP - s)
    N3 = 0.25_DP * (1.0_DP + r) * (1.0_DP + s)
    N4 = 0.25_DP * (1.0_DP - r) * (1.0_DP + s)

    N = 0.0_DP
    N(1,1) = N1; N(1,3) = N2; N(1,5) = N3; N(1,7) = N4
    N(2,2) = N1; N(2,4) = N2; N(2,6) = N3; N(2,8) = N4
END FUNCTION MontaNtri

! Calcula o vetor de forças nodais equivalentes Fb devido a forças de corpo do elemento Q4, usando integração de Gauss

FUNCTION MontaFBetri(ele, IJ, X, Y, te, b) RESULT(Fb)
    USE precision_module, ONLY : DP

    INTEGER, INTENT(IN) :: ele
    INTEGER, INTENT(IN) :: IJ(:,:)
    REAL(KIND=DP), INTENT(IN) :: X(4), Y(4)
    REAL(KIND=DP), INTENT(IN) :: te
    REAL(KIND=DP), INTENT(IN) :: b(2)

    REAL(KIND=DP) :: Fb(8)

    INTEGER, PARAMETER :: NUM_GAUSS_PTS_Q4 = 4
    REAL(KIND=DP), PARAMETER :: GPS_COORD_Q4 = 1.0_DP / SQRT(3.0_DP)

    REAL(KIND=DP), DIMENSION(NUM_GAUSS_PTS_Q4) :: r_pts, s_pts
    REAL(KIND=DP), DIMENSION(NUM_GAUSS_PTS_Q4) :: w_pts_q4

    REAL(KIND=DP) :: N_matrix(2,8)
    REAL(KIND=DP) :: dNrs(2,4), J(2,2), dJ

    INTEGER :: i_pt
    INTEGER :: element_type

    ! element_type = IJ(ele, 1)

    ! IF (element_type == 2) THEN
        Fb = 0.0_DP

        r_pts(1) = -GPS_COORD_Q4
        s_pts(1) = -GPS_COORD_Q4
        w_pts_q4(1) = 1.0_DP

        r_pts(2) = GPS_COORD_Q4
        s_pts(2) = -GPS_COORD_Q4
        w_pts_q4(2) = 1.0_DP

        r_pts(3) = GPS_COORD_Q4
        s_pts(3) = GPS_COORD_Q4
        w_pts_q4(3) = 1.0_DP

        r_pts(4) = -GPS_COORD_Q4
        s_pts(4) = GPS_COORD_Q4
        w_pts_q4(4) = 1.0_DP

        DO i_pt = 1, NUM_GAUSS_PTS_Q4

            dNrs = MontadNtri(r_pts(i_pt), s_pts(i_pt))
            J = MontaJtri(dNrs, X, Y)
            dJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)

            IF (ABS(dJ) < 1.0E-12_DP) THEN
                WRITE(*,*) 'ERRO em MontaFBetri: Elemento #', ele, ': Determinante Jacobiano muito pequeno no ponto de Gauss (', r_pts(i_pt), ',', s_pts(i_pt), ')'
                STOP 'Jacobiano singular em MontaFBetri'
            END IF

            N_matrix = MontaNtri(r_pts(i_pt), s_pts(i_pt))

            Fb = Fb + MATMUL(TRANSPOSE(N_matrix), b) * ABS(dJ) * te * w_pts_q4(i_pt)

        END DO
    ! END IF

END FUNCTION MontaFBetri

! Calcula o vetor global de forças de corpo F geradas no elemento Q4, com base nos dados de entrada

FUNCTION ForcaGlobaltri(nn, ESP, XY, IJ, nfb, FB_input) RESULT(F)
    USE precision_module, ONLY : DP
    INTEGER, INTENT(IN) :: nn
    REAL(KIND=DP), INTENT(IN) :: ESP(:)
    REAL(KIND=DP), INTENT(IN) :: XY(:,:)
    INTEGER, INTENT(IN) :: IJ(:,:)
    INTEGER, INTENT(IN) :: nfb
    REAL(KIND=DP), INTENT(IN) :: FB_input(:,:)

    REAL(KIND=DP) :: F(2*nn)
    REAL(KIND=DP) :: b(2)
    INTEGER :: i, j
    INTEGER :: ele, dir
    REAL(KIND=DP) :: val
    REAL(KIND=DP) :: X(4), Y(4)
    REAL(KIND=DP) :: te
    REAL(KIND=DP) :: Fb(8)
    INTEGER :: ge(8)
    INTEGER :: element_type

    F = 0.0_DP

    DO i = 1, nfb
        ele = INT(FB_input(i,1))
        dir = INT(FB_input(i,2))
        val = FB_input(i,3)

        element_type = IJ(ele, 1)

        IF (element_type == 2) THEN
            b = 0.0_DP
            b(dir) = val

            CALL MontaXYtri(ele, IJ, XY, X, Y)
            te = ESP(ele)

            Fb = MontaFBetri(ele, IJ, X,Y,te,b)

            ge = Montaglstri(ele,IJ)

            DO j = 1, 8
                F(ge(j)) = F(ge(j)) + Fb(j)
            END DO
        END IF
    END DO

END FUNCTION ForcaGlobaltri

END MODULE Corpotri_module