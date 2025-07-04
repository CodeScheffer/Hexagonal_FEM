MODULE Corpo_module
    USE precision_module, ONLY : DP
    USE Material_Module, ONLY : MontaCEPT
    USE Elemento_Module, ONLY : MontaXY, MontaJ, MontadN, CorrigedN, MontaB, Montagls

    IMPLICIT NONE

    REAL(KIND=DP), PARAMETER :: IPI = ACOS(-1.0_DP)

CONTAINS

! Monta as funçõe de interpolação

FUNCTION MontaN(r, s) RESULT(N)

    REAL(KIND=DP), INTENT(IN) :: r, s

    REAL(KIND=DP) :: N(2,12)

    REAL(KIND=DP) :: N1, N2, N3, N4, N5, N6

    REAL(KIND=DP) :: den_tipo1
    REAL(KIND=DP) :: sqrt3_dp

    sqrt3_dp = SQRT(3.0_DP)

    den_tipo1 = 2.0_DP * (-r**2 - s**2 + 3.0_DP)

    IF (ABS(den_tipo1) < 1.0E-12_DP) THEN
        WRITE(*,*) 'ERRO: Denominador zero em MontaN!'
        STOP 'Denominador singular'
    END IF

    N1 = (1.0_DP/den_tipo1) * (1.0_DP - 2.0_DP*s/sqrt3_dp) * (1.0_DP + r - s/sqrt3_dp) * (1.0_DP + r + s/sqrt3_dp) * (1.0_DP + 2.0_DP*s/sqrt3_dp)
    N2 = (1.0_DP/den_tipo1) * (1.0_DP + r - s/sqrt3_dp) * (1.0_DP + r + s/sqrt3_dp) * (1.0_DP + 2.0_DP*s/sqrt3_dp) * (1.0_DP - r + s/sqrt3_dp)
    N3 = (1.0_DP/den_tipo1) * (1.0_DP + r + s/sqrt3_dp) * (1.0_DP + 2.0_DP*s/sqrt3_dp) * (1.0_DP - r + s/sqrt3_dp) * (1.0_DP - r - s/sqrt3_dp)
    N4 = (1.0_DP/den_tipo1) * (1.0_DP + 2.0_DP*s/sqrt3_dp) * (1.0_DP - r + s/sqrt3_dp) * (1.0_DP - r - s/sqrt3_dp) * (1.0_DP - 2.0_DP*s/sqrt3_dp)
    N5 = (1.0_DP/den_tipo1) * (1.0_DP - r + s/sqrt3_dp) * (1.0_DP - r - s/sqrt3_dp) * (1.0_DP - 2.0_DP*s/sqrt3_dp) * (1.0_DP + r - s/sqrt3_dp)
    N6 = (1.0_DP/den_tipo1) * (1.0_DP - r - s/sqrt3_dp) * (1.0_DP - 2.0_DP*s/sqrt3_dp) * (1.0_DP + r - s/sqrt3_dp) * (1.0_DP + r + s/sqrt3_dp)


    N = 0.0_DP

    N(1,1) = N1
    N(2,2) = N1

    N(1,3) = N2
    N(2,4) = N2

    N(1,5) = N3
    N(2,6) = N3

    N(1,7) = N4
    N(2,8) = N4

    N(1,9) = N5
    N(2,10) = N5

    N(1,11) = N6
    N(2,12) = N6

END FUNCTION MontaN

! Testa para ver se o somatório delas em um pondo é 1, para ver se tem algum erro com elas

SUBROUTINE TesteFuncaoForma()
    USE precision_module, ONLY: DP
    IMPLICIT NONE
    REAL(KIND=DP) :: r_test = 0.5_DP, s_test = 0.5_DP
    REAL(KIND=DP) :: N(2,12)
    INTEGER :: i
    
    N = MontaN(r_test, s_test)
    
    WRITE(*,*) 'TESTE DAS FUNÇÕES DE FORMA EM (0.5, 0.5)'
    WRITE(*,*) 'Soma de todas N_i:', SUM(N(1,:)), ' (deve ser ≈1.0)'
END SUBROUTINE TesteFuncaoForma

! Calcula o vetor de forças nodais equivalentes Fb devido a forças de corpo do elemento hexagonal, usando integração de Gauss

FUNCTION MontaFBe(ele, IJ, X_coords, Y_coords, te, b) RESULT(Fb)
    USE precision_module, ONLY : DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ele
    INTEGER, INTENT(IN) :: IJ(:,:)
    REAL(KIND=DP), INTENT(IN) :: X_coords(6), Y_coords(6), te
    REAL(KIND=DP), INTENT(IN) :: b(2)

    REAL(KIND=DP) :: Fb(12)
    REAL(KIND=DP) :: x_gauss(7), y_gauss(7), w_gauss(7)
    REAL(KIND=DP) :: N_matrix(2,12)
    REAL(KIND=DP) :: area_hex
    INTEGER :: i, k

    Fb = 0.0_DP

    IF (IJ(ele, 1) == 0) THEN

        ! Coordenadas dos pontos de Gauss
        y_gauss = (/ 0.0_DP, 0.64807407_DP, 0.64807407_DP, 0.0_DP, -0.64807407_DP, -0.64807407_DP, 0.0_DP /)
        x_gauss = (/ 0.748331477_DP, 0.374165739_DP, -0.374165739_DP, -0.748331477_DP, -0.374165739_DP, 0.374165739_DP, 0.0_DP /)
        w_gauss = (/ 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.255952381_DP /)

        ! Calcula área real do hexágono
        area_hex = 0.0_DP
        DO k = 1, 6
            area_hex = area_hex + (X_coords(k) * Y_coords(MOD(k,6)+1) - X_coords(MOD(k,6)+1) * Y_coords(k))
        END DO
        area_hex = 0.5_DP * ABS(area_hex)

        ! Loop de integração nos 7 pontos
        DO i = 1, 7
            N_matrix = MontaN(x_gauss(i), y_gauss(i))
            Fb = Fb + MATMUL(TRANSPOSE(N_matrix), b) * w_gauss(i) * te * area_hex
        END DO

    END IF

END FUNCTION MontaFBe

! Calcula o vetor global de forças de corpo F geradas no elementos hexagonal, com o dados de entrada

FUNCTION ForcaGlobal(nn, ESP, XY, IJ, nfb, FB_input) RESULT(F)
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
    REAL(KIND=DP) :: X(6), Y(6)
    REAL(KIND=DP) :: te
    REAL(KIND=DP) :: Fb(12)
    INTEGER :: ge(12)
    INTEGER :: element_type

    F = 0.0_DP

    DO i = 1, nfb
        ele = INT(FB_input(i,1))
        dir = INT(FB_input(i,2))
        val = FB_input(i,3)

        element_type = IJ(ele, 1)

        IF (element_type == 0) THEN
            b = 0.0_DP
            b(dir) = val

            CALL MontaXY(ele, IJ, XY, X, Y)
            te = ESP(ele)

            Fb = MontaFBe(ele, IJ, X, Y, te, b)

            ge = Montagls(ele,IJ)

            DO j = 1, 12
                F(ge(j)) = F(ge(j)) + Fb(j)
            END DO
        END IF
    END DO

END FUNCTION ForcaGlobal

END MODULE Corpo_module