MODULE TensaoDef
    USE precision_module, ONLY : DP
    USE Material_Module
	USE Elemento_module

    IMPLICIT NONE

CONTAINS

! Função que calcula a tensão em um ponto (r,s) de um elemento hexagonal

FUNCTION CalculaTensaoElemento(r,s,X,Y,E,nu,Ue) RESULT(sigma)
    USE precision_module, ONLY : DP
    USE Material_Module, ONLY : MontaCEPT
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: r, s
    REAL(KIND=DP), INTENT(IN) :: X(6), Y(6)
    REAL(KIND=DP), INTENT(IN) :: E, nu
    REAL(KIND=DP), INTENT(IN) :: Ue(12)

    REAL(KIND=DP) :: Cv(3,3)
    REAL(KIND=DP) :: dNrs(2,6)
    REAL(KIND=DP) :: J(2,2)
    REAL(KIND=DP) :: dNxy(2,6)
    REAL(KIND=DP) :: B(3,12)
    REAL(KIND=DP) :: sigma(3)

    Cv = MontaCEPT(E,nu)

    CALL MontadN(r,s, dNrs)

    J = MontaJ(dNrs,X,Y)

    dNxy = CorrigedN(dNrs,J)

    B = MontaB(dNxy)

    ! Tensão no ponto: sigma = Cv * B * Ue
    sigma = MATMUL(Cv, MATMUL(B, Ue))

END FUNCTION CalculaTensaoElemento
    
! Função que calcula a deformação em um ponto (r,s) de um elemento hexagonal
	
FUNCTION CalculaDeformacaoElemento(r,s,X,Y,Ue) RESULT(epsilon)
    USE precision_module, ONLY : DP
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: r, s
    REAL(KIND=DP), INTENT(IN) :: X(6), Y(6)
    REAL(KIND=DP), INTENT(IN) :: Ue(12)

    REAL(KIND=DP) :: dNrs(2,6)
    REAL(KIND=DP) :: J(2,2)
    REAL(KIND=DP) :: dNxy(2,6)
    REAL(KIND=DP) :: B(3,12)
    REAL(KIND=DP) :: epsilon(3) ! Vetor de deformações (epsilon_x, epsilon_y, gamma_xy)

    CALL MontadN(r,s, dNrs)

    J = MontaJ(dNrs,X,Y)

    dNxy = CorrigedN(dNrs,J)

    B = MontaB(dNxy)

    ! Deformação no ponto: epsilon = B * Ue
    epsilon = MATMUL(B, Ue)

END FUNCTION CalculaDeformacaoElemento

! Utiliza a função CalculaTensaoElemento para calcular a tensão nos ponto de Gauss e no centro do elemento

SUBROUTINE CalculaTensao_Hexa(ne, MAT, IJ, XY, U, &
                                  Msigma_centro_out, Msigma_pontos_gauss_hexa_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne
        REAL(KIND=DP), INTENT(IN) :: MAT(:,:), XY(:,:), U(:)
        INTEGER, INTENT(IN) :: IJ(:,:)
        REAL(KIND=DP), INTENT(INOUT) :: Msigma_centro_out(ne, 3)
        REAL(KIND=DP), INTENT(INOUT) :: Msigma_pontos_gauss_hexa_out(ne, 7, 3)

        INTEGER :: ele, gp_idx
        REAL(KIND=DP) :: E_val, nu_val
        REAL(KIND=DP) :: X_coords_hexa(6), Y_coords_hexa(6)
        INTEGER :: gle_hexa(12)
        REAL(KIND=DP) :: Ue_local_hexa(12)
        REAL(KIND=DP) :: sigma_local(3)
        REAL(KIND=DP) :: sigma_ponto_gauss(3)
        
        INTEGER, PARAMETER :: NUM_GAUSS_PTS_HEXA = 7
        REAL(KIND=DP), PARAMETER :: PI_VAL = ACOS(-1.0_DP)
        REAL(KIND=DP) :: r_hexa_gps(NUM_GAUSS_PTS_HEXA)
        REAL(KIND=DP) :: s_hexa_gps(NUM_GAUSS_PTS_HEXA)
        REAL(KIND=DP), PARAMETER :: RING_RADIUS_HEXA = 0.748331477354788_DP
        INTEGER :: j_angle_idx
        INTEGER :: current_gp_idx

        r_hexa_gps(1) = 0.0_DP
        s_hexa_gps(1) = 0.0_DP 

        DO j_angle_idx = 0, 5 
            current_gp_idx = j_angle_idx + 2 
            r_hexa_gps(current_gp_idx) = RING_RADIUS_HEXA * COS(REAL(j_angle_idx, KIND=DP) * PI_VAL / 3.0_DP)
            s_hexa_gps(current_gp_idx) = RING_RADIUS_HEXA * SIN(REAL(j_angle_idx, KIND=DP) * PI_VAL / 3.0_DP)
        END DO
        
        DO ele = 1, ne

            IF (IJ(ele, 1) == 0) THEN
                E_val = MAT(ele,1)
                nu_val = MAT(ele,2)

                CALL MontaXY(ele,IJ,XY,X_coords_hexa,Y_coords_hexa)
                gle_hexa = Montagls(ele,IJ)
                Ue_local_hexa = U(gle_hexa)

                sigma_local = CalculaTensaoElemento(0.0_DP, 0.0_DP, X_coords_hexa, Y_coords_hexa, &
                                                     E_val, nu_val, Ue_local_hexa)
                Msigma_centro_out(ele,1) = sigma_local(1)
                Msigma_centro_out(ele,2) = sigma_local(2)
                Msigma_centro_out(ele,3) = sigma_local(3)

                DO gp_idx = 1, NUM_GAUSS_PTS_HEXA
                    sigma_ponto_gauss = CalculaTensaoElemento(r_hexa_gps(gp_idx), s_hexa_gps(gp_idx), &
                                                                X_coords_hexa, Y_coords_hexa, E_val, nu_val, Ue_local_hexa)
                    Msigma_pontos_gauss_hexa_out(ele, gp_idx, 1) = sigma_ponto_gauss(1)
                    Msigma_pontos_gauss_hexa_out(ele, gp_idx, 2) = sigma_ponto_gauss(2)
                    Msigma_pontos_gauss_hexa_out(ele, gp_idx, 3) = sigma_ponto_gauss(3)
                END DO
            END IF
        END DO
END SUBROUTINE CalculaTensao_Hexa

! Utiliza a função CalculaDeformacaoElemento para calcular a deformacao nos ponto de Gauss e no centro do elemento

SUBROUTINE CalculaDeformacao_Hexa(ne, IJ, XY, U, &
                                      Mepsilon_centro_out, Mepsilon_pontos_gauss_hexa_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne
        INTEGER, INTENT(IN) :: IJ(:,:)
        REAL(KIND=DP), INTENT(IN) :: XY(:,:)
        REAL(KIND=DP), INTENT(IN) :: U(:)
        REAL(KIND=DP), INTENT(INOUT) :: Mepsilon_centro_out(ne, 3)
        REAL(KIND=DP), INTENT(INOUT) :: Mepsilon_pontos_gauss_hexa_out(ne, 7, 3)

        INTEGER :: ele, gp_idx
        REAL(KIND=DP) :: X_coords_hexa(6), Y_coords_hexa(6)
        INTEGER :: gle_hexa(12)
        REAL(KIND=DP) :: Ue_local_hexa(12)
        REAL(KIND=DP) :: epsilon_local(3)
        REAL(KIND=DP) :: epsilon_ponto_gauss(3)

        INTEGER, PARAMETER :: NUM_GAUSS_PTS_HEXA = 7
        REAL(KIND=DP), PARAMETER :: PI_VAL = ACOS(-1.0_DP)
        REAL(KIND=DP) :: r_hexa_gps(NUM_GAUSS_PTS_HEXA)
        REAL(KIND=DP) :: s_hexa_gps(NUM_GAUSS_PTS_HEXA)
        REAL(KIND=DP), PARAMETER :: RING_RADIUS_HEXA = 0.748331477354788_DP
        INTEGER :: j_angle_idx
        INTEGER :: current_gp_idx

        r_hexa_gps(1) = 0.0_DP
        s_hexa_gps(1) = 0.0_DP

        DO j_angle_idx = 0, 5
            current_gp_idx = j_angle_idx + 2
            r_hexa_gps(current_gp_idx) = RING_RADIUS_HEXA * COS(REAL(j_angle_idx, KIND=DP) * PI_VAL / 3.0_DP)
            s_hexa_gps(current_gp_idx) = RING_RADIUS_HEXA * SIN(REAL(j_angle_idx, KIND=DP) * PI_VAL / 3.0_DP)
        END DO
        
        DO ele = 1, ne

            IF (IJ(ele, 1) == 0) THEN
                CALL MontaXY(ele,IJ,XY,X_coords_hexa,Y_coords_hexa)
                gle_hexa = Montagls(ele,IJ)
                Ue_local_hexa = U(gle_hexa)

                epsilon_local = CalculaDeformacaoElemento(0.0_DP, 0.0_DP, X_coords_hexa, Y_coords_hexa, Ue_local_hexa)
                Mepsilon_centro_out(ele,1) = epsilon_local(1)
                Mepsilon_centro_out(ele,2) = epsilon_local(2)
                Mepsilon_centro_out(ele,3) = epsilon_local(3)

                DO gp_idx = 1, NUM_GAUSS_PTS_HEXA
                    epsilon_ponto_gauss = CalculaDeformacaoElemento(r_hexa_gps(gp_idx), s_hexa_gps(gp_idx), &
                                                                     X_coords_hexa, Y_coords_hexa, Ue_local_hexa)
                    Mepsilon_pontos_gauss_hexa_out(ele, gp_idx, 1) = epsilon_ponto_gauss(1)
                    Mepsilon_pontos_gauss_hexa_out(ele, gp_idx, 2) = epsilon_ponto_gauss(2)
                    Mepsilon_pontos_gauss_hexa_out(ele, gp_idx, 3) = epsilon_ponto_gauss(3)
                END DO
            END IF
        END DO
END SUBROUTINE CalculaDeformacao_Hexa

END MODULE TensaoDef