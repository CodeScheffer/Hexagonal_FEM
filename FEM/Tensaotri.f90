MODULE TensaoDeftri
    USE precision_module, ONLY : DP
    USE Material_Module
	USE Elementotri_module

    IMPLICIT NONE

CONTAINS

! Calcula a tencao em um ponto (r, s) do elemento triangular

FUNCTION CalculaTensaoElementotri(r,s,X,Y,E,nu,Ue) RESULT(sigma)
        USE precision_module, ONLY : DP
        USE Material_Module, ONLY : MontaCEPT
        IMPLICIT NONE
        REAL(KIND=DP), INTENT(IN) :: r, s
        REAL(KIND=DP), INTENT(IN) :: X(4), Y(4)
        REAL(KIND=DP), INTENT(IN) :: E, nu
        REAL(KIND=DP), INTENT(IN) :: Ue(8)

        REAL(KIND=DP) :: Cv(3,3)
        REAL(KIND=DP) :: dNrs(2,4)
        REAL(KIND=DP) :: J(2,2)
        REAL(KIND=DP) :: dNxy(2,4)
        REAL(KIND=DP) :: B(3,8)
        REAL(KIND=DP) :: sigma(3)

        Cv = MontaCEPT(E,nu)

        dNrs = MontadNtri(r,s)

        J = MontaJtri(dNrs,X,Y)

        dNxy = CorrigedNtri(dNrs,J)

        B = MontaBtri(dNxy)

        ! Tensão no ponto: sigma = Cv * B * Ue
        sigma = MATMUL(Cv, MATMUL(B, Ue))

END FUNCTION CalculaTensaoElementotri

! Calcula a deformacao em um ponto (r, s) do elemento triangular

FUNCTION CalculaDeformacaoElementotri(r,s,X,Y,Ue) RESULT(epsilon)
        USE precision_module, ONLY : DP
        IMPLICIT NONE
        REAL(KIND=DP), INTENT(IN) :: r, s
        REAL(KIND=DP), INTENT(IN) :: X(4), Y(4)
        REAL(KIND=DP), INTENT(IN) :: Ue(8)

        REAL(KIND=DP) :: dNrs(2,4)
        REAL(KIND=DP) :: J(2,2)
        REAL(KIND=DP) :: dNxy(2,4)
        REAL(KIND=DP) :: B(3,8)
        REAL(KIND=DP) :: epsilon(3) ! Vetor de deformações (epsilon_x, epsilon_y, gamma_xy)

        dNrs = MontadNtri(r,s)

        J = MontaJtri(dNrs,X,Y)

        dNxy = CorrigedNtri(dNrs,J)

        B = MontaBtri(dNxy)

        ! Deformação no ponto: epsilon = B * Ue
        epsilon = MATMUL(B, Ue)

END FUNCTION CalculaDeformacaoElementotri

! Calcula a tencao no centro e nos pontos de Gauss

SUBROUTINE CalculaTensao_Tri(ne, MAT, IJ, XY, U, &
                                 Msigma_centro_out, Msigma_pontos_gauss_tri_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne
        REAL(KIND=DP), INTENT(IN) :: MAT(:,:), XY(:,:), U(:)
        INTEGER, INTENT(IN) :: IJ(:,:)
        REAL(KIND=DP), INTENT(INOUT) :: Msigma_centro_out(ne, 3)
        REAL(KIND=DP), INTENT(INOUT) :: Msigma_pontos_gauss_tri_out(ne, 4, 3)

        INTEGER :: ele, gp_idx
        REAL(KIND=DP) :: E_val, nu_val
        REAL(KIND=DP) :: X_coords_q4(4), Y_coords_q4(4)
        INTEGER :: gle_q4(8)
        REAL(KIND=DP) :: Ue_local_q4(8)
        REAL(KIND=DP) :: sigma_local(3)
        REAL(KIND=DP) :: sigma_ponto_gauss(3)
        
        INTEGER, PARAMETER :: NUM_PONTOS_GAUSS_Q4 = 4
        REAL(KIND=DP), PARAMETER :: GPS_COORD_Q4_LOCAL = 1.0_DP / SQRT(3.0_DP)
        REAL(KIND=DP) :: r_gauss_q4(NUM_PONTOS_GAUSS_Q4), s_gauss_q4(NUM_PONTOS_GAUSS_Q4)

        r_gauss_q4(1) = -GPS_COORD_Q4_LOCAL; s_gauss_q4(1) = -GPS_COORD_Q4_LOCAL
        r_gauss_q4(2) = GPS_COORD_Q4_LOCAL;  s_gauss_q4(2) = -GPS_COORD_Q4_LOCAL
        r_gauss_q4(3) = GPS_COORD_Q4_LOCAL;  s_gauss_q4(3) = GPS_COORD_Q4_LOCAL
        r_gauss_q4(4) = -GPS_COORD_Q4_LOCAL; s_gauss_q4(4) = GPS_COORD_Q4_LOCAL

        DO ele = 1, ne
            IF (IJ(ele, 1) == 2) THEN
                E_val = MAT(ele,1)
                nu_val = MAT(ele,2)

                CALL MontaXYtri(ele,IJ,XY,X_coords_q4,Y_coords_q4)
                gle_q4 = Montaglstri(ele,IJ)
                Ue_local_q4 = U(gle_q4)

                sigma_local = CalculaTensaoElementotri(0.0_DP, 0.0_DP, X_coords_q4, Y_coords_q4, &
                                                         E_val, nu_val, Ue_local_q4)
                Msigma_centro_out(ele,1) = sigma_local(1)
                Msigma_centro_out(ele,2) = sigma_local(2)
                Msigma_centro_out(ele,3) = sigma_local(3)

                DO gp_idx = 1, NUM_PONTOS_GAUSS_Q4
                    sigma_ponto_gauss = CalculaTensaoElementotri(r_gauss_q4(gp_idx), s_gauss_q4(gp_idx), &
                                                                  X_coords_q4, Y_coords_q4, E_val, nu_val, Ue_local_q4)
                    Msigma_pontos_gauss_tri_out(ele, gp_idx, 1) = sigma_ponto_gauss(1)
                    Msigma_pontos_gauss_tri_out(ele, gp_idx, 2) = sigma_ponto_gauss(2)
                    Msigma_pontos_gauss_tri_out(ele, gp_idx, 3) = sigma_ponto_gauss(3)
                END DO
            END IF
        END DO
END SUBROUTINE CalculaTensao_Tri

! Calcula a deformacao no centro e nos ponto de Gauss

SUBROUTINE CalculaDeformacao_Tri(ne, IJ, XY, U, &
                                     Mepsilon_centro_out, Mepsilon_pontos_gauss_tri_out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne
        INTEGER, INTENT(IN) :: IJ(:,:)
        REAL(KIND=DP), INTENT(IN) :: XY(:,:)
        REAL(KIND=DP), INTENT(IN) :: U(:)
        REAL(KIND=DP), INTENT(INOUT) :: Mepsilon_centro_out(ne, 3)
        REAL(KIND=DP), INTENT(INOUT) :: Mepsilon_pontos_gauss_tri_out(ne, 4, 3)

        INTEGER :: ele, gp_idx
        REAL(KIND=DP) :: X_coords_q4(4), Y_coords_q4(4)
        INTEGER :: gle_q4(8)
        REAL(KIND=DP) :: Ue_local_q4(8)
        REAL(KIND=DP) :: epsilon_local(3)
        REAL(KIND=DP) :: epsilon_ponto_gauss(3)

        INTEGER, PARAMETER :: NUM_PONTOS_GAUSS_Q4 = 4
        REAL(KIND=DP), PARAMETER :: GPS_COORD_Q4_LOCAL = 1.0_DP / SQRT(3.0_DP)
        REAL(KIND=DP) :: r_gauss_q4(NUM_PONTOS_GAUSS_Q4), s_gauss_q4(NUM_PONTOS_GAUSS_Q4)

        r_gauss_q4(1) = -GPS_COORD_Q4_LOCAL; s_gauss_q4(1) = -GPS_COORD_Q4_LOCAL
        r_gauss_q4(2) = GPS_COORD_Q4_LOCAL;  s_gauss_q4(2) = -GPS_COORD_Q4_LOCAL
        r_gauss_q4(3) = GPS_COORD_Q4_LOCAL;  s_gauss_q4(3) = GPS_COORD_Q4_LOCAL
        r_gauss_q4(4) = -GPS_COORD_Q4_LOCAL; s_gauss_q4(4) = GPS_COORD_Q4_LOCAL

        DO ele = 1, ne
            IF (IJ(ele, 1) == 2) THEN
                CALL MontaXYtri(ele,IJ,XY,X_coords_q4,Y_coords_q4)
                gle_q4 = Montaglstri(ele,IJ)
                Ue_local_q4 = U(gle_q4)

                epsilon_local = CalculaDeformacaoElementotri(0.0_DP, 0.0_DP, X_coords_q4, Y_coords_q4, Ue_local_q4)
                Mepsilon_centro_out(ele,1) = epsilon_local(1)
                Mepsilon_centro_out(ele,2) = epsilon_local(2)
                Mepsilon_centro_out(ele,3) = epsilon_local(3)

                DO gp_idx = 1, NUM_PONTOS_GAUSS_Q4
                    epsilon_ponto_gauss = CalculaDeformacaoElementotri(r_gauss_q4(gp_idx), s_gauss_q4(gp_idx), &
                                                                        X_coords_q4, Y_coords_q4, Ue_local_q4)
                    Mepsilon_pontos_gauss_tri_out(ele, gp_idx, 1) = epsilon_ponto_gauss(1)
                    Mepsilon_pontos_gauss_tri_out(ele, gp_idx, 2) = epsilon_ponto_gauss(2)
                    Mepsilon_pontos_gauss_tri_out(ele, gp_idx, 3) = epsilon_ponto_gauss(3)
                END DO
            END IF
        END DO
END SUBROUTINE CalculaDeformacao_Tri

END MODULE TensaoDeftri