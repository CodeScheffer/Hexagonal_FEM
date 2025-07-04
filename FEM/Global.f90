MODULE Global_Module
    USE precision_module, ONLY : DP
    USE Elemento_Module
    USE Elementotri_Module
    USE Material_Module, ONLY : MontaCEPT
	USE TensaoDef
	USE TensaoDeftri

    IMPLICIT NONE

CONTAINS

! Monta a matriz de rigidez global para o elemento hexagonal

FUNCTION RigidezGlobal_Hexa(nn, ne, MAT, ESP, XY, IJ) RESULT(K_global_hexa)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nn, ne
        REAL(KIND=DP), INTENT(IN) :: MAT(:,:), ESP(:), XY(:,:)
        INTEGER, INTENT(IN) :: IJ(:,:)

        REAL(KIND=DP) :: K_global_hexa(2*nn, 2*nn)
        REAL(KIND=DP) :: Ke_local(12,12)
        REAL(KIND=DP) :: X_coords_hexa(6), Y_coords_hexa(6)
        INTEGER :: gle(12)
        INTEGER :: i, j, ele                     
        REAL(KIND=DP) :: E_val, nu_val, te_val
        INTEGER :: num_local_dofs_hexa

        K_global_hexa = 0.0_DP
        num_local_dofs_hexa = 12

        DO ele = 1, ne
            ! Verifica se o elemento é do tipo Hexagonal (0)
            IF (IJ(ele, 1) == 0) THEN
                E_val = MAT(ele,1)
                nu_val = MAT(ele,2)
                te_val = ESP(ele)

                ! Monta as coordenadas locais para o elemento hexagonal
                CALL MontaXY(ele, IJ, XY, X_coords_hexa, Y_coords_hexa)
                ! Monta a matriz de rigidez local para o elemento hexagonal
                CALL MontaKe(X_coords_hexa, Y_coords_hexa, E_val, nu_val, te_val, Ke_local)
                ! Obtém os graus de liberdade globais para o elemento hexagonal
                gle = Montagls(ele, IJ) 

                ! Mapeia Ke_local para K_global_hexa
                DO i = 1, num_local_dofs_hexa
                    DO j = 1, num_local_dofs_hexa
                        K_global_hexa(gle(i), gle(j)) = K_global_hexa(gle(i), gle(j)) + Ke_local(i,j)
                    END DO
                END DO
            END IF
        END DO

END FUNCTION RigidezGlobal_Hexa

!Monta a matriz de rigidez global para o elemento Q4/Q4 Degenerado

FUNCTION RigidezGlobal_Tri(nn, ne, MAT, ESP, XY, IJ) RESULT(K_global_tri)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nn, ne
        REAL(KIND=DP), INTENT(IN) :: MAT(:,:), ESP(:), XY(:,:)
        INTEGER, INTENT(IN) :: IJ(:,:)

        REAL(KIND=DP) :: K_global_tri(2*nn, 2*nn)
        REAL(KIND=DP) :: Ke_local(8,8)
        REAL(KIND=DP) :: X_coords_q4(4), Y_coords_q4(4)
        INTEGER :: gle(8)
        INTEGER :: i, j, ele                     
        REAL(KIND=DP) :: E_val, nu_val, te_val
        INTEGER :: num_local_dofs_tri

        K_global_tri = 0.0_DP
        num_local_dofs_tri = 8

        DO ele = 1, ne
            IF (IJ(ele, 1) == 2) THEN
                E_val = MAT(ele,1)
                nu_val = MAT(ele,2)
                te_val = ESP(ele)

                CALL MontaXYtri(ele, IJ, XY, X_coords_q4, Y_coords_q4)
                CALL MontaKetri(ele, IJ, X_coords_q4, Y_coords_q4, E_val, nu_val, te_val, Ke_local)
                gle = Montaglstri(ele, IJ) 

                DO i = 1, num_local_dofs_tri
                    DO j = 1, num_local_dofs_tri
                        K_global_tri(gle(i), gle(j)) = K_global_tri(gle(i), gle(j)) + Ke_local(i,j)
                    END DO
                END DO
            END IF
        END DO

END FUNCTION RigidezGlobal_Tri

! Soma a matriz global de cada tipo de elemento para gerar a matriz de rigidez global

FUNCTION RigidezGlobal(nn, ne, MAT, ESP, XY, IJ) RESULT(K_total)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nn, ne
        REAL(KIND=DP), INTENT(IN) :: MAT(:,:), ESP(:), XY(:,:)
        INTEGER, INTENT(IN) :: IJ(:,:)

        REAL(KIND=DP) :: K_total(2*nn, 2*nn)
        REAL(KIND=DP) :: K_hexa_contrib(2*nn, 2*nn)
        REAL(KIND=DP) :: K_tri_contrib(2*nn, 2*nn)
        INTEGER :: ele

        K_total = 0.0_DP

        K_hexa_contrib = RigidezGlobal_Hexa(nn, ne, MAT, ESP, XY, IJ)

        K_tri_contrib = RigidezGlobal_Tri(nn, ne, MAT, ESP, XY, IJ)

        K_total = K_hexa_contrib + K_tri_contrib

        DO ele = 1, ne
            SELECT CASE (IJ(ele, 1))
            CASE (0, 2)
            CASE DEFAULT
                WRITE(*,*) 'ERRO em RigidezGlobal: Tipo de elemento desconhecido para o elemento ', ele, ': ', IJ(ele, 1)
                STOP 'Tipo de elemento inválido na montagem global.'
            END SELECT
        END DO

END FUNCTION RigidezGlobal
        
! Monta o vetor de forças concentradas (nodais) global
		
FUNCTION ForcaCglobal(nn, nf, FC_input) RESULT(FC_global)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nn, nf
    REAL(KIND=DP), INTENT(IN) :: FC_input(:,:)

    REAL(KIND=DP) :: FC_global(2*nn)

    INTEGER :: i, no, dir
    REAL(KIND=DP) :: valor_forca

    FC_global = 0.0_DP

    DO i = 1, nf
        no      = INT(FC_input(i,1)) ! Nó onde a força é aplicada
        dir     = INT(FC_input(i,2)) ! Direção (1 para X, 2 para Y)
        valor_forca = FC_input(i,3)  ! Valor da força

        FC_global(2*(no-1) + dir) = FC_global(2*(no-1) + dir) + valor_forca
    END DO

END FUNCTION ForcaCglobal

! Calcula as deformações nos elementos da malha, no centro e nos pontos de Gauss, para elementos triangulares e hexaédricos

SUBROUTINE CalculaDeformacao(ne, IJ, XY, U, Mepsilon_centro, Mepsilon_pontos_gauss_tri, &
                                  Mepsilon_pontos_gauss_hexa)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne
        INTEGER, INTENT(IN) :: IJ(:,:)
        REAL(KIND=DP), INTENT(IN) :: XY(:,:)
        REAL(KIND=DP), INTENT(IN) :: U(:)

        REAL(KIND=DP), INTENT(OUT) :: Mepsilon_centro(ne, 3)
        REAL(KIND=DP), INTENT(OUT) :: Mepsilon_pontos_gauss_tri(ne, 4, 3)
        REAL(KIND=DP), INTENT(OUT) :: Mepsilon_pontos_gauss_hexa(ne, 7, 3) 

        Mepsilon_centro = 0.0_DP
        Mepsilon_pontos_gauss_tri = 0.0_DP
        Mepsilon_pontos_gauss_hexa = 0.0_DP

        CALL CalculaDeformacao_Hexa(ne, IJ, XY, U, Mepsilon_centro, Mepsilon_pontos_gauss_hexa)

        CALL CalculaDeformacao_Tri(ne, IJ, XY, U, Mepsilon_centro, Mepsilon_pontos_gauss_tri)
        
END SUBROUTINE CalculaDeformacao

! Calcula as tensões nos elementos da malha, no centro e nos pontos de Gauss, para elementos triangulares e hexaédricos

SUBROUTINE CalculaTensao(ne, MAT, IJ, XY, U, Msigma_centro, Msigma_pontos_gauss_tri, &
                             Msigma_pontos_gauss_hexa)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne
        REAL(KIND=DP), INTENT(IN) :: MAT(:,:), XY(:,:), U(:)
        INTEGER, INTENT(IN) :: IJ(:,:)

        REAL(KIND=DP), INTENT(OUT) :: Msigma_centro(ne, 3)
        REAL(KIND=DP), INTENT(OUT) :: Msigma_pontos_gauss_tri(ne, 4, 3)
        REAL(KIND=DP), INTENT(OUT) :: Msigma_pontos_gauss_hexa(ne, 7, 3) 

        Msigma_centro = 0.0_DP
        Msigma_pontos_gauss_tri = 0.0_DP
        Msigma_pontos_gauss_hexa = 0.0_DP

        CALL CalculaTensao_Hexa(ne, MAT, IJ, XY, U, Msigma_centro, Msigma_pontos_gauss_hexa)
        
        CALL CalculaTensao_Tri(ne, MAT, IJ, XY, U, Msigma_centro, Msigma_pontos_gauss_tri)
        
END SUBROUTINE CalculaTensao
	
END MODULE Global_Module
