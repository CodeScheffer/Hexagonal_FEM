PROGRAM Bilinear
    USE precision_module, ONLY : DP
    USE Dados_Module
    USE Material_Module
    USE Elemento_Module, ONLY : MontaXY, Montagls, MontaKe, TesteDerivadas
    USE Elementotri_Module, ONLY : MontaXYtri, Montaglstri, MontaKetri
    USE Global_Module, ONLY : RigidezGlobal, ForcaCglobal, CalculaTensao, CalculaDeformacao 
    USE Apoios_Module
    USE Contorno_Functions_Module
    USE Contorno_Functions_Moduletri
	USE Corpo_module, ONLY : MontaFBe, ForcaGlobal, TesteFuncaoForma
	USE Corpotri_module, ONLY : MontaFBetri, ForcaGlobaltri
	USE TensaoDef, ONLY : CalculaTensaoElemento, CalculaDeformacaoElemento
	USE TensaoDeftri, ONLY : CalculaTensaoElementotri, CalculaDeformacaoElementotri

    IMPLICIT NONE

    INTEGER :: nn, ne, nf, np, na, nfb
    REAL(KIND=DP), ALLOCATABLE :: XY(:,:), MAT(:,:), ESP(:)
    INTEGER, ALLOCATABLE :: IJ(:,:)
    REAL(KIND=DP), ALLOCATABLE :: FC(:,:), P(:,:), AP(:,:), FB(:,:)

    REAL(KIND=DP), ALLOCATABLE :: Msigma_centro(:,:)
    REAL(KIND=DP), ALLOCATABLE :: Mepsilon_centro(:,:)
    REAL(KIND=DP), ALLOCATABLE :: Msigma_pontos_gauss_tri(:,:,:) 
    REAL(KIND=DP), ALLOCATABLE :: Mepsilon_pontos_gauss_tri(:,:,:) 
    
    REAL(KIND=DP), ALLOCATABLE :: Msigma_pontos_gauss_hexa(:,:,:) 
    REAL(KIND=DP), ALLOCATABLE :: Mepsilon_pontos_gauss_hexa(:,:,:) 

    REAL(KIND=DP), ALLOCATABLE :: K(:,:), F(:), FCg(:), U(:) 
	REAL(KIND=DP), ALLOCATABLE :: FTRI(:), FHEX(:), FBODY(:)
	REAL(KIND=DP), ALLOCATABLE :: FDH(:), FDT(:), FDISTRIB(:)
	
    WRITE(*,*) "Iniciando programa Bilinear."

    CALL Dados(nn, XY, ne, IJ, MAT, ESP, nf, FC, np, P, na, AP, nfb, FB)

	! Aloca a matriz de rigidez
    ALLOCATE(K(2*nn, 2*nn))
	! Aloca o vetor de forças 
    ALLOCATE(F(2*nn))
	! Aloca o vetor de forças concentrada globais
    ALLOCATE(FCg(2*nn))
	! Aloca o vetor de deslocamentos
    ALLOCATE(U(2*nn))
	! Aloca o vetor de forças de corpo globais do elemento Q4
	ALLOCATE(FTRI(2*nn))
	! Aloca o vetor de forças de corpo dos elementos Hexagonais 
    ALLOCATE(FHEX(2*nn))
	! Aloca o vetor global de forças de corpo
    ALLOCATE(FBODY(2*nn))
	! Aloca o vetor de forças distribuidas
	ALLOCATE(FDISTRIB(2*nn))
	! Aloca o vetor de forças distribuidas
    ALLOCATE(FDT(2*nn))
	! Aloca o vetor de forças distribuidas
    ALLOCATE(FDH(2*nn))

	! Aloca o vetor de tensões no centro dos elementos
    ALLOCATE(Msigma_centro(ne, 3))
	! Aloca o vetor de deformação do centro do elemento
    ALLOCATE(Mepsilon_centro(ne, 3))
    ALLOCATE(Msigma_pontos_gauss_tri(ne, 4, 3))
    ALLOCATE(Mepsilon_pontos_gauss_tri(ne, 4, 3))
    ALLOCATE(Msigma_pontos_gauss_hexa(ne, 7, 3)) 
    ALLOCATE(Mepsilon_pontos_gauss_hexa(ne, 7, 3)) 

    K = 0.0_DP
    F = 0.0_DP
    FCg = 0.0_DP
	FTRI = 0.0_DP
    FHEX = 0.0_DP
    FBODY = 0.0_DP
	FDH = 0.0_DP
    FDT = 0.0_DP
    FDISTRIB = 0.0_DP
    U = 0.0_DP
    Msigma_centro = 0.0_DP
    Mepsilon_centro = 0.0_DP
    Msigma_pontos_gauss_tri = 0.0_DP
    Mepsilon_pontos_gauss_tri = 0.0_DP
    Msigma_pontos_gauss_hexa = 0.0_DP
    Mepsilon_pontos_gauss_hexa = 0.0_DP


    K = RigidezGlobal(nn, ne, MAT, ESP, XY, IJ)

    FCg = ForcaCglobal(nn, nf, FC)

    FTRI = ForcaGlobaltri(nn, ESP, XY, IJ, nfb, FB)
	
	FHEX = ForcaGlobal(nn, ESP, XY, IJ, nfb, FB)
	
	FBODY = FHEX + FTRI

	FDH = ForcatGlobal(nn, ESP, XY, IJ, np, P)
	
	FDT = ForcatGlobaltri(nn, ESP, XY, IJ, np, P)

	FDISTRIB = FDH + FDT

    F = FBODY + FCg + FDISTRIB

    CALL TesteFuncaoForma()
    CALL TesteDerivadas()

    CALL AplicaCCH(nn, na, AP, K, F)

    CALL solve_linear_system(K, F, U, 2*nn)

    CALL imprime_resultados(U, nn)
    
    CALL CalculaTensao(ne, MAT, IJ, XY, U, Msigma_centro, &
                       Msigma_pontos_gauss_tri, Msigma_pontos_gauss_hexa)

    CALL CalculaDeformacao(ne, IJ, XY, U, Mepsilon_centro, &
                           Mepsilon_pontos_gauss_tri, Mepsilon_pontos_gauss_hexa)

    CALL imprime_tensoes_centro(Msigma_centro, Mepsilon_centro, ne)

    CALL imprime_tensoes_deformacoes_pontos_gauss(ne, IJ, &
                                                  Msigma_pontos_gauss_tri, Mepsilon_pontos_gauss_tri, &
                                                  Msigma_pontos_gauss_hexa, Mepsilon_pontos_gauss_hexa)
    
    CALL exporta_resultados_para_python(nn, XY, ne, IJ, U, Msigma_centro, &
	Mepsilon_centro, Msigma_pontos_gauss_tri, &
	Mepsilon_pontos_gauss_tri, Msigma_pontos_gauss_hexa, &
	Mepsilon_pontos_gauss_hexa)

    IF (ALLOCATED(XY)) DEALLOCATE(XY)
    IF (ALLOCATED(MAT)) DEALLOCATE(MAT)
    IF (ALLOCATED(ESP)) DEALLOCATE(ESP)
    IF (ALLOCATED(IJ)) DEALLOCATE(IJ)
    IF (nf > 0 .AND. ALLOCATED(FC)) DEALLOCATE(FC)
    IF (np > 0 .AND. ALLOCATED(P)) DEALLOCATE(P)
    IF (na > 0 .AND. ALLOCATED(AP)) DEALLOCATE(AP)
    IF (nfb > 0 .AND. ALLOCATED(FB)) DEALLOCATE(FB)
    
    IF (ALLOCATED(K)) DEALLOCATE(K)
    IF (ALLOCATED(F)) DEALLOCATE(F)
    IF (ALLOCATED(FCg)) DEALLOCATE(FCg)
    IF (ALLOCATED(U)) DEALLOCATE(U)
    
    IF (ALLOCATED(Msigma_centro)) DEALLOCATE(Msigma_centro)
    IF (ALLOCATED(Mepsilon_centro)) DEALLOCATE(Mepsilon_centro)
    IF (ALLOCATED(Msigma_pontos_gauss_tri)) DEALLOCATE(Msigma_pontos_gauss_tri)
    IF (ALLOCATED(Mepsilon_pontos_gauss_tri)) DEALLOCATE(Mepsilon_pontos_gauss_tri)
    IF (ALLOCATED(Msigma_pontos_gauss_hexa)) DEALLOCATE(Msigma_pontos_gauss_hexa)
    IF (ALLOCATED(Mepsilon_pontos_gauss_hexa)) DEALLOCATE(Mepsilon_pontos_gauss_hexa)

    WRITE(*,*) "Programa Bilinear finalizado e memória desalocada."

CONTAINS

	! Resolve o sistema linear, utilizando decomposição LU com pivoteamento parcial

    SUBROUTINE solve_linear_system(K_matrix, F_vector, U_solution, n_dofs)
        USE precision_module, ONLY : DP
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n_dofs
        REAL(KIND=DP), INTENT(INOUT) :: K_matrix(n_dofs, n_dofs)
        REAL(KIND=DP), INTENT(IN)    :: F_vector(n_dofs)
        REAL(KIND=DP), INTENT(OUT)   :: U_solution(n_dofs)

        INTEGER :: i, j, k, pivot_row
        REAL(KIND=DP) :: factor, max_pivot, temp_val
        INTEGER, ALLOCATABLE :: P_permutations(:)
        REAL(KIND=DP), ALLOCATABLE :: B_permuted(:)

        ALLOCATE(P_permutations(n_dofs))
        ALLOCATE(B_permuted(n_dofs))

        DO i = 1, n_dofs
            P_permutations(i) = i
            B_permuted(i) = F_vector(i)
        END DO

        WRITE(*,*) "Iniciando a Decomposição LU com Pivoteamento Parcial. n_dofs=", n_dofs

        DO k = 1, n_dofs - 1
            max_pivot = ABS(K_matrix(k, k))
            pivot_row = k

            DO i = k + 1, n_dofs
                IF (ABS(K_matrix(i, k)) > max_pivot) THEN
                    max_pivot = ABS(K_matrix(i, k))
                    pivot_row = i
                END IF
            END DO

            IF (max_pivot < 1.0E-18_DP) THEN
                WRITE(*,*) "ERRO: Pivô muito pequeno na coluna ", k, " (max_pivot = ", max_pivot, "). Sistema singular ou mal-condicionado."
                STOP "Falha na Decomposição LU (pivô muito pequeno)."
            END IF

            IF (pivot_row /= k) THEN
                WRITE(*,*) "Trocando linhas ", k, " e ", pivot_row, " para pivoteamento."
                DO j = 1, n_dofs
                    temp_val = K_matrix(k, j)
                    K_matrix(k, j) = K_matrix(pivot_row, j)
                    K_matrix(pivot_row, j) = temp_val
                END DO
                temp_val = P_permutations(k)
                P_permutations(k) = P_permutations(pivot_row)
                P_permutations(pivot_row) = NINT(temp_val) 
            END IF

            DO i = k + 1, n_dofs
                factor = K_matrix(i, k) / K_matrix(k, k)
                K_matrix(i, k) = factor 

                DO j = k + 1, n_dofs
                    K_matrix(i, j) = K_matrix(i, j) - factor * K_matrix(k, j)
                END DO
            END DO
        END DO

        IF (ABS(K_matrix(n_dofs, n_dofs)) < 1.0E-18_DP) THEN
            WRITE(*,*) "ERRO: Último pivô (diagonal) muito pequeno em K_matrix(", n_dofs, ",", n_dofs, "). Sistema singular ou mal-condicionado."
            STOP "Falha na Decomposição LU (último pivô muito pequeno)."
        END IF

        WRITE(*,*) "Decomposição LU concluída."

        DO i = 1, n_dofs
            B_permuted(i) = F_vector(P_permutations(i))
        END DO

        WRITE(*,*) "Iniciando Forward Substitution."
        DO i = 1, n_dofs
            temp_val = B_permuted(i)
            DO j = 1, i - 1
                temp_val = temp_val - K_matrix(i, j) * B_permuted(j)
            END DO
            B_permuted(i) = temp_val 
        END DO
        WRITE(*,*) "Forward Substitution concluída."

        WRITE(*,*) "Iniciando Backward Substitution."
        DO i = n_dofs, 1, -1
            temp_val = B_permuted(i)
            DO j = i + 1, n_dofs
                temp_val = temp_val - K_matrix(i, j) * U_solution(j)
            END DO
            U_solution(i) = temp_val / K_matrix(i, i) 
        END DO
        WRITE(*,*) "Backward Substitution concluída."

        IF (ALLOCATED(P_permutations)) DEALLOCATE(P_permutations)
        IF (ALLOCATED(B_permuted)) DEALLOCATE(B_permuted)

        WRITE(*,*) "Sistema resolvido com Decomposição LU."

    END SUBROUTINE solve_linear_system

! Imprime resultados dos deslocamentos no terminal

    SUBROUTINE imprime_resultados(U_in, nn_nodes)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nn_nodes
        REAL(KIND=DP), INTENT(IN) :: U_in(2 * nn_nodes)
        INTEGER :: i

        WRITE(*,*) " "
        WRITE(*,*) "===== RESULTADOS DOS DESLOCAMENTOS ====="
        DO i = 1, nn_nodes
            WRITE(*,*) "Nó ", i, ": Dx = ", U_in(2*i-1), " Dy = ", U_in(2*i)
        END DO
        WRITE(*,*) "========================================="
        WRITE(*,*) " "
    END SUBROUTINE imprime_resultados
    
! Imprime resultados no centro do elemento para tensões
	
    SUBROUTINE imprime_tensoes_centro(Msigma_centro_in, Mepsilon_centro_in, ne_elements)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne_elements
        REAL(KIND=DP), INTENT(IN) :: Msigma_centro_in(ne_elements, 3)
        REAL(KIND=DP), INTENT(IN) :: Mepsilon_centro_in(ne_elements, 3)

        INTEGER :: i

        WRITE(*,*) " "
        WRITE(*,*) "===== RESULTADOS DAS TENSÕES NO CENTRO DOS ELEMENTOS ====="
        WRITE(*,*) "Elemento | Sigma_XX (N/m^2) | Sigma_YY (N/m^2) | Sigma_XY (N/m^2)"
        WRITE(*,*) "----------------------------------------------------------------"
        DO i = 1, ne_elements
            WRITE(*, '(I5, 3(1X, F20.10))') i, Msigma_centro_in(i,1), Msigma_centro_in(i,2), Msigma_centro_in(i,3)
        END DO
        WRITE(*,*) "=============================================================="
        WRITE(*,*) " "

        WRITE(*,*) "===== RESULTADOS DAS DEFORMAÇÕES NO CENTRO DOS ELEMENTOS ====="
        WRITE(*,*) "Elemento | Epsilon_XX | Epsilon_YY | Epsilon_XY (Gamma_XY)"
        WRITE(*,*) "----------------------------------------------------------------"
        DO i = 1, ne_elements
            WRITE(*, '(I5, 3(1X, F20.10))') i, Mepsilon_centro_in(i,1), Mepsilon_centro_in(i,2), Mepsilon_centro_in(i,3)
        END DO
        WRITE(*,*) "=============================================================="
        WRITE(*,*) " "
    END SUBROUTINE imprime_tensoes_centro
    
! Imprime reultados nos pontos de Gauss para tensõe
	
    SUBROUTINE imprime_tensoes_deformacoes_pontos_gauss(ne_elements, IJ_in, &
                                                        Msigma_tri_in, Mepsilon_tri_in, &
                                                        Msigma_hexa_in, Mepsilon_hexa_in)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ne_elements
        INTEGER, INTENT(IN) :: IJ_in(:, :) 
        REAL(KIND=DP), INTENT(IN) :: Msigma_tri_in(ne_elements, 4, 3)
        REAL(KIND=DP), INTENT(IN) :: Mepsilon_tri_in(ne_elements, 4, 3)
        REAL(KIND=DP), INTENT(IN) :: Msigma_hexa_in(ne_elements, 7, 3) 
        REAL(KIND=DP), INTENT(IN) :: Mepsilon_hexa_in(ne_elements, 7, 3) 

        INTEGER :: i_ele, i_gp, element_type

        WRITE(*,*) " "
        WRITE(*,*) "===== TENSÕES E DEFORMAÇÕES NOS PONTOS DE GAUSS ====="

        DO i_ele = 1, ne_elements
            element_type = IJ_in(i_ele, 1)

            SELECT CASE (element_type)
            CASE (2) 
                WRITE(*,*) "----------------------------------------------------------------"
                WRITE(*,*) "Elemento Q4 No. ", i_ele
                DO i_gp = 1, 4
                    WRITE(*, '(A,I1,A,3F22.10,A,3F22.10)') "  Ponto Gauss ", i_gp, ": Sigma = (", &
                                     Msigma_tri_in(i_ele, i_gp, 1), Msigma_tri_in(i_ele, i_gp, 2), Msigma_tri_in(i_ele, i_gp, 3), &
                                     ") Epsilon = (", Mepsilon_tri_in(i_ele, i_gp, 1), Mepsilon_tri_in(i_ele, i_gp, 2), Mepsilon_tri_in(i_ele, i_gp, 3), ")"
                END DO
            CASE (0) 
                WRITE(*,*) "----------------------------------------------------------------"
                WRITE(*,*) "Elemento Hexagonal No. ", i_ele
                DO i_gp = 1, 7 
                    WRITE(*, '(A,I1,A,3F22.10,A,3F22.10)') "  Ponto Gauss ", i_gp, ": Sigma = (", &
                                     Msigma_hexa_in(i_ele, i_gp, 1), Msigma_hexa_in(i_ele, i_gp, 2), Msigma_hexa_in(i_ele, i_gp, 3), &
                                     ") Epsilon = (", Mepsilon_hexa_in(i_ele, i_gp, 1), Mepsilon_hexa_in(i_ele, i_gp, 2), Mepsilon_hexa_in(i_ele, i_gp, 3), ")"
                END DO
            CASE DEFAULT
                WRITE(*,*) "Tipo de elemento desconhecido para impressão de pontos de Gauss no elemento ", i_ele, ": ", element_type
                STOP "Tipo de elemento inválido na exportacao de conectividade."
            END SELECT
        END DO
        WRITE(*,*) "======================================================"
        WRITE(*,*) " "
    END SUBROUTINE imprime_tensoes_deformacoes_pontos_gauss
    
! Exporta .pos no diretório, para pós procesamento
	
    SUBROUTINE exporta_resultados_para_python(nn, XY, ne, IJ, U, &
                                              Msigma_centro, Mepsilon_centro, &
                                              Msigma_pontos_gauss_tri, Mepsilon_pontos_gauss_tri, &
                                              Msigma_pontos_gauss_hexa, Mepsilon_pontos_gauss_hexa)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nn, ne
        REAL(KIND=DP), INTENT(IN) :: XY(nn, 2)
        INTEGER, INTENT(IN) :: IJ(:, :) 
        REAL(KIND=DP), INTENT(IN) :: U(2*nn)
        REAL(KIND=DP), INTENT(IN) :: Msigma_centro(ne, 3)
        REAL(KIND=DP), INTENT(IN) :: Mepsilon_centro(ne, 3)
        REAL(KIND=DP), INTENT(IN) :: Msigma_pontos_gauss_tri(ne, 4, 3)
        REAL(KIND=DP), INTENT(IN) :: Mepsilon_pontos_gauss_tri(ne, 4, 3)
        REAL(KIND=DP), INTENT(IN) :: Msigma_pontos_gauss_hexa(ne, 7, 3) 
        REAL(KIND=DP), INTENT(IN) :: Mepsilon_pontos_gauss_hexa(ne, 7, 3) 

        INTEGER :: unit_num, i, j, k, iostat_val, element_type, num_nodes_in_elem
        CHARACTER(LEN=100) :: filename

        filename = "resultados.pos"
        OPEN(NEWUNIT=unit_num, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE', IOSTAT=iostat_val)

        IF (iostat_val /= 0) THEN
            WRITE(*,*) "ERRO: Nao foi possivel abrir o arquivo '", TRIM(filename), "' (Codigo IOSTAT: ", iostat_val, ")."
            WRITE(*,*) "Verifique as permissoes de escrita na pasta de destino."
            RETURN
        END IF

        WRITE(unit_num, *) "# Resultados da Simulacao FEM"
        WRITE(unit_num, *) "# Nos: ", nn
        WRITE(unit_num, *) "# Elementos: ", ne
        WRITE(unit_num, *) ""

        WRITE(unit_num, *) "# Coordenadas dos Nos (ID, X, Y)"
        DO i = 1, nn
            WRITE(unit_num, '(I5, 2(1X, F18.10))') i, XY(i,1), XY(i,2)
        END DO
        WRITE(unit_num, *) ""

        WRITE(unit_num, *) "# Conectividade dos Elementos (ID, Tipo_Elemento, No1, No2, ...)"
        DO i = 1, ne
            element_type = IJ(i, 1)

            SELECT CASE (element_type)
            CASE (2) 
                num_nodes_in_elem = 4
                WRITE(unit_num, '(I5, 1X, I2, *(1X, I5))') i, element_type, (IJ(i,j), j=2, 2+num_nodes_in_elem-1)
            CASE (0) 
                num_nodes_in_elem = 6
                WRITE(unit_num, '(I5, 1X, I2, *(1X, I5))') i, element_type, (IJ(i,j), j=2, 2+num_nodes_in_elem-1)
            CASE DEFAULT
                WRITE(*,*) "ERRO: Tipo de elemento desconhecido em exporta_resultados_para_python para elemento ", i, ": ", element_type
                STOP "Tipo de elemento invalido na exportacao de conectividade."
            END SELECT
        END DO
        WRITE(unit_num, *) ""

        WRITE(unit_num, *) "# Deslocamentos Nodais (ID do No, Dx, Dy)"
        DO i = 1, nn
            WRITE(unit_num, '(I5, 2(1X, F18.10))') i, U(2*i-1), U(2*i)
        END DO
        WRITE(unit_num, *) ""

        WRITE(unit_num, *) "# Tensoes nos Centros dos Elementos (ID do Elemento, Sigma_XX, Sigma_YY, Sigma_XY)"
        DO i = 1, ne
            WRITE(unit_num, '(I5, 3(1X, F30.20))') i, Msigma_centro(i,1), Msigma_centro(i,2), Msigma_centro(i,3)
        END DO
        WRITE(unit_num, *) ""

        WRITE(unit_num, *) "# Deformacoes nos Centros dos Elementos (ID do Elemento, Epsilon_XX, Epsilon_YY, Epsilon_XY)"
        DO i = 1, ne
            WRITE(unit_num, '(I5, 3(1X, F30.20))') i, Mepsilon_centro(i,1), Mepsilon_centro(i,2), Mepsilon_centro(i,3)
        END DO
        WRITE(unit_num, *) ""

        WRITE(unit_num, *) "# Tensoes e Deformacoes nos Pontos de Gauss"
        DO i = 1, ne
            element_type = IJ(i, 1)

            SELECT CASE (element_type)
            CASE (2) 
                WRITE(unit_num, '(A, I5, A)') "# Elemento Q4 ID: ", i, " (4 Pontos de Gauss)"
                DO j = 1, 4
                    WRITE(unit_num, '(A, I1, A, 3F18.10, A, 3F18.10)') "  GP ", j, ": Sigma=(", &
                                     Msigma_pontos_gauss_tri(i, j, 1), Msigma_pontos_gauss_tri(i, j, 2), Msigma_pontos_gauss_tri(i, j, 3), &
                                     ") Epsilon=(", Mepsilon_pontos_gauss_tri(i, j, 1), Mepsilon_pontos_gauss_tri(i, j, 2), Mepsilon_pontos_gauss_tri(i, j, 3), ")"
                END DO
            CASE (0) 
                WRITE(unit_num, '(A, I5, A)') "# Elemento Hexa ID: ", i, " (7 Pontos de Gauss)"
                DO j = 1, 7 
                    WRITE(unit_num, '(A, I1, A, 3F18.10, A, 3F18.10)') "  GP ", j, ": Sigma=(", &
                                     Msigma_pontos_gauss_hexa(i, j, 1), Msigma_pontos_gauss_hexa(i, j, 2), Msigma_pontos_gauss_hexa(i, j, 3), &
                                     ") Epsilon=(", Mepsilon_pontos_gauss_hexa(i, j, 1), Mepsilon_pontos_gauss_hexa(i, j, 2), Mepsilon_pontos_gauss_hexa(i, j, 3), ")"
                END DO
            CASE DEFAULT
            END SELECT
        END DO
        WRITE(unit_num, *) ""

        CLOSE(UNIT=unit_num)
        WRITE(*,*) "Resultados exportados para ", TRIM(filename)

    END SUBROUTINE exporta_resultados_para_python

END PROGRAM Bilinear
