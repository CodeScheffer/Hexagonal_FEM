MODULE Contorno_Functions_Module
    USE precision_module, ONLY : DP
    USE Elemento_Module, ONLY : Montagls

    IMPLICIT NONE

TYPE face_data_type
    REAL(KIND=DP) :: dJ
    REAL(KIND=DP), DIMENSION(2) :: n_vec, v_vec
END TYPE face_data_type

CONTAINS

! Retorna os dois nós que formam a aresta (face) do elemento hexagonal

FUNCTION FACES(face_in) RESULT(p_nodes)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: face_in

    INTEGER, DIMENSION(2) :: p_nodes

    IF (face_in == 1) THEN
        p_nodes(1) = 1
        p_nodes(2) = 2
    ELSE IF (face_in == 2) THEN
        p_nodes(1) = 2
        p_nodes(2) = 3
    ELSE IF (face_in == 3) THEN
        p_nodes(1) = 3
        p_nodes(2) = 4
    ELSE IF (face_in == 4) THEN
        p_nodes(1) = 4
        p_nodes(2) = 5
    ELSE IF (face_in == 5) THEN
        p_nodes(1) = 5
        p_nodes(2) = 6
    ELSE IF (face_in == 6) THEN
        p_nodes(1) = 6
        p_nodes(2) = 1
    ELSE
        WRITE(*,*) 'ERRO em FACES:: face ', face_in, ' não existe'
        STOP 'Erro na função FACES'
    END IF
END FUNCTION FACES

! Calcula a matriz de interpolação N numa face específica do elemento hexagonal, em função da coordenada local

FUNCTION MontaN_Face(zeta_in, p1_in, p2_in, face_id_in) RESULT(N_mat_face)

    USE precision_module, ONLY : DP
    USE Corpo_module, ONLY : MontaN

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: zeta_in
    INTEGER, INTENT(IN) :: p1_in, p2_in
    INTEGER, INTENT(IN) :: face_id_in

    REAL(KIND=DP), DIMENSION(2, 12) :: N_mat_face
    REAL(KIND=DP) :: r_coord, s_coord
    
    REAL(KIND=DP), DIMENSION(2, 12) :: N_matrix_from_MontaN 

    REAL(KIND=DP), DIMENSION(6) :: N_scalar_values 

    INTEGER :: k_node

    N_mat_face = 0.0_DP

    SELECT CASE (face_id_in)
        CASE (1)
            r_coord = zeta_in
            s_coord = -1.0_DP
        CASE (2)
            r_coord = 1.0_DP
            s_coord = zeta_in
        CASE (3)
            r_coord = -zeta_in
            s_coord = 1.0_DP
        CASE (4)
            r_coord = -1.0_DP
            s_coord = zeta_in
        CASE (5)
            r_coord = zeta_in
            s_coord = -1.0_DP
        CASE (6)
            r_coord = 1.0_DP
            s_coord = -zeta_in
        CASE DEFAULT
            WRITE(*,*) 'ERRO em MontaN_Face: face_id_in ', face_id_in, ' não é válida.'
            STOP 'Erro de mapeamento em MontaN_Face'
    END SELECT

    N_matrix_from_MontaN = MontaN(r_coord, s_coord) 

    DO k_node = 1, 6

        N_scalar_values(k_node) = N_matrix_from_MontaN(1, 2 * k_node - 1) 
    END DO

    DO k_node = 1, 6
        N_mat_face(1, 2 * k_node - 1) = N_scalar_values(k_node) 
        N_mat_face(2, 2 * k_node)     = N_scalar_values(k_node) 
    END DO

END FUNCTION MontaN_Face

! Calcula o comprimento, vetor tangente e vetor normal de uma face (aresta) do elemento

FUNCTION Facesnv(p1_in, p2_in, e_in, IJ_in, XY_in) RESULT(face_data)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p1_in, p2_in, e_in

    INTEGER, DIMENSION(:,:), INTENT(IN) :: IJ_in

    REAL(KIND=DP), DIMENSION(:,:), INTENT(IN) :: XY_in

    TYPE(face_data_type) :: face_data

    INTEGER :: n1, n2
    REAL(KIND=DP) :: x1, y1, x2, y2
    REAL(KIND=DP) :: deltaX, deltaY, L

    n1 = IJ_in(e_in, p1_in + 1)
    n2 = IJ_in(e_in, p2_in + 1)

    x1 = XY_in(n1, 1)
    y1 = XY_in(n1, 2)
    x2 = XY_in(n2, 1)
    y2 = XY_in(n2, 2)

    deltaX = x2 - x1
    deltaY = y2 - y1

    L = SQRT(deltaX**2 + deltaY**2)

    face_data%dJ = L / 2.0_DP

    IF (L > 0.0_DP) THEN
        face_data%v_vec(1) = (1.0_DP / L) * deltaX
        face_data%v_vec(2) = (1.0_DP / L) * deltaY
    ELSE
        face_data%v_vec = 0.0_DP
    END IF

    IF (L > 0.0_DP) THEN
        face_data%n_vec(1) = (1.0_DP / L) * deltaY
        face_data%n_vec(2) = (1.0_DP / L) * (-deltaX)
    ELSE
        face_data%n_vec = 0.0_DP
    END IF
END FUNCTION Facesnv

! Calcula o vetor de forças equivalentes nodais F devido a um carregamento distribuído linear aplicado na face do elemento hexagonal, usando integração de Gauss

FUNCTION Ftblinear(dJ_in, vetor_dir_local_in, p1_in, p2_in, face_id_in, carregamento_in, te_in) RESULT(F_vec)

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: dJ_in, carregamento_in, te_in
    REAL(KIND=DP), DIMENSION(2), INTENT(IN) :: vetor_dir_local_in
    INTEGER, INTENT(IN) :: p1_in, p2_in, face_id_in

    REAL(KIND=DP), DIMENSION(12) :: F_vec

    INTEGER, PARAMETER :: NUM_GAUSS_PTS_1D = 3
    REAL(KIND=DP), DIMENSION(NUM_GAUSS_PTS_1D) :: pg, wg
    
    REAL(KIND=DP) :: zeta_gauss
    REAL(KIND=DP), DIMENSION(2, 12) :: N_gauss
    INTEGER :: i

    pg(1) = -SQRT(3.0_DP / 5.0_DP)
    pg(2) = 0.0_DP
    pg(3) = SQRT(3.0_DP / 5.0_DP)
    
    wg(1) = 5.0_DP / 9.0_DP
    wg(2) = 8.0_DP / 9.0_DP
    wg(3) = 5.0_DP / 9.0_DP

    F_vec = 0.0_DP

    DO i = 1, NUM_GAUSS_PTS_1D
        zeta_gauss = pg(i)

        N_gauss = MontaN_Face(zeta_gauss, p1_in, p2_in, face_id_in)

        F_vec = F_vec + carregamento_in * MATMUL(TRANSPOSE(N_gauss), vetor_dir_local_in) * dJ_in * te_in * wg(i)
    END DO
    
END FUNCTION Ftblinear

! Calcula o vetor global de forças equivalentes Fglobal geradas pelo carregamentos distribuídos aplicados nas faces de elementos hexagonais

FUNCTION ForcatGlobal(nn_in, ESP_in, XY_in, IJ_in, np_in, P_in) RESULT(F_global)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nn_in, np_in
    REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: ESP_in
    REAL(KIND=DP), DIMENSION(:,:), INTENT(IN) :: XY_in
    INTEGER, DIMENSION(:,:), INTENT(IN) :: IJ_in
    REAL(KIND=DP), DIMENSION(:,:), INTENT(IN) :: P_in

    REAL(KIND=DP), DIMENSION(2 * nn_in) :: F_global

    INTEGER :: i, ele, face, dir
    REAL(KIND=DP) :: val, te

    INTEGER, DIMENSION(2) :: p_nodes_face
    REAL(KIND=DP), DIMENSION(2) :: vetor_dir_final
    REAL(KIND=DP), DIMENSION(12) :: F_local_elementar
    INTEGER, DIMENSION(12) :: global_dofs
    TYPE(face_data_type) :: current_face_data
    INTEGER :: j, p1_local, p2_local

    F_global = 0.0_DP

    DO i = 1, np_in
        ele = NINT(P_in(i,1))
        face = NINT(P_in(i,2))
        dir = NINT(P_in(i,3))
        val = P_in(i,4)

        te = ESP_in(ele)

        p_nodes_face = FACES(face)
        p1_local = p_nodes_face(1)
        p2_local = p_nodes_face(2)

        current_face_data = Facesnv(p1_local, p2_local, ele, IJ_in, XY_in)

        IF (dir == 1) THEN
            vetor_dir_final = current_face_data%n_vec
        ELSE
            vetor_dir_final = current_face_data%v_vec
        END IF

        F_local_elementar = Ftblinear(current_face_data%dJ, vetor_dir_final, &
                                       p1_local, p2_local, face, val, te)

        global_dofs = Montagls(ele, IJ_in)

        DO j = 1, 12
            F_global(global_dofs(j)) = F_global(global_dofs(j)) + F_local_elementar(j)
        END DO
    END DO
END FUNCTION ForcatGlobal

END MODULE Contorno_Functions_Module