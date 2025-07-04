MODULE Contorno_Functions_Moduletri
    USE precision_module, ONLY : DP
    USE Elementotri_Module, ONLY : Montaglstri 
    IMPLICIT NONE

    TYPE face_data_type
        REAL(KIND=DP) :: dJ
        REAL(KIND=DP), DIMENSION(2) :: n_vec, v_vec
    END TYPE face_data_type

CONTAINS

! Retorna os dois nós que formam a aresta (face) do elemento Q4

FUNCTION FACES_TRI(face_in) RESULT(p_nodes)

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
        p_nodes(2) = 1
    ELSE
        WRITE(*,*) 'ERRO em FACES_TRI: face ', face_in, ' não existe para elemento de 4 nós.'
        STOP 'Erro na função FACES_TRI'
    END IF
END FUNCTION FACES_TRI

! Calcula a matriz de interpolação N numa face específica do elemento Q4, em função da coordenada local

FUNCTION MontaN_Facetri(zeta_in, p1_in, p2_in) RESULT(N_mat_face)

    USE precision_module, ONLY : DP
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: zeta_in
    INTEGER, INTENT(IN) :: p1_in, p2_in

    REAL(KIND=DP), DIMENSION(2, 8) :: N_mat_face 

    REAL(KIND=DP) :: N1, N2

    N_mat_face = 0.0_DP

    N1 = (1.0_DP / 2.0_DP) * (1.0_DP - zeta_in)
    N2 = (1.0_DP / 2.0_DP) * (1.0_DP + zeta_in)

    N_mat_face(1, 2 * p1_in - 1) = N1
    N_mat_face(2, 2 * p1_in)     = N1

    N_mat_face(1, 2 * p2_in - 1) = N2
    N_mat_face(2, 2 * p2_in)     = N2

END FUNCTION MontaN_Facetri

! Calcula o comprimento, vetor tangente e vetor normal de uma face (aresta) do elemento

FUNCTION Facesnvtri(p1_in, p2_in, e_in, IJ_in, XY_in) RESULT(face_data)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p1_in, p2_in, e_in
    INTEGER, DIMENSION(:,:), INTENT(IN) :: IJ_in
    REAL(KIND=DP), DIMENSION(:,:), INTENT(IN) :: XY_in

    TYPE(face_data_type) :: face_data

    INTEGER :: n1_global, n2_global
    REAL(KIND=DP) :: x1, y1, x2, y2
    REAL(KIND=DP) :: deltaX, deltaY, L

    n1_global = IJ_in(e_in, p1_in + 1) 

    n2_global = IJ_in(e_in, p2_in + 1)

    x1 = XY_in(n1_global, 1)
    y1 = XY_in(n1_global, 2)
    x2 = XY_in(n2_global, 1)
    y2 = XY_in(n2_global, 2)

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
END FUNCTION Facesnvtri

! Calcula o vetor de forças equivalentes nodais F devido a um carregamento distribuído linear aplicado na face do elemento Q4, usando integração de Gauss

FUNCTION Ftblineartri(dJ_in, vetor_dir_local_in, p1_in, p2_in, carregamento_in, te_in) RESULT(F_vec)

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: dJ_in, carregamento_in, te_in
    REAL(KIND=DP), DIMENSION(2), INTENT(IN) :: vetor_dir_local_in
    INTEGER, INTENT(IN) :: p1_in, p2_in

    REAL(KIND=DP), DIMENSION(8) :: F_vec 

    REAL(KIND=DP), DIMENSION(2) :: pg
    REAL(KIND=DP) :: zeta_gauss
    REAL(KIND=DP), DIMENSION(2, 8) :: N_gauss

    INTEGER :: i

    pg(1) = -1.0_DP / SQRT(3.0_DP)
    pg(2) = 1.0_DP / SQRT(3.0_DP)

    F_vec = 0.0_DP

    DO i = 1, 2
        zeta_gauss = pg(i)

        N_gauss = MontaN_Facetri(zeta_gauss, p1_in, p2_in)

        F_vec = F_vec + carregamento_in * MATMUL(TRANSPOSE(N_gauss), vetor_dir_local_in) * dJ_in * te_in
    END DO
END FUNCTION Ftblineartri

! Calcula o vetor global de forças equivalentes Fglobal geradas pelo carregamentos distribuídos aplicados nas faces de elementos hexagonais

FUNCTION ForcatGlobaltri(nn_in, ESP_in, XY_in, IJ_in, np_in, P_in) RESULT(F_global)

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

    REAL(KIND=DP), DIMENSION(8) :: F_local_elementar 

    INTEGER, DIMENSION(8) :: global_dofs 

    TYPE(face_data_type) :: current_face_data

    INTEGER :: j, p1_local, p2_local

    F_global = 0.0_DP

    DO i = 1, np_in
        ele = NINT(P_in(i,1))
        face = NINT(P_in(i,2))
        dir = NINT(P_in(i,3))
        val = P_in(i,4)

        te = ESP_in(ele)

        p_nodes_face = FACES_TRI(face) 
        p1_local = p_nodes_face(1)
        p2_local = p_nodes_face(2)

        current_face_data = Facesnvtri(p1_local, p2_local, ele, IJ_in, XY_in)

        IF (dir == 1) THEN
            vetor_dir_final = current_face_data%n_vec
        ELSE
            vetor_dir_final = current_face_data%v_vec
        END IF

        F_local_elementar = Ftblineartri(current_face_data%dJ, vetor_dir_final, &
                                         p1_local, p2_local, val, te)

        global_dofs = Montaglstri(ele, IJ_in) 

        DO j = 1, 8 
            F_global(global_dofs(j)) = F_global(global_dofs(j)) + F_local_elementar(j)
        END DO
    END DO
END FUNCTION ForcatGlobaltri

END MODULE Contorno_Functions_Moduletri