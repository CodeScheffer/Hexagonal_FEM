MODULE Elementotri_Module
    USE precision_module, ONLY : DP
    USE Material_Module, ONLY : MontaCEPT
    IMPLICIT NONE
CONTAINS

! Derivadas dass funções de interpolação do Q4

FUNCTION MontadNtri(r, s) RESULT(dNrs)
        REAL(KIND=DP) :: r, s
        REAL(KIND=DP) :: dNrs(2,4)

        dNrs(1,1) = 0.25_DP * (s - 1.0_DP)
        dNrs(1,2) = 0.25_DP * (1.0_DP - s)
        dNrs(1,3) = 0.25_DP * (1.0_DP + s)
        dNrs(1,4) = -0.25_DP * (1.0_DP + s)

        dNrs(2,1) = 0.25_DP * (r - 1.0_DP)
        dNrs(2,2) = -0.25_DP * (1.0_DP + r)
        dNrs(2,3) = 0.25_DP * (1.0_DP + r)
        dNrs(2,4) = 0.25_DP * (1.0_DP - r)
END FUNCTION MontadNtri

! Monta a matriz jacobiana

FUNCTION MontaJtri(dNrs, X, Y) RESULT(J)
        REAL(KIND=DP), INTENT(IN) :: dNrs(2,4), X(4), Y(4)
        REAL(KIND=DP) :: J(2,2)
        INTEGER :: i

        J = 0.0_DP
        DO i = 1, 4
            J(1,1) = J(1,1) + dNrs(1,i) * X(i)
            J(1,2) = J(1,2) + dNrs(1,i) * Y(i)
            J(2,1) = J(2,1) + dNrs(2,i) * X(i)
            J(2,2) = J(2,2) + dNrs(2,i) * Y(i)
        END DO
END FUNCTION MontaJtri

! Posições dos nós do elemento Q4

SUBROUTINE MontaXYtri(e, IJ_global, XY, X, Y)
        INTEGER, INTENT(IN) :: e
        INTEGER, INTENT(IN) :: IJ_global(:,:)
        REAL(KIND=DP), INTENT(IN) :: XY(:,:)
        REAL(KIND=DP), INTENT(OUT) :: X(4), Y(4)
        INTEGER :: n, no

        DO n = 1, 4
            no = IJ_global(e, n + 1)
            X(n) = XY(no,1)
            Y(n) = XY(no,2)
        END DO
END SUBROUTINE MontaXYtri

! Corrige sitema de coordenadas

FUNCTION CorrigedNtri(dNrs, J) RESULT(dNxy)
        REAL(KIND=DP), INTENT(IN) :: dNrs(2,4), J(2,2)
        REAL(KIND=DP) :: dNxy(2,4), invJ(2,2), detJ

        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        IF (ABS(detJ) < 1.0E-12_DP) THEN
            WRITE(*,*) 'ERRO: Determinante Jacobiano muito pequeno ou zero em CorrigedNtri!'
            STOP 'Jacobiano singular'
        END IF

        invJ(1,1) = J(2,2)/detJ
        invJ(1,2) = -J(1,2)/detJ
        invJ(2,1) = -J(2,1)/detJ
        invJ(2,2) = J(1,1)/detJ

        dNxy = MATMUL(invJ, dNrs)
END FUNCTION CorrigedNtri

! Monta a matriz B

FUNCTION MontaBtri(dNxy) RESULT(B)
        REAL(KIND=DP), INTENT(IN) :: dNxy(2,4)
        REAL(KIND=DP) :: B(3,8)
        INTEGER :: i, c

        B = 0.0_DP
        c = 1
        DO i = 1, 4
            B(1,c) = dNxy(1,i)
            B(2,c+1) = dNxy(2,i)
            B(3,c) = dNxy(2,i)
            B(3,c+1) = dNxy(1,i)
            c = c + 2
        END DO
END FUNCTION MontaBtri

! Graus de liberdade

FUNCTION Montaglstri(e, IJ) RESULT(gls)
        INTEGER, INTENT(IN) :: e
        INTEGER, INTENT(IN) :: IJ(:,:)
        INTEGER :: gls(8), i, j, c, no

        c = 1
        DO i = 1, 4
            no = IJ(e,i+1)
            DO j = 1, 2
                gls(c) = 2*(no - 1) + j
                c = c + 1
            END DO
        END DO
END FUNCTION Montaglstri

! Monta a matriz de rigidez do elemento Q4

SUBROUTINE MontaKetri(ele, IJ, X_coords, Y_coords, E_val, nu_val, te_val, Ke_matrix)

        USE Material_Module, ONLY : MontaCEPT
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ele
        INTEGER, INTENT(IN) :: IJ(:,:)
        REAL(KIND=DP), INTENT(IN) :: X_coords(4), Y_coords(4), E_val, nu_val, te_val
        REAL(KIND=DP), INTENT(OUT) :: Ke_matrix(8,8)

        INTEGER, PARAMETER :: NUM_GAUSS_PTS_Q4 = 4
        REAL(KIND=DP), PARAMETER :: GPS_COORD_Q4 = 1.0_DP / SQRT(3.0_DP)

        REAL(KIND=DP), DIMENSION(NUM_GAUSS_PTS_Q4) :: r_pts, s_pts
        REAL(KIND=DP), DIMENSION(NUM_GAUSS_PTS_Q4) :: w_pts_q4

        REAL(KIND=DP) :: dNrs_local(2,4), J_local(2,2), dNxy_local(2,4)
        REAL(KIND=DP) :: B_local(3,8), D_local(3,3)
        REAL(KIND=DP) :: detJ_local
        INTEGER :: i_pt

        Ke_matrix = 0.0_DP

        IF (IJ(ele, 1) == 2) THEN
		
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

                dNrs_local = MontadNtri(r_pts(i_pt), s_pts(i_pt))

                J_local = MontaJtri(dNrs_local, X_coords, Y_coords)

                detJ_local = J_local(1,1)*J_local(2,2) - J_local(1,2)*J_local(2,1)

				 WRITE(*,*) 'Elemento: ', ele, ', Ponto de Gauss: ', i_pt, ', detJ = ', detJ_local

                IF (ABS(detJ_local) < 1.0E-12_DP) THEN
                    WRITE(*,*) 'ERRO em MontaKetri: Determinante Jacobiano muito pequeno ou zero no ponto de Gauss (', r_pts(i_pt), ',', s_pts(i_pt), ')'
                    STOP 'Jacobiano singular em MontaKetri'
                END IF

                dNxy_local = CorrigedNtri(dNrs_local, J_local)

                B_local = MontaBtri(dNxy_local)

                D_local = MontaCEPT(E_val, nu_val)

                Ke_matrix = Ke_matrix + MATMUL(TRANSPOSE(B_local), MATMUL(D_local, B_local)) * &
                                        detJ_local * w_pts_q4(i_pt) * te_val

            END DO
        END IF

    END SUBROUTINE MontaKetri

END MODULE Elementotri_Module
