MODULE Dados_Module
    USE precision_module, ONLY : DP
    IMPLICIT NONE
CONTAINS
    SUBROUTINE Dados(nn_out, XY_out, ne_out, IJ_out, MAT_out, ESP_out, &
                     nf_out, FC_out, np_out, P_out, na_out, AP_out, nfb_out, FB_out)
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: nn_out, ne_out, nf_out, np_out, na_out, nfb_out
        REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: XY_out, MAT_out, FC_out, P_out, AP_out
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: IJ_out
        REAL(KIND=DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ESP_out
		REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: FB_out
		INTEGER :: i_elem


        ! Número de nós
        nn_out = 6

        ! Número de elementos
        ne_out = 4

        ! Número de forças concentradas
        nf_out = 2

		! Forcas distribuida
        np_out = 0
		
		! Npumero de Apoios
        na_out = 5
		
		! Forcas de corpo
		nfb_out = ne_out

        ALLOCATE(XY_out(nn_out, 2))
        ALLOCATE(IJ_out(ne_out, 4))
        ALLOCATE(MAT_out(ne_out, 2))
        ALLOCATE(ESP_out(ne_out))

        IF (nf_out > 0) THEN
            ALLOCATE(FC_out(nf_out, 3))
        ELSE
            ALLOCATE(FC_out(0,0))
        END IF

        IF (np_out > 0) THEN
            ALLOCATE(P_out(np_out, 4))
        ELSE
            ALLOCATE(P_out(0,0))
        END IF

        IF (na_out > 0) THEN
            ALLOCATE(AP_out(na_out, 3))
        ELSE
            ALLOCATE(AP_out(0,0))
        END IF
		
		IF (nfb_out > 0) THEN
            ALLOCATE(FB_out(nfb_out, 3))
        ELSE
            ALLOCATE(FB_out(0,0))
        END IF
		
        ! Coordenadas XY
        XY_out(1, 1) = 0.0_DP; XY_out(1, 2) = 0.0_DP
        XY_out(2, 1) = 0.5_DP; XY_out(2, 2) = 0.0_DP
        XY_out(3, 1) = 1.0_DP; XY_out(3, 2) = 0.0_DP
        XY_out(4, 1) = 0.0_DP; XY_out(4, 2) = 0.1_DP
        XY_out(5, 1) = 0.5_DP; XY_out(5, 2) = 0.1_DP
        XY_out(6, 1) = 1.0_DP; XY_out(6, 2) = 0.1_DP

        ! Conectividades IJ
        IJ_out(1, 1) = 2; IJ_out(1, 2) = 1; IJ_out(1, 3) = 2; IJ_out(1, 4) = 5; IJ_out(1, 5) = 5
        IJ_out(2, 1) = 2; IJ_out(2, 2) = 2; IJ_out(2, 3) = 3; IJ_out(2, 4) = 6; IJ_out(2, 5) = 6
        IJ_out(3, 1) = 2; IJ_out(3, 2) = 1; IJ_out(3, 3) = 5; IJ_out(3, 4) = 4; IJ_out(3, 5) = 4
        IJ_out(4, 1) = 2; IJ_out(4, 2) = 2; IJ_out(4, 3) = 6; IJ_out(4, 4) = 5; IJ_out(4, 5) = 5

        ! Material MAT
        MAT_out(1, 1) = 1.0E+9_DP; MAT_out(1, 2) = 0.3_DP
        MAT_out(2, 1) = 1.0E+9_DP; MAT_out(2, 2) = 0.3_DP
        MAT_out(3, 1) = 1.0E+9_DP; MAT_out(3, 2) = 0.3_DP
        MAT_out(4, 1) = 1.0E+9_DP; MAT_out(4, 2) = 0.3_DP

        ! Espessuras ESP
        ESP_out(1) = 0.01_DP
        ESP_out(2) = 0.01_DP
        ESP_out(3) = 0.01_DP
        ESP_out(4) = 0.01_DP

        ! Nó 3 (força em X), Nó 6 (força em X)
        FC_out(1, 1) = 3.0_DP; FC_out(1, 2) = 1.0_DP; FC_out(1, 3) = 0.5_DP ! Nó 3, direcao 1 (X), valor 1N
        FC_out(2, 1) = 6.0_DP; FC_out(2, 2) = 1.0_DP; FC_out(2, 3) = 0.5_DP ! Nó 6, direcao 1 (X), valor 1N

        ! Apoios AP
        ! Nó 1: Fixo em X e Y (Dx=0, Dy=0)
        AP_out(1, 1) = 1.0_DP; AP_out(1, 2) = 1.0_DP; AP_out(1, 3) = 0.0_DP ! Nó 1, direcao 1 (X), valor 0
        AP_out(2, 1) = 1.0_DP; AP_out(2, 2) = 2.0_DP; AP_out(2, 3) = 0.0_DP ! Nó 1, direcao 2 (Y), valor 0
        ! Nó 2: Rolamento (Fixo em Y)
        AP_out(3, 1) = 2.0_DP; AP_out(3, 2) = 2.0_DP; AP_out(3, 3) = 0.0_DP ! Nó 2, direcao 2 (Y), valor 0
        ! Nó 3: Rolamento (Fixo em Y)
        AP_out(4, 1) = 3.0_DP; AP_out(4, 2) = 2.0_DP; AP_out(4, 3) = 0.0_DP ! Nó 3, direcao 2 (Y), valor 0
		AP_out(5, 1) = 4.0_DP; AP_out(5, 2) = 1.0_DP; AP_out(5, 3) = 0.0_DP ! Nó 4, direcao 1 (X), valor 0

        DO i_elem = 1, ne_out
            FB_out(i_elem, 1) = REAL(i_elem, KIND=DP) ! ID do elemento
            FB_out(i_elem, 2) = 2.0_DP                 ! Direção Y
            FB_out(i_elem, 3) = -9.81_DP
        END DO
		
    END SUBROUTINE Dados
END MODULE Dados_Module