MODULE Dados_Module
    USE precision_module, ONLY : DP
    IMPLICIT NONE
CONTAINS
    SUBROUTINE Dados(nn_out, XY_out, ne_out, IJ_out, MAT_out, ESP_out, &
                     nf_out, FC_out, np_out, P_out, na_out, AP_out, nfb_out, FB_out)
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: nn_out, ne_out, nf_out, np_out, na_out, nfb_out
        REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: XY_out, MAT_out, FC_out, AP_out
        REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: P_out
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: IJ_out
        REAL(KIND=DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ESP_out
	REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: FB_out
	INTEGER :: i_elem

        ! Número de nós
        nn_out = 10

        ! Número de elementos
        ne_out = 5

        ! Número de forças concentradas
        nf_out = 3

        ! Forças distribuidas
        np_out = 0

        ! Número de apoios
        na_out = 7

		! número de forças de corpo
		nfb_out = ne_out

        ALLOCATE(XY_out(nn_out, 2))
        ALLOCATE(IJ_out(ne_out, 7))
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

         ! Coordenadas XY nó 1
         XY_out(1, 1) = 0.0_DP
         XY_out(1, 2) = 0.0_DP
		! Coordenada XY nó 2...		
         XY_out(2, 1) = 0.25_DP
         XY_out(2, 2) = 0.0_DP
		 
         XY_out(3, 1) = 0.75_DP
         XY_out(3, 2) = 0.0_DP
		 
         XY_out(4, 1) = 1.0_DP
         XY_out(4, 2) = 0.0_DP
		 
         XY_out(5, 1) = 0.0_DP
         XY_out(5, 2) = 0.05_DP
		 
         XY_out(6, 1) = 1.0_DP
         XY_out(6, 2) = 0.05_DP
		 
         XY_out(7, 1) = 0.0_DP
         XY_out(7, 2) = 0.1_DP
		 
         XY_out(8, 1) = 0.25_DP
         XY_out(8, 2) = 0.1_DP
		 
         XY_out(9, 1) = 0.75_DP
         XY_out(9, 2) = 0.1_DP
		 
         XY_out(10, 1) = 1.0_DP
         XY_out(10, 2) = 0.1_DP
		
		! Elemento 1, (Tipo 2 == Q4), (Tipo = Primeira coluna)
         IJ_out(1, 1) = 2
         IJ_out(1, 2) = 1
         IJ_out(1, 3) = 2
         IJ_out(1, 4) = 5
         IJ_out(1, 5) = 5

         ! Elemento 2 (6 nós: 3, 4, 5, 6, 21, 20), (tipo 0 == Hexagonal (tipo = primeira coluna)
         IJ_out(2, 1) = 0
         IJ_out(2, 2) = 2
         IJ_out(2, 3) = 3
         IJ_out(2, 4) = 6
         IJ_out(2, 5) = 9
         IJ_out(2, 6) = 8
         IJ_out(2, 7) = 5

         IJ_out(3, 1) = 2
         IJ_out(3, 2) = 3
         IJ_out(3, 3) = 4
         IJ_out(3, 4) = 6
         IJ_out(3, 5) = 6

         IJ_out(4, 1) = 2
         IJ_out(4, 2) = 5
         IJ_out(4, 3) = 8
         IJ_out(4, 4) = 7
         IJ_out(4, 5) = 7

         IJ_out(5, 1) = 2
         IJ_out(5, 2) = 9
         IJ_out(5, 3) = 6
         IJ_out(5, 4) = 10
         IJ_out(5, 5) = 10
		 
		 ! Material do elemento (E, nu)
         ! Todos os 25 elementos tem o mesmo material
         DO i_elem = 1, 5
             MAT_out(i_elem, 1) = 1.0E+9_DP
             MAT_out(i_elem, 2) = 0.3_DP
         END DO
		 
		 ! Espessura do elemento
         ! Todos os 25 elementos têm a mesma espessura
         DO i_elem = 1, 5
             ESP_out(i_elem) = 0.01_DP
         END DO
		 
		 ! Forças Concentradas (Nó, Direção, Magnitude)
         ! Carregamento Uniformemente Distribuído de 1000 N Total sobre os 3 nós
         FC_out(1, 1) = 10.0_DP
         FC_out(1, 2) = 1.0_DP
         FC_out(1, 3) = 0.25_DP

         FC_out(2, 1) = 6.0_DP
         FC_out(2, 2) = 1.0_DP
         FC_out(2, 3) = 0.5_DP

         FC_out(3, 1) = 4.0_DP
         FC_out(3, 2) = 1.0_DP
         FC_out(3, 3) = 0.25_DP
		 
		 !P_out(1, 1) = 24.0_DP
		 !P_out(1, 2) = 2.0_DP
		 !P_out(1, 3) = 2.0_DP
		 !P_out(1, 4) = 500.0_DP
		 
		 !P_out(2, 1) = 22.0_DP
		 !P_out(2, 2) = 1.0_DP
		 !P_out(2, 3) = 2.0_DP
		 !P_out(2, 4) = 500.0_DP
		
		! Apoios
         AP_out(1, 1) = 1.0_DP   ! Nó 1
         AP_out(1, 2) = 1.0_DP   ! Direção X
         AP_out(1, 3) = 0.0_DP   ! Valor da restrição (deslocamento zero)

         AP_out(2, 1) = 5.0_DP
         AP_out(2, 2) = 1.0_DP
         AP_out(2, 3) = 0.0_DP

         AP_out(3, 1) = 7.0_DP
         AP_out(3, 2) = 1.0_DP
         AP_out(3, 3) = 0.0_DP

         AP_out(4, 1) = 1.0_DP
         AP_out(4, 2) = 2.0_DP
         AP_out(4, 3) = 0.0_DP

         AP_out(5, 1) = 2.0_DP
         AP_out(5, 2) = 2.0_DP
         AP_out(5, 3) = 0.0_DP

         AP_out(6, 1) = 3.0_DP
         AP_out(6, 2) = 2.0_DP
         AP_out(6, 3) = 0.0_DP

         AP_out(7, 1) = 4.0_DP
         AP_out(7, 2) = 2.0_DP
         AP_out(7, 3) = 0.0_DP
		
        ! Indice do elemento, direção (1=X, 2=Y), valor da força de corpo
        DO i_elem = 1, ne_out
            FB_out(i_elem, 1) = REAL(i_elem, KIND=DP) ! ID do elemento
            FB_out(i_elem, 2) = 2.0_DP                 ! Direção Y
            FB_out(i_elem, 3) = 0.0_DP
        END DO

    END SUBROUTINE Dados
END MODULE Dados_Module
