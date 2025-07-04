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
        nn_out = 63

        ! Número de elementos
        ne_out = 25

        ! Número de forças concentradas
        nf_out = 3

        ! Forças distribuidas
        np_out = 0

        ! Número de apoios
        na_out = 8

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

         ! Coordenadas XY
         XY_out(1, 1) = 0.0_DP
         XY_out(1, 2) = 0.0_DP
         XY_out(2, 1) = 0.58_DP
         XY_out(2, 2) = 0.0_DP
         XY_out(3, 1) = 1.73_DP
         XY_out(3, 2) = 0.0_DP
         XY_out(4, 1) = 2.39_DP
         XY_out(4, 2) = 0.0_DP
         XY_out(5, 1) = 3.28_DP
         XY_out(5, 2) = 0.0_DP
         XY_out(6, 1) = 4.05_DP
         XY_out(6, 2) = 0.0_DP
         XY_out(7, 1) = 5.2_DP
         XY_out(7, 2) = 0.0_DP
         XY_out(8, 1) = 5.86_DP
         XY_out(8, 2) = 0.0_DP
         XY_out(9, 1) = 6.75_DP
         XY_out(9, 2) = 0.0_DP
         XY_out(10, 1) = 7.52_DP
         XY_out(10, 2) = 0.0_DP
         XY_out(11, 1) = 8.67_DP
         XY_out(11, 2) = 0.0_DP
         XY_out(12, 1) = 9.33_DP
         XY_out(12, 2) = 0.0_DP
         XY_out(13, 1) = 10.22_DP
         XY_out(13, 2) = 0.0_DP
         XY_out(14, 1) = 10.91_DP
         XY_out(14, 2) = 0.0_DP
         XY_out(15, 1) = 12.06_DP
         XY_out(15, 2) = 0.0_DP
         XY_out(16, 1) = 12.72_DP
         XY_out(16, 2) = 0.0_DP
         XY_out(17, 1) = 13.61_DP
         XY_out(17, 2) = 0.0_DP
         XY_out(18, 1) = 14.3_DP
         XY_out(18, 2) = 0.0_DP
         XY_out(19, 1) = 0.0_DP
         XY_out(19, 2) = 1.0_DP
         XY_out(20, 1) = 2.31_DP
         XY_out(20, 2) = 1.0_DP
         XY_out(21, 1) = 3.46_DP
         XY_out(21, 2) = 1.0_DP
         XY_out(22, 1) = 5.77_DP
         XY_out(22, 2) = 1.0_DP
         XY_out(23, 1) = 6.92_DP
         XY_out(23, 2) = 1.0_DP
         XY_out(24, 1) = 9.23_DP
         XY_out(24, 2) = 1.0_DP
         XY_out(25, 1) = 10.31_DP
         XY_out(25, 2) = 1.0_DP
         XY_out(26, 1) = 12.62_DP
         XY_out(26, 2) = 1.0_DP
         XY_out(27, 1) = 13.7_DP
         XY_out(27, 2) = 1.0_DP
         XY_out(28, 1) = 0.58_DP
         XY_out(28, 2) = 2.0_DP
         XY_out(29, 1) = 1.73_DP
         XY_out(29, 2) = 2.0_DP
         XY_out(30, 1) = 4.04_DP
         XY_out(30, 2) = 2.0_DP
         XY_out(31, 1) = 5.19_DP
         XY_out(31, 2) = 2.0_DP
         XY_out(32, 1) = 7.5_DP
         XY_out(32, 2) = 2.0_DP
         XY_out(33, 1) = 8.65_DP
         XY_out(33, 2) = 2.0_DP
         XY_out(34, 1) = 10.89_DP
         XY_out(34, 2) = 2.0_DP
         XY_out(35, 1) = 12.04_DP
         XY_out(35, 2) = 2.0_DP
         XY_out(36, 1) = 14.28_DP
         XY_out(36, 2) = 2.0_DP
         XY_out(37, 1) = 0.0_DP
         XY_out(37, 2) = 3.0_DP
         XY_out(38, 1) = 2.31_DP
         XY_out(38, 2) = 3.0_DP
         XY_out(39, 1) = 3.46_DP
         XY_out(39, 2) = 3.0_DP
         XY_out(40, 1) = 5.77_DP
         XY_out(40, 2) = 3.0_DP
         XY_out(41, 1) = 6.92_DP
         XY_out(41, 2) = 3.0_DP
         XY_out(42, 1) = 9.23_DP
         XY_out(42, 2) = 3.0_DP
         XY_out(43, 1) = 10.31_DP
         XY_out(43, 2) = 3.0_DP
         XY_out(44, 1) = 12.62_DP
         XY_out(44, 2) = 3.0_DP
         XY_out(45, 1) = 13.7_DP
         XY_out(45, 2) = 3.0_DP
         XY_out(46, 1) = 0.0_DP
         XY_out(46, 2) = 4.0_DP
         XY_out(47, 1) = 0.58_DP
         XY_out(47, 2) = 4.0_DP
         XY_out(48, 1) = 1.73_DP
         XY_out(48, 2) = 4.0_DP
         XY_out(49, 1) = 2.39_DP
         XY_out(49, 2) = 4.0_DP
         XY_out(50, 1) = 3.28_DP
         XY_out(50, 2) = 4.0_DP
         XY_out(51, 1) = 4.05_DP
         XY_out(51, 2) = 4.0_DP
         XY_out(52, 1) = 5.2_DP
         XY_out(52, 2) = 4.0_DP
         XY_out(53, 1) = 5.86_DP
         XY_out(53, 2) = 4.0_DP
         XY_out(54, 1) = 6.75_DP
         XY_out(54, 2) = 4.0_DP
         XY_out(55, 1) = 7.52_DP
         XY_out(55, 2) = 4.0_DP
         XY_out(56, 1) = 8.67_DP
         XY_out(56, 2) = 4.0_DP
         XY_out(57, 1) = 9.33_DP
         XY_out(57, 2) = 4.0_DP
         XY_out(58, 1) = 10.22_DP
         XY_out(58, 2) = 4.0_DP
         XY_out(59, 1) = 10.91_DP
         XY_out(59, 2) = 4.0_DP
         XY_out(60, 1) = 12.06_DP
         XY_out(60, 2) = 4.0_DP
         XY_out(61, 1) = 12.72_DP
         XY_out(61, 2) = 4.0_DP
         XY_out(62, 1) = 13.61_DP
         XY_out(62, 2) = 4.0_DP
		 XY_out(63, 1) = 14.3_DP
         XY_out(63, 2) = 4.0_DP
		
		 ! Elemento 1
         IJ_out(1, 1) = 0   ! tipo Hexagonal
         IJ_out(1, 2) = 19
         IJ_out(1, 3) = 2
         IJ_out(1, 4) = 3
         IJ_out(1, 5) = 20
         IJ_out(1, 6) = 29
         IJ_out(1, 7) = 28

         ! Elemento 2
         IJ_out(2, 1) = 0
         IJ_out(2, 2) = 3
         IJ_out(2, 3) = 4
         IJ_out(2, 4) = 5
         IJ_out(2, 5) = 6
         IJ_out(2, 6) = 21
         IJ_out(2, 7) = 20

         ! Elemento 3
         IJ_out(3, 1) = 0
         IJ_out(3, 2) = 21
         IJ_out(3, 3) = 6
         IJ_out(3, 4) = 7
         IJ_out(3, 5) = 22
         IJ_out(3, 6) = 31
         IJ_out(3, 7) = 30

         ! Elemento 4
         IJ_out(4, 1) = 0
         IJ_out(4, 2) = 7
         IJ_out(4, 3) = 8
         IJ_out(4, 4) = 9
         IJ_out(4, 5) = 10
         IJ_out(4, 6) = 23
         IJ_out(4, 7) = 22

         ! Elemento 5
         IJ_out(5, 1) = 0
         IJ_out(5, 2) = 23
         IJ_out(5, 3) = 10
         IJ_out(5, 4) = 11
         IJ_out(5, 5) = 24
         IJ_out(5, 6) = 33
         IJ_out(5, 7) = 32

         ! Elemento 6
         IJ_out(6, 1) = 0
         IJ_out(6, 2) = 11
         IJ_out(6, 3) = 12
         IJ_out(6, 4) = 13
         IJ_out(6, 5) = 14
         IJ_out(6, 6) = 25
         IJ_out(6, 7) = 24

         ! Elemento 7
         IJ_out(7, 1) = 0
         IJ_out(7, 2) = 25
         IJ_out(7, 3) = 14
         IJ_out(7, 4) = 15
         IJ_out(7, 5) = 26
         IJ_out(7, 6) = 35
         IJ_out(7, 7) = 34

         ! Elemento 8
         IJ_out(8, 1) = 0
         IJ_out(8, 2) = 15
         IJ_out(8, 3) = 16
         IJ_out(8, 4) = 17
         IJ_out(8, 5) = 18
         IJ_out(8, 6) = 27
         IJ_out(8, 7) = 26

         ! Elemento 9
         IJ_out(9, 1) = 0
         IJ_out(9, 2) = 37
         IJ_out(9, 3) = 28
         IJ_out(9, 4) = 29
         IJ_out(9, 5) = 38
         IJ_out(9, 6) = 48
         IJ_out(9, 7) = 47

         ! Elemento 10
         IJ_out(10, 1) = 0
         IJ_out(10, 2) = 29
         IJ_out(10, 3) = 20
         IJ_out(10, 4) = 21
         IJ_out(10, 5) = 30
         IJ_out(10, 6) = 39
         IJ_out(10, 7) = 38

         ! Elemento 11
         IJ_out(11, 1) = 0
         IJ_out(11, 2) = 39
         IJ_out(11, 3) = 30
         IJ_out(11, 4) = 31
         IJ_out(11, 5) = 40
         IJ_out(11, 6) = 52
         IJ_out(11, 7) = 51

         ! Elemento 12
         IJ_out(12, 1) = 0
         IJ_out(12, 2) = 31
         IJ_out(12, 3) = 22
         IJ_out(12, 4) = 23
         IJ_out(12, 5) = 32
         IJ_out(12, 6) = 41
         IJ_out(12, 7) = 40

         ! Elemento 13
         IJ_out(13, 1) = 0
         IJ_out(13, 2) = 41
         IJ_out(13, 3) = 32
         IJ_out(13, 4) = 33
         IJ_out(13, 5) = 42
         IJ_out(13, 6) = 56
         IJ_out(13, 7) = 55

         ! Elemento 14
         IJ_out(14, 1) = 0
         IJ_out(14, 2) = 33
         IJ_out(14, 3) = 24
         IJ_out(14, 4) = 25
         IJ_out(14, 5) = 34
         IJ_out(14, 6) = 43
         IJ_out(14, 7) = 42

         ! Elemento 15
         IJ_out(15, 1) = 0
         IJ_out(15, 2) = 43
         IJ_out(15, 3) = 34
         IJ_out(15, 4) = 35
         IJ_out(15, 5) = 44
         IJ_out(15, 6) = 60
         IJ_out(15, 7) = 59

         ! Elemento 16
         IJ_out(16, 1) = 0
         IJ_out(16, 2) = 35
         IJ_out(16, 3) = 26
         IJ_out(16, 4) = 27
         IJ_out(16, 5) = 36
         IJ_out(16, 6) = 45
         IJ_out(16, 7) = 44

         ! Elemento 17
         IJ_out(17, 1) = 0
         IJ_out(17, 2) = 48
         IJ_out(17, 3) = 38
         IJ_out(17, 4) = 39
         IJ_out(17, 5) = 51
         IJ_out(17, 6) = 50
         IJ_out(17, 7) = 49

         ! Elemento 18
         IJ_out(18, 1) = 0
         IJ_out(18, 2) = 52
         IJ_out(18, 3) = 40
         IJ_out(18, 4) = 41
         IJ_out(18, 5) = 55
         IJ_out(18, 6) = 54
         IJ_out(18, 7) = 53

         ! Elemento 19
         IJ_out(19, 1) = 0
         IJ_out(19, 2) = 56
         IJ_out(19, 3) = 42
         IJ_out(19, 4) = 43
         IJ_out(19, 5) = 59
         IJ_out(19, 6) = 58
         IJ_out(19, 7) = 57

         ! Elemento 20
         IJ_out(20, 1) = 0
         IJ_out(20, 2) = 60
         IJ_out(20, 3) = 44
         IJ_out(20, 4) = 45
         IJ_out(20, 5) = 63
         IJ_out(20, 6) = 62
         IJ_out(20, 7) = 61

         ! Elemento 21
         IJ_out(21, 1) = 2 ! Triangular Degenerado (Q4)
         IJ_out(21, 2) = 1
         IJ_out(21, 3) = 2
         IJ_out(21, 4) = 19
         IJ_out(21, 5) = 19
         IJ_out(21, 6) = 0
         IJ_out(21, 7) = 0

         ! Elemento 22
         IJ_out(22, 1) = 2
         IJ_out(22, 2) = 18
         IJ_out(22, 3) = 36
         IJ_out(22, 4) = 27
         IJ_out(22, 5) = 27
         IJ_out(22, 6) = 0
         IJ_out(22, 7) = 0

         ! Elemento 23
         IJ_out(23, 1) = 2
         IJ_out(23, 2) = 19
         IJ_out(23, 3) = 28
         IJ_out(23, 4) = 37
         IJ_out(23, 5) = 37
         IJ_out(23, 6) = 0
         IJ_out(23, 7) = 0

         ! Elemento 24
         IJ_out(24, 1) = 2
         IJ_out(24, 2) = 45
         IJ_out(24, 3) = 36
         IJ_out(24, 4) = 63
         IJ_out(24, 5) = 63
         IJ_out(24, 6) = 0
         IJ_out(24, 7) = 0

         ! Elemento 25
         IJ_out(25, 1) = 2
         IJ_out(25, 2) = 37
         IJ_out(25, 3) = 47
         IJ_out(25, 4) = 46
         IJ_out(25, 5) = 46
         IJ_out(25, 6) = 0
         IJ_out(25, 7) = 0
		 
		 ! Material do elemento (E, nu)
         ! Todos os 25 elementos tem o mesmo material
         DO i_elem = 1, 25
             MAT_out(i_elem, 1) = 1.0E+9_DP
             MAT_out(i_elem, 2) = 0.0_DP
         END DO
		 
		 ! Espessura do elemento
         ! Todos os 25 elementos têm a mesma espessura
         DO i_elem = 1, 25
             ESP_out(i_elem) = 0.01_DP
         END DO
		 
		 ! Forças Concentradas (Nó, Direção, Magnitude)
         ! Carregamento Uniformemente Distribuído de 1000 N Total sobre os 3 nós
         FC_out(1, 1) = 63.0_DP
         FC_out(1, 2) = 1.0_DP
         FC_out(1, 3) = 20.0_DP

         FC_out(2, 1) = 36.0_DP
         FC_out(2, 2) = 1.0_DP
         FC_out(2, 3) = 40.0_DP

         FC_out(3, 1) = 18.0_DP
         FC_out(3, 2) = 1.0_DP
         FC_out(3, 3) = 20.0_DP
		 
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

         AP_out(2, 1) = 1.0_DP
         AP_out(2, 2) = 2.0_DP
         AP_out(2, 3) = 0.0_DP

         AP_out(3, 1) = 19.0_DP
         AP_out(3, 2) = 1.0_DP
         AP_out(3, 3) = 0.0_DP

         AP_out(4, 1) = 19.0_DP
         AP_out(4, 2) = 2.0_DP
         AP_out(4, 3) = 0.0_DP

         AP_out(5, 1) = 37.0_DP
         AP_out(5, 2) = 1.0_DP
         AP_out(5, 3) = 0.0_DP

         AP_out(6, 1) = 37.0_DP
         AP_out(6, 2) = 2.0_DP
         AP_out(6, 3) = 0.0_DP

         AP_out(7, 1) = 46.0_DP
         AP_out(7, 2) = 1.0_DP
         AP_out(7, 3) = 0.0_DP

         AP_out(8, 1) = 46.0_DP
         AP_out(8, 2) = 2.0_DP
         AP_out(8, 3) = 0.0_DP
		
        ! Indice do elemento, direção (1=X, 2=Y), valor da força de corpo
        DO i_elem = 1, ne_out
            FB_out(i_elem, 1) = REAL(i_elem, KIND=DP) ! ID do elemento
            FB_out(i_elem, 2) = 2.0_DP                 ! Direção Y
            FB_out(i_elem, 3) = -9.81_DP
        END DO

    END SUBROUTINE Dados
END MODULE Dados_Module
