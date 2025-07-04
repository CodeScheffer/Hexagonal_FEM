MODULE Elemento_Module
    USE precision_module, ONLY : DP
    USE Material_Module, ONLY : MontaCEPT
    IMPLICIT NONE
    
    REAL(KIND=DP), PARAMETER :: PI = ACOS(-1.0_DP)
    
CONTAINS

! Derivadas das funções de interpolação

SUBROUTINE MontadN(r, s, dNrs)
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: r, s
    REAL(KIND=DP), DIMENSION(2,6), INTENT(OUT) :: dNrs

   ! dN1/dr
          dNrs(1,1) = (1.0_DP/81.0_DP)*(2*sqrt(3.0_DP)*s - 3)*(2*sqrt(3.0_DP)*s+ 3)*(-r*(3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s + 3) + 9*(r + 1)*(r**2 + s**2 - 3))/(r**2 + s**2 - 3)**2
    ! dN2/dr
          dNrs(1,2) = (1.0_DP/162.0_DP)*(2*sqrt(3.0_DP)*s + 3)*(2*r*(-3*r +sqrt(3.0_DP)*s + 3)*(3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s+ 3) + 3*(r**2 + s**2 - 3)*(-(-3*r + sqrt(3.0_DP)*s + 3)*(3*r -sqrt(3.0_DP)*s + 3) - (-3*r + sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s + 3) + (3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s + 3)))/(r**2 + s**2 - 3)**2
    ! dN3/dr
          dNrs(1,3) = (1.0_DP/162.0_DP)*(2*sqrt(3.0_DP)*s + 3)*(-2*r*(-3*r + sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s - 3)*(3*r + sqrt(3.0_DP)*s+ 3) + 3*(r**2 + s**2 - 3)*((-3*r + sqrt(3.0_DP)*s + 3)*(3*r +sqrt(3.0_DP)*s - 3) + (-3*r + sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s + 3) - (3*r + sqrt(3.0_DP)*s - 3)*(3*r + sqrt(3.0_DP)*s + 3)))/(r**2 + s**2 - 3)**2
    ! dN4/dr
          dNrs(1,4) = (1.0_DP/81.0_DP)*(2*sqrt(3.0_DP)*s - 3)*(2*sqrt(3.0_DP)*s+ 3)*(r*(-3*r + sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s - 3) + 9*(r - 1)*(r**2 + s**2 - 3))/(r**2 + s**2 - 3)**2
    ! dN5/dr
          dNrs(1,5) = (1.0_DP/162.0_DP)*(2*sqrt(3.0_DP)*s - 3)*(2*r*(-3*r +sqrt(3.0_DP)*s + 3)*(3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s- 3) + 3*(r**2 + s**2 - 3)*(-(-3*r + sqrt(3.0_DP)*s + 3)*(3*r -sqrt(3.0_DP)*s + 3) - (-3*r + sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s - 3) + (3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s - 3)))/(r**2 + s**2 - 3)**2
    ! dN6/dr
          dNrs(1,6) = (1.0_DP/162.0_DP)*(2*sqrt(3.0_DP)*s - 3)*(2*r*(3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s - 3)*(3*r + sqrt(3.0_DP)*s + 3) - 3*(r**2 + s**2 - 3)*((3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s - 3) + (3*r - sqrt(3.0_DP)*s + 3)*(3*r + sqrt(3.0_DP)*s +3) + (3*r + sqrt(3.0_DP)*s - 3)*(3*r + sqrt(3.0_DP)*s + 3)))/(r**2+ s**2 - 3)**2
    ! dN1/ds
          dNrs(2,1) = (2.0_DP/9.0_DP)*s*(6*r**4 + 12*r**3 - 4*r**2*s**2 - 6*r**2 - 27*r - 2*s**4 + 12*s**2 - 18)/(r**4 + 2*r**2*s**2 - 6*r**2+ s**4 - 6*s**2 + 9)
    ! dN2/ds
          dNrs(2,2) = (1.0_DP/18.0_DP)*(6*sqrt(3.0_DP)*r**5 - 12*r**4*s + 3*sqrt(3.0_DP)*r**4 - 12*sqrt(3.0_DP)*r**3*s**2 - 48*r**3*s - 30*sqrt(3.0_DP)*r**3 + 8*r**2*s**3 + 6*sqrt(3.0_DP)*r**2*s**2 + 12*r**2*s- 18*sqrt(3.0_DP)*r**2 - 2*sqrt(3.0_DP)*r*s**4 + 30*sqrt(3.0_DP)*r*s**2 + 108*r*s + 36*sqrt(3.0_DP)*r + 4*s**5 + 3*sqrt(3.0_DP)*s**4 -24*s**3 - 18*sqrt(3.0_DP)*s**2 + 36*s + 27*sqrt(3.0_DP))/(r**4 + 2*r**2*s**2 - 6*r**2 + s**4 - 6*s**2 + 9)
    ! dN3/ds
          dNrs(2,3) = (1.0_DP/18.0_DP)*(-6*sqrt(3.0_DP)*r**5 - 12*r**4*s + 3*sqrt(3.0_DP)*r**4 + 12*sqrt(3.0_DP)*r**3*s**2 + 48*r**3*s + 30*sqrt(3.0_DP)*r**3 + 8*r**2*s**3 + 6*sqrt(3.0_DP)*r**2*s**2 + 12*r**2*s - 18*sqrt(3.0_DP)*r**2 + 2*sqrt(3.0_DP)*r*s**4 - 30*sqrt(3.0_DP)*r*s**2 - 108*r*s - 36*sqrt(3.0_DP)*r + 4*s**5 + 3*sqrt(3.0_DP)*s**4 -24*s**3 - 18*sqrt(3.0_DP)*s**2 + 36*s + 27*sqrt(3.0_DP))/(r**4 + 2*r**2*s**2 - 6*r**2 + s**4 - 6*s**2 + 9)
    ! dN4/ds
          dNrs(2,4) = (2.0_DP/9.0_DP)*s*(6*r**4 - 12*r**3 - 4*r**2*s**2 - 6*r**2 + 27*r - 2*s**4 + 12*s**2 - 18)/(r**4 + 2*r**2*s**2 - 6*r**2+ s**4 - 6*s**2 + 9)
    ! dN5/ds
          dNrs(2,5) = (1.0_DP/18.0_DP)*(6*sqrt(3.0_DP)*r**5 - 12*r**4*s - 3*sqrt(3.0_DP)*r**4 - 12*sqrt(3.0_DP)*r**3*s**2 + 48*r**3*s - 30*sqrt(3.0_DP)*r**3 + 8*r**2*s**3 - 6*sqrt(3.0_DP)*r**2*s**2 + 12*r**2*s+ 18*sqrt(3.0_DP)*r**2 - 2*sqrt(3.0_DP)*r*s**4 + 30*sqrt(3.0_DP)*r*s**2 - 108*r*s + 36*sqrt(3.0_DP)*r + 4*s**5 - 3*sqrt(3.0_DP)*s**4 -24*s**3 + 18*sqrt(3.0_DP)*s**2 + 36*s - 27*sqrt(3.0_DP))/(r**4 + 2*r**2*s**2 - 6*r**2 + s**4 - 6*s**2 + 9)
    ! dN6/ds
          dNrs(2,6) = (1.0_DP/18.0_DP)*(-6*sqrt(3.0_DP)*r**5 - 12*r**4*s - 3*sqrt(3.0_DP)*r**4 + 12*sqrt(3.0_DP)*r**3*s**2 - 48*r**3*s + 30*sqrt(3.0_DP)*r**3 + 8*r**2*s**3 - 6*sqrt(3.0_DP)*r**2*s**2 + 12*r**2*s+ 18*sqrt(3.0_DP)*r**2 + 2*sqrt(3.0_DP)*r*s**4 - 30*sqrt(3.0_DP)*r*s**2 + 108*r*s - 36*sqrt(3.0_DP)*r + 4*s**5 - 3*sqrt(3.0_DP)*s**4 -24*s**3 + 18*sqrt(3.0_DP)*s**2 + 36*s - 27*sqrt(3.0_DP))/(r**4 + 2*r**2*s**2 - 6*r**2 + s**4 - 6*s**2 + 9)

END SUBROUTINE MontadN

! Monta o Jacobiano

FUNCTION MontaJ(dNrs, X, Y) RESULT(J)

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: dNrs(2,6), X(6), Y(6)

    REAL(KIND=DP) :: J(2,2)

    INTEGER :: i

    J = 0.0_DP

    DO i = 1, 6
        J(1,1) = J(1,1) + dNrs(1,i) * X(i)
        J(1,2) = J(1,2) + dNrs(1,i) * Y(i)
        J(2,1) = J(2,1) + dNrs(2,i) * X(i)
        J(2,2) = J(2,2) + dNrs(2,i) * Y(i)
    END DO

END FUNCTION MontaJ

! Corrige para o sistema global X e Y

FUNCTION CorrigedN(dNrs, J) RESULT(dNxy)
        REAL(KIND=DP), INTENT(IN) :: dNrs(2,6), J(2,2)
        REAL(KIND=DP) :: dNxy(2,6), invJ(2,2), detJ

		! Determinante Jacobiano
        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        IF (ABS(detJ) < 1.0E-12_DP) THEN
            WRITE(*,*) 'ERRO: Determinante Jacobiano muito pequeno ou zero em CorrigedNtri!'
            STOP 'Jacobiano singular'
        END IF

		!Inverte a matriz do jacobiano
        invJ(1,1) = J(2,2)/detJ
        invJ(1,2) = -J(1,2)/detJ
        invJ(2,1) = -J(2,1)/detJ
        invJ(2,2) = J(1,1)/detJ

        dNxy = MATMUL(invJ, dNrs)
		
END FUNCTION CorrigedN

! Monta a matriz B, que contem as derivadas das funções de interpolação no sistema global de coordenada

FUNCTION MontaB(dNxy) RESULT(B)

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: dNxy(2,6)

    REAL(KIND=DP) :: B(3,12)

    INTEGER :: i, c

    B = 0.0_DP

    c = 1

    DO i = 1, 6
        B(1,c) = dNxy(1,i)

        B(2,c+1) = dNxy(2,i)

        B(3,c) = dNxy(2,i)

        B(3,c+1) = dNxy(1,i)

        c = c + 2
    END DO

END FUNCTION MontaB

! Testa se o somatorio das derivadas em alguns pontos é 0, devido ao tamanho das equações

SUBROUTINE TesteDerivadas()
    USE precision_module, ONLY: DP
    IMPLICIT NONE

    REAL(KIND=DP) :: r, s
    REAL(KIND=DP) :: dNrs(2,6)
    INTEGER :: i

    REAL(KIND=DP), DIMENSION(10,2) :: pontos_teste

    ! Definição dos pontos (r, s) para teste
    pontos_teste = RESHAPE( (/ &
        0.0_DP,  0.0_DP,  &
        0.5_DP,  0.0_DP,  &
       -0.5_DP,  0.0_DP,  &
        0.0_DP,  0.5_DP,  &
        0.0_DP, -0.5_DP,  &
        0.5_DP,  0.5_DP,  &
       -0.5_DP,  0.5_DP,  &
        0.5_DP, -0.5_DP,  &
       -0.5_DP, -0.5_DP,  &
        0.9_DP,  0.3_DP   /), (/10,2/))

    WRITE(*,*) '=============================================='
    WRITE(*,*) 'Teste das Derivadas das Funções de Forma'
    WRITE(*,*) '=============================================='

    DO i = 1, 10
        r = pontos_teste(i,1)
        s = pontos_teste(i,2)

        CALL MontadN(r, s, dNrs)

        WRITE(*,*) 'Ponto (r,s) = (', r, ',', s, ')'
        WRITE(*,'(A, F14.10)') ' Soma das derivadas em r: ', SUM(dNrs(1,:))
        WRITE(*,'(A, F14.10)') ' Soma das derivadas em s: ', SUM(dNrs(2,:))
        WRITE(*,*) '----------------------------------------------'
    END DO

END SUBROUTINE TesteDerivadas

! extrai as coordenada XY dos nós do elemento hexagonal

SUBROUTINE MontaXY(e, IJ, XY, X, Y)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: e
    INTEGER, INTENT(IN) :: IJ(:,:)
    REAL(KIND=DP), INTENT(IN) :: XY(:,:)
    REAL(KIND=DP), INTENT(OUT) :: X(6), Y(6)

    INTEGER :: n, no

    DO n = 1, 6
        no = IJ(e,n+1)
        X(n) = XY(no,1)
        Y(n) = XY(no,2)
    END DO
END SUBROUTINE MontaXY

! Monta os graus de liberdades de cada nó

FUNCTION Montagls(e, IJ) RESULT(gls)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: e
    INTEGER, INTENT(IN) :: IJ(:,:)

    INTEGER :: gls(12)

    INTEGER :: i, j, c, no

    c = 1

    DO i = 1, 6
        no = IJ(e,i+1)

        DO j = 1, 2
            gls(c) = 2*(no - 1) + j
            c = c + 1
        END DO
    END DO

END FUNCTION Montagls

! Monta a matriz de rigidez do elemento hexagonal

SUBROUTINE MontaKe(X_coords, Y_coords, E_val, nu_val, te_val, Ke_matrix)
    USE Material_Module, ONLY : MontaCEPT
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: X_coords(6), Y_coords(6)
    REAL(KIND=DP), INTENT(IN) :: E_val, nu_val, te_val
    REAL(KIND=DP), INTENT(OUT) :: Ke_matrix(12,12)

    ! Dados dos pontos de integração
    REAL(KIND=DP), PARAMETER :: r_gauss(7) = (/ 0.748331477_DP, 0.374165739_DP, -0.374165739_DP, -0.748331477_DP, -0.374165739_DP, 0.374165739_DP, 0.0_DP /)
    REAL(KIND=DP), PARAMETER :: s_gauss(7) = (/ 0.0_DP, 0.64807407_DP, 0.64807407_DP, 0.0_DP, -0.64807407_DP, -0.64807407_DP, 0.0_DP /)
    REAL(KIND=DP), PARAMETER :: w_gauss(7) = (/ 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.124007937_DP, 0.255952381_DP /)

    REAL(KIND=DP) :: r, s
    REAL(KIND=DP) :: dNrs(2,6), dNxy(2,6), J(2,2)
    REAL(KIND=DP) :: B(3,12), D(3,3)
    REAL(KIND=DP) :: area_hex
    INTEGER :: i, k

    Ke_matrix = 0.0_DP
    D = MontaCEPT(E_val, nu_val)

    area_hex = 0.0_DP
    DO k = 1, 6
        area_hex = area_hex + (X_coords(k) * Y_coords(MOD(k,6)+1) - X_coords(MOD(k,6)+1) * Y_coords(k))
    END DO
    area_hex = 0.5_DP * ABS(area_hex)

    ! Loop nos pontos de Gauss
    DO i = 1, 7
        r = r_gauss(i)
        s = s_gauss(i)

        CALL MontadN(r, s, dNrs)
        J = MontaJ(dNrs, X_coords, Y_coords)

        dNxy = CorrigedN(dNrs, J)
        B = MontaB(dNxy)

        Ke_matrix = Ke_matrix + MATMUL(TRANSPOSE(B), MATMUL(D, B)) * area_hex * w_gauss(i) * te_val
    END DO
END SUBROUTINE MontaKe



!!!!!!!!!!!!!!!!!!!!!!!!!!! Monte carlo

! SUBROUTINE MontaKe(X_coords, Y_coords, E_val, nu_val, te_val, Ke_matrix)
    ! USE Material_Module, ONLY : MontaCEPT
    ! IMPLICIT NONE

    ! ! Entradas
    ! REAL(KIND=DP), INTENT(IN) :: X_coords(6), Y_coords(6)
    ! REAL(KIND=DP), INTENT(IN) :: E_val, nu_val, te_val
    ! REAL(KIND=DP), INTENT(OUT) :: Ke_matrix(12,12)

    ! ! Configuração Monte Carlo
    ! INTEGER, PARAMETER :: NPTS = 100  ! <<< Número de pontos a sserem utilizados
    ! REAL(KIND=DP) :: r, s
    ! REAL(KIND=DP) :: dNrs_local(2,6), J_local(2,2), dNxy_local(2,6)
    ! REAL(KIND=DP) :: B_local(3,12), D_local(3,3)
    ! REAL(KIND=DP) :: detJ_local, peso
    ! INTEGER :: j, aceitos

    ! REAL(KIND=DP) :: area_hex

    ! ! Inicializa
    ! Ke_matrix = 0.0_DP
    ! D_local = MontaCEPT(E_val, nu_val)

    ! ! Área do hexágono regular
    ! area_hex = (3.0_DP * SQRT(3.0_DP)) / 2.0_DP

    ! aceitos = 0

    ! DO j = 1, NPTS
        ! ! Gera ponto aleatório no retângulo [-1,1]x[-sqrt(3),sqrt(3)]
        ! r = 2.0_DP * (RAND() - 0.5_DP)
        ! s = (2.0_DP * (RAND() - 0.5_DP)) * SQRT(3.0_DP)

        ! ! Teste se está dentro do hexágono regular centrado na origem
        ! IF ( (ABS(r) <= 1.0_DP) .AND. (ABS(s) <= SQRT(3.0_DP)) ) THEN
            ! IF (ABS(s) <= SQRT(3.0_DP) * (1.0_DP - ABS(r))) THEN
                ! aceitos = aceitos + 1

                ! CALL MontadN(r, s, dNrs_local)
                ! J_local = MontaJ(dNrs_local, X_coords, Y_coords)
                ! detJ_local = J_local(1,1)*J_local(2,2) - J_local(1,2)*J_local(2,1)

                ! IF (ABS(detJ_local) < 1.0E-12_DP) CYCLE

                ! dNxy_local = CorrigedN(dNrs_local, J_local)
                ! B_local = MontaB(dNxy_local)

                ! Ke_matrix = Ke_matrix + MATMUL(TRANSPOSE(B_local), MATMUL(D_local, B_local)) * detJ_local
            ! END IF
        ! END IF
    ! END DO

    ! ! Peso final (média vezes área)
    ! peso = (area_hex / REAL(aceitos, DP))

    ! Ke_matrix = Ke_matrix * peso * te_val

    ! WRITE(*,*) 'Monte Carlo: Pontos aceitos =', aceitos, ' de ', NPTS

! END SUBROUTINE MontaKe

END MODULE Elemento_Module