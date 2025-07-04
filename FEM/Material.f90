MODULE Material_Module
    USE precision_module, ONLY : DP
    IMPLICIT NONE
CONTAINS

! Calcula a matriz constitutiva para tensões planas 2D para material isotrópico

    FUNCTION MontaCEPT(E, nu) RESULT(Cv)
        IMPLICIT NONE

        REAL(KIND=DP), INTENT(IN) :: E, nu

        REAL(KIND=DP) :: c
        REAL(KIND=DP), DIMENSION(3, 3) :: Cv

        c = E / (1.0_DP - (nu**2))

        Cv(1, 1) = c
        Cv(1, 2) = c * nu
        Cv(1, 3) = 0.0_DP

        Cv(2, 1) = c * nu
        Cv(2, 2) = c
        Cv(2, 3) = 0.0_DP

        Cv(3, 1) = 0.0_DP
        Cv(3, 2) = 0.0_DP
        Cv(3, 3) = c * ((1.0_DP - nu) / 2.0_DP)

    END FUNCTION MontaCEPT
END MODULE Material_Module