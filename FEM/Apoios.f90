MODULE Apoios_Module
    USE precision_module, ONLY : DP
    IMPLICIT NONE
CONTAINS

! Aplica as condições de contorno essenciais

    SUBROUTINE AplicaCCH(nn, na, AP, K, F)
        INTEGER, INTENT(IN) :: nn, na
        REAL(KIND=DP), INTENT(INOUT) :: K(2*nn, 2*nn), F(2*nn)
        REAL(KIND=DP), INTENT(IN) :: AP(na,3)
        INTEGER :: l, no, gll, glg, i

        DO l = 1, na
            no = NINT(AP(l,1))
            gll = NINT(AP(l,2))
            glg = 2*(no-1) + gll

            DO i = 1, 2*nn
                K(glg,i) = 0.0_DP
                K(i,glg) = 0.0_DP
            END DO

            K(glg,glg) = 1.0_DP
            F(glg) = 0.0_DP
        END DO
    END SUBROUTINE AplicaCCH
END MODULE Apoios_Module