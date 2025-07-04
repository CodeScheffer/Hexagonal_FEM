MODULE precision_module
    IMPLICIT NONE
    ! Define um parâmetro de dupla precisão
    INTEGER, PARAMETER :: DP = KIND(0.0D0) ! KIND(0.0D0) garante precisão dupla
	REAL(KIND=DP), PARAMETER :: PI = ACOS(-1.0_DP)
END MODULE precision_module