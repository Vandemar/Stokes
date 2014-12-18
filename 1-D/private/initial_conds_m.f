      MODULE INITIAL_CONDS_M
!
!     This module contains the global parameters for the implementation
!     of various finite difference methods. This includes global 
!     variables as well as routines to initialize/setup the initial
!     conditions array.  Currently this is only for 1-D implementations
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!=======================================================================
!
      TYPE INIT_CONDS
        REAL (KIND = KIND(0.D0)), ALLOCATABLE :: V0(:)
        INTEGER :: IC_TYPE
        CONTAINS
        PROCEDURE :: INIT_SQ_WAVE
        PROCEDURE :: INIT_SEMI_CIRCLE
        PROCEDURE :: INIT_GAUSS_PULSE
        PROCEDURE :: DELETE_ICS
      END TYPE INIT_CONDS
!
!=======================================================================
      CONTAINS
!=======================================================================
!
      SUBROUTINE INIT_SQ_WAVE (IC_OBJ, X, M, H)
!
        CLASS (INIT_CONDS) :: IC_OBJ
        INTEGER :: M
        REAL (KIND = KIND(0.D0)) :: H
!
!       declare local variables
!
        INTEGER :: I
        REAL (KIND = KIND(0.D0)) :: X(:)
!
!-----------------------------------------------------------------------
!
        DO I = 1, SIZE(X)
          IF ((X(I) .GE. 0.25D0) .AND. (X(I) .LE. 0.75D0)) THEN
            IC_OBJ%V0(I) = 1.D0
          ELSE
            IC_OBJ%V0(I) = 0.D0
          END IF
        END DO
!
      RETURN
      END SUBROUTINE INIT_SQ_WAVE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE INIT_SEMI_CIRCLE (IC_OBJ, X, M, H)
!
        CLASS (INIT_CONDS) :: IC_OBJ
        INTEGER :: M
        REAL (KIND = KIND(0.D0)) :: H
!
!       declare local variables
!
        INTEGER :: I
        REAL (KIND = KIND(0.D0)) :: X(:)
!
!-----------------------------------------------------------------------
!
        DO I = 1, SIZE(X)
          IC_OBJ%V0(I) = SQRT(0.25D0 - (X(I) - 0.5D0)**2)
        END DO
!
      RETURN
      END SUBROUTINE INIT_SEMI_CIRCLE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE INIT_GAUSS_PULSE (IC_OBJ, X, M, H)
!
        CLASS (INIT_CONDS) :: IC_OBJ
        INTEGER :: M
        REAL (KIND = KIND(0.D0)) :: H
!
!       declare local variables
!
        INTEGER :: I, J
        REAL (KIND = KIND(0.D0)) :: X(:)
!
!-----------------------------------------------------------------------
!
        IC_OBJ%V0 = 0.D0
        DO I = 1, SIZE(X)
          IC_OBJ%V0(I) = EXP(-256.D0 * (X(I) - 0.5D0)**2)
          IF (IC_OBJ%V0(I) .EQ. 1.D0) THEN
            DO J = I, SIZE(X) 
              IC_OBJ%V0(J+1) = IC_OBJ%V0(M + 1 - J)
            END DO
            EXIT
          END IF
        END DO
!
!         debug with cosine wave
!
!         IC_OBJ%V0 = 0.D0
!         DO I = 1, SIZE(X)
!           IC_OBJ%V0(I) = COS(2.D0 * 4.D0 * ATAN(1.D0) * X(I))
!         END DO
!
      RETURN
      END SUBROUTINE INIT_GAUSS_PULSE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE DELETE_ICS (IC_OBJ)
!
        CLASS (INIT_CONDS) :: IC_OBJ
!
!-----------------------------------------------------------------------
!
        IF (ALLOCATED(IC_OBJ%V0)) DEALLOCATE (IC_OBJ%V0)
!
      RETURN
      END SUBROUTINE DELETE_ICS
!
!-----------------------------------------------------------------------
      END MODULE INITIAL_CONDS_M
