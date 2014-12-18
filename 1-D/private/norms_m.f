      MODULE NORMS_M
!
!     This module contains the norm derived type and any associated
!     type components and/or type bound procedures. The norm type bound
!     procedures are discrete norms.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!=======================================================================
!
      TYPE NORMS
        CONTAINS
        PROCEDURE :: NORM
        PROCEDURE :: TWO_NORM
        PROCEDURE :: ONE_NORM
        PROCEDURE :: MAX_NORM
      END TYPE NORMS
!
!-----------------------------------------------------------------------
!
      PRIVATE :: TWO_NORM, ONE_NORM, MAX_NORM
!
!=======================================================================
      CONTAINS
!=======================================================================
!
      SUBROUTINE NORM (NORM_OBJ, ARRAY, H, NORM_ID, NORM_VAL)
!
        CLASS (NORMS) :: NORM_OBJ
        REAL (KIND = KIND(0.D0)), INTENT (IN) :: ARRAY(:)
        REAL (KIND = KIND(0.D0)), INTENT (OUT) :: NORM_VAL
        REAL (KIND = KIND(0.D0)) :: H
        INTEGER :: NORM_ID
!
!-----------------------------------------------------------------------
!
        SELECT CASE (NORM_ID)
        CASE (1)
          CALL ONE_NORM (NORM_OBJ, ARRAY, H, NORM_VAL)
        CASE (2)
          CALL TWO_NORM (NORM_OBJ, ARRAY, H, NORM_VAL)
        CASE (3)
          CALL MAX_NORM (NORM_OBJ, ARRAY, NORM_VAL)
        END SELECT
!
      RETURN
      END SUBROUTINE NORM
!
!-----------------------------------------------------------------------
!
      SUBROUTINE TWO_NORM (NORM_OBJ, ARRAY, H, NORM_VAL)
!
!     declare dummy variables
!
      CLASS (NORMS) :: NORM_OBJ
      REAL (KIND = KIND(0.D0)), INTENT (IN) :: ARRAY(:)
      REAL (KIND = KIND(0.D0)), INTENT (OUT) :: NORM_VAL
      REAL (KIND = KIND(0.D0)) :: H
!
!     declare local variables
!
      REAL (KIND = KIND(0.D0)), ALLOCATABLE :: ABS_ARRAY(:) 
      REAL (KIND = KIND(0.D0)) :: FAC, M
      INTEGER :: I
!
!-----------------------------------------------------------------------
!
      IF (.NOT. ALLOCATED(ABS_ARRAY)) ALLOCATE (ABS_ARRAY(SIZE(ARRAY)))
!
      ABS_ARRAY = ABS(ARRAY)
!
      NORM_VAL = 0.D0
!
      DO I = 1, SIZE(ABS_ARRAY)
        NORM_VAL = NORM_VAL + ABS_ARRAY(I)**2
      END DO
!
      M = 1.D0 / H
      FAC = 1.D0 / DSQRT(M)
!
      NORM_VAL = DSQRT(H * NORM_VAL)
!
      DEALLOCATE (ABS_ARRAY)
!
      RETURN
      END SUBROUTINE TWO_NORM
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ONE_NORM (NORM_OBJ, ARRAY, H, NORM_VAL)
!
!     declare dummy variables
!
      CLASS (NORMS) :: NORM_OBJ
      REAL (KIND = KIND(0.D0)), INTENT (IN) :: ARRAY(:)
      REAL (KIND = KIND(0.D0)), INTENT (OUT) :: NORM_VAL
      REAL (KIND = KIND(0.D0)) :: H
!
!     declare local variables
!
      REAL (KIND = KIND(0.D0)), ALLOCATABLE :: ABS_ARRAY(:)
      REAL (KIND = KIND(0.D0)) :: FAC, M
      INTEGER :: I
!
!-----------------------------------------------------------------------
!
      IF (.NOT. ALLOCATED(ABS_ARRAY)) ALLOCATE(ABS_ARRAY(SIZE(ARRAY)))
!
      ABS_ARRAY = ABS(ARRAY)
!
      NORM_VAL = 0.D0
!
      DO I = 1, SIZE(ABS_ARRAY)
        NORM_VAL = NORM_VAL + ABS_ARRAY(I)
      END DO
!
      M = 1.D0 / H
      FAC = 1.D0 / M
!
      NORM_VAL = FAC * NORM_VAL
!
      DEALLOCATE (ABS_ARRAY)
!
      RETURN
      END SUBROUTINE ONE_NORM
!
!-----------------------------------------------------------------------
!
      SUBROUTINE MAX_NORM (NORM_OBJ, ARRAY, NORM_VAL)
!
!     declare dummy variables
!
      CLASS (NORMS) :: NORM_OBJ
      REAL (KIND = KIND(0.D0)), INTENT (IN) :: ARRAY(:)
      REAL (KIND = KIND(0.D0)), INTENT (OUT) ::  NORM_VAL
!
!     declare local variables
!
      REAL (KIND = KIND(0.D0)), ALLOCATABLE :: ABS_ARRAY(:)
!
!-----------------------------------------------------------------------
!
      IF (.NOT. ALLOCATED(ABS_ARRAY)) ALLOCATE(ABS_ARRAY(SIZE(ARRAY)))
!
      ABS_ARRAY = ABS(ARRAY)
!
      NORM_VAL = 0.D0
!
      NORM_VAL = MAXVAL(ABS_ARRAY)
!
      DEALLOCATE (ABS_ARRAY)
!
      RETURN
      END SUBROUTINE MAX_NORM
!
!-----------------------------------------------------------------------
      END MODULE NORMS_M
