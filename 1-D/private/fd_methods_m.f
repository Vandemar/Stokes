      MODULE FD_METHODS_M
!
!     This module contains the solution update schemes for various 
!     finite difference methods.  Currently, these are only 1-D
!     implementations. there are also derived types pertaining
!     to global parameters, etc. and various TBPs
!
!-----------------------------------------------------------------------
!
      USE INITIAL_CONDS_M
      USE NORMS_M
      IMPLICIT NONE
!
!=======================================================================
!
      TYPE ANALYSIS_CONTROL
        INTEGER :: M_CASE
      END TYPE ANALYSIS_CONTROL
!
!-----------------------------------------------------------------------
!
      TYPE GLOBAL_PARAMS
        REAL (KIND = KIND(0.D0)) :: H, DT, T, T_END, SIG, A
        INTEGER :: NDIM, M, MAX_ITER
        CHARACTER :: BCS
      END TYPE GLOBAL_PARAMS
!
!-----------------------------------------------------------------------
!
      TYPE EXACT_SOL
        REAL (KIND = KIND(0.D0)), ALLOCATABLE :: U(:), X(:)
        REAL (KIND = KIND(0.D0)) :: H, A, T, DT
        INTEGER :: M
        TYPE (INIT_CONDS) :: ICS
        CONTAINS
        PROCEDURE :: GET_EXACT_DATA
        PROCEDURE :: EXACT_SOL_UPDATE
        PROCEDURE :: EXACT_SOL_OUTPUT
        PROCEDURE :: DELETE_EXACT_SOL
      END TYPE EXACT_SOL
!
!-----------------------------------------------------------------------
!
        TYPE :: FINITE_DIFF
        REAL (KIND = KIND(0.D0)), ALLOCATABLE :: U_B(:), U_E(:),        &
     &                                           U_HALF(:), X(:),       &
     &                                           NORMS(:)
        INTEGER :: FD_ID, ITER
        TYPE (GLOBAL_PARAMS) :: GP
        TYPE (INIT_CONDS) :: ICS
        CONTAINS
        PROCEDURE :: GET_FD_DATA
        PROCEDURE :: GET_NORMS
        PROCEDURE :: FD_OUTPUT
        PROCEDURE :: DELETE
      END TYPE FINITE_DIFF
!
!-----------------------------------------------------------------------
!
!     extensions of type FINITE_DIFF
!
      TYPE, EXTENDS (FINITE_DIFF) :: UPWIND
        CONTAINS
        PROCEDURE :: FD_UPWIND
      END TYPE UPWIND
!
!-----------------------------------------------------------------------
!
      TYPE, EXTENDS (FINITE_DIFF) :: LAX_FRIEDRICHS
        CONTAINS
        PROCEDURE :: FD_LAX_FRIEDRICHS
      END TYPE LAX_FRIEDRICHS
!
!-----------------------------------------------------------------------
!
      TYPE, EXTENDS (FINITE_DIFF) :: LAX_WENDROFF
        CONTAINS
        PROCEDURE :: FD_LAX_WENDROFF
      END TYPE LAX_WENDROFF
!
!=======================================================================
      CONTAINS
!=======================================================================
!
      SUBROUTINE FD_UPWIND (FD)
!
        CLASS (UPWIND) :: FD
!
!       declare local variables
!
        INTEGER :: I
!
!-----------------------------------------------------------------------
!
        ASSOCIATE (U_B => FD%U_B, U_E => FD%U_E, V0 => FD%ICS%V0,       &
     &             SIG => FD%GP%SIG)
!
!       upwind scheme
!
        FORALL (I = 2 : SIZE(U_B))
          U_E(I) = U_B(I) + SIG * (U_B(I - 1) - U_B(I))
        END FORALL
!
!       handle periodic boundary conditions
!
        U_E(1) = U_E(SIZE(U_B))
        U_B = U_E
!
        IF (FD%GP%T + FD%GP%DT .GT. FD%GP%T_END) THEN
          FD%GP%DT = 9.D0 - FD%GP%T
          SIG = FD%GP%DT * FD%GP%A / FD%GP%H
        END IF
!
        FD%GP%T = FD%ITER * FD%GP%DT
!
        END ASSOCIATE
!
      RETURN
      END SUBROUTINE FD_UPWIND
!
!-----------------------------------------------------------------------
!
      SUBROUTINE FD_LAX_FRIEDRICHS (FD)
!
        CLASS (LAX_FRIEDRICHS) :: FD
!
!       declare local variables
!
        INTEGER :: I
!
!-----------------------------------------------------------------------
!
        ASSOCIATE (U_B => FD%U_B, U_E => FD%U_E, V0 => FD%ICS%V0,       &
     &             SIG => FD%GP%SIG)
!
!       lax-friedrichs scheme
!
        FORALL (I = 2 : SIZE(U_B)-1)
         U_E(I) = 0.5D0 * (U_B(I-1) + U_B(I+1)) + 0.5D0 * SIG *         &
     &            (U_B(I-1) - U_B(I+1)) 
        END FORALL
!
!       handle periodic boundary conditions
!
        U_E(SIZE(U_B)) = 0.5D0 * (U_B(SIZE(U_B)-1) + U_B(1)) + 0.5D0 *  &
     &                 SIG * (U_B(SIZE(U_B)-1) - U_B(1))
        U_E(1) = U_E(SIZE(U_B))
        U_B = U_E
!
        IF (FD%GP%T + FD%GP%DT .GT. FD%GP%T_END) THEN
          FD%GP%DT = 9.D0 - FD%GP%T
          SIG = 0.9D0
        END IF
!
        FD%GP%T = FD%ITER * FD%GP%DT
!
        END ASSOCIATE
!
      RETURN
      END SUBROUTINE FD_LAX_FRIEDRICHS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE FD_LAX_WENDROFF (FD)
!
        CLASS (LAX_WENDROFF) :: FD
!
!       declare local variables
!
        INTEGER :: I
!
!-----------------------------------------------------------------------
!
        ASSOCIATE (U_B => FD%U_B, U_E => FD%U_E, V0 => FD%ICS%V0,       &
     &             U_HALF => FD%U_HALF, SIG => FD%GP%SIG)
!
        U_HALF = 0.D0
         U_E = 0.D0
!
          FORALL (I = 1 : (SIZE(U_B) - 1))
            U_HALF(I) = 0.5D0 * (U_B(I) + U_B(I+1)) + SIG * 0.5D0 *     &
     &                  (U_B(I) - U_B(I+1))
          END FORALL
          U_HALF(SIZE(U_HALF)) = U_HALF(1)
!
          FORALL (I = 2 : (SIZE(U_B)))
            U_E(I) = U_B(I) + SIG * (U_HALF(I-1) - U_HALF(I))
          END FORALL
          U_E(1) = U_E(SIZE(U_B))
          U_B = U_E
!
        IF (FD%GP%T + FD%GP%DT .GT. FD%GP%T_END) THEN
          FD%GP%DT = 9.D0 - FD%GP%T
          SIG = FD%GP%DT * FD%GP%A / FD%GP%H
        END IF
!
        FD%GP%T = FD%ITER * FD%GP%DT
!
        END ASSOCIATE
!
      RETURN
      END SUBROUTINE FD_LAX_WENDROFF
!
!-----------------------------------------------------------------------
!
      SUBROUTINE DELETE (FD)
!
        CLASS (FINITE_DIFF), INTENT (INOUT) :: FD
!
!-----------------------------------------------------------------------
!
!       deallocate the base type's allocatable arrays
!
        IF (ALLOCATED(FD%U_B)) DEALLOCATE(FD%U_B)
        IF (ALLOCATED(FD%U_E)) DEALLOCATE(FD%U_E)
        IF (ALLOCATED(FD%U_HALF)) DEALLOCATE(FD%U_HALF)
        IF (ALLOCATED(FD%X)) DEALLOCATE(FD%X)
        IF (ALLOCATED(FD%NORMS)) DEALLOCATE(FD%NORMS)
!
!       call the delete TBP for the initial conditions derived type
!       this deallocates the initial conditions array
!
        CALL FD%ICS%DELETE_ICS
!
      RETURN
      END SUBROUTINE DELETE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_FD_DATA (FD, M_CASE)
!
        CLASS (FINITE_DIFF) :: FD
        INTEGER :: M_CASE
!
!       declare local variables
!
        INTEGER :: I
!
!-----------------------------------------------------------------------
!
        ASSOCIATE (H => FD%GP%H, DT => FD%GP%DT, T => FD%GP%T,          &
     &             T_END => FD%GP%T_END, SIG => FD%GP%SIG,              &
     &             NDIM => FD%GP%NDIM, M => FD%GP%M, A => FD%GP%A,      &
     &             V0 => FD%ICS%V0, IC_TYPE => FD%ICS%IC_TYPE)
!
!       set problem data.  set M actually to M+1 to account for indexing
!       from 1 instead of from 0
!
!       set the finite difference method id for purposes of naming the
!       file output
!
        SELECT TYPE (FD)
          TYPE IS (UPWIND)
            FD%FD_ID = 1
          TYPE IS (LAX_FRIEDRICHS)
            FD%FD_ID = 2
          TYPE IS (LAX_WENDROFF)
            FD%FD_ID = 3
        CLASS DEFAULT
        END SELECT
!
!       set the H value manually.  as of now the code does not loop over
!       the various H values per the assignment.
!
        SELECT CASE (M_CASE)
          CASE (1)
            H = 1.D0 / 64.D0
            M = 64
          CASE (2) 
            H = 1.D0 / 128.D0
            M = 128
          CASE (3)
            H = 1.D0 / 256.D0
            M = 256
          CASE (4)
            H = 1.D0 / 512.D0
            M = 512
          CASE (5)
            H = 1.D0 / 1024.D0
            M = 1024
          CASE (6)
            H = 1.D0 / 2048.D0
            M = 2048
          CASE (7)
            H = 1.D0 / 4096.D0
            M = 4096
          CASE (8)
            H = 1.D0 / 8192.D0
            M = 8192
          CASE (9)
            H = 1.D0 / 16384.D0
            M = 16384
          CASE (10)
            H = 1.D0 / 32768.D0
            M = 32768
          CASE (11)
            H = 1.D0 / 65536.D0
            M = 65536
        END SELECT
!
!       assign the rest of the global parameters
!
        T = 0.D0
        T_END = 9.D0
        SIG = 0.9D0
        A = 3.D0
        DT = SIG * H / A
        NDIM = 1
!
!       set the maximum time step integer for the time step loop
!       in the main program.  
!
!        FD%GP%MAX_ITER = NINT (9.D0 / DT) + 1 
        FD%GP%MAX_ITER = INT(9.D0 / DT)
        FD%ITER = 0
!
        print *, 'the max iter is: ', FD%GP%MAX_ITER
        print *, 'the max time is: ', FD%GP%MAX_ITER * DT
!
!       allocate space for the allocatable arrays in each derived type
!
        IF (.NOT. ALLOCATED(FD%U_B)) ALLOCATE(FD%U_B(M + 1))
        IF (.NOT. ALLOCATED(FD%U_E)) ALLOCATE(FD%U_E(M + 1))
        IF (.NOT. ALLOCATED(FD%U_HALF)) ALLOCATE(FD%U_HALF(M + 1))
        IF (.NOT. ALLOCATED(FD%ICS%V0)) ALLOCATE(FD%ICS%V0(M + 1))
        IF (.NOT. ALLOCATED(FD%X)) ALLOCATE(FD%X(M + 1))
!
!       initialize the allocatable arrays
!
        FD%U_B = 0.D0
        FD%U_E = 0.D0
        FD%U_HALF = 0.D0
        FD%ICS%V0 = 0.D0
        FD%X = 0.D0
!
!       populate the space vectore
!
        FORALL (I = 1 : M + 1)
          FD%X(I) = H * (I - 1.D0) 
        END FORALL
!
!       call routines to populate the initial conditions vector
!
        IF (IC_TYPE .EQ. 1) CALL FD%ICS%INIT_SQ_WAVE (FD%X, M, H)
        IF (IC_TYPE .EQ. 2) CALL FD%ICS%INIT_SEMI_CIRCLE (FD%X, M, H)
        IF (IC_TYPE .EQ. 3) CALL FD%ICS%INIT_GAUSS_PULSE (FD%X, M, H)
!
        END ASSOCIATE
!
      RETURN
      END SUBROUTINE GET_FD_DATA
!
!-----------------------------------------------------------------------
!
      SUBROUTINE EXACT_SOL_UPDATE (SOL)
!
        CLASS (EXACT_SOL) :: SOL
!
!       declare local variables
!
        INTEGER :: I
        REAL (KIND = KIND(0.D0)), ALLOCATABLE :: XT(:)
!
!-----------------------------------------------------------------------
!
!       note: this gives the same result as realizing that at T = 9
!       the exact solution is the same as at T = 0.
!
        IF (.NOT. ALLOCATED(XT)) ALLOCATE(XT(SIZE(SOL%X)))
        DO I = 1, SIZE(XT)
          XT(I) = SOL%X(I) - SOL%A * SOL%DT
          IF ((XT(I) .LT. 0.D0) .OR. (XT(I) .GT. 1.D0)) THEN
            XT(I) = XT(I) - FLOOR(XT(I))
          END IF
        END DO
!
        IF (SOL%ICS%IC_TYPE .EQ. 1) CALL SOL%ICS%INIT_SQ_WAVE (XT,      &
     &                                   SOL%M, SOL%H)
        IF (SOL%ICS%IC_TYPE .EQ. 2) CALL SOL%ICS%INIT_SEMI_CIRCLE (XT,  &
     &                                   SOL%M, SOL%H)
        IF (SOL%ICS%IC_TYPE .EQ. 3) CALL SOL%ICS%INIT_GAUSS_PULSE (XT,  &
     &                                   SOL%M, SOL%H) 
!
      SOL%X = XT
      SOL%U = SOL%ICS%V0
!
      IF (ALLOCATED(XT)) DEALLOCATE (XT)
      RETURN
      END SUBROUTINE EXACT_SOL_UPDATE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_EXACT_DATA (SOL, M_CASE)
!
        CLASS (EXACT_SOL) :: SOL
        INTEGER :: M_CASE
!
!       declare local variables
!
        INTEGER :: I
!
!-----------------------------------------------------------------------
!
        ASSOCIATE (IC_TYPE => SOL%ICS%IC_TYPE, M => SOL%M, H => SOL%H)
!
        SOL%A = 3.D0
        SOL%T = 0.D0
        SELECT CASE (M_CASE)
          CASE (1)
            H = 1.D0 / 64.D0
            M = 64
          CASE (2)
            H = 1.D0 / 128.D0
            M = 128
          CASE (3)
            H = 1.D0 / 256.D0
            M = 256
          CASE (4)
            H = 1.D0 / 512.D0
            M = 512
          CASE (5)
            H = 1.D0 / 1024.D0
            M = 1024
          CASE (6)
            H = 1.D0 / 2048.D0
            M = 2048
          CASE (7)
            H = 1.D0 / 4096.D0
            M = 4096
          CASE (8)
            H = 1.D0 / 8192.D0
            M = 8192
          CASE (9)
            H = 1.D0 / 16384.D0
            M = 16384
          CASE (10)
            H = 1.D0 / 32768.D0
            M = 32768
          CASE (11)
            H = 1.D0 / 65536.D0
            M = 65536
        END SELECT
!
        IF (.NOT. ALLOCATED(SOL%X)) ALLOCATE(SOL%X(M + 1))
        IF (.NOT. ALLOCATED(SOL%U)) ALLOCATE(SOL%U(M + 1))
        IF (.NOT. ALLOCATED(SOL%ICS%V0)) ALLOCATE(SOL%ICS%V0(M + 1))
!
!       initialize allocatable arrays
!
        SOL%X = 0.D0
        SOL%U = 0.D0
        SOL%ICS%V0 = 0.D0
!
!       populate the space vectore
!
        FORALL (I = 1 : M + 1)
          SOL%X(I) = H * (I - 1.D0)
        END FORALL            
!
!       call routine to populate the initial conditions vector
!
        IF (IC_TYPE .EQ. 1) CALL SOL%ICS%INIT_SQ_WAVE (SOL%X, M, H)
        IF (IC_TYPE .EQ. 2) CALL SOL%ICS%INIT_SEMI_CIRCLE (SOL%X, M, H)
        IF (IC_TYPE .EQ. 3) CALL SOL%ICS%INIT_GAUSS_PULSE (SOL%X, M, H)
!
        END ASSOCIATE
!
      RETURN
      END SUBROUTINE GET_EXACT_DATA
!
!-----------------------------------------------------------------------
!
      SUBROUTINE DELETE_EXACT_SOL (SOL)
!
        CLASS (EXACT_SOL) :: SOL
!
!-----------------------------------------------------------------------
!
        IF (ALLOCATED(SOL%U)) DEALLOCATE(SOL%U)
        IF (ALLOCATED(SOL%X)) DEALLOCATE(SOL%X)
        CALL SOL%ICS%DELETE_ICS
!
      RETURN
      END SUBROUTINE DELETE_EXACT_SOL
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GET_NORMS (FD, SOL)
!
        CLASS (FINITE_DIFF), INTENT (INOUT) :: FD
        TYPE (EXACT_SOL), INTENT (IN) :: SOL
!
!       declare local variables
!
        TYPE (NORMS) :: NORM_OBJ
        REAL (KIND = KIND(0.D0)), ALLOCATABLE :: ARRAY(:)
        REAL (KIND = KIND(0.D0)) :: NORM_VAL
!
!-----------------------------------------------------------------------
!
        IF (.NOT. ALLOCATED(FD%NORMS)) ALLOCATE(FD%NORMS(3))
        IF (.NOT. ALLOCATED(ARRAY)) ALLOCATE(ARRAY(SIZE(FD%U_E)))
        FD%NORMS = 0.D0
        ARRAY = 0.D0
        ARRAY = SOL%U - FD%U_B
        CALL NORM_OBJ%NORM(ARRAY, FD%GP%H, 1, NORM_VAL)
        FD%NORMS(1) = NORM_VAL
        CALL NORM_OBJ%NORM(ARRAY, FD%GP%H, 2, NORM_VAL)
        FD%NORMS(2) = NORM_VAL
        CALL NORM_OBJ%NORM(ARRAY, FD%GP%H, 3, NORM_VAL)
        FD%NORMS(3) = NORM_VAL
!
!        IF (FD%ICS%IC_TYPE .EQ. 3) THEN
!          print *, 'norms: ', fd%gp%h, fd%norms(2)
!        END IF
!
        IF (ALLOCATED(ARRAY)) DEALLOCATE(ARRAY)
      RETURN
      END SUBROUTINE GET_NORMS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE FD_OUTPUT (FD, A_CONTROL, SOL)
!
        CLASS (FINITE_DIFF) :: FD
        TYPE (ANALYSIS_CONTROL) :: A_CONTROL
        TYPE (EXACT_SOL) :: SOL
!
!       declare local variables
!
        INTEGER :: I
!
!-----------------------------------------------------------------------
!
!       this output routine outputs data specific to the computing
!       assignment number 1.  The output is the spatial discretization
!       vector, the approximate solution and the exact solutino at the
!       required discretizations per the assignment.  The output is also
!       the one norm, two norm, and max norm, in that order for all
!       analysis cases
!
!       output for upwind
!
        IF (FD%FD_ID .EQ. 1) THEN
          IF (A_CONTROL%M_CASE .EQ. 1) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u1.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm1.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u1.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm1.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u1.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm1.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I)
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 2) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/upwind1_norm2.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/upwind2_norm2.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/upwind3_norm2.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 3) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u3.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm3.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u3.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm3.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u3.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm3.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 4) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/upwind1_norm4.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/upwind2_norm4.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/upwind3_norm4.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 5) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u5.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm5.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u5.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm5.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u5.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm5.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 6) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/upwind1_norm6.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/upwind2_norm6.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/upwind3_norm6.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 7) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u7.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm7.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u7.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm7.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u7.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm7.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 8) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u8.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm8.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u8.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm8.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u8.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm8.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 9) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u9.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm9.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u9.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm9.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u9.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm9.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 10) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u10.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm10.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u10.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm10.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u10.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm10.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 11) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/upwind1_u11.txt")
                OPEN (unit = 4, file = "norms/upwind1_norm11.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/upwind2_u11.txt")
                OPEN (unit = 4, file = "norms/upwind2_norm11.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/upwind3_u11.txt")
                OPEN (unit = 4, file = "norms/upwind3_norm11.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
        END IF
!
!       output for lax-friedrichs
!
        IF (FD%FD_ID .EQ. 2) THEN
          IF (A_CONTROL%M_CASE .EQ. 1) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u1.txt")
                OPEN (unit = 4, file = "norms/lf1_norm1.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u1.txt")
                OPEN (unit = 4, file = "norms/lf2_norm1.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u1.txt")
                OPEN (unit = 4, file = "norms/lf3_norm1.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 2) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/lf1_norm2.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/lf2_norm2.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/lf3_norm2.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 3) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u3.txt")
                OPEN (unit = 4, file = "norms/lf1_norm3.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u3.txt")
                OPEN (unit = 4, file = "norms/lf2_norm3.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u3.txt")
                OPEN (unit = 4, file = "norms/lf3_norm3.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 4) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/lf1_norm4.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/lf2_norm4.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/lf3_norm4.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 5) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u5.txt")
                OPEN (unit = 4, file = "norms/lf1_norm5.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u5.txt")
                OPEN (unit = 4, file = "norms/lf2_norm5.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u5.txt")
                OPEN (unit = 4, file = "norms/lf3_norm5.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 6) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/lf1_norm6.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/lf2_norm6.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/lf3_norm6.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 7) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u7.txt")
                OPEN (unit = 4, file = "norms/lf1_norm7.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u7.txt")
                OPEN (unit = 4, file = "norms/lf2_norm7.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u7.txt")
                OPEN (unit = 4, file = "norms/lf3_norm7.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 8) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u8.txt")
                OPEN (unit = 4, file = "norms/lf1_norm8.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u8.txt")
                OPEN (unit = 4, file = "norms/lf2_norm8.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u8.txt")
                OPEN (unit = 4, file = "norms/lf3_norm8.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 9) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u9.txt")
                OPEN (unit = 4, file = "norms/lf1_norm9.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u9.txt")
                OPEN (unit = 4, file = "norms/lf2_norm9.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u9.txt")
                OPEN (unit = 4, file = "norms/lf3_norm9.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 10) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u10.txt")
                OPEN (unit = 4, file = "norms/lf1_norm10.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u10.txt")
                OPEN (unit = 4, file = "norms/lf2_norm10.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u10.txt")
                OPEN (unit = 4, file = "norms/lf3_norm10.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 11) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lf1_u11.txt")
                OPEN (unit = 4, file = "norms/lf1_norm11.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lf2_u11.txt")
                OPEN (unit = 4, file = "norms/lf2_norm11.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lf3_u11.txt")
                OPEN (unit = 4, file = "norms/lf3_norm11.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
        END IF
!
!       output for lax-wendroff
!
        IF (FD%FD_ID .EQ. 3) THEN
          IF (A_CONTROL%M_CASE .EQ. 1) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u1.txt")
                OPEN (unit = 4, file = "norms/lw1_norm1.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u1.txt")
                OPEN (unit = 4, file = "norms/lw2_norm1.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u1.txt")
                OPEN (unit = 4, file = "norms/lw3_norm1.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 2) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/lw1_norm2.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/lw2_norm2.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/lw3_norm2.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 3) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u3.txt")
                OPEN (unit = 4, file = "norms/lw1_norm3.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u3.txt")
                OPEN (unit = 4, file = "norms/lw2_norm3.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u3.txt")
                OPEN (unit = 4, file = "norms/lw3_norm3.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 4) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/lw1_norm4.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/lw2_norm4.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/lw3_norm4.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 5) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u5.txt")
                OPEN (unit = 4, file = "norms/lw1_norm5.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u5.txt")
                OPEN (unit = 4, file = "norms/lw2_norm5.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u5.txt")
                OPEN (unit = 4, file = "norms/lw3_norm5.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 6) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 4, file = "norms/lw1_norm6.txt")
              CASE (2)
                OPEN (unit = 4, file = "norms/lw2_norm6.txt")
              CASE (3)
                OPEN (unit = 4, file = "norms/lw3_norm6.txt")
            END SELECT
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 7) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u7.txt")
                OPEN (unit = 4, file = "norms/lw1_norm7.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u7.txt")
                OPEN (unit = 4, file = "norms/lw2_norm7.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u7.txt")
                OPEN (unit = 4, file = "norms/lw3_norm7.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 8) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u8.txt")
                OPEN (unit = 4, file = "norms/lw1_norm8.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u8.txt")
                OPEN (unit = 4, file = "norms/lw2_norm8.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u8.txt")
                OPEN (unit = 4, file = "norms/lw3_norm8.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 9) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u9.txt")
                OPEN (unit = 4, file = "norms/lw1_norm9.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u9.txt")
                OPEN (unit = 4, file = "norms/lw2_norm9.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u9.txt")
                OPEN (unit = 4, file = "norms/lw3_norm9.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 10) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u10.txt")
                OPEN (unit = 4, file = "norms/lw1_norm10.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u10.txt")
                OPEN (unit = 4, file = "norms/lw2_norm10.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u10.txt")
                OPEN (unit = 4, file = "norms/lw3_norm10.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
!
          IF (A_CONTROL%M_CASE .EQ. 11) THEN
            SELECT CASE (FD%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file = "output/lw1_u11.txt")
                OPEN (unit = 4, file = "norms/lw1_norm11.txt")
              CASE (2)
                OPEN (unit = 2, file = "output/lw2_u11.txt")
                OPEN (unit = 4, file = "norms/lw2_norm11.txt")
              CASE (3)
                OPEN (unit = 2, file = "output/lw3_u11.txt")
                OPEN (unit = 4, file = "norms/lw3_norm11.txt")
            END SELECT
            DO I = 1, SIZE(FD%U_E)
              WRITE(2,*) FD%X(I),' ',FD%U_E(I),' ',SOL%U(I) 
            END DO
            CLOSE (2)
            WRITE(4,*) FD%GP%H, ' ',FD%NORMS(1), ' ',FD%NORMS(2), ' ',  &
     &                 FD%NORMS(3)
            CLOSE (4)
          END IF
        END IF
!
      RETURN
      END SUBROUTINE FD_OUTPUT
!
!-----------------------------------------------------------------------
!
      SUBROUTINE EXACT_SOL_OUTPUT (SOL, A_CONTROL)
!
        CLASS (EXACT_SOL) :: SOL
        TYPE (ANALYSIS_CONTROL) :: A_CONTROL
!
!       declare local variables
!
        INTEGER :: I
!
!-----------------------------------------------------------------------
!
        IF (A_CONTROL%M_CASE .EQ. 1) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file="output/exact_plots/exact1_u1.txt")
              CASE (2)
                OPEN (unit = 2, file="output/exact_plots/exact2_u1.txt")
              CASE (3)
                OPEN (unit = 2, file="output/exact_plots/exact3_u1.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
        IF (A_CONTROL%M_CASE .EQ. 3) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file="output/exact_plots/exact1_u3.txt")
              CASE (2)
                OPEN (unit = 2, file="output/exact_plots/exact2_u3.txt")
              CASE (3)
                OPEN (unit = 2, file="output/exact_plots/exact3_u3.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
        IF (A_CONTROL%M_CASE .EQ. 5) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file="output/exact_plots/exact1_u5.txt")
              CASE (2)
                OPEN (unit = 2, file="output/exact_plots/exact2_u5.txt")
              CASE (3)
                OPEN (unit = 2, file="output/exact_plots/exact3_u5.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
        IF (A_CONTROL%M_CASE .EQ. 7) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file="output/exact_plots/exact1_u7.txt")
              CASE (2)
                OPEN (unit = 2, file="output/exact_plots/exact2_u7.txt")
              CASE (3)
                OPEN (unit = 2, file="output/exact_plots/exact3_u7.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
        IF (A_CONTROL%M_CASE .EQ. 8) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file="output/exact_plots/exact1_u8.txt")
              CASE (2)
                OPEN (unit = 2, file="output/exact_plots/exact2_u8.txt")
              CASE (3)
                OPEN (unit = 2, file="output/exact_plots/exact3_u8.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
        IF (A_CONTROL%M_CASE .EQ. 9) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2, file="output/exact_plots/exact1_u9.txt")
              CASE (2)
                OPEN (unit = 2, file="output/exact_plots/exact2_u9.txt")
              CASE (3)
                OPEN (unit = 2, file="output/exact_plots/exact3_u9.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
        IF (A_CONTROL%M_CASE .EQ. 10) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2,file="output/exact_plots/exact1_u10.txt")
              CASE (2)
                OPEN (unit = 2,file="output/exact_plots/exact2_u10.txt")
              CASE (3)
                OPEN (unit = 2,file="output/exact_plots/exact3_u10.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
        IF (A_CONTROL%M_CASE .EQ. 11) THEN
            SELECT CASE (SOL%ICS%IC_TYPE)
              CASE (1)
                OPEN (unit = 2,file="output/exact_plots/exact1_u11.txt")
              CASE (2)
                OPEN (unit = 2,file="output/exact_plots/exact2_u11.txt")
              CASE (3)
                OPEN (unit = 2,file="output/exact_plots/exact3_u11.txt")
            END SELECT
          DO I = 1, SIZE(SOL%U)
            WRITE(2,*) SOL%X(I),' ',SOL%U(I) 
          END DO
          CLOSE (2)
        END IF
!
      RETURN
      END SUBROUTINE EXACT_SOL_OUTPUT
!
!-----------------------------------------------------------------------
!
      END MODULE FD_METHODS_M
