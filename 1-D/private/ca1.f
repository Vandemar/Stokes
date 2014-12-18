      PROGRAM CA1
!
!       steven wopschall
!       mat228a
!       computing assignment number 1
!       10/21/14
!
!       this program runs the various numerical methods for computing
!       assignment number 1.
!
!-----------------------------------------------------------------------
!
        USE FD_METHODS_M
        USE INITIAL_CONDS_M 
!
!       declare the derived types associated with the three finite 
!       difference methods and the exact solution
!
        TYPE (UPWIND) :: UPW
        TYPE (LAX_FRIEDRICHS) :: LF
        TYPE (LAX_WENDROFF) :: LW
        TYPE (EXACT_SOL) :: SOL
!
!       declare the analysis control derived type
!
        TYPE (ANALYSIS_CONTROL) :: A_CONTROL
!
!       declare any local variables
!
        INTEGER :: I, J, K
!
!       declare the maximum time step integer that acts as the upper
!       bound on the time step loop
!
        INTEGER :: MAX_ITER
!
!       enter the initial condition loop.  this loops over the three
!       initial condition cases for the three finite difference schemes
!
        DO I = 1, 3
!
          print *, 'we are at the ith IC iteration: ', i
!
!         assign the initial conditions type to the IC loop index
!
          UPW%ICS%IC_TYPE = I
          LF%ICS%IC_TYPE = I
          LW%ICS%IC_TYPE = I
          SOL%ICS%IC_TYPE = I
!
!         enter the loop over the space discretization
!
          DO J = 1, 11
!
!           set the analysis control variable, M_CASE, for the number of grid
!           points.  note that the grid points run from 0 to M, but we will
!           index from 1.  There will be an M+1 adjustment where necessary
!
            A_CONTROL%M_CASE = J
!
            print *, 'we are at the jth iteration: ', j
!
!           call the type bound procedures that get the finite difference
!           method's data.  the actual problem parameters are hard coded
!           into this type bound procedure on the base type, which all
!           extended types inherit
!
            CALL UPW%GET_FD_DATA (A_CONTROL%M_CASE)
            CALL LF%GET_FD_DATA (A_CONTROL%M_CASE)
            CALL LW%GET_FD_DATA (A_CONTROL%M_CASE)
            CALL SOL%GET_EXACT_DATA (A_CONTROL%M_CASE)
!
!           initialize the beginning step solution to the initial conditions
!
            UPW%U_B = UPW%ICS%V0
            LF%U_B = LF%ICS%V0
            LW%U_B = LW%ICS%V0
!
!           enter the time step iteration loop.  note: the maximum 
!           time step value is the same for all methods
!
            MAX_ITER = UPW%GP%MAX_ITER 
!
            DO K = 1, MAX_ITER
!
              UPW%ITER = K
              LF%ITER = K
              LW%ITER = K
!
              IF ((UPW%GP%T .GT. UPW%GP%T_END) .OR.                     &
     &           (UPW%GP%DT .EQ. 0.D0)) EXIT
              IF ((LF%GP%T .GT. LF%GP%T_END) .OR.                       &
     &           (LF%GP%DT .EQ. 0.D0)) EXIT
              IF ((LW%GP%T .GT. LW%GP%T_END) .OR.                       &
     &           (LW%GP%DT .EQ. 0.D0)) EXIT
!
              CALL UPW%FD_UPWIND
              CALL LF%FD_LAX_FRIEDRICHS
              CALL LW%FD_LAX_WENDROFF
!
              SOL%DT = UPW%GP%DT
              SOL%T = UPW%GP%T
              CALL SOL%EXACT_SOL_UPDATE
!
!           end the time step loop
!
            END DO 
!
!           get the norms associated with the last time step
!
            CALL UPW%GET_NORMS (SOL)
            CALL LF%GET_NORMS (SOL)
            CALL LW%GET_NORMS (SOL)
!
!           call output routines
!
            CALL UPW%FD_OUTPUT (A_CONTROL, SOL)
            CALL LF%FD_OUTPUT (A_CONTROL, SOL)
            CALL LW%FD_OUTPUT (A_CONTROL, SOL)
            CALL SOL%EXACT_SOL_OUTPUT(A_CONTROL)
!
!           call the delete routines to get ready to reallocate the
!           various allocatable component arrays for the next
!           discretization loop
!
            CALL UPW%DELETE
            CALL LF%DELETE
            CALL LW%DELETE
            CALL SOL%DELETE_EXACT_SOL
!
!         end the discretization loop
!
          END DO
!
!       end the initial conditions loop
!  
        END DO
!
!-----------------------------------------------------------------------
!
      END PROGRAM CA1
