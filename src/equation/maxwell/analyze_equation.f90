!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/flexi)
!
! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "defines.h"

!==================================================================================================================================
!> Contains analyze routines specific to the Maxwell's equation
!==================================================================================================================================
MODULE MOD_AnalyzeEquation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyzeEquation
  MODULE PROCEDURE InitAnalyzeEquation
END INTERFACE

INTERFACE AnalyzeEquation
  MODULE PROCEDURE AnalyzeEquation
END INTERFACE

INTERFACE FinalizeAnalyzeEquation
  MODULE PROCEDURE FinalizeAnalyzeEquation
END INTERFACE


PUBLIC :: AnalyzeEquation
PUBLIC :: InitAnalyzeEquation
PUBLIC :: FinalizeAnalyzeEquation
PUBLIC :: DefineParametersAnalyzeEquation
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE DefineParametersAnalyzeEquation

!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Analyze_Vars
!USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Get the various analysis/output variables 

! Initialize eval routines

END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
!USE MOD_Globals
!USE MOD_PreProc
!USE MOD_Analyze_Vars
!USE MOD_AnalyzeEquation_Vars
!USE MOD_Restart_Vars,       ONLY: RestartTime
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!CHARACTER(LEN=40)               :: formatStr
!==================================================================================================================================
! Attention: during the initialization phase no face data / gradients available!
!IF(ABS(Time-RestartTime) .GT. 1.E-12) THEN
!END IF

!IF(doCalcMeanFlux)THEN
!  CALL CalcMeanFlux(Time,MeanFlux)
!  IF(MPIroot) THEN
!    WRITE(formatStr,'(A5,I1,A8)')'(A14,',PP_nVar,'ES16.7)'
!    DO i=1,nBCs
!      IF(Boundarytype(i,BC_TYPE) .EQ. 1) CYCLE
!      WRITE(UNIT_StdOut,*)'MeanFlux ',TRIM(BoundaryName(i)),' : '
!      WRITE(UNIT_StdOut,formatStr)'              ',MeanFlux(:,i)
!    END DO
!  END IF
!END IF  !(doCalcBodyforces)

END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
!USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
