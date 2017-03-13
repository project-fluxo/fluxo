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
!> Contains analyze routines specific to the linear scalar advection equation
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


PUBLIC:: AnalyzeEquation
PUBLIC:: InitAnalyzeEquation
PUBLIC:: FinalizeAnalyzeEquation
PUBLIC:: DefineParametersAnalyzeEquation
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for analyze Linadv 
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("AnalyzeEquation")
CALL prms%CreateLogicalOption('CalcMinMax'   , "Set true to compute minimum and maximum of ssolution"         , '.FALSE.')
END SUBROUTINE DefineParametersAnalyzeEquation


!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_AnalyzeEquation_Vars,ONLY:doCalcMinMax
USE MOD_Analyze_Vars,       ONLY:doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_ReadInTools,        ONLY: GETLOGICAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Get the various analysis/output variables 
doCalcMinMax  = GETLOGICAL('CalcMinMax','.FALSE.')

! Initialize eval routines
IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doCalcMinMax)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"SolutionMinimum"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"SolutionMaximum"'
  END IF  !doCalcMinMax
END IF !MPIroot.AND.doAnalyzeToFile

END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> Calculates equation related analyze variables
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_AnalyzeEquation_Vars,ONLY:doCalcMinMax
USE MOD_Analyze_Vars,       ONLY:doAnalyzeToFile,A2F_iVar,A2F_Data
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Umin,Umax
CHARACTER(LEN=40)               :: formatStr
!==================================================================================================================================
! Attention: during the initialization phase no face data / gradients available!
!IF(ABS(Time-RestartTime) .GT. 1.E-12) THEN
!END IF
IF(doCalcMinMax)THEN
  CALL CalcMinMax(Umin,Umax)
  IF(MPIroot)THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',1,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' Mininum    : ',Umin
    WRITE(UNIT_StdOut,formatStr)' Maximum    : ',Umax
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Umin
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Umax
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF !doCalcMinMax

END SUBROUTINE AnalyzeEquation




!==================================================================================================================================
!> Calculates minimum and Maximum of solution
!==================================================================================================================================
SUBROUTINE CalcMinMax(Umin,Umax)
! MODULES
USE MOD_Globals
USE MOD_DG_Vars            ,ONLY:U
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Umin
REAL,INTENT(OUT)                :: Umax
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
Umin=MINVAL(U)
Umax=MAXVAL(U)
#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Umin,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,Umax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(Umin,0           ,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(Umax,0           ,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
END SUBROUTINE CalcMinMax 


!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
!USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!==================================================================================================================================
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
