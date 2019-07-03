!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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
!> Subroutines defining one specific testcase with all necessary variables
!==================================================================================================================================
MODULE MOD_Testcase_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitAnalyzeTestCase
  MODULE PROCEDURE InitAnalyzeTestCase
END INTERFACE

INTERFACE AnalyzeTestCase
  MODULE PROCEDURE AnalyzeTestCase
END INTERFACE


PUBLIC:: DefineParametersAnalyzeTestcase 
PUBLIC:: InitAnalyzeTestCase
PUBLIC:: AnalyzeTestCase

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeTestcase()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("AnalyzeTestcase")
CALL prms%CreateLogicalOption('doTCanalyze', &
     "switch off/on TC_analyze" , '.FALSE.')
END SUBROUTINE DefineParametersAnalyzeTestcase


!==================================================================================================================================
!> Initialize Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE InitAnalyzeTestcase()
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars,     ONLY:doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_Testcase_Vars,    ONLY:doTCanalyze
USE MOD_ReadInTools,      ONLY: GETLOGICAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
doTCanalyze = GETLOGICAL('doTCanalyze','.FALSE.')

!prepare AnalyzeToFile
IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doTCanalyze)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"TCanalyze1"'
  END IF  !doTCanalyze
END IF !MPIroot & doAnalyzeToFile
END SUBROUTINE InitAnalyzeTestcase

!==================================================================================================================================
!> Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE AnalyzeTestcase(Time)
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars,     ONLY:doAnalyzeToFile,A2F_iVar,A2F_Data
USE MOD_Testcase_Vars,    ONLY:doTCanalyze
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: time         !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tmp
!==================================================================================================================================
tmp=-999.
IF(MPIroot) THEN
  IF(doTCanalyze)THEN
    WRITE(UNIT_StdOut,*)'Default testcase analyze :', tmp
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=tmp
    END IF
  END IF  !doTCanalyze
END IF !MPIroot & doAnalyzeToFile
END SUBROUTINE AnalyzeTestcase


END MODULE MOD_Testcase_Analyze
