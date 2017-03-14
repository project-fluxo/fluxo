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
!> This Testcase initializes a background state called Ueq, which can be taken from an exactfunction (testcase_exactfunc) or read
!> from a mesh file that contains the MHD equilibrium solution. 
!> The time-derivative of the state Ut_Eq=Ut(Eq) is evaluated and introduced as a source, such that the new time derivative becomes
!> Ut_new = Ut - Ut_Eq(U_eq)
!> basically forcing  U_eq to be a steady state
!==================================================================================================================================
MODULE MOD_Testcase
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitTestcase
  MODULE PROCEDURE InitTestcase
END INTERFACE

INTERFACE FinalizeTestcase
  MODULE PROCEDURE FinalizeTestcase
END INTERFACE


PUBLIC:: InitTestcase
PUBLIC:: FinalizeTestcase
PUBLIC:: DefineParametersTestcase 

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersTestcase()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Testcase")
CALL prms%CreateStringOption('WhichTestcase', &
     "name of the testcase that is supposed to be used in simulation. Check if its = current testcase." &
     ,"mhd_equilibrium")
CALL prms%CreateIntOption('EquilibriumStateIni', &
     "=-1: Default: no equilibrium state used. U_t(U_eq)=0. Sanity check if code is compiled with this testcase."//&
     "= 0 : Use U_eq=exactFunc(IniExactFunc) for equilibrium state "//&
     "> 0 : Use U_eq=exactFunc(EquilibriumStateIni) for equilibrium state "//&
     "=-2 : Read U_eq from MeshFile"//&
     "=-3 : Read U_eq from the solution of a stateFile" &
     ,"-1")
CALL prms%CreateStringOption('EquilibriumStateFile', &
     "name of statefile to used for eq. solution (only if EquilibriumStateIni=-3)")
CALL prms%CreateLogicalOption('EquilibriumDivBcorr', &
     "switch to compute B from a vector potential instead of using B directly" , '.FALSE.')
CALL prms%CreateLogicalOption('CalcErrorToEquilibrium', &
     "switch for TC_analyze: compute difference of |U-Ueq|" , '.FALSE.')
CALL prms%CreateLogicalOption('CalcDeltaBEnergy', &
     "switch for TC_analyze: compute Energy of 1/(2mu0)|B-Beq|^2" , '.FALSE.')
CALL prms%CreateIntOption('EquilibriumDisturbFunc', &
     "=0: Default: disturbance case number = iniExactFunc"// &
     ">0: specific disturbance function added to initial state." &
     ,"0")
END SUBROUTINE DefineParametersTestcase

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE InitTestcase()
! MODULES
USE MOD_Globals
USE MOD_Testcase_Vars
USE MOD_EquilibriumState   ,ONLY: InitEquilibriumState
USE MOD_ReadInTools        ,ONLY: GETINT,GETLOGICAL,GETSTR
USE MOD_Restart_Vars       ,ONLY: DoRestart
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE "MHD EQUILIBRIUM"...'
WhichTestcase=GETSTR('WhichTestcase','mhd_equilibrium')
IF(INDEX(TRIM(WhichTestcase),'mhd_equilibrium').EQ.0)THEN
  CALL abort(__STAMP__,&
       "compiled with testcase mhd_equilibrium, but specified testcase is : "//TRIM(whichTestcase))
END IF !check whichtestcase
doTCsource=.TRUE.
EquilibriumStateIni=GETINT('EquilibriumStateIni','-1')
IF(EquilibriumStateIni.EQ.-3) THEN
  EquilibriumStateFile=GETSTR('EquilibriumStateFile')
END IF
EquilibriumDivBcorr=GETLOGICAL('EquilibriumDivBcorr','.FALSE.')
!======EQUILIBRIUM STUFF====
CALL InitEquilibriumState()
!TESTCASE ANALYZE
doCalcErrorToEquilibrium = GETLOGICAL('CalcErrorToEquilibrium','.FALSE.')
doCalcDeltaBEnergy       = GETLOGICAL('CalcDeltaBEnergy','.FALSE.')
IF(.NOT.doRestart)THEN
  EquilibriumDisturbFunc=GETINT('EquilibriumDisturbFunc','0')
END IF

SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitTestcase


!=================================================================================================================================
!> Finalize Variables 
!=================================================================================================================================
SUBROUTINE FinalizeTestcase()
! MODULES
USE MOD_Testcase_Vars
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
SDEALLOCATE(Ueq)
SDEALLOCATE(Ueq_BC)
SDEALLOCATE(Uteq)
END SUBROUTINE FinalizeTestcase

END MODULE MOD_Testcase
