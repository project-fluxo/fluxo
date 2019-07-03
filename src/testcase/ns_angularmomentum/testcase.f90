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
!> Module defining one specific testcase with all necessary variables
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

INTERFACE DefineParametersTestcase
  MODULE PROCEDURE DefineParametersTestcase
END INTERFACE

INTERFACE FinalizeTestcase
  MODULE PROCEDURE FinalizeTestcase
END INTERFACE

PUBLIC:: InitTestcase
PUBLIC:: DefineParametersTestcase
PUBLIC:: FinalizeTestcase

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
     "name of the testcase that is supposed to be used in simulation. Check if its = current testcase.","ns_angularmomentum")

END SUBROUTINE DefineParametersTestcase


!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE InitTestcase()
! MODULES
USE MOD_Globals
USE MOD_Testcase_Vars
USE MOD_ReadInTools,      ONLY: GETSTR
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE "NS_ANGULARMOMENTUM"...'
WhichTestcase=GETSTR('WhichTestcase','ns_angularmomentum') !not mandatory here

IF(INDEX(TRIM(WhichTestcase),'ns_angularmomentum').EQ.0)THEN
  CALL abort(__STAMP__,&
       "compiled with testcase ns_angularmomentum, but specified testcase is : "//TRIM(whichTestcase))
END IF


TestcaseInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitTestcase

!==================================================================================================================================
!> Finalize testcase data 
!==================================================================================================================================
SUBROUTINE FinalizeTestcase()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeTestcase

END MODULE MOD_Testcase
