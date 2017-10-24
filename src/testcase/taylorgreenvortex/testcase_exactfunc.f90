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
MODULE MOD_Testcase_ExactFunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE TestcaseExactFunc
  MODULE PROCEDURE TestcaseExactFunc
END INTERFACE

PUBLIC:: TestcaseExactFunc

CONTAINS

!==================================================================================================================================
!> Specifies all the initial conditions. For TGV the ExactFunc is already available in equation.f90
!==================================================================================================================================
SUBROUTINE TestcaseExactFunc(ExactFunction,t,x,resu,resu_t,resu_tt) 
! MODULES
USE MOD_Globals,      ONLY: Abort
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction    !< determines the exact function
REAL,INTENT(IN)                 :: t                !< current simulation time
REAL,INTENT(IN)                 :: x(3)             !< position in physical coordinates
REAL,INTENT(OUT)                :: resu(PP_nVar)    !< exact fuction evaluated at tIn, returning state in conservative variables
REAL,INTENT(OUT)                :: resu_t(PP_nVar)  !< first time deriv of exact fuction
REAL,INTENT(OUT)                :: resu_tt(PP_nVar) !< second time deriv of exact fuction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL abort(__STAMP__,'Exactfunction not specified!')
Resu=-1.
Resu_t =-1.
Resu_tt=-1.
END SUBROUTINE TestcaseExactFunc

END MODULE MOD_Testcase_ExactFunc
