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

!==================================================================================================================================
!> Contains the variables of your testcase that should be globally accessible!
!==================================================================================================================================
MODULE MOD_TestCase_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL            :: doTCPreTimeStep=.FALSE.    !< compute something before the timestep
LOGICAL            :: doTCSource     =.FALSE.    !< compute source terms for testcase
LOGICAL            :: doTCanalyze                !< switch on analyze routines
CHARACTER(LEN=255) :: whichTestcase              !< input variable, to be able to check if wanted testcase was compiled
REAL               :: last_Ekin_comp             !< used in analyze to compute -dEkin/dt
REAL               :: last_Emag_comp             !< used in analyze to compute -dEmag/dt
LOGICAL            :: TestcaseInitIsDone=.FALSE. !< Switch to check TestcaseInit status
!==================================================================================================================================

END MODULE MOD_TestCase_Vars
