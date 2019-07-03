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
!=================================================================================================================================
#include "defines.h"

!=================================================================================================================================
!> Module to handle commandline arguments 
!=================================================================================================================================
MODULE MOD_Commandline_Arguments
IMPLICIT NONE

! Global variables for command line argument parsing
INTEGER                              :: nArgs              ! number of command line argumens

INTERFACE ParseCommandlineArguments
  MODULE PROCEDURE ParseCommandlineArguments
END INTERFACE ParseCommandlineArguments

!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Reads all commandline arguments
!==================================================================================================================================
SUBROUTINE ParseCommandlineArguments()
! MODULES
USE MOD_Globals
USE MOD_StringTools     ,ONLY: STRICMP
USE MOD_Restart_Vars    ,ONLY: RestartFile
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iArg                             
CHARACTER(LEN=255)      :: tmp
LOGICAL,ALLOCATABLE     :: alreadyRead(:)
!==================================================================================================================================
! Get number of command line arguments
nArgs=COMMAND_ARGUMENT_COUNT()
ALLOCATE(alreadyRead(nArgs))
alreadyRead = .FALSE.

! Get keyword arguments (arbitrary order)
doGenerateUnittestReferenceData = .FALSE.
doPrintHelp = 0
DO iArg = 1, nArgs
  CALL GET_COMMAND_ARGUMENT(iArg,tmp)
  IF (STRICMP(tmp, "--generateUnittestReferenceData")) THEN
    doGenerateUnittestReferenceData = .TRUE.
    alreadyRead(iArg) = .TRUE.
  END IF
  IF (STRICMP(tmp, "--help").OR.STRICMP(tmp,"-h")) THEN
    doPrintHelp = 1
    alreadyRead(iArg) = .TRUE.
  END IF
  IF (STRICMP(tmp, "--markdown")) THEN
    doPrintHelp = 2
    alreadyRead(iArg) = .TRUE.
  END IF
END DO ! iArg = 1, nArgs

! Get parameter file
ParameterFile = ""
DO iArg = 1,nArgs
  IF (.NOT.alreadyRead(iArg)) THEN
    CALL GET_COMMAND_ARGUMENT(iArg,ParameterFile)
    alreadyRead(iArg) = .TRUE.
    EXIT
  END IF
END DO

! Get restart file 
RestartFile = ""
DO iArg = 1,nArgs
  IF (.NOT.alreadyRead(iArg)) THEN
    CALL GET_COMMAND_ARGUMENT(iArg,RestartFile)
    alreadyRead(iArg) = .TRUE.
    EXIT
  END IF
END DO

END SUBROUTINE ParseCommandlineArguments

END MODULE MOD_Commandline_Arguments
