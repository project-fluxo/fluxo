!==================================================================================================================================
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

!==================================================================================================================================
!> Contains global variables used by the restart module
!==================================================================================================================================
MODULE MOD_Restart_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER            :: nVar_Restart                  !< number of variables in restart file
INTEGER            :: N_Restart = 0                 !< polynomial degree of restart solution
INTEGER            :: nElems_Restart                !< number of elements in restart file
LOGICAl            :: RestartInitIsDone   = .FALSE. !< flag if restart routines are finished
LOGICAl            :: DoRestart           = .FALSE. !< flag whether a restart should actually be performed
LOGICAL            :: InterpolateSolution = .FALSE. !< flag whether restart solution should be interpolated
                                                    !< if node type or polynomial degree are different
CHARACTER(LEN=300) :: RestartFile                   !< name of restart file
CHARACTER(LEN=255) :: NodeType_Restart              !< node type of restart file
REAL               :: RestartTime                   !< time at which computation is resumed
!==================================================================================================================================
END MODULE MOD_Restart_Vars
