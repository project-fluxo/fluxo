!==================================================================================================================================
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2017 - 2019 Marvin Bohm
! Copyright (c) 2020 - 2020 Andrés Rueda
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
!> Contains global variables used by the positivitiy preservation modules.
!==================================================================================================================================

MODULE MOD_PositivityPreservation_Vars
!===================================================================================================================================
! Contains global variables used by the Positivity Preservation modul§
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: PositivityPreservationInitIsDone = .FALSE.
LOGICAL           :: PositivityPreservationInitFirst  = .TRUE.
real              :: PositivityEpsilon
real              :: PositivityDelta
real, allocatable :: Jac(:,:,:,:)  ! det of jacobian.
real, allocatable :: vol      (:)  ! Volume of cell.
!===================================================================================================================================

END MODULE MOD_PositivityPreservation_Vars
