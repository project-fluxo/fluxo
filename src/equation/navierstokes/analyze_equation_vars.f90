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

!==================================================================================================================================
!> Contains global variables used by the equation specific analyze modules.
!==================================================================================================================================
MODULE MOD_AnalyzeEquation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
! precomputed variables 

! Analyze features
LOGICAL              :: doCalcBulkVelocity
LOGICAL              :: doCalcWallVelocity
LOGICAL              :: doCalcBodyForces
LOGICAL              :: doCalcEntropy

! Variables for the specific analyze routines
! WallVelocity
REAL,ALLOCATABLE     :: meanV(:),minV(:),maxV(:)    ! size nBCs
REAL                 :: Entropy
!==================================================================================================================================
END MODULE MOD_AnalyzeEquation_Vars
