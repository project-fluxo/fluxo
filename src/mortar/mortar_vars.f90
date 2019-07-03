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
!> Variables used for mortars: mortar interpolation and projection matrices
!==================================================================================================================================
MODULE MOD_Mortar_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!> 1D-Mortar Operator: interpolation full interval 0: [-1,1] to left interval 1: [-1,0] and right intervall 2: [0,1]
REAL,ALLOCATABLE,TARGET :: M_0_1(:,:),M_0_2(:,:)
!> 1D-Mortar Operator: projection left interval 1: [-1,0] and right intervall 2: [0,1] to full intervall 0: [-1,1]
REAL,ALLOCATABLE,TARGET :: M_1_0(:,:),M_2_0(:,:)
LOGICAL                 :: MortarInitIsDone=.FALSE. !< marks whether mortar init routines are complete
!==================================================================================================================================
END MODULE MOD_Mortar_Vars
