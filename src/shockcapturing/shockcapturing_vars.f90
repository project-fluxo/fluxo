!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
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
#include "defines.h"

!==================================================================================================================================
!> Contains global variables used by the shock capturing modules.
!==================================================================================================================================
MODULE MOD_ShockCapturing_Vars
!===================================================================================================================================
! MODULES
use MOD_Indicators, only: Indicator_PerssonPeraire
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: ShockCapturingInitIsDone=.FALSE. !< marks whether the shock capturing init routine is complete
integer           :: ShockIndicator                   !< Identifier for the shock indicator
character(len=255):: ShockIndicatorQuantity           !< Shock indicator quantity to be used
! The shock indicator itself:
type(Indicator_PerssonPeraire) :: Shock_Indicator

!===================================================================================================================================

END MODULE MOD_ShockCapturing_Vars
