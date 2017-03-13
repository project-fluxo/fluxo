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
!> Contains parameters for the linear scalar advection equation
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: doCalcSource            !< logical to define if a source term (e.g. exactfunc) is added
REAL                :: AdvVel(3)               !< Advection velocity
REAL                :: DiffC                   !< Diffusion constant
REAL                :: IniWavenumber(3)        !< wavenumbers in 3 directions (sinus periodic with exactfunc=6)
INTEGER             :: IniExactFunc            !< Exact Function for initialization
INTEGER,PARAMETER   :: nAuxVar=1               !< Auxiliary variable for DiscType=3
INTEGER,ALLOCATABLE :: nBCByType(:)            !< Number of sides for each boundary
INTEGER,ALLOCATABLE :: BCSideID(:,:)           !< SideIDs for BC types

CHARACTER(LEN=255),PARAMETER :: StrVarNames(1)=(/'Solution'/) !< Variable names for output

LOGICAL           :: EquationInitIsDone=.FALSE. !< Init switch  
#if (PP_DiscType==2)
!procedure pointers for split form DG
INTEGER             :: WhichVolumeFlux          !< for split-form DG, two-point average flux
PROCEDURE(),POINTER :: VolumeFluxAverageVec     !< procedure pointer to two-point average flux
#endif /*PP_DiscType==2*/

!==================================================================================================================================
END MODULE MOD_Equation_Vars
