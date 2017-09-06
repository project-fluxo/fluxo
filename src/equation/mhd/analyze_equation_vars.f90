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
!> Contains parameters for the analyze equation module
!==================================================================================================================================
MODULE MOD_AnalyzeEquation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
! Analyze features
LOGICAL              :: doCalcDivergence   !<handle for analyze: compute divergence of B
LOGICAL              :: doCalcBulk         !<handle for analyze: compute Bulk integral of all variables
LOGICAL              :: doCalcEnergy       !<handle for analyze: compute kinetic and magnetic energy
LOGICAL              :: doCalcEntropy      !<handle for analyze: compute entropy
REAL                 :: Energy(2)          !< store kinetic and magnetic energy to compute growth rate to 
                                           !< last analyze step 
REAL                 :: Entropy            !< store entropy to compute change to last step 

!==================================================================================================================================
END MODULE MOD_AnalyzeEquation_Vars
