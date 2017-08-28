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

!=================================================================================================================================
!> Contains the variables of your testcase that should be globally accessible!
!=================================================================================================================================
MODULE MOD_TestCase_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!---------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!---------------------------------------------------------------------------------------------------------------------------------
LOGICAL                  :: doTCPreTimeStep=.FALSE.  !< compute something before the timestep
LOGICAL                  :: doTCSource               !< compute source terms for testcase
CHARACTER(LEN=255)       :: whichTestcase            !< input variable, to be able to check if wanted testcase was compiled
LOGICAL                  :: EvalEquilibrium          !< switch for TC_exactfunc
LOGICAL                  :: doCalcErrorToEquilibrium !< switch for TC_analyze: compute difference of |U-Ueq|
LOGICAL                  :: doCalcDeltaBEnergy       !< switch for TC_analyze: compute Energy of 1/(2mu0) |B-Beq|^2
REAL                     :: deltaB_Energy            !< stored to compute growth rate 
LOGICAL                  :: doCalcAngularMomentum    !< switch for TC analyze: compute total Angular Momentum 
REAL                     :: RotationCenter(3)        !< center around which the angular momentum will be computed
INTEGER                  :: EquilibriumStateIni      !< =-1: Default: no equilibrium state used. U_t(U_eq)=0. 
                                                     !<      Sanity check if code is compiled with this testcase
                                                     !< = 0 : Use U_eq=exactFunc(IniExactFunc) for equilibrium state
                                                     !< > 0 : Use U_eq=exactFunc(EquilibriumStateIni) for equilibrium state
                                                     !< =-2 : Read U_eq from MeshFile
                                                     !< =-3 : Read U_eq from the solution of a stateFile

LOGICAL                  :: EquilibriumDivBcorr     !< switch to compute B from a vector potential instead of using B directly
INTEGER                  :: EquilibriumDisturbFunc  !< =0: Default: use same number as iniExactfunc for disturbance added to 
                                                    !< initialization state, else specific function
CHARACTER(LEN=255)       :: EquilibriumStateFile    !< if equilibrium is read from a state



LOGICAL                  :: EqBCexists               !< BCtype =21/29 exists: equilibrium state dirichlet BC
REAL,ALLOCATABLE         :: InputEq(:,:,:,:,:)       !< equilibrium data from Mesh
REAL,ALLOCATABLE         :: Ueq(:,:,:,:,:)           !< full equilibrium solution
REAL,ALLOCATABLE         :: Uteq(:,:,:,:,:)          !< DG time derivative of Ueq
LOGICAL                  :: TestcaseInitIsDone=.FALSE. !< Switch to check TestcaseInit status
!=================================================================================================================================

END MODULE MOD_TestCase_Vars
