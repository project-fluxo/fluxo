!==================================================================================================================================
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
!> Module for the shock capturing routines
!==================================================================================================================================
MODULE MOD_ShockCapturing
#if SHOCKCAPTURE
! MODULES
IMPLICIT NONE
PRIVATE
! ----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersShockCapturing
   MODULE PROCEDURE DefineParametersShockCapturing
END INTERFACE

INTERFACE InitShockCapturing
   MODULE PROCEDURE InitShockCapturing
END INTERFACE

INTERFACE FinalizeShockCapturing
   MODULE PROCEDURE FinalizeShockCapturing
END INTERFACE

PUBLIC :: DefineParametersShockCapturing
PUBLIC :: InitShockCapturing
PUBLIC :: FinalizeShockCapturing
public :: InitShockCapturingAfterAdapt
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersShockCapturing()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
#if SHOCK_NFVSE
use MOD_NFVSE       ,only: DefineParametersNFVSE
#endif /*SHOCK_NFVSE*/
use MOD_Equation_Vars, only: IndicatorQuantityNames, nIndVar
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
integer            :: i
character(len=255) :: IndicatorQuantities,fmt
!==================================================================================================================================
CALL prms%SetSection("ShockCapturing")

CALL prms%CreateIntOption(     "ShockIndicator",  " Specifies the quantity to be used as shock-indicator:\n"//&
                                              "   1: Persson-Peraire indicator\n"&
                                             ,"1")

write(fmt,'(A,I0,A)') '(',nIndVar,'A)'
write(IndicatorQuantities,fmt) ('  * '//trim(IndicatorQuantityNames(i))//'\n', i=1, nIndVar)

CALL prms%CreateStringOption("ShockIndicatorQuantity",  " Specifies the quantity to be used for the shock indicator. One of the following:\n"//&
                                              trim(IndicatorQuantities)&
                                             ,"DensityTimesPressure")

#if SHOCK_NFVSE
call DefineParametersNFVSE()
#endif /*SHOCK_NFVSE*/

END SUBROUTINE DefineParametersShockCapturing

SUBROUTINE InitShockCapturing()
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ShockCapturing_Vars
USE MOD_ReadInTools
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone
#if SHOCK_NFVSE
use MOD_NFVSE             , only: InitNFVSE
#endif /*SHOCK_NFVSE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!============================================================================================================================

IF (ShockCapturingInitIsDone.OR.(.NOT.InterpolationInitIsDone)) THEN
  SWRITE(*,*) "InitShockCapturing not ready to be called or already called."
  RETURN
END IF
IF (PP_N.LT.2) THEN
  CALL abort(__STAMP__,'Polynomial Degree too small for Shock Capturing!',999,999.)
  RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SHOCKCAPTURING...'

ShockIndicator         = GETINT('ShockIndicator','1')
ShockIndicatorQuantity = GETSTR('ShockIndicatorQuantity','DensityTimesPressure')

#if SHOCK_NFVSE
call InitNFVSE()
#endif /*SHOCK_NFVSE*/

if (ShockIndicator==1) then
  call Shock_Indicator % construct(ShockIndicatorQuantity)
else
  CALL abort(__STAMP__,'Shock indicator not defined!',999,999.)
  RETURN
end if

ShockCapturingInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SHOCKCAPTURING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitShockCapturing
!===================================================================================================================================
!> Reinitializes all variables that need reinitialization after the h-adaptation
!===================================================================================================================================
SUBROUTINE InitShockCapturingAfterAdapt(ChangeElem,nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld)
  use MOD_NFVSE    , only: InitNFVSEAfterAdaptation
  use MOD_Mesh_Vars, only: nElems
  implicit none
  !-arguments-----------------------------------
  integer, intent(in) :: ChangeElem(8,nElems)
  integer, intent(in) :: nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld
  !---------------------------------------------
  
#if SHOCK_NFVSE
  call InitNFVSEAfterAdaptation(ChangeElem,nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld)
#endif /*SHOCK_NFVSE*/
END SUBROUTINE InitShockCapturingAfterAdapt
  
SUBROUTINE FinalizeShockCapturing()
!============================================================================================================================
!> Deallocate all global shock capturing variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ShockCapturing_Vars
#if SHOCK_NFVSE
use MOD_NFVSE             , only: FinalizeNFVSE
#endif /*SHOCK_NFVSE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!============================================================================================================================
IF (.NOT.ShockCapturingInitIsDone) THEN
  WRITE(UNIT_stdOut,*) "InitShockCapturing was not called before."
  RETURN
END IF
ShockCapturingInitIsDone = .FALSE.
call Shock_Indicator % destruct
#if SHOCK_NFVSE
call FinalizeNFVSE()
#endif /*SHOCK_NFVSE*/

END SUBROUTINE FinalizeShockCapturing

#endif /*SHOCKCAPTURE*/
END MODULE MOD_ShockCapturing
