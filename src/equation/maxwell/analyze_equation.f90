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
#include "defines.h"

!==================================================================================================================================
!> Contains analyze routines specific to the Maxwell's equation
!==================================================================================================================================
MODULE MOD_AnalyzeEquation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyzeEquation
  MODULE PROCEDURE InitAnalyzeEquation
END INTERFACE

INTERFACE AnalyzeEquation
  MODULE PROCEDURE AnalyzeEquation
END INTERFACE

INTERFACE FinalizeAnalyzeEquation
  MODULE PROCEDURE FinalizeAnalyzeEquation
END INTERFACE


PUBLIC :: AnalyzeEquation
PUBLIC :: InitAnalyzeEquation
PUBLIC :: FinalizeAnalyzeEquation
PUBLIC :: DefineParametersAnalyzeEquation
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%CreateLogicalOption('CalcEnergy', "Set true to compute the integrated electric and magnetic energy"&
           , '.FALSE.')
END SUBROUTINE DefineParametersAnalyzeEquation

!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars,       ONLY:doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_AnalyzeEquation_Vars
USE MOD_ReadInTools,        ONLY: GETLOGICAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Get the various analysis/output variables 
doCalcEnergy           = GETLOGICAL('CalcEnergy','.FALSE.')
IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doCalcEnergy)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"ElectricEnergy"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"MagneticEnergy"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"TotalEnergy"'
  END IF !doCalcEnergy
END IF !MPIroot & doAnalyzeToFile

! Initialize eval routines

END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,        ONLY:doAnalyzeToFile,A2F_iVar,A2F_data
USE MOD_AnalyzeEquation_Vars,ONLY:doCalcEnergy
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=40)    :: formatStr
REAL                 :: Energy(3) 
!==================================================================================================================================
! Calculate divergence 

IF(doCalcEnergy)THEN
  !store last energies
  CALL CalcPotentialEnergy(Energy)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A28,',3,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' elec./magn./total Energy : ',Energy(1:3)
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(1) !kineticEnergy
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(2) !magneticEnergy
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(3) !TotalEnergy
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF !doCalcEnergy

END SUBROUTINE AnalyzeEquation

!==================================================================================================================================
!> Calculates electric and magnetic Energy over whole domain (normalized with volume)
!==================================================================================================================================
SUBROUTINE CalcPotentialEnergy(Energy)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Equation_Vars,      ONLY: smu0,eps0
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Energy(3) !< kinetic ,magnetic and total energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k
REAL                            :: IntegrationWeight
!==================================================================================================================================
Energy=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem)
    Energy(1)  = Energy(1)+SUM(U(1:3,i,j,k,iElem)**2)*IntegrationWeight
    Energy(2)  = Energy(2)+SUM(U(4:6,i,j,k,iElem)**2)*IntegrationWeight
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Energy,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(Energy         ,0  ,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
Energy(1)=0.5*eps0*Energy(1)/vol
Energy(2)=0.5*smu0*Energy(2)/vol
Energy(3)=Energy(1)+Energy(2)


END SUBROUTINE CalcPotentialEnergy

!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
!USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
