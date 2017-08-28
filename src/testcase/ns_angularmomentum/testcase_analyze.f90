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
#include "defines.h"


!==================================================================================================================================
!> Subroutines defining one specific testcase with all necessary variables
!==================================================================================================================================
MODULE MOD_Testcase_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitAnalyzeTestCase
  MODULE PROCEDURE InitAnalyzeTestCase
END INTERFACE

INTERFACE AnalyzeTestCase
  MODULE PROCEDURE AnalyzeTestCase
END INTERFACE


PUBLIC:: DefineParametersAnalyzeTestcase 
PUBLIC:: InitAnalyzeTestCase
PUBLIC:: AnalyzeTestCase

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeTestcase()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("AnalyzeTestcase")
CALL prms%CreateLogicalOption('doTCanalyze', &
     "switch off/on TC_analyze" , '.FALSE.')
CALL prms%CreateRealArrayOption('TC_RotationCenter'   , "center around which the angular momentum will be computed.","0.,0.,0.")
END SUBROUTINE DefineParametersAnalyzeTestcase


!==================================================================================================================================
!> Initialize Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE InitAnalyzeTestcase()
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars,     ONLY:doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_Testcase_Vars,    ONLY:doTCanalyze,RotationCenter
USE MOD_ReadInTools,      ONLY: GETLOGICAL,GETREALARRAY
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

doTCanalyze = GETLOGICAL('doTCanalyze','.TRUE.')

RotationCenter=GETREALARRAY('TC_RotationCenter',3,'0.0,0.0,0.0')
!prepare AnalyzeToFile
IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doTCanalyze)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentumX"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentumY"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentumZ"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentum_t_X"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentum_t_Y"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentum_t_Z"'
  END IF  !doTCanalyze
END IF !MPIroot & doAnalyzeToFile
END SUBROUTINE InitAnalyzeTestcase

!==================================================================================================================================
!> Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE AnalyzeTestcase(Time)
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars,     ONLY:doAnalyzeToFile,A2F_iVar,A2F_Data
USE MOD_Testcase_Vars,    ONLY:doTCanalyze
USE MOD_DG_Vars,          ONLY: U,Ut
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: time         !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: AngMom(3),AngMom_t(3)
!==================================================================================================================================
IF(doTCanalyze)THEN
  CALL CalcAngMom(AngMom,U)
  CALL CalcAngMom(AngMom_t,Ut)
  IF(MPIroot) THEN
    WRITE(UNIT_StdOut,'(A, E21.13)')' Total Angular Momentum     : ', SQRT(SUM(AngMom**2))
    WRITE(UNIT_StdOut,'(A,3E21.13)')' Angular Momentum (X,Y,Z)   : ', AngMom(:)
    WRITE(UNIT_StdOut,'(A,3E21.13)')' Angular Momentum_t (X,Y,Z) : ', AngMom_t(:) 
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=AngMom(:)
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=AngMom_t(:)
    END IF !doAnalyzeToFile
  END IF !MPIroot 
END IF  !doTCanalyze
END SUBROUTINE AnalyzeTestcase

!==================================================================================================================================
!> Calculate the angular Momentum 
!==================================================================================================================================
SUBROUTINE CalcAngMom(AngMom,U_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems,Elem_xGP
USE MOD_Testcase_Vars,      ONLY: RotationCenter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: U_in(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: AngMom(3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k
REAL                            :: IntegrationWeight,r(3)
!==================================================================================================================================
AngMom=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem)
    r=Elem_xGP(:,i,j,k,iElem)-RotationCenter(:)
    AngMom(:)  = AngMom(:)+CROSS(r(:),U_in(2:4,i,j,k,iElem))*IntegrationWeight
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,AngMom,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(AngMom         ,0  ,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
END SUBROUTINE CalcAngMom

END MODULE MOD_Testcase_Analyze
