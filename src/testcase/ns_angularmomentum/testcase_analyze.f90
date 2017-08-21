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


PUBLIC:: InitAnalyzeTestCase
PUBLIC:: AnalyzeTestCase

CONTAINS


!==================================================================================================================================
!> Initialize Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE InitAnalyzeTestcase()
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars,     ONLY:doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_Testcase_Vars,    ONLY:doTCanalyze
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
!prepare AnalyzeToFile
IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doTCanalyze)THEN
WRITE(*,*) 'INIT ANALYZE TESTCASE'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentumX"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentumY"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"AngularMomentumZ"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"DiffAngularMomentumX"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"DiffAngularMomentumY"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"DiffAngularMomentumZ"'
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
USE MOD_Testcase_Vars,    ONLY:doTCanalyze,AngMomInit
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: time         !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: AngMom(3)
!==================================================================================================================================
IF(doTCanalyze)THEN
  CALL CalcAngMom(AngMom)
  IF(MPIroot) THEN
    IF(time.LT.1e-08)THEN
      AngMomInit=AngMom
    END IF
    WRITE(UNIT_StdOut,'(A, E21.13)')'Total Angular Momentum   : ', SQRT(SUM(AngMom**2))
    WRITE(UNIT_StdOut,'(A,3E21.13)')'Angular Momentum (X,Y,Z) : ', AngMom(:)
    WRITE(UNIT_StdOut,'(A,3E21.13)')'Difference to init       : ', AngMom(:)-AngMomInit(:)
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=AngMom(:)
      A2F_iVar=A2F_iVar+3
      A2F_Data(A2F_iVar-2:A2F_iVar)=AngMom(:)-AngMomInit(:)
    END IF
  END IF !MPIroot & doAnalyzeToFile
END IF  !doTCanalyze
END SUBROUTINE AnalyzeTestcase

!==================================================================================================================================
!> Calculate the angular Momentum 
!==================================================================================================================================
SUBROUTINE CalcAngMom(AngMom)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems,Elem_xGP
USE MOD_DG_Vars,            ONLY: U
USE MOD_Testcase_Vars,      ONLY: RotationCenter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
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
    AngMom(:)  = AngMom(:)+CROSS(r(:),U(2:4,i,j,k,iElem))*IntegrationWeight
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
