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

INTERFACE CalcAngMom
  MODULE PROCEDURE CalcAngMom
END INTERFACE


PUBLIC:: DefineParametersAnalyzeTestcase 
PUBLIC:: InitAnalyzeTestCase
PUBLIC:: AnalyzeTestCase
PUBLIC:: CalcAngMom 

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
CALL prms%CreateLogicalOption('CalcErrorToEquilibrium', &
     "switch for TC_analyze: compute difference of |U-Ueq|" , '.FALSE.')
CALL prms%CreateLogicalOption('CalcDeltaBEnergy', &
     "switch for TC_analyze: compute Energy of 1/(2mu0)|B-Beq|^2" , '.FALSE.')
CALL prms%CreateLogicalOption('CalcAngularMomentum', &
     "switch for TC_analyze: compute total Angular momentum" , '.FALSE.')
CALL prms%CreateRealArrayOption('TC_RotationCenter', &
     "center around which the angular momentum will be computed.","0.,0.,0.")
END SUBROUTINE DefineParametersAnalyzeTestcase

!==================================================================================================================================
!> Initialize Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE InitAnalyzeTestcase()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,      ONLY: GETLOGICAL,GETREALARRAY
USE MOD_Analyze_Vars,     ONLY:doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_Testcase_Vars,    ONLY:doCalcErrorToEquilibrium,doCalcDeltaBEnergy
USE MOD_Testcase_Vars,    ONLY:doCalcAngularMomentum,RotationCenter
USE MOD_Equation_Vars,    ONLY:StrVarNames
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER    :: iVar
!==================================================================================================================================
!TESTCASE ANALYZE
doCalcErrorToEquilibrium = GETLOGICAL('CalcErrorToEquilibrium','.FALSE.')
doCalcDeltaBEnergy       = GETLOGICAL('CalcDeltaBEnergy','.FALSE.')
doCalcAngularMomentum    = GETLOGICAL('CalcAngularMomentum','.FALSE.')
RotationCenter=GETREALARRAY('TC_RotationCenter',3,'0.0,0.0,0.0')
!prepare AnalyzeToFile
IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doCalcErrorToEquilibrium)THEN
    DO iVar=1,PP_nVar
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"L2_eq_'//TRIM(StrVarNames(iVar))//'"'
    END DO
    DO iVar=1,PP_nVar
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"Linf_eq_'//TRIM(StrVarNames(iVar))//'"'
    END DO
  END IF !CalcErrorToEquilibrium
  IF(doCalcDeltaBEnergy)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"DeltaB_Energy"'
  END IF !CalcErrorToEquilibrium
  IF(doCalcAngularMomentum)THEN
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
  END IF !CalcErrorToEquilibrium
END IF !MPIroot & doAnalyzeToFile
END SUBROUTINE InitAnalyzeTestcase

!==================================================================================================================================
!> Testcase specific analyze routines
!==================================================================================================================================
SUBROUTINE AnalyzeTestcase(Time)
! MODULES
USE MOD_Globals
USE MOD_Analyze_Vars,     ONLY:doAnalyzeToFile,A2F_iVar,A2F_Data,Analyze_dt
USE MOD_Testcase_Vars,    ONLY:doCalcErrorToEquilibrium,doCalcDeltaBEnergy,deltaB_Energy
USE MOD_Testcase_Vars,    ONLY:doCalcAngularMomentum
USE MOD_Restart_Vars,     ONLY:RestartTime
USE MOD_DG_Vars,          ONLY: U,Ut
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: time         !< current simulation time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)               :: formatStr
REAL                            :: L2_eq_Error(PP_nVar),Linf_eq_Error(PP_nVar)
REAL                            :: tmp,AngMom(3),AngMom_t(3)
!==================================================================================================================================
IF(doCalcErrorToEquilibrium)THEN
  CALL CalcErrorToEquilibrium(L2_eq_Error,Linf_eq_Error)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',PP_nVar,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' L2_eq      : ',L2_eq_Error
    WRITE(UNIT_StdOut,formatStr)' Linf_eq    : ',Linf_eq_Error
    IF(doAnalyzeToFile)THEN
      A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=L2_eq_Error(:)
      A2F_iVar=A2F_iVar+PP_nVar
      A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=Linf_eq_Error(:)
      A2F_iVar=A2F_iVar+PP_nVar
    END IF !doAnalyzeToFile  
  END IF !MPIroot
END IF
IF(doCalcDeltaBEnergy)THEN
  tmp=deltaB_Energy
  CALL CalcDeltaBEnergy(deltaB_Energy)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A15,',1,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' deltaBEnergy: ',deltaB_Energy
    IF((time-RestartTime).GE.Analyze_dt)THEN
      tmp=LOG(deltaB_Energy /tmp)/Analyze_dt
      WRITE(UNIT_StdOut,formatStr)' growth rate : ',tmp
    END IF !time>Analyze_dt
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=deltaB_Energy
    END IF !doAnalyzeToFile  
  END IF !MPIroot
END IF
IF(doCalcAngularMomentum)THEN
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
END IF !doCalcAngularMomentum
END SUBROUTINE AnalyzeTestcase


!==================================================================================================================================
!> Calculates L_infinfity and L_2 norms of U-Ueq 
!==================================================================================================================================
SUBROUTINE CalcErrorToEquilibrium(L2_eq_Error,Linf_eq_Error)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars       ,ONLY: sJ,nElems
USE MOD_DG_Vars         ,ONLY: U
USE MOD_Testcase_Vars   ,ONLY: Ueq
USE MOD_ChangeBasis     ,ONLY: ChangeBasis3D
USE MOD_Analyze_Vars    ,ONLY: NAnalyze,Vdm_GaussN_NAnalyze
USE MOD_Analyze_Vars    ,ONLY: wGPVolAnalyze,Vol
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: L2_eq_Error(PP_nVar),Linf_eq_Error(PP_nVar) !< L2 and Linf difference of U-Ueq
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,k,l,m
REAL                            :: U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: Ueq_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                            :: IntegrationWeight
!==================================================================================================================================
! Calculate error norms
Linf_eq_Error(:)=-1.E10
L2_eq_Error(:)=0.
! Interpolate values of Error-Grid from GP's
DO iElem=1,nElems
   ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the jacobian ;-)
   J_N(1,:,:,:)=1./sJ(:,:,:,iElem)
   CALL ChangeBasis3D(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_NAnalyze(1:1,:,:,:))
   ! Interpolate the solution to the analyze grid
   CALL ChangeBasis3D(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Ueq(1:PP_nVar,:,:,:,iElem),Ueq_NAnalyze(1:PP_nVar,:,:,:))
   CALL ChangeBasis3D(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,U(1:PP_nVar,:,:,:,iElem),U_NAnalyze(1:PP_nVar,:,:,:))
   DO m=0,NAnalyze
     DO l=0,NAnalyze
       DO k=0,NAnalyze
         Linf_eq_Error = MAX(Linf_eq_Error,ABS(U_NAnalyze(:,k,l,m) - Ueq_NAnalyze(:,k,l,m)))
         IntegrationWeight = wGPVolAnalyze(k,l,m)*J_NAnalyze(1,k,l,m)
         ! To sum over the elements, We compute here the square of the L_2 error
         L2_eq_Error = L2_eq_Error+((U_NAnalyze(:,k,l,m) - Ueq_NAnalyze(:,k,l,m))**2)*IntegrationWeight
       END DO ! k
     END DO ! l
   END DO ! m
END DO ! iElem=1,nElems

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,L2_eq_Error  ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,Linf_eq_Error,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(L2_eq_Error  ,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(Linf_eq_Error,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L2_eq_Error = SQRT(L2_eq_Error/Vol)

END SUBROUTINE CalcErrorToEquilibrium

!==================================================================================================================================
!> Calculates  magnetic Energy of B-Beq over whole domain (normalized with volume)
!==================================================================================================================================
SUBROUTINE CalcDeltaBEnergy(dbEnergy)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Testcase_Vars,      ONLY: Ueq
USE MOD_Equation_Vars,      ONLY: s2mu_0
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: dbEnergy !< mean magnetic energy  1/(2mu0)*|B-Beq|^2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k
REAL                            :: IntegrationWeight
!==================================================================================================================================
dbEnergy=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem)
    dbEnergy  = dbEnergy+SUM((U(6:8,i,j,k,iElem)-Ueq(6:8,i,j,k,iElem))**2)*IntegrationWeight
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,dbEnergy,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(dbEnergy         ,0  ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
dbEnergy= s2mu_0*dbEnergy/vol
END SUBROUTINE CalcDeltaBEnergy

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
