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
!> Contains analyze routines specific to the magnetohydrodynamic equations
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


PUBLIC:: AnalyzeEquation
PUBLIC:: InitAnalyzeEquation
PUBLIC:: FinalizeAnalyzeEquation
PUBLIC:: DefineParametersAnalyzeEquation
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for analyze MHD 
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
CALL prms%SetSection("AnalyzeEquation")
CALL prms%CreateLogicalOption('CalcDivergence'   , "Set true to compute the current divergence of the magnetic field" &
           , '.FALSE.')
CALL prms%CreateLogicalOption('CalcBulk',"Set true to compute the integrated mean value of each state variable over whole domain"&
           , '.FALSE.')
CALL prms%CreateLogicalOption('CalcEnergy', "Set true to compute the integrated kinetic and full and disturbed magnetic energy"&
           , '.FALSE.')
CALL prms%CreateLogicalOption('CalcEntropy', "Set true to compute the integrated entropy"&
           , '.FALSE.')
END SUBROUTINE DefineParametersAnalyzeEquation


!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars,      ONLY:StrVarNames
USE MOD_Analyze_Vars,       ONLY:doAnalyzeToFile,A2F_iVar,A2F_VarNames
USE MOD_AnalyzeEquation_Vars
USE MOD_ReadInTools,        ONLY: GETLOGICAL,GETINT,GETINTARRAY,GETREAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iVar

!==================================================================================================================================
! Get the various analysis/output variables 
doCalcDivergence       = GETLOGICAL('CalcDivergence','.FALSE.')
doCalcBulk             = GETLOGICAL('CalcBulk'      ,'.FALSE.')
doCalcEnergy           = GETLOGICAL('CalcEnergy'    ,'.FALSE.')
doCalcEntropy          = GETLOGICAL('CalcEntropy'   ,'.FALSE.')
IF(doCalcEnergy)THEN
END IF

IF(MPIroot.AND.doAnalyzeToFile) THEN
  IF(doCalcDivergence)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"L2_divB"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"Linf_divB"'
  END IF  !doCalcDivergence
  IF(doCalcBulk)THEN
    DO iVar=1,PP_nVar
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"Bulk_'//TRIM(StrVarNames(iVar))//'"'
    END DO !iVar 
  END IF ! doCalcBulk
  IF(doCalcEnergy)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"KineticEnergy"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"MagneticEnergy"'
  END IF !doCalcEnergy
  IF(doCalcEntropy)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"Entropy"'
  END IF !doCalcEntropy
END IF !MPIroot & doAnalyzeToFile


END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> execute the analyze steps 
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,        ONLY:doAnalyzeToFile,A2F_iVar,A2F_data
USE MOD_Analyze_Vars,        ONLY:Analyze_dt
USE MOD_Restart_Vars,        ONLY:RestartTime
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=40)    :: formatStr
REAL                 :: L2_divB,Linf_divB
REAL                 :: bulk(1:PP_nVar)
REAL                 :: tmp(2) 
!==================================================================================================================================
! Calculate divergence 
IF(doCalcDivergence)THEN
  CALL CalcDivergence(L2_divB,Linf_divB)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A21,',2,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' L2,Linf divB : ',L2_divB,Linf_divB
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=L2_divB
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Linf_divB
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF  !(doCalcDivergence)

IF(doCalcBulk)THEN
  CALL CalcBulk(bulk)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A21,',PP_nVar,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' Bulk : ',bulk(:)
    IF(doAnalyzeToFile)THEN
      A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=bulk(:)
      A2F_iVar=A2F_iVar+PP_nVar
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF

IF(doCalcEnergy)THEN
  !store last energies
  tmp(1:2)=Energy(1:2)
  CALL CalcEnergy(Energy)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A21,',3,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' kin./magn. Energy : ',Energy(1:2)
    IF((time-RestartTime).GE.Analyze_dt)THEN
      tmp(1:2)=LOG(Energy(1:2) /tmp(1:2))/Analyze_dt
      WRITE(UNIT_StdOut,formatStr)' growth rates    : ',tmp(1:2)
    END IF !time>Analyze_dt
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(1) !kineticEnergy
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(2) !magneticEnergy
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF !doCalcEnergy

IF(doCalcEntropy)THEN
  tmp(1)=Entropy
  CALL CalcEntropy(Entropy)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A)')'(A21,ES21.12)'
    WRITE(UNIT_StdOut,formatStr)  ' Entropy      : ',Entropy
    IF((time-RestartTime).GE.Analyze_dt)THEN
      WRITE(UNIT_StdOut,formatStr)' dEntropy/dt  : ',(Entropy-tmp(1))/Analyze_dt
    END IF !time>Analyze_dt
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Entropy
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF !doCalcEnergy
END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Calculates the divergence of the magnetic field
!==================================================================================================================================
SUBROUTINE CalcDivergence(L2_divB,Linf_divB)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_Mesh_Vars          ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde   ! metrics
USE MOD_DG_Vars            ,ONLY:U,D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: L2_divB
REAL,INTENT(OUT)                :: Linf_divB
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k,l
REAL,DIMENSION(1:3)             :: gradB_xi,gradB_eta,gradB_zeta  
REAL                            :: divB_loc
!==================================================================================================================================
! Needed for the computation of the forcing term for the channel flow
L2_divB=0.
Linf_divB=0.
DO iElem=1,nElems
  ! Compute the gradient in the reference system
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    gradB_xi  = 0.
    gradB_eta = 0.
    gradB_zeta= 0.
    DO l=0,N
      gradB_xi  (1:3)= gradB_xi  (1:3) + D(i,l) * U(6:8,l,j,k,iElem)
      gradB_eta (1:3)= gradB_eta (1:3) + D(j,l) * U(6:8,i,l,k,iElem)
      gradB_zeta(1:3)= gradB_zeta(1:3) + D(k,l) * U(6:8,i,j,l,iElem)
    END DO ! l 
    divB_loc = sJ(i,j,k,iElem) * (                                   &   
               Metrics_fTilde(1,i,j,k,iElem) * gradB_xi  (1) + & !gradUx(1 
               Metrics_gTilde(1,i,j,k,iElem) * gradB_eta (1) + & 
               Metrics_hTilde(1,i,j,k,iElem) * gradB_zeta(1) + &
               Metrics_fTilde(2,i,j,k,iElem) * gradB_xi  (2) + & !gradUy(2 
               Metrics_gTilde(2,i,j,k,iElem) * gradB_eta (2) + & 
               Metrics_hTilde(2,i,j,k,iElem) * gradB_zeta(2) + &   
               Metrics_fTilde(3,i,j,k,iElem) * gradB_xi  (3) + & !gradUz(3 
               Metrics_gTilde(3,i,j,k,iElem) * gradB_eta (3) + & 
               Metrics_hTilde(3,i,j,k,iElem) * gradB_zeta(3)   ) 
    L2_divB   = L2_divB+(divB_loc)**2*wGPVol(i,j,k)/sJ(i,j,k,iElem)
    Linf_divB = MAX(Linf_divB,ABS(divB_loc))
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,L2_divB,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(L2_divB        ,0  ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Linf_divB,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(Linf_divB     ,0  ,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF
#endif

L2_divB=SQRT(L2_divB/Vol)

END SUBROUTINE CalcDivergence 


!==================================================================================================================================
!> Calculates bulk velocities over whole domain
!==================================================================================================================================
SUBROUTINE CalcBulk(Bulk)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Bulk(1:PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k
#if MPI
REAL                            :: box(1:PP_nVar)
#endif
!==================================================================================================================================
Bulk=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    Bulk =Bulk+U(:,i,j,k,iElem)*wGPVol(i,j,k)/sJ(i,j,k,iElem)
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
Box=Bulk
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,box,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  Bulk=Box
ELSE
  CALL MPI_REDUCE(Box         ,0  ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif

Bulk    =Bulk/Vol

END SUBROUTINE CalcBulk


!==================================================================================================================================
!> Calculates kinetic and magnetic Energy over whole domain (normalized with volume)
!==================================================================================================================================
SUBROUTINE CalcEnergy(Energy)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Equation_Vars,      ONLY: s2mu_0
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Energy(2) !< kinetic and magnetic energy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k
REAL                            :: IntegrationWeight
!==================================================================================================================================
Energy=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem)
    Energy(1)  = Energy(1)+SUM(U(2:4,i,j,k,iElem)**2)/U(1,i,j,k,iElem)*IntegrationWeight
    Energy(2)  = Energy(2)+SUM(U(6:8,i,j,k,iElem)**2)*IntegrationWeight
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Energy,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(Energy         ,0  ,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
Energy(1)=0.5*Energy(1)/vol
Energy(2)= s2mu_0*Energy(2)/vol


END SUBROUTINE CalcEnergy
!==================================================================================================================================
!> Calculates  Entropy over whole domain Entropy=rho*s/(kappa-1), s=ln(p rho^(-kappa))=ln(p)-kappa*ln(rho) 
!==================================================================================================================================
SUBROUTINE CalcEntropy(Entropy)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Equation_Vars,      ONLY: kappa,sKappaM1,ConsToPrim
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Entropy !< Entropy
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k
REAL                            :: ent_loc,prim(PP_nVar)
!==================================================================================================================================
Entropy=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim(prim,U(:,i,j,k,iElem))
    ent_loc=prim(1)*(LOG(prim(5))-kappa*LOG(prim(1)))
    Entropy  = Entropy+ent_loc*wGPVol(i,j,k)/sJ(i,j,k,iElem)
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Entropy,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(Entropy     ,0      ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
Entropy=Entropy*sKappaM1


END SUBROUTINE CalcEntropy


!==================================================================================================================================
!> Finalizes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
