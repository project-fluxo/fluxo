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
CALL prms%CreateLogicalOption('CalcBulk',"Set true to compute the integrated mean value of U and Ut over whole domain"&
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
  IF(doCalcBulk)THEN
    DO iVar=1,PP_nVar
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"Bulk_'//TRIM(StrVarNames(iVar))//'"'
    END DO !iVar 
    DO iVar=1,PP_nVar
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"Bulk_t_'//TRIM(StrVarNames(iVar))//'"'
    END DO !iVar 
  END IF ! doCalcBulk
  IF(doCalcDivergence)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"maxDivB"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"max|[Bn]|"'
  END IF  !doCalcDivergence
  IF(doCalcEnergy)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"KineticEnergy"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"MagneticEnergy"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"TotalEnergy"'
#ifdef PP_GLM
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"PsiEnergy"'
#endif /*PP_GLM*/
  END IF !doCalcEnergy
  IF(doCalcEntropy)THEN
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"Entropy"'
    A2F_iVar=A2F_iVar+1
    A2F_VarNames(A2F_iVar)='"dSdU_Ut"'
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
REAL                 :: maxJumpB,maxDivB
REAL                 :: bulk(1:PP_nVar),bulk_t(1:PP_nVar)
REAL                 :: tmp(2),dSdU_Ut 
!==================================================================================================================================
! Calculate divergence 
IF(doCalcBulk)THEN
  CALL CalcBulk(bulk,bulk_t)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',PP_nVar,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)'Bulk       : ',bulk(:)
    WRITE(UNIT_StdOut,formatStr)'Bulk_t     : ',bulk(:)
    IF(doAnalyzeToFile)THEN
      A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=bulk(:)
      A2F_iVar=A2F_iVar+PP_nVar
      A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=bulk_t(:)
      A2F_iVar=A2F_iVar+PP_nVar
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF

IF(doCalcDivergence)THEN
  CALL CalcDivergence(maxDivB,maxJumpB)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A21,',2,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)'maxdivB,maxJumpB  : ',maxDivB,maxJumpB
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=maxDivB
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=maxJumpB
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF  !(doCalcDivergence)

IF(doCalcEnergy)THEN
  !store last energies
  tmp(1:2)=Energy(1:2)
  CALL CalcEnergy(Energy)
  IF(MPIroot) THEN
#ifdef PP_GLM
    WRITE(formatStr,'(A5,I1,A7)')'(A21,',4,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)  'Ekin/mag/tot/psi  : ',Energy(1:4)
#else
    WRITE(formatStr,'(A5,I1,A7)')'(A21,',3,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)  'Ekin/mag/total    : ',Energy(1:3)
#endif 
    IF((time-RestartTime).GE.Analyze_dt)THEN
      tmp(1:2)=LOG(Energy(1:2) /tmp(1:2))/Analyze_dt
      WRITE(UNIT_StdOut,formatStr)'  ...growth rates : ',tmp(1:2)
    END IF !time>Analyze_dt
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(1) !kineticEnergy
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(2) !magneticEnergy
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(3) !totalEnergy
#ifdef PP_GLM
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Energy(4) !PsiEnergy
#endif 
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF !doCalcEnergy

IF(doCalcEntropy)THEN
  tmp(1)=Entropy
  CALL CalcEntropy(Entropy,dSdU_Ut)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A)')'(A21,ES21.12)'
    WRITE(UNIT_StdOut,formatStr)  'Entropy           : ',Entropy
    IF((time-RestartTime).GE.Analyze_dt)THEN
      WRITE(UNIT_StdOut,formatStr)'  ...dEntropy/dt  : ',(Entropy-tmp(1))/Analyze_dt
    END IF !time>Analyze_dt
    WRITE(UNIT_StdOut,formatStr)  'dSdU*Ut           : ',dSdU_Ut
    IF(doAnalyzeToFile)THEN
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=Entropy
      A2F_iVar=A2F_iVar+1
      A2F_Data(A2F_iVar)=dSdU_Ut
    END IF !doAnalyzeToFile
  END IF !MPIroot
END IF !doCalcEntropy
END SUBROUTINE AnalyzeEquation


!==================================================================================================================================
!> Calculates the maximum of the discrete divergence of the magnetic field and the maximum Jumps of B*n
!==================================================================================================================================
SUBROUTINE CalcDivergence(maxDivB,maxJumpB)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars      ,ONLY: nElems
USE MOD_Mesh_Vars      ,ONLY: sJ,metrics_ftilde,metrics_gtilde,metrics_htilde,NormVec
USE MOD_Mesh_Vars      ,ONLY: ElemToSide,firstInnerSide,LastMPISide_MINE
USE MOD_DG_Vars        ,ONLY: D,U,U_master,U_slave 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: maxDivB
REAL,INTENT(OUT)                :: maxJumpB
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: i,j,k,l,iElem,SideID
REAL    :: divB_loc,Btilde(3,0:PP_N,0:PP_N,0:PP_N)
#if MPI
REAL                            :: box(2)
#endif
!==================================================================================================================================
  maxDivB=-1.0e20
  maxJumpB=-1.0e20
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Btilde(1,i,j,k)=SUM(Metrics_ftilde(:,i,j,k,iElem)*U(6:8,i,j,k,iElem))
      Btilde(2,i,j,k)=SUM(Metrics_gtilde(:,i,j,k,iElem)*U(6:8,i,j,k,iElem))
      Btilde(3,i,j,k)=SUM(Metrics_htilde(:,i,j,k,iElem)*U(6:8,i,j,k,iElem))
    END DO; END DO; END DO ! i,j,k
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      divB_loc=0.
      DO l=0,PP_N
        divB_loc=divB_loc + D(i,l)*Btilde(1,l,j,k) + &
                            D(j,l)*Btilde(2,i,l,k) + &
                            D(k,l)*Btilde(3,i,j,l)
      END DO ! l
      divB_loc = sJ(i,j,k,iElem) * divB_loc 
      maxDivB=MAX(maxDivB,ABS(divB_loc))
  
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
  DO SideID=firstInnerSide,LastMPISide_MINE
    maxJumpB=MAX(maxJumpB,  &
                 MAXVAL(ABS( NormVec(1,:,:,SideID)*(U_slave(6,:,:,SideID)-U_master(6,:,:,SideID)) &
                            +NormVec(2,:,:,SideID)*(U_slave(7,:,:,SideID)-U_master(7,:,:,SideID)) &
                            +NormVec(3,:,:,SideID)*(U_slave(8,:,:,SideID)-U_master(8,:,:,SideID)) )) )
  END DO !SideID
#if MPI
  Box=(/maxDivB,maxJumpB/)
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,box,2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
    maxDivB =Box(1)
    maxJumpB=Box(2)
  ELSE
    CALL MPI_REDUCE(Box         ,0  ,2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  END IF
#endif

END SUBROUTINE CalcDivergence 


!==================================================================================================================================
!> Calculates bulk velocities over whole domain
!==================================================================================================================================
SUBROUTINE CalcBulk(Bulk,Bulk_t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U,Ut
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Bulk(1:PP_nVar)
REAL,INTENT(OUT)                :: Bulk_t(1:PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: iElem,i,j,k
#if MPI
REAL                            :: box(1:2*PP_nVar)
#endif
!==================================================================================================================================
Bulk=0.
Bulk_t=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    Bulk   = Bulk  + U(:,i,j,k,iElem)*wGPVol(i,j,k)/sJ(i,j,k,iElem)
    Bulk_t = Bulk_t+Ut(:,i,j,k,iElem)*wGPVol(i,j,k)/sJ(i,j,k,iElem)
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
Box=(/Bulk,Bulk_t/)
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,box,2*PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  Bulk=Box(1:PP_nVar)
  Bulk_t=Box(PP_nVar+1:2*PP_nVar)
ELSE
  CALL MPI_REDUCE(Box         ,0  ,2*PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif

Bulk    = Bulk  /Vol
Bulk_t  = Bulk_t/Vol

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
REAL,INTENT(OUT)                :: Energy(4) !< kinetic and magnetic energy / total and psi divcorr
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
    Energy(3)  = Energy(3)+U(5,i,j,k,iElem)*IntegrationWeight
#ifdef PP_GLM
    Energy(4)  = Energy(4)+U(PP_nVar,i,j,k,iElem)*IntegrationWeight
#endif
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Energy,4,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(Energy         ,0  ,4,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
Energy(1)=0.5*Energy(1)/vol
Energy(2)= s2mu_0*Energy(2)/vol
Energy(4)= 0.5*Energy(4)/vol
END SUBROUTINE CalcEnergy


!==================================================================================================================================
!> Calculates  Entropy over whole domain Entropy=rho*s/(kappa-1), s=ln(p rho^(-kappa))=ln(p)-kappa*ln(rho) 
!==================================================================================================================================
SUBROUTINE CalcEntropy(Entropy,dSdU_Ut)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U,Ut
USE MOD_Equation_Vars,      ONLY: kappa,sKappaM1,ConsToPrim,ConsToEntropy
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Entropy !< Entropy
REAL,INTENT(OUT)             :: dSdU_Ut !< dS/dU * U_t
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: iElem,i,j,k
REAL                         :: ent_loc,prim(PP_nVar),dSdU(PP_nVar)
#if MPI
REAL                         :: box(2)
#endif 
!==================================================================================================================================
Entropy=0.
dSdU_Ut=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim(prim,U(:,i,j,k,iElem))
    ent_loc  = -prim(1)*(LOG(prim(5))-kappa*LOG(prim(1)))
    Entropy  = Entropy+ent_loc*wGPVol(i,j,k)/sJ(i,j,k,iElem)
    dSdU     = ConsToEntropy(U(:,i,j,k,iElem))
    ent_loc  = SUM(dSdU(:)*Ut(:,i,j,k,iElem))
    dSdU_Ut  = dSdU_Ut+ent_loc*wGPVol(i,j,k)/sJ(i,j,k,iElem)
  END DO; END DO; END DO !i,j,k
END DO ! iElem

#if MPI
box(1) = Entropy
box(2) = dSdU_Ut
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,box,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  Entropy = box(1)
  dSdU_Ut = box(2)
ELSE
  CALL MPI_REDUCE(box         ,0  ,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
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
