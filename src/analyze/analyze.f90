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
!> Basic routines performing an analysis of the solution valid for all equation systems and testcases. 
!==================================================================================================================================
MODULE MOD_Analyze
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyze
  MODULE PROCEDURE InitAnalyze
END INTERFACE

INTERFACE Analyze
  MODULE PROCEDURE Analyze
END INTERFACE

INTERFACE FinalizeAnalyze
  MODULE PROCEDURE FinalizeAnalyze
END INTERFACE


PUBLIC:: DefineParametersAnalyze
PUBLIC:: InitAnalyze
PUBLIC:: Analyze
PUBLIC:: FinalizeAnalyze
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyze()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
USE MOD_AnalyzeEquation ,ONLY: DefineParametersAnalyzeEquation
USE MOD_Testcase_Analyze,ONLY: DefineParametersAnalyzeTestcase
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Analyze")
CALL prms%CreateLogicalOption('CalcErrorNorms' , "Set true to compute L2 and LInf error norms at analyze step.",&
                                                 '.TRUE.')
CALL prms%CreateLogicalOption('AnalyzeToFile',   "Set true to output result of error norms to a file (CalcErrorNorms=T)",&
                                                 '.FALSE.')
CALL prms%CreateRealOption(   'Analyze_dt',      "Specifies time intervall at which analysis routines are called.",&
                                                 '0.')
CALL prms%CreateIntOption(    'nWriteData' ,     "Intervall as multiple of Analyze_dt at which HDF5 files "//&
                                                 "(e.g. State,TimeAvg,Fluc) are written.",&
                                                 '1')
CALL prms%CreateIntOption(    'AnalyzeExactFunc',"possible to chosse another Exact for Analyis purposes ")
CALL prms%CreateLogicalOption('CalcMeanFlux',"Computes the integral of the flux over a boundary", '.FALSE.')
CALL prms%CreateIntOption(    'NAnalyze'   ,     "Polynomial degree at which analysis is performed (e.g. for L2 errors). "//&
                                                 "Default: 2*N.")
CALL DefineParametersAnalyzeEquation()
CALL DefineParametersAnalyzeTestcase()
END SUBROUTINE DefineParametersAnalyze

!==================================================================================================================================
!> Initializes variables necessary for analyze subroutines
!> - provides basic quantities like global domain volume, surface area of boundary conditions 
!>   or precomputed surface and volume integration weights
!> - initializes other specific analysis and benchmarking routines
!==================================================================================================================================
SUBROUTINE InitAnalyze()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_Equation_Vars,      ONLY: IniExactFunc
USE MOD_ReadInTools,        ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_StringTools,        ONLY: INTTOSTR
USE MOD_Interpolation_Vars, ONLY: xGP,wGP,wBary,InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: nBCs,SurfElem,nSides,AnalyzeSide,sJ,nElems
USE MOD_Mesh_Vars,          ONLY: BoundaryType,BoundaryName
USE MOD_Equation_Vars,      ONLY: StrVarNames
USE MOD_AnalyzeEquation,    ONLY: InitAnalyzeEquation
USE MOD_Testcase_Analyze,   ONLY: InitAnalyzeTestcase
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iVar,i,j,k,iSurf,iBC,iElem,iSide
!==================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.AnalyzeInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,'InitAnalyse not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ANALYZE...'

! Get the various analysis/output variables 
doCalcErrorNorms    =GETLOGICAL('CalcErrorNorms','.TRUE.')
IF(doCalcErrorNorms)THEN
  AnalyzeExactFunc =GETINT('AnalyzeExactFunc',INTTOSTR(IniExactFunc)) 
END IF
doCalcMeanFlux   =GETLOGICAL('CalcMeanFlux','.FALSE.')

Analyze_dt       =GETREAL('Analyze_dt','0.')
nWriteData       =GETINT('nWriteData','1')
NAnalyze         =GETINT('NAnalyze'   ,INTTOSTR(2*(PP_N+1)))

WriteData_dt = Analyze_dt*nWriteData

doAnalyzeToFile  =GETLOGICAL('AnalyzeToFile','.TRUE.')

IF(MPIroot.AND.doAnalyzeToFile)THEN
  ALLOCATE(A2F_VarNames(200)) !max of 100 entries for output
  ALLOCATE(A2F_Data(200))
  A2F_Data(:)=0.
  A2F_VarNames(:)=''
  A2F_iVar=0
  A2F_iVar=A2F_iVar+7
  A2F_VarNames(1)='"t_sim"' 
  A2F_VarNames(2)='"timesteps"' 
  A2F_VarNames(3)='"t_CPU"' 
  A2F_VarNames(4)='"PID"' 
  A2F_VarNames(5)='"DOF"' 
  A2F_VarNames(6)='"nElems"' 
  A2F_VarNames(7)='"nRanks"' 
  IF(doCalcErrorNorms)THEN
    DO iVar=1,PP_nVar
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"L2_'//TRIM(StrVarNames(iVar))//'"'
    END DO
    DO iVar=1,PP_nVar
      A2F_iVar=A2F_iVar+1
      A2F_VarNames(A2F_iVar)='"Linf_'//TRIM(StrVarNames(iVar))//'"'
    END DO
  END IF !doCalcErrorNorms
  IF(doCalcMeanFlux)THEN
    DO iBC=1,nBCs
      IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE
      DO iVar=1,PP_nVar
        A2F_iVar=A2F_iVar+1
        A2F_VarNames(A2F_iVar)='"MeanFlux_'//TRIM(BoundaryName(iBC))//'_'//TRIM(StrVarNames(iVar))//'"'
      END DO !iVar
    END DO
  END IF  !doCalcMeanFlux
END IF !MPIroot & doAnalyzeToFile
  

! precompute integration weights 
ALLOCATE(wGPSurf(0:PP_N,0:PP_N),wGPVol(0:PP_N,0:PP_N,0:PP_N))
DO j=0,PP_N; DO i=0,PP_N
  wGPSurf(i,j)  = wGP(i)*wGP(j) 
END DO; END DO
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  wGPVol(i,j,k) = wGP(i)*wGP(j)*wGP(k) 
END DO; END DO; END DO

! precompute volume of the domain
ALLOCATE(ElemVol(nElems))
ElemVol=0.
DO iElem=1,nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    ElemVol(iElem)=ElemVol(iElem)+wGPVol(i,j,k)/sJ(i,j,k,iElem)
  END DO; END DO; END DO !i,j,k
END DO ! iElem
Vol=SUM(ElemVol)


! compute surface of each boundary
ALLOCATE(Surf(nBCs))
Surf=0.
DO iSide=1,nSides
  iSurf=AnalyzeSide(iSide)
  IF(iSurf.EQ.0) CYCLE
  DO j=0,PP_N; DO i=0,PP_N
    Surf(iSurf)=Surf(iSurf)+wGPSurf(i,j)*SurfElem(i,j,iSide)
  END DO; END DO
END DO
#if MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Vol ,1   ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Surf,nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif /*MPI*/
IF(MPIROOT)THEN
DO iBC=1,nBCs
  IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE
  WRITE(UNIT_StdOut,'(A,A33,A3,E22.15)') ' |                     Surface of | ', &
                                          TRIM(BoundaryName(iBC)),' | ',Surf(iBC)
END DO
WRITE(UNIT_StdOut,'(A,A33,A3,E22.15)')   ' |                   total Volume | ', '',' | ',vol
END IF !MPIROOT

! Initialize eval routines
CALL InitAnalyzeBasis(PP_N,NAnalyze,xGP,wBary)


CALL InitAnalyzeEquation()


CALL InitAnalyzeTestcase()

IF(MPIroot.AND.doAnalyzeToFile)THEN
  A2F_nVars=A2F_iVar
  IF(A2F_nVars.GT.200) STOP 'max 200 entries for AnalyzeToFileOutput reached!'
END IF

AnalyzeInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitAnalyze



!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!> - Builds Vandermonde to interpolate the solution onto a Gauss-Lobatto mesh at a higher polynomial degree
!> - Precomputes volume interpolation weights
!==================================================================================================================================
SUBROUTINE InitAnalyzeBasis(N_in,Nanalyze_in,xGP,wBary)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars, ONLY: wGPVolAnalyze,Vdm_GaussN_NAnalyze
USE MOD_Basis,        ONLY: InitializeVandermonde
USE MOD_Interpolation,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars,ONLY: NodeTypeGL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: N_in                  !< input polynomial degree
INTEGER,INTENT(IN)               :: Nanalyze_in           !< polynomial degree of analysis polynomial
REAL,INTENT(IN)                  :: xGP(0:N_in)           !< interpolation points
REAL,INTENT(IN)                  :: wBary(0:N_in)         !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: XiAnalyze(0:Nanalyze_in)
REAL                             :: wAnalyze(0:NAnalyze_in)  ! GL integration weights used for the analyze
INTEGER                          :: i,j,k
!==================================================================================================================================
ALLOCATE(wGPVolAnalyze(0:Nanalyze_in,0:Nanalyze_in,0:Nanalyze_in),Vdm_GaussN_NAnalyze(0:NAnalyze_in,0:N_in))
CALL GetNodesAndWeights(NAnalyze_in,NodeTypeGL,XiAnalyze,wAnalyze)
CALL InitializeVandermonde(N_in,NAnalyze_in,wBary,xGP,XiAnalyze,Vdm_GaussN_NAnalyze)

DO k=0,Nanalyze_in; DO j=0,Nanalyze_in; DO i=0,Nanalyze_in
  wGPVolAnalyze(i,j,k) = wAnalyze(i)*wAnalyze(j)*wAnalyze(k)
END DO; END DO; END DO

END SUBROUTINE InitAnalyzeBasis



!==================================================================================================================================
!> Controls analysis routines and is called at analyze time levels
!> - calls generic error norm computation
!> - calls equation system specific analysis and testcase specific analysis
!==================================================================================================================================
SUBROUTINE Analyze(Time,iter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars
USE MOD_AnalyzeEquation,    ONLY: AnalyzeEquation
USE MOD_Testcase_Analyze,   ONLY: AnalyzeTestcase
USE MOD_Output_Vars,        ONLY: ProjectName
USE MOD_Mesh_Vars,          ONLY: nGlobalElems
USE MOD_Mesh_Vars,          ONLY: nBCs,BoundaryType,BoundaryName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time
INTEGER(KIND=8),INTENT(IN)      :: iter
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=40)               :: formatStr
INTEGER                         :: iBC
REAL                            :: CalcTime
REAL                            :: L_Inf_Error(PP_nVar),L_2_Error(PP_nVar)
REAL                            :: MeanFlux(PP_nVar,nBCs)
!==================================================================================================================================
! Graphical output
CalcTime=FLUXOTIME()
SWRITE(UNIT_StdOut,'(A14,ES16.7)')' Sim time   : ',Time

IF(MPIroot.AND.doAnalyzeToFile) THEN
  A2F_iVar=0  !reset counter
  A2F_iVar=A2F_iVar+7
  !A2F_Data(1:4) Time,iter,tCPU,PID  set in AnalyzeToFile
  A2F_Data(5)=REAL(nGlobalElems*(PP_N+1)**3) !DOF
  A2F_Data(6)=REAL(nGlobalElems)
  A2F_Data(7)=REAL(nProcessors)
END IF !MPIroot & AnalyzeToFile
! Calculate error norms
IF(doCalcErrorNorms)THEN
  CALL CalcErrorNorms(Time,L_2_Error,L_Inf_Error)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A7)')'(A14,',PP_nVar,'ES16.7)'
    WRITE(UNIT_StdOut,formatStr)' L_2        : ',L_2_Error
    WRITE(UNIT_StdOut,formatStr)' L_inf      : ',L_Inf_Error
    IF(doAnalyzeToFile)THEN
      A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=L_2_Error(:)
      A2F_iVar=A2F_iVar+PP_nVar
      A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=L_Inf_Error(:)
      A2F_iVar=A2F_iVar+PP_nVar
    END IF !doAnalyzeToFile  
  END IF !MPIroot
END IF  ! ErrorNorms
IF(doCalcMeanFlux)THEN
  CALL CalcMeanFlux(Time,MeanFlux)
  IF(MPIroot) THEN
    WRITE(formatStr,'(A5,I1,A9)')'(A14,',PP_nVar,'ES22.15)'
    DO iBC=1,nBCs
      IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE
      WRITE(UNIT_StdOut,*)'MeanFlux ',TRIM(BoundaryName(iBC)),' : '
      WRITE(UNIT_StdOut,formatStr)'              ',MeanFlux(:,iBC)
      IF(doAnalyzeToFile)THEN
        A2F_Data(A2F_iVar+1:A2F_iVar+PP_nVar)=MeanFlux(:,iBC)
        A2F_iVar=A2F_iVar+PP_nVar
      END IF !doAnalyzeToFile
    END DO
  END IF
END IF  !(doCalcMeanFlux)

! Attention: during the initialization phase no face data / gradients available!
!END IF

CALL AnalyzeEquation(Time)

CALL AnalyzeTestcase(Time)

IF(MPIroot.AND.doAnalyzeToFile) CALL AnalyzeToFile(Time,CalcTime,iter)

IF(MPIroot .AND. (Time.GT.0.)) THEN
  WRITE(UNIT_StdOut,'(132("."))')
  WRITE(UNIT_stdOut,'(A,A,A,F8.2,A)') ' FLUXO RUNNING ',TRIM(ProjectName),'... [',CalcTime-StartTime,' sec ]'
  WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_StdOut,*)
END IF

END SUBROUTINE Analyze



!==================================================================================================================================
!> Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!==================================================================================================================================
SUBROUTINE CalcErrorNorms(Time,L_2_Error,L_Inf_Error)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: Elem_xGP,sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Equation,           ONLY: ExactFunc
USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
USE MOD_Analyze_Vars,       ONLY: AnalyzeExactFunc
USE MOD_Analyze_Vars,       ONLY: NAnalyze,Vdm_GaussN_NAnalyze
USE MOD_Analyze_Vars,       ONLY: wGPVolAnalyze,Vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time                   !< current simulation time
REAL,INTENT(OUT)                :: L_2_Error(  PP_nVar)   !< L2 error of the solution
REAL,INTENT(OUT)                :: L_Inf_Error(PP_nVar)   !< LInf error of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem,k,l,m
REAL                            :: U_exact(PP_nVar)
REAL                            :: U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: Coords_NAnalyze(3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                            :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                            :: IntegrationWeight
!==================================================================================================================================
! Calculate error norms
L_Inf_Error(:)=-1.E10
L_2_Error(:)=0.
! Interpolate values of Error-Grid from GP's
DO iElem=1,nElems
  ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
  CALL ChangeBasis3D(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Elem_xGP(1:3,:,:,:,iElem),Coords_NAnalyze(1:3,:,:,:))
  ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  CALL ChangeBasis3D(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_NAnalyze(1:1,:,:,:))
  ! Interpolate the solution to the analyze grid
  CALL ChangeBasis3D(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,U(1:PP_nVar,:,:,:,iElem),U_NAnalyze(1:PP_nVar,:,:,:))
  DO m=0,NAnalyze
    DO l=0,NAnalyze
      DO k=0,NAnalyze
        CALL ExactFunc(AnalyzeExactFunc,time,Coords_NAnalyze(1:3,k,l,m),U_exact)
        L_Inf_Error = MAX(L_Inf_Error,abs(U_NAnalyze(:,k,l,m) - U_exact))
        IntegrationWeight = wGPVolAnalyze(k,l,m)*J_NAnalyze(1,k,l,m)
        ! To sum over the elements, We compute here the square of the L_2 error
         L_2_Error = L_2_Error+(U_NAnalyze(:,k,l,m) - U_exact)*(U_NAnalyze(:,k,l,m) - U_exact)*IntegrationWeight
      END DO ! k
    END DO ! l
  END DO ! m
END DO ! iElem=1,nElems

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,L_2_Error  ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,L_Inf_Error,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(L_2_Error  ,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(L_Inf_Error,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/Vol)

END SUBROUTINE CalcErrorNorms

!==================================================================================================================================
!> Integrate state variable U over a set of given analyzeSide (meanFlux)
!==================================================================================================================================
SUBROUTINE CalcMeanFlux(Time,MeanFlux)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Analyze_Vars,      ONLY: wGPSurf
USE MOD_DG_Vars,           ONLY: Flux
USE MOD_Mesh_Vars,         ONLY: nSides,nBCs,AnalyzeSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)               :: MeanFlux(PP_nVar,nBCs) 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSurf,iSide,i,j
!==================================================================================================================================
! Compute one DG stage for the current time to have the right fluxes
!CALL DGTimeDerivative(Time) !already called in timedisc
MeanFlux=0.

DO iSide=1,nSides
  iSurf=AnalyzeSide(iSide)
  IF(iSurf.EQ.0) CYCLE
  DO j=0,PP_N; DO i=0,PP_N
    MeanFlux(:,iSurf)=MeanFlux(:,iSurf)+Flux(:,i,j,iSide)*wGPSurf(i,j)
  END DO; END DO !i,j
END DO !iSide

#if MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,MeanFlux ,PP_nVar*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(MeanFlux    ,0        ,PP_nVar*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
!DO iBC=1,nBCs
!  IF(Boundarytype(iBC,BC_TYPE) .EQ. 1) CYCLE
!  MeanFlux(:,iBC)=MeanFlux(:,iBC)/SurfBC(iBC)
!END DO

END SUBROUTINE CalcMeanFlux


!==================================================================================================================================
!> Writes all analyze data to out file
!==================================================================================================================================
SUBROUTINE AnalyzeToFile(Time,CalcTime,iter)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars ,ONLY:IterRestart,CalcTimeRestart,A2F_VarNames,A2F_nVars,A2F_Data
USE MOD_Restart_Vars ,ONLY:RestartTime
USE MOD_Output_Vars  ,ONLY:ProjectName
USE MOD_TimeDisc_Vars,ONLY: PID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: TIME                         !< physical time
REAL,INTENT(IN)                :: CalcTime                     !< computational time
INTEGER(KIND=8),INTENT(IN)     :: iter                         !< number of timesteps
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                           :: DummyIter,Dummytime,DummyPID ! Dummy values for file handling
INTEGER                        :: openStat,stat                ! File IO status
INTEGER                        :: ioUnit 
INTEGER                        :: iVar 
CHARACTER(LEN=1)               :: DummyFirst 
CHARACTER(LEN=50)              :: formatStr                    ! format string for the output and Tecplot header
CHARACTER(LEN=300)             :: Filename                     ! Output filename,
LOGICAL                        :: fileExists                   ! Error handler for file
!==================================================================================================================================
Filename = 'out.'//TRIM(ProjectName)//'.dat'
! Check for file
INQUIRE(FILE = Filename, EXIST = fileExists)
!! File processing starts here open old and extratct information or create new file.
ioUnit=GETFREEUNIT()
IF(fileExists .AND. Time .NE. 0.0)THEN ! File exists and append data
  OPEN(UNIT     = ioUnit     , &
       FILE     = Filename   , &
       FORM     = 'FORMATTED', &
       STATUS   = 'OLD'      , &
       POSITION = 'APPEND'   , &
       RECL     = 50000      , &
       IOSTAT = openStat           )
  IF(openStat.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR: cannot open '//TRIM(Filename))
  END IF
  ! If we have a restart we need to find the position from where to move on. 
  ! Read the values from the previous analyse interval, get the CPUtime
  IF(Time .EQ. RestartTime) THEN
    WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'Searching for time stamp...'

    Dummytime = 0.0
    REWIND(ioUnit)
    READ(ioUnit,*)
    READ(ioUnit,*)
    READ(ioUnit,*)
    READ(ioUnit,*)
    ! Loop until we have found the position
    stat=0
    DO WHILE ((ABS(Dummytime-RestartTime).GT.1.0E-8*RestartTime) .AND. (stat.eq.0))
      READ(ioUnit,'(A1)',IOSTAT=stat) DummyFirst
      IF(INDEX(' "VZT',DummyFirst) .NE.0) THEN
        CYCLE
      ELSE
        BACKSPACE(ioUnit)
      END IF
      READ(ioUnit,*,IOSTAT=stat) Dummytime,DummyIter,CalcTimeRestart,DummyPID
    END DO
    IterRestart = INT(DummyIter)
    PID         = DummyPID
    WRITE(UNIT_stdOut,*)' Found. time: ' ,Dummytime,' with ', IterRestart,' timesteps and CPU time: ', &
              CalcTimeRestart
    WRITE(ioUnit,'(A)')'VARIABLES ='
    ! Variable names 
    DO iVar=1,A2F_nVars-1
      WRITE(ioUnit,'(A,1X)',ADVANCE="no")TRIM(A2F_VarNames(iVar))
    END DO
    WRITE(ioUnit,'(A)')TRIM(A2F_VarNames(A2F_nVars))
    WRITE(ioUnit,'(A,E23.14E5,A1)') 'ZONE T="Analysis,'//TRIM(ProjectName)//' restarttime=',Time,'"'
  END IF
ELSE ! No restart create new file
  OPEN(UNIT   = ioUnit       ,&
       FILE   = Filename     ,&
       STATUS = 'Unknown'    ,&
       ACCESS = 'SEQUENTIAL' ,&
       IOSTAT = openStat                 )
  IF (openStat.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR: cannot open '//TRIM(Filename))
  END IF
  ! Greate a new file with the Tecplot header etc. 
  WRITE(ioUnit,'(A)')'TITLE="Analysis,'//TRIM(ProjectName)//'"'
  WRITE(ioUnit,'(A)')'VARIABLES ='
  ! Variable names 
  DO iVar=1,A2F_nVars-1
    WRITE(ioUnit,'(A,1X)',ADVANCE="no")TRIM(A2F_VarNames(iVar))
  END DO
  WRITE(ioUnit,'(A)')TRIM(A2F_VarNames(A2F_nVars))
  WRITE(ioUnit,'(A)') 'ZONE T="Analysis,'//TRIM(ProjectName)//'"'
END IF
! Fill Time,timesteps and tCPU
A2F_Data(1)=Time  !Simulation time
A2F_Data(2)=REAL(iter+iterRestart) !timesteps
IF (iter+iterRestart .EQ. 0) THEN
  A2F_Data(3)=0.
ELSE
  A2F_Data(3)=CalcTime-StartTime+CalcTimeRestart !tCPU
END IF
A2F_Data(4)=PID
! Create format string for the variable output
WRITE(formatStr,'(A10,I3,A14)')'(E23.14E5,',A2F_nVars-1,'(1X,E23.14E5))'
WRITE(ioUnit,formatStr) A2F_Data(1:A2F_nVars)
CLOSE(ioUnit) ! outputfile
END SUBROUTINE AnalyzeToFile

!==================================================================================================================================
!> Finalizes variables necessary for analyze subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyze()
! MODULES
USE MOD_AnalyzeEquation,   ONLY: FinalizeAnalyzeEquation
USE MOD_Analyze_Vars,      ONLY: AnalyzeInitIsDone,wGPSurf,wGPVol,Surf,wGPVolAnalyze,Vdm_GaussN_NAnalyze,ElemVol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL FinalizeAnalyzeEquation()
SDEALLOCATE(Vdm_GaussN_NAnalyze)
SDEALLOCATE(wGPVolAnalyze)
SDEALLOCATE(Surf)
SDEALLOCATE(wGPVol)
SDEALLOCATE(wGPSurf)
SDEALLOCATE(ElemVol)
AnalyzeInitIsDone = .FALSE.
END SUBROUTINE FinalizeAnalyze


END MODULE MOD_Analyze
