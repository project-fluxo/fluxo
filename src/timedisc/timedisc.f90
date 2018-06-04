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
!> Module for the global time stepping temporal discretization
!==================================================================================================================================
MODULE MOD_TimeDisc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitTimeDisc
  MODULE PROCEDURE InitTimeDisc
END INTERFACE

INTERFACE TimeDisc
  MODULE PROCEDURE TimeDisc
END INTERFACE

INTERFACE FinalizeTimeDisc
  MODULE PROCEDURE FinalizeTimeDisc
END INTERFACE

PUBLIC:: InitTimeDisc
PUBLIC:: TimeDisc
PUBLIC:: FinalizeTimeDisc
PUBLIC:: DefineParametersTimeDisc
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersTimeDisc()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("TimeDisc")
CALL prms%CreateStringOption('TimeDiscMethod', "Specifies the type of time-discretization to be used, \ne.g. the name of"//&
                                               " a specific Runge-Kutta scheme. Possible values:\n"//&
                                               "  * standardrk3-3\n  * carpenterrk4-5\n  * niegemannrk4-14\n"//&
                                               "  * toulorgerk4-8c\n  * toulorgerk3-7c\n  * toulorgerk4-8f\n"//&
                                               "  * ketchesonrk4-20\n  * ketchesonrk4-18", value='CarpenterRK4-5')
CALL prms%CreateRealOption(  'TEnd',           "End time of the simulation (mandatory).")
CALL prms%CreateRealOption(  'CFLScale',       "Scaling factor for the theoretical CFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateRealOption(  'DFLScale',       "Scaling factor for the theoretical DFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateIntOption(   'maxIter',        "Stop simulation when specified number of timesteps has been performed.", value='-1')
CALL prms%CreateIntOption(   'NCalcTimeStepMax',"Compute dt at least after every Nth timestep.", value='1')
END SUBROUTINE DefineParametersTimeDisc

!==================================================================================================================================
!> Get information for end time and max time steps from ini file
!==================================================================================================================================
SUBROUTINE InitTimeDisc()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars
USE MOD_ReadInTools    ,ONLY:GETREAL,GETINT,GETSTR
USE MOD_StringTools    ,ONLY:LowCase,StripSpaces
USE MOD_Mesh_Vars      ,ONLY:nElems
USE MOD_IO_HDF5        ,ONLY:AddToElemData
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255):: TimeDiscMethod
!==================================================================================================================================
TimeDiscMethod = GETSTR('TimeDiscMethod','Carpenter RK4-5')
CALL StripSpaces(TimeDiscMethod)
CALL LowCase(TimeDiscMethod)

CALL SetTimeDiscCoefs(TimeDiscMethod)
SELECT CASE(TimeDiscType)
CASE('LSERKW2')
  TimeStep=>TimeStepByLSERKW2
CASE('LSERKK3')
  TimeStep=>TimeStepByLSERKK3
END SELECT

IF(TimeDiscInitIsDone)THEN
   SWRITE(*,*) "InitTimeDisc already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TIMEDISC...'

! Read the end time TEnd from ini file
TEnd     = GETREAL('TEnd')
! Read the normalized CFL number
CFLScale = GETREAL('CFLScale')
#if PARABOLIC
! Read the normalized DFL number
DFLScale = GETREAL('DFLScale')
#endif /*PARABOLIC*/
CALL fillCFL_DFL(PP_N)
! Set timestep to a large number
dt=HUGE(1.)
! Read max number of iterations to perform
maxIter = GETINT('maxIter','-1')
nCalcTimeStepMax = GETINT('nCalcTimeStepMax','1')

SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: '//TRIM(TimeDiscName)
ALLOCATE(dtElem(nElems))
dtElem=0.

CALL AddToElemData('dt',dtElem)

TimediscInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TIMEDISC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitTimeDisc



!==================================================================================================================================
!> GTS Temporal discretization
!==================================================================================================================================
SUBROUTINE TimeDisc()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_TimeDisc_Vars       ,ONLY: TEnd,t,dt,tAnalyze,ViscousTimeStep,maxIter,Timestep,nRKStages,nCalcTimeStepMax,PID
USE MOD_Analyze_Vars        ,ONLY: Analyze_dt,WriteData_dt,tWriteData,nWriteData
USE MOD_Analyze             ,ONLY: Analyze
USE MOD_Testcase_vars       ,ONLY: doTCpreTimeStep
USE MOD_Testcase_Pre        ,ONLY: CalcPreTimeStep
USE MOD_Restart_Vars        ,ONLY: DoRestart,RestartTime
USE MOD_CalcTimeStep        ,ONLY: CalcTimeStep
USE MOD_Output              ,ONLY: Visualize,PrintStatusLine
USE MOD_HDF5_Output         ,ONLY: WriteState
USE MOD_Mesh_Vars           ,ONLY: nGlobalElems
USE MOD_DG                  ,ONLY: DGTimeDerivative
USE MOD_DG_Vars             ,ONLY: U
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: dt_Min,dt_MinOld,dtAnalyze,dtEnd,tStart
INTEGER(KIND=8)              :: iter,iter_loc
REAL                         :: CalcTimeStart,CalcTimeEnd
INTEGER                      :: TimeArray(8)              !< Array for system time
INTEGER                      :: errType,nCalcTimestep,writeCounter
LOGICAL                      :: doAnalyze,doFinalize
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')

! write number of grid cells and dofs only once per computation
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#GridCells : ',REAL(nGlobalElems)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs      : ',REAL(nGlobalElems)*REAL((PP_N+1)**3)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#Procs     : ',REAL(nProcessors)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs/Proc : ',REAL(nGlobalElems*(PP_N+1)**3/nProcessors)

IF(.NOT.DoRestart)THEN
  t=0.
  SWRITE(UNIT_StdOut,*)'WRITING INITIAL SOLUTION:'
ELSE
  t=RestartTime
  SWRITE(UNIT_StdOut,*)'REWRITING SOLUTION:'
END IF

! Determine next write time, since it will be written into output file
tWriteData=MIN(t+WriteData_dt,tEnd)
tAnalyze=MIN(t+Analyze_dt,tEnd)

! --- Perform some preparational steps ---

! Do first RK stage of first timestep to fill gradients
dt_Min=CALCTIMESTEP(errType)
CALL DGTimeDerivative(t)


! Write the state at time=0, i.e. the initial condition
CALL WriteState(OutputTime=t, FutureTime=tWriteData,isErrorFile=.FALSE.)

CALL Visualize(t,U)

! No computation needed if tEnd=tStart!
IF((t.GE.tEnd).OR.maxIter.EQ.0) RETURN

iter=0
iter_loc=0
writeCounter=0
doAnalyze=.FALSE.
doFinalize=.FALSE.
! compute initial timestep
dt=CALCTIMESTEP(errType)
nCalcTimestep=0
dt_MinOld=-999.
IF(errType.NE.0) CALL abort(__STAMP__,&
  'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)

! Run initial analyze
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,*) 'Analyze of initial solution:' 
CALL Analyze(t,iter)

IF(MPIroot)THEN
  WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_StdOut,'(A,ES16.7)')'Initial Timestep  : ', dt
  IF(ViscousTimeStep) WRITE(UNIT_StdOut,'(A)')' Viscous timestep dominates! '
  WRITE(UNIT_StdOut,*)'CALCULATION RUNNING...'
END IF ! MPIroot


! Run computation
tStart = t
CalcTimeStart=FLUXOTIME()
DO
  IF(nCalcTimestepMax.EQ.1)THEN
    dt_Min=CALCTIMESTEP(errType)
  ELSE
    ! be careful, this is using an estimator, when to recompute the timestep
    IF(nCalcTimestep.LT.1)THEN
      dt_Min=CALCTIMESTEP(errType)
      nCalcTimestep=nCalcTimeStepMax
      dt_MinOld=dt_Min
    END IF
    nCalcTimestep=nCalcTimestep-1
  END IF !nCalcTimeStepMax <>1
  IF(errType.NE.0)THEN !error in time step computation
    CALL WriteState(OutputTime=t, FutureTime=tWriteData,isErrorFile=.TRUE.)
    CALL abort(__STAMP__,&
   'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)
  END IF

  dt=dt_Min
  dtAnalyze=HUGE(1.)
  IF(tAnalyze-t.LE.dt*(1.+1.E-4))THEN
    dtAnalyze=tAnalyze-t; doAnalyze=.TRUE.
  END IF
  dt=MIN(dt,dtAnalyze)

  dtEnd=HUGE(1.)
  IF(tEnd-t    .LE.dt*(1.+1.E-4))THEN
    dtEnd=tEnd-t;         doAnalyze=.TRUE.; doFinalize=.TRUE.
  END IF
  dt=MIN(dt,dtEnd)

  IF(doTCpreTimeStep) CALL CalcPreTimeStep(t,dt)


  !CALL PrintStatusLine(t,dt,tStart,tEnd)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Perform Timestep using a global time stepping routine, attention: only RK3 has time dependent BC
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL TimeStep(t)

  iter=iter+1
  iter_loc=iter_loc+1
  t=t+dt
  IF(iter.EQ.maxIter)THEN
    tEnd=t; tAnalyze=t; tWriteData=t
    doAnalyze=.TRUE.; doFinalize=.TRUE.
  END IF

  ! Analyze and output now
  IF(doAnalyze) THEN
    ! Call DG operator to fill face data, fluxes, gradients for analyze
    CALL DGTimeDerivative(t)

    CalcTimeEnd=FLUXOTIME()

    IF(MPIroot)THEN
      ! Get calculation time per DOF
      PID=(CalcTimeEnd-CalcTimeStart)*REAL(nProcessors)/(REAL(nGlobalElems)*REAL((PP_N+1)**3)*REAL(iter_loc))
      CALL DATE_AND_TIME(values=TimeArray) ! get System time
      WRITE(UNIT_StdOut,'(132("-"))')
      WRITE(UNIT_stdOut,'(A,I2.2,A1,I2.2,A1,I4.4,A1,I2.2,A1,I2.2,A1,I2.2)') &
        ' Sys date   :    ',timeArray(3),'.',timeArray(2),'.',timeArray(1),' ',timeArray(5),':',timeArray(6),':',timeArray(7)
      WRITE(UNIT_stdOut,'(A,ES12.5,A,I4)')' CALCULATION TIME PER TSTEP/DOF: [',PID,' sec ], nRKstages:',nRKstages
      WRITE(UNIT_StdOut,'(A,ES16.7)')' Timestep   : ',dt_Min
      IF(ViscousTimeStep) WRITE(UNIT_StdOut,'(A)')' Viscous timestep dominates! '
      WRITE(UNIT_stdOut,'(A,ES16.7)')'#Timesteps  : ',REAL(iter)
    END IF !MPIroot

    ! Visualize data and write solution
    writeCounter=writeCounter+1
    IF((writeCounter.EQ.nWriteData).OR.doFinalize)THEN
      ! Visualize data
      CALL Visualize(t,U)
      ! Write state to file
      CALL WriteState(OutputTime=t,FutureTime=tWriteData,isErrorFile=.FALSE.)
      writeCounter=0
      tWriteData=MIN(tAnalyze+WriteData_dt,tEnd)
    END IF

    ! do analysis
    CALL Analyze(t,iter)
    iter_loc=0
    CalcTimeStart=FLUXOTIME()
    tAnalyze=  MIN(tAnalyze+Analyze_dt,  tEnd)
    doAnalyze=.FALSE.
  END IF

  IF(doFinalize) EXIT
END DO

END SUBROUTINE TimeDisc


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE TimeStepByLSERKW2(t)
! MODULES
USE MOD_PreProc
USE MOD_Vector
USE MOD_DG           ,ONLY: DGTimeDerivative
USE MOD_DG_Vars      ,ONLY: U,Ut,nTotalU
USE MOD_TimeDisc_Vars,ONLY: dt,RKA,RKb,RKc,nRKStages,CurrentStage
USE MOD_Mesh_Vars    ,ONLY: nElems
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: t                                     !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: Ut_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) ! temporal variable for Ut
REAL     :: tStage,b_dt(1:nRKStages)
INTEGER  :: iStage
!===================================================================================================================================
! Premultiply with dt
b_dt=RKb*dt

! First evaluation of DG operator already done in timedisc
CurrentStage=1
tStage=t
CALL DGTimeDerivative(tStage)
CALL VCopy(nTotalU,Ut_temp,Ut)               !Ut_temp = Ut
CALL VAXPBY(nTotalU,U,Ut,ConstIn=b_dt(1))    !U       = U + Ut*b_dt(1)


! Following steps
DO iStage=2,nRKStages
  CurrentStage=iStage
  tStage=t+dt*RKc(iStage)
  CALL DGTimeDerivative(tStage)
  CALL VAXPBY(nTotalU,Ut_temp,Ut,ConstOut=-RKA(iStage)) !Ut_temp = Ut - Ut_temp*RKA(iStage)
  CALL VAXPBY(nTotalU,U,Ut_temp,ConstIn =b_dt(iStage))  !U       = U + Ut_temp*b_dt(iStage)

END DO
CurrentStage=1

END SUBROUTINE TimeStepByLSERKW2


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration:  3 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE TimeStepByLSERKK3(t)
! MODULES
USE MOD_PreProc
USE MOD_Vector
USE MOD_DG           ,ONLY: DGTimeDerivative
USE MOD_DG_Vars      ,ONLY: U,Ut,nTotalU
USE MOD_TimeDisc_Vars,ONLY: dt,RKdelta,RKg1,RKg2,RKg3,RKb,RKc,nRKStages,CurrentStage
USE MOD_Mesh_Vars    ,ONLY: nElems
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: t                                     !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: S2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL     :: UPrev(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL     :: tStage,b_dt(1:nRKStages)
INTEGER  :: iStage
!===================================================================================================================================

! Premultiply with dt
b_dt=RKb*dt

! Nomenclature:
! S1 == U, S2 == S2, S3 == UPrev

CurrentStage=1
tStage=t
CALL VCopy(nTotalU,Uprev,U)                    !Uprev=U
CALL VCopy(nTotalU,S2,U)                       !S2=U
CALL DGTimeDerivative(t)
CALL VAXPBY(nTotalU,U,Ut,ConstIn=b_dt(1))      !U      = U + Ut*b_dt(1)

DO iStage=2,nRKStages
  CurrentStage=iStage
  tStage=t+dt*RKc(iStage)
  CALL DGTimeDerivative(tStage)
  CALL VAXPBY(nTotalU,S2,U,ConstIn=RKdelta(iStage))                !S2 = S2 + U*RKdelta(iStage)
  CALL VAXPBY(nTotalU,U,S2,ConstOut=RKg1(iStage),ConstIn=RKg2(iStage)) !U = RKg1(iStage)*U + RKg2(iStage)*S2
  CALL VAXPBY(nTotalU,U,Uprev,ConstIn=RKg3(iStage))                !U = U + RKg3(ek)*Uprev
  CALL VAXPBY(nTotalU,U,Ut,ConstIn=b_dt(iStage))                   !U = U + Ut*b_dt(iStage)
END DO
CurrentStage=1

END SUBROUTINE TimeStepByLSERKK3


!===================================================================================================================================
!> Scaling of the CFL number, from paper GASSNER, KOPRIVA, "A comparision of the Gauss and Gauss-Lobatto
!> Discontinuous Galerkin Spectral Element Method for Wave Propagation Problems" .
!> For N=1-10 input CFLscale can now be (nearly) 1. and will be scaled adequately depending on
!> polynomial degree N, NodeType and TimeDisc method.
!===================================================================================================================================
SUBROUTINE FillCFL_DFL(Nin)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars,ONLY:CFLScale,CFLScaleAlpha
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY:DFLScale,DFLScaleAlpha,RelativeDFL
#endif /*PARABOLIC*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nin !< input polynomial degree for 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: alpha
!===================================================================================================================================
! CFL in DG depends on the polynomial degree
! Runge-Kutta methods
alpha    = CFLScaleAlpha(MIN(15,Nin))
CFLScale = CFLScale*alpha
IF((Nin.GT.15).OR.(CFLScale.GT.alpha))THEN
  SWRITE(UNIT_StdOut,'(132("!"))')
  SWRITE(UNIT_StdOut,'(A)')'Warning: The chosen CFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_StdOut,'(132("!"))')
END IF
!scale with 2N+1
CFLScale = CFLScale/(2.*Nin+1.)
SWRITE(UNIT_stdOut,'(A,ES16.7)') '   CFL:',CFLScale

#if PARABOLIC
!########################### DFL ########################################
! DFL in DG depends on the polynomial degree
! since DFl is only on real axis, stability numbers are defined for RK3 and then scaled for RK4

alpha = DFLScaleAlpha(MIN(10,Nin))*RelativeDFL
DFLScale=DFLScale*alpha
IF((Nin.GT.10).OR.(DFLScale.GT.alpha))THEN
  SWRITE(UNIT_StdOut,'(132("!"))')
  SWRITE(UNIT_StdOut,'(A)')'Warning: The chosen DFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_StdOut,'(132("!"))')
END IF
DFLScale = DFLScale/(2.*Nin+1.)**2
SWRITE(UNIT_stdOut,'(A,ES16.7)') '   DFL:',DFLScale
#endif /*PARABOLIC*/
END SUBROUTINE fillCFL_DFL

!==================================================================================================================================
!> Finalizes variables necessary for timedisc subroutines
!==================================================================================================================================
SUBROUTINE FinalizeTimeDisc()
! MODULES
USE MOD_TimeDisc_Vars
IMPLICIT NONE
!==================================================================================================================================
TimeDiscInitIsDone = .FALSE.
SDEALLOCATE(dtElem)
SDEALLOCATE(RKA)
SDEALLOCATE(RKb)
SDEALLOCATE(RKc)
SDEALLOCATE(RKg1)
SDEALLOCATE(RKg2)
SDEALLOCATE(RKg3)
SDEALLOCATE(RKdelta)
NULLIFY(TimeStep)
END SUBROUTINE FinalizeTimeDisc


END MODULE MOD_TimeDisc
