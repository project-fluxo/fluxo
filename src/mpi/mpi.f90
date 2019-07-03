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
!> \brief Contains the routines that set up communicators and control non-blocking communication
!==================================================================================================================================
MODULE MOD_MPI
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitMPI
  MODULE PROCEDURE InitMPI
END INTERFACE

PUBLIC::InitMPI

#if MPI
INTERFACE InitMPIvars
  MODULE PROCEDURE InitMPIvars
END INTERFACE

!INTERFACE StartReceiveMPIData
!  MODULE PROCEDURE StartReceiveMPIData
!END INTERFACE
!
!INTERFACE StartSendMPIData
!  MODULE PROCEDURE StartSendMPIData
!END INTERFACE

!INTERFACE FinishExchangeMPIData
!  MODULE PROCEDURE FinishExchangeMPIData
!END INTERFACE

INTERFACE FinalizeMPI
  MODULE PROCEDURE FinalizeMPI
END INTERFACE

PUBLIC::InitMPIvars
PUBLIC::StartReceiveMPIData
PUBLIC::StartSendMPIData
PUBLIC::FinishExchangeMPIData
PUBLIC::FinalizeMPI
#endif
!==================================================================================================================================

PUBLIC::DefineParametersMPI
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersMPI()
! MODULES
USE MOD_ReadInTools,              ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("MPI")
CALL prms%CreateIntOption('GroupSize', "Define size of MPI subgroups, used to e.g. perform grouped IO, where group master\n"//&
                                       "collects and outputs data.",&
                                       '0')
END SUBROUTINE DefineParametersMPI


!==================================================================================================================================
!> Basic mpi initialization. Calls initialization routine of the mpi library and sets myRank, nProcessors and MPIRoot. If the code 
!> is not compiled with mpi, InitMPI sets standard values for these variables.
!==================================================================================================================================
SUBROUTINE InitMPI()
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if MPI
CALL MPI_INIT(iError)
IF(iError .NE. 0) &
  CALL Abort(__STAMP__,'Error in MPI_INIT',iError)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank     , iError)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nProcessors, iError)
IF(iError .NE. 0) &
  CALL Abort(__STAMP__,'Could not get rank and number of processors',iError)
MPIRoot=(myRank .EQ. 0)
#else  /*MPI*/
myRank      = 0
myLocalRank = 0
nProcessors = 1
MPIRoot     =.TRUE.
MPILocalRoot=.TRUE.
#endif  /*MPI*/

! At this point the initialization is not completed. We first have to create a new MPI communicator. 
END SUBROUTINE InitMPI



#if MPI
!==================================================================================================================================
!> Initialize derived mpi variables used for communication
!==================================================================================================================================
SUBROUTINE InitMPIVars()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_ReadinTools,             ONLY: GETINT
USE MOD_Interpolation_Vars,      ONLY: InterpolationInitIsDone
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: color,groupsize
!==================================================================================================================================
IF(.NOT.InterpolationInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,'InitMPITypes called before InitInterpolation')
END IF

ALLOCATE(MPIRequest_U(nNbProcs,2)    )
ALLOCATE(MPIRequest_Flux(nNbProcs,2) )
MPIRequest_U      = MPI_REQUEST_NULL
MPIRequest_Flux   = MPI_REQUEST_NULL

#if PARABOLIC
ALLOCATE(MPIRequest_Lifting(nNbProcs,3,2))
MPIRequest_Lifting = MPI_REQUEST_NULL
#endif /*PARABOLIC*/

DataSizeSide      =PP_nVar*(PP_N+1)**2

! split communicator into smaller groups (e.g. for local nodes)
GroupSize=GETINT('GroupSize','0')
IF(GroupSize.LT.1)THEN ! group procs by node
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myRank,myRank,MPI_COMM_NODE,iError)
ELSE ! use groupsize
  color=myRank/GroupSize
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,myRank,MPI_COMM_NODE,iError)
END IF
CALL MPI_COMM_RANK(MPI_COMM_NODE,myLocalRank,iError)
CALL MPI_COMM_SIZE(MPI_COMM_NODE,nLocalProcs,iError)
MPILocalRoot=(myLocalRank .EQ. 0)

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS=MPI_COMM_NULL
MPI_COMM_WORKERS=MPI_COMM_NULL
myLeaderRank=-1
myWorkerRank=-1
IF(myLocalRank.EQ.0)THEN
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,0,myRank,MPI_COMM_LEADERS,iError)
  CALL MPI_COMM_RANK( MPI_COMM_LEADERS,myLeaderRank,iError)
  CALL MPI_COMM_SIZE( MPI_COMM_LEADERS,nLeaderProcs,iError)
  nWorkerProcs=nProcessors-nLeaderProcs
ELSE
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,1,myRank,MPI_COMM_WORKERS,iError)
  CALL MPI_COMM_RANK( MPI_COMM_WORKERS,myWorkerRank,iError)
  CALL MPI_COMM_SIZE( MPI_COMM_WORKERS,nWorkerProcs,iError)
  nLeaderProcs=nProcessors-nWorkerProcs
END IF
END SUBROUTINE InitMPIvars



!==================================================================================================================================
!> Subroutine that controls the receive operations for the face data that has to be exchanged between processors.
!>  (1:PP_nVar,0:PP_N,0:PP_N,LowerBound:upperBound) to an array (1:DataSize,LowerBound:upperBound)
!==================================================================================================================================
SUBROUTINE StartReceiveMPIData(FaceData,DataSize,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: SendID                                   !< defines the send / receive direction -> 1=send MINE 
                                                                        !< / receive YOUR, 2=send YOUR / receive MINE 
INTEGER,INTENT(IN)          :: DataSize                                 !< size of one entry in array (e.g. one side: 
                                                                        !< nVar*(N+1)**2
INTEGER,INTENT(IN)          :: LowerBound                               !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)          :: UpperBound                               !< upper side index for last dimension of FaceData
INTEGER,INTENT(OUT)         :: MPIRequest(nNbProcs)                     !< communication handles
REAL,INTENT(OUT)            :: FaceData(DataSize,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iNBProc
!==================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =DataSize*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(FaceData(:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,MPIRequest(iNbProc),iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPIData



!==================================================================================================================================
!> Subroutine that performs the send operations for the face data that has to be exchanged between processors.
!> this routine has no interface, since the FaceData is cast from an array
!>  (1:PP_nVar,0:PP_N,0:PP_N,LowerBound:upperBound) to an array (1:DataSize,LowerBound:upperBound)
!==================================================================================================================================
SUBROUTINE StartSendMPIData(FaceData,DataSize,LowerBound,UpperBound,MPIRequest,SendID)
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: SendID                                   !< defines the send / receive direction -> 1=send MINE 
                                                                        !< / receive YOUR, 2=send YOUR / receive MINE 
INTEGER,INTENT(IN)          :: DataSize                                 !< size of one entry in array (e.g. one side:
                                                                        !< nVar*(N+1)*(N+1))
INTEGER,INTENT(IN)          :: LowerBound                               !< lower side index for last dimension of FaceData
INTEGER,INTENT(IN)          :: UpperBound                               !< upper side index for last dimension of FaceData
INTEGER,INTENT(OUT)         :: MPIRequest(nNbProcs)                     !< communication handles
REAL,INTENT(IN)             :: FaceData(DataSize,LowerBound:UpperBound) !< the complete face data (for inner, BC and MPI sides).
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iNBProc
!==================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =DataSize*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(FaceData(:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,MPIRequest(iNbProc),iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPIData


!==================================================================================================================================
!> We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
!==================================================================================================================================
SUBROUTINE FinishExchangeMPIData(nRequests,MPIRequest)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nRequests             !< size of the handles 
INTEGER,INTENT(INOUT)       :: MPIRequest(nRequests) !< communication handles
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL MPI_WaitAll(nRequests,MPIRequest,MPI_STATUSES_IGNORE,iError)
END SUBROUTINE FinishExchangeMPIData

!==================================================================================================================================
!> Deallocate MPI arrays
!==================================================================================================================================
SUBROUTINE FinalizeMPI()
! MODULES
USE MOD_MPI_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(MPIRequest_U)
SDEALLOCATE(MPIRequest_Flux)
SDEALLOCATE(nMPISides_send)
SDEALLOCATE(OffsetMPISides_send)
SDEALLOCATE(nMPISides_rec)
SDEALLOCATE(OffsetMPISides_rec)
#if PARABOLIC
SDEALLOCATE(MPIRequest_Lifting)
#endif
END SUBROUTINE FinalizeMPI

#endif /*MPI*/

END MODULE MOD_MPI
