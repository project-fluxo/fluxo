!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
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
!> Provides routines for visualization and ASCII output of time series. Initializes some variables used for HDF5 output.
!==================================================================================================================================
MODULE MOD_Output
! MODULES
USE MOD_ReadInTools
!USE ISO_C_BINDING
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

! Output format for state visualization
INTEGER,PARAMETER :: OUTPUTFORMAT_NONE         = 0
INTEGER,PARAMETER :: OUTPUTFORMAT_PARAVIEW_3D_SINGLE     = 1
INTEGER,PARAMETER :: OUTPUTFORMAT_PARAVIEW_3D_MULTI      = 2
INTEGER,PARAMETER :: OUTPUTFORMAT_PARAVIEW_2D_SINGLE     = 3
INTEGER,PARAMETER :: OUTPUTFORMAT_PARAVIEW_2D_MULTI      = 4

! Output format for ASCII data files
INTEGER,PARAMETER :: ASCIIOUTPUTFORMAT_CSV     = 0
INTEGER,PARAMETER :: ASCIIOUTPUTFORMAT_TECPLOT = 1

INTERFACE DefineParametersOutput
  MODULE PROCEDURE DefineParametersOutput
END INTERFACE

INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

INTERFACE PrintStatusLine
  MODULE PROCEDURE PrintStatusLine
END INTERFACE

INTERFACE Visualize
  MODULE PROCEDURE Visualize
END INTERFACE

INTERFACE VisualizeAny
  MODULE PROCEDURE VisualizeAny
END INTERFACE

INTERFACE InitOutputToFile
  MODULE PROCEDURE InitOutputToFile
END INTERFACE

INTERFACE OutputToFile
  MODULE PROCEDURE OutputToFile
END INTERFACE

INTERFACE FinalizeOutput
  MODULE PROCEDURE FinalizeOutput
END INTERFACE

PUBLIC:: InitOutput
PUBLIC:: PrintStatusLine
PUBLIC:: Visualize
PUBLIC:: VisualizeAny
PUBLIC:: InitOutputToFile
PUBLIC:: OutputToFile
PUBLIC:: FinalizeOutput
!==================================================================================================================================

PUBLIC::DefineParametersOutput
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersOutput()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Output")
CALL prms%CreateIntOption(          'NVisu',       "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateStringOption(       'ProjectName', "Name of the current simulation (mandatory).")
CALL prms%CreateLogicalOption(      'Logging',     "Write log files containing debug output.", '.FALSE.')
CALL prms%CreateLogicalOption(      'ErrorFiles',  "Write error files containing error output.", '.TRUE.')
CALL prms%CreateLogicalOption( 'PrimVisuDefault',  "Visualize the primitive variables by default.", '.FALSE.')
CALL prms%CreateIntOption('OutputFormat',"File format for visualization: 0: None, 1: ParaView Single vtu File,"//&
                                          "2: ParaView vtu files per proc+pvtu link file,"//&
                                          "3: 2D (zeta=-1 of each element) ParaView, single vtu file,"//&
                                          "4: 2D (zeta=-1 of each element) ParaView, vtu files per proc+pvtu link file", '1')
CALL prms%CreateIntOption('ASCIIOutputFormat',"File format for ASCII files, e.g. body forces: 0: CSV, 1: Tecplot." , '0')
CALL prms%CreateLogicalOption(      'doPrintStatusLine','Print: percentage of time, ...', '.FALSE.')
CALL prms%CreateLogicalOption(      'ColoredOutput','Colorize stdout', '.FALSE.')
CALL prms%CreateRealArrayOption('VisuBoundingBox'    , "Bounding box to reduce visualization output (multiple are possible!)"//&
                                                       "input (/xmin,xmax,ymin,ymax,zmin,zmax/)",&
                                                multiple=.TRUE.)
END SUBROUTINE DefineParametersOutput

!==================================================================================================================================
!> Initialize all output variables.
!==================================================================================================================================
SUBROUTINE InitOutput()
! MODULES
USE MOD_Globals
USE MOD_Preproc 
USE MOD_Output_Vars
USE MOD_Equation_Vars     ,ONLY:StrVarNames
USE MOD_ReadInTools       ,ONLY:GETSTR,GETLOGICAL,GETINT
USE MOD_StringTools       ,ONLY:INTTOSTR
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars,ONLY:InterpolationInitIsDone,NodeTypeVISU,NodeType
!USE ISO_C_BINDING,         ONLY: C_NULL_CHAR
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iBox,OpenStat, nVars
CHARACTER(LEN=8)               :: StrDate
CHARACTER(LEN=10)              :: StrTime
!==================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.OutputInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'InitOutput not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

NVisu=GETINT('NVisu',INTTOSTR(PP_N))

! Gauss/Gl -> Visu : computation -> visualization
ALLOCATE(Vdm_GaussN_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVISU,Vdm_GaussN_NVisu)

! Name for all output files
ProjectName=GETSTR('ProjectName')
Logging    =GETLOGICAL('Logging')
ErrorFiles =GETLOGICAL('ErrorFiles')
PrimVisuDefault=GETLOGICAL('PrimVisuDefault')

doPrintStatusLine=GETLOGICAL("doPrintStatusLine",".FALSE.")

WRITE(ErrorFileName,'(A,A8,I6.6,A4)')TRIM(ProjectName),'_ERRORS_',myRank,'.out'

! Get output format for state visualization
OutputFormat = GETINT('OutputFormat','0')
IF(OutputFormat.NE.0)THEN
  nBoundingBoxes=COUNTOPTION('VisuBoundingBox')
  IF(nBoundingBoxes.GT.0)THEN 
    ALLOCATE(VisuBoundingBox(6,nBoundingBoxes))
    DO iBox=1,nBoundingBoxes
      VisuBoundingBox(:,iBox)  = GETREALARRAY('VisuBoundingBox',6)
    END DO 
  END IF !nBoundingBoxes>0
END IF !outputFormat /= 0
! Get output format for ASCII data files
ASCIIOutputFormat = GETINT('ASCIIOutputFormat','0')

! Open file for logging
IF(Logging)THEN
  WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
  OPEN(UNIT=UNIT_logOut,  &
       FILE=LogFile,      &
       STATUS='REPLACE',  &
       ACTION='WRITE',    &
       IOSTAT=OpenStat)
  CALL DATE_AND_TIME(StrDate,StrTime)
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,'(132("#"))')
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                      StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
END IF  ! Logging

! Set the default number of output variables and allocate
! -------------------------------------------------------
#if defined (linearscalaradvection)
nOutVars=2
#elif (defined (mhd) | defined (navierstokes))
nOutVars=PP_nVar
#else
!default
nOutVars=PP_nVar
#endif
! If shock-capturing is activated output extra quantities
#if SHOCK_NFVSE
nOutvars = nOutvars + 1
#if NFVSE_CORR
nOutvars = nOutvars + 1
#endif /*NFVSE_CORR*/
#endif /*SHOCK_NFVSE*/
allocate(strvarnames_tmp(nOutVars))

! Set the default names
! ---------------------
#if defined (linearscalaradvection)
strvarnames_tmp(1)=StrVarnames(1)
strvarnames_tmp(2)='ExactSolution'
#endif /*linearscalaradvection*/
nVars = PP_nVar

#if SHOCK_NFVSE
nVars = nVars+1
strvarnames_tmp(nVars) = 'BlendingFunction'
#if NFVSE_CORR
nVars = nVars+1
strvarnames_tmp(nVars) = 'BlendingFunction_old'
#endif /*NFVSE_CORR*/
#endif /*SHOCK_NFVSE*/


OutputInitIsDone =.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOutput

SUBROUTINE PrintStatusLine(t,dt,tStart,tEnd) 
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars , ONLY: doPrintStatusLine
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN) :: t,dt,tStart,tEnd
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: percent,time_remaining,mins,secs,hours
!==================================================================================================================================
IF(.NOT.doPrintStatusLine) RETURN

IF(MPIroot)THEN
#ifdef INTEL
  OPEN(UNIT_stdOut,CARRIAGECONTROL='fortran')
#endif
  percent = (t-tStart) / (tend-tStart) 
  CALL CPU_TIME(time_remaining)
  IF (percent.GT.0.0) time_remaining = time_remaining/percent - time_remaining
  percent = percent*100.
  secs = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  mins = MOD(time_remaining,60.)
  time_remaining = time_remaining / 60
  hours = MOD(time_remaining,24.)
  WRITE(UNIT_stdOut,'(A,E10.4,A,E10.4,A,F6.2,A,I4,A1,I0.2,A1,I0.2,A1)',ADVANCE='NO') 'Time = ', t, &
      ' dt = ', dt, '  ', percent, '% complete, est. Time Remaining = ',INT(hours),':',INT(mins),':',INT(secs), ACHAR(13)
#ifdef INTEL
  CLOSE(UNIT_stdOut)
#endif
END IF
END SUBROUTINE PrintStatusLine

!==================================================================================================================================
!> Supersample DG dataset at (equidistant) visualization points and output to file.
!> Currently only Paraview binary format is supported.
!==================================================================================================================================
SUBROUTINE Visualize(OutputTime,Uin,FileTypeStrIn,PrimVisuOpt,StrVarNames_opt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars  ,ONLY:Elem_xGP,nElems
USE MOD_Equation_Vars,ONLY:StrVarNames
#if defined(linearscalaradvection)
USE MOD_Equation      ,ONLY: ExactFunc
USE MOD_Analyze_Vars  ,ONLY: AnalyzeExactFunc
#elif (defined(mhd) || defined(navierstokes))
USE MOD_Equation_Vars ,ONLY:StrVarnamesPrim,ConsToPrim
#endif /*defined(mhd)*/
#if SHOCK_NFVSE
use MOD_NFVSE_Vars ,only: alpha
#if NFVSE_CORR
use MOD_NFVSE_Vars ,only: alpha_old
#endif /*NFVSE_CORR*/
#endif /*SHOCK_NFVSE*/
USE MOD_Output_Vars,ONLY:OutputFormat
USE MOD_Mesh_Vars  ,ONLY:Elem_xGP,nElems
USE MOD_Output_Vars,ONLY:NVisu,Vdm_GaussN_NVisu, strvarnames_tmp, nOutVars, PrimVisuDefault
USE MOD_ChangeBasis,ONLY:ChangeBasis3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)               :: OutputTime               !< simulation time of output
REAL,INTENT(IN)               :: Uin(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution vector to be visualized
CHARACTER(LEN=255),OPTIONAL,INTENT(IN) :: FileTypeStrIn
LOGICAL,OPTIONAL,INTENT(IN) :: PrimVisuOpt
CHARACTER(LEN=255),OPTIONAL,INTENT(IN) :: StrVarNames_Opt(PP_nVar)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem
REAL,ALLOCATABLE              :: Coords_NVisu(:,:,:,:,:) 
REAL,ALLOCATABLE              :: U_NVisu(:,:,:,:,:)
CHARACTER(LEN=255)            :: FileTypeStr
LOGICAL                       :: PrimVisu
#if defined (linearscalaradvection)
REAL                          :: Uex(1,0:PP_N,0:PP_N,0:PP_N)
#elif (defined (mhd) || defined (navierstokes))
REAL                          :: cons(PP_nVar)
#endif
INTEGER                       :: i,j,k
integer                       :: nVars
!==================================================================================================================================
IF(outputFormat.LE.0) RETURN

IF(PRESENT(FileTypeStrIn))THEN
  FileTypeStr=FileTypeStrIn
ELSE
  FileTypeStr='Solution'
END IF !PRESENT(FileTypeStrIn)
IF(PRESENT(PrimVisuOpt))THEN
  PrimVisu=PrimVisuOpt
ELSE
  PrimVisu=PrimVisuDefault
END IF
! Specify output names
#if (defined (mhd) | defined (navierstokes))
IF(PrimVisu)THEN
  strvarnames_tmp(1:PP_nVar)=StrVarnamesPrim
ELSE
  IF(PRESENT(StrVarNames_opt))THEN
    strvarnames_tmp(1:PP_nVar)=StrVarnames_opt
  ELSE
    strvarnames_tmp(1:PP_nVar)=StrVarnames
  END IF
END IF
#else
!default
IF(PRESENT(StrVarNames_opt))THEN
  strvarnames_tmp(1:PP_nVar)=StrVarnames_opt
ELSE
  strvarnames_tmp(1:PP_nVar)=StrVarnames
END IF
#endif

ALLOCATE(U_NVisu(1:nOutVars,0:NVisu,0:NVisu,0:NVisu,1:nElems))
U_NVisu = 0.
ALLOCATE(Coords_NVisu(1:3,0:NVisu,0:NVisu,0:NVisu,1:nElems))

DO iElem=1,nElems
    ! Create coordinates of visualization points (TODO: can't this be done only once at the beginning?)
    CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
    ! Interpolate solution onto visu grid
    CALL ChangeBasis3D(PP_nVar,PP_N,NVisu,Vdm_GaussN_NVisu,Uin(1:PP_nVar,:,:,:,iElem),U_NVisu(1:PP_nVar,:,:,:,iElem))
#if defined (linearscalaradvection)
  DO k=0,Nvisu ; DO j=0,Nvisu ; DO i=0,Nvisu
     CALL ExactFunc(AnalyzeExactFunc,OutputTime,Coords_NVisu(1:3,i,j,k,iElem),U_nVisu(2:2,i,j,k,iElem))
  END DO ; END DO ; END DO
#elif (defined (mhd) || defined (navierstokes))
  IF(PrimVisu)THEN
    DO k=0,Nvisu ; DO j=0,Nvisu ; DO i=0,Nvisu
      cons=U_NVisu(1:PP_nVar,i,j,k,iElem)
      CALL ConsToPrim(U_NVisu(1:PP_nVar,i,j,k,iElem),cons)
    END DO ; END DO ; END DO
  END IF !PrimVisu
#endif /*linadv,navierstokes,mhd*/
  nVars = PP_nVar
#if SHOCK_NFVSE
  nVars = nVars+1
  U_NVisu(nVars,:,:,:,iElem) = alpha(iElem)
#if NFVSE_CORR
  nVars = nVars+1
  U_NVisu(nVars ,:,:,:,iElem) = alpha_old(iElem)
#endif /*NFVSE_CORR*/
#endif /*SHOCK_NFVSE*/
END DO !iElem
CALL VisualizeAny(OutputTime,nOutvars,Nvisu,.FALSE.,Coords_Nvisu,U_Nvisu,FileTypeStr,strvarnames_tmp)
DEALLOCATE(U_NVisu)
DEALLOCATE(Coords_NVisu) 
END SUBROUTINE Visualize

!==================================================================================================================================
!> Use any input to write to paraview
!==================================================================================================================================
SUBROUTINE VisualizeAny(OutputTime,nOutVars,NIn,On_xGP,CoordsIn,Uin,FileTypeStrIn,VarNamesIn)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars,ONLY:Projectname,OutputFormat
USE MOD_Output_Vars,ONLY:nBoundingBoxes,VisuBoundingBox
USE MOD_Mesh_Vars  ,ONLY:nElems
USE MOD_Output_Vars,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_ChangeBasis,ONLY:ChangeBasis3D
!USE MOD_VTK        ,ONLY:WriteDataToVTK3D
USE MOD_VTK        ,ONLY:WriteDataToVTK
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nOutVars,Nin 
LOGICAL,INTENT(IN)            :: On_xGP                    !< if input is on xGP points             
REAL,INTENT(IN)               :: OutputTime               !< simulation time of output
REAL,INTENT(IN)               :: CoordsIn(3,0:Nin,0:Nin,0:Nin,1:nElems) !< solution vector to be visualized
REAL,INTENT(IN)               :: Uin(nOutVars,0:Nin,0:Nin,0:Nin,1:nElems) !< solution vector to be visualized
CHARACTER(LEN=255),INTENT(IN) :: FileTypeStrIn
CHARACTER(LEN=255),INTENT(IN) :: VarNamesIn(nOutVars)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: i,j,k,iElem,jElem,iBox,nVisuElems
LOGICAL,ALLOCATABLE           :: ElemInBoxes(:)
INTEGER,ALLOCATABLE           :: VisuElem(:)
REAL,ALLOCATABLE              :: Coords_NVisu(:,:,:,:,:) 
REAL,ALLOCATABLE              :: U_NVisu(:,:,:,:,:)
CHARACTER(LEN=255)            :: FileString
!==================================================================================================================================
IF(outputFormat.LE.0) RETURN
IF(nBoundingBoxes.GT.0)THEN
  ALLOCATE(ElemInBoxes(nElems))
  ElemInBoxes=.FALSE.
  jElem=0
  DO iElem=1,nElems 
    BBLOOP: DO ibox=1,nBoundingBoxes
      DO k=0,Nin ; DO j=0,Nin; DO i=0,Nin
        IF(ElemInBoxes(iElem)) EXIT BBLOOP !exit if element already in a box
        ElemInBoxes(iElem)=(((CoordsIn(1,i,j,k,iElem)).GT.visuBoundingBox(1,iBox)).AND. &
                            ((CoordsIn(1,i,j,k,iElem)).LT.visuBoundingBox(2,iBox)).AND. &
                            ((CoordsIn(2,i,j,k,iElem)).GT.visuBoundingBox(3,iBox)).AND. &
                            ((CoordsIn(2,i,j,k,iElem)).LT.visuBoundingBox(4,iBox)).AND. &
                            ((CoordsIn(3,i,j,k,iElem)).GT.visuBoundingBox(5,iBox)).AND. &
                            ((CoordsIn(3,i,j,k,iElem)).LT.visuBoundingBox(6,iBox)) )
      END DO; END DO; END DO !i,j,k
    END DO BBLOOP
  END DO !iElem
  nVisuElems=COUNT(ElemInBoxes)
  ALLOCATE(VisuElem(nVisuElems))
  jElem=0
  DO iElem=1,nElems
    IF(ElemInBoxes(iElem))THEN
      jElem=jElem+1 
      VisuElem(jElem)=iElem
    END IF
  END DO
  DEALLOCATE(ElemInBoxes)
ELSE
  nVisuElems=nElems
  ALLOCATE(VisuElem(nVisuElems))
  DO iElem=1,nElems
    VisuElem(iElem)=iElem
  END DO
END IF
ALLOCATE(Coords_NVisu(1:3,0:NVisu,0:NVisu,0:NVisu,1:nVisuElems))
ALLOCATE(U_NVisu(1:nOutVars,0:NVisu,0:NVisu,0:NVisu,1:nVisuElems))
IF(on_xGP)THEN
  U_NVisu = 0.
  jElem=0
  DO iElem=1,nVisuElems
    ! Create coordinates of visualization points
    CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_NVisu,CoordsIn(1:3,:,:,:,VisuElem(iElem)),Coords_NVisu(1:3,:,:,:,iElem))
    ! Interpolate solution onto visu grid
    CALL ChangeBasis3D(nOutVars,PP_N,NVisu,Vdm_GaussN_NVisu,Uin(1:nOutVars,:,:,:,VisuElem(iElem)),U_NVisu(1:nOutVars,:,:,:,iElem))
  END DO !iElem
ELSE
  !use input
  Coords_NVisu = CoordsIn(1:3,:,:,:,VisuElem(:))
  U_NVisu      = Uin(:,:,:,:,VisuElem(:))
END IF !on_xGP

  
! Visualize data
SELECT CASE(OutputFormat)
CASE(OUTPUTFORMAT_PARAVIEW_3D_SINGLE)
  !old routines
  !FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_singleOLD_'//TRIM(FileTypeStrIn),OutputTime))//'.vtu'
  !CALL WriteDataToVTK3D(        NVisu,nVisuElems,nOutVars,VarNamesIn,Coords_NVisu(1:3,:,:,:,:), &
  !                            U_NVisu,TRIM(FileString))
  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileTypeStrIn),OutputTime)) !without .vtu
  CALL WriteDataToVTK(3,3,(/NVisu,Nvisu,Nvisu/),nVisuElems,nOutVars,VarNamesIn,Coords_NVisu(1:3,:,:,:,:), &
                              U_Nvisu(:,:,:,:,:),TRIM(FileString),MPI_SingleFile=.TRUE.)
CASE(OUTPUTFORMAT_PARAVIEW_3D_MULTI)
  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileTypeStrIn),OutputTime)) !without .vtu
  CALL WriteDataToVTK(3,3,(/NVisu,Nvisu,Nvisu/),nVisuElems,nOutVars,VarNamesIn,Coords_NVisu(1:3,:,:,:,:), &
                              U_Nvisu(:,:,:,:,:),TRIM(FileString),MPI_SingleFile=.FALSE.)
CASE(OUTPUTFORMAT_PARAVIEW_2D_SINGLE)
  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileTypeStrIn)//'2D',OutputTime)) !without .vtu
  CALL WriteDataToVTK(2,3,(/NVisu,Nvisu/),nVisuElems,nOutVars,VarNamesIn,Coords_NVisu(1:3,:,:,0,:), &
                              U_Nvisu(:,:,:,0,:),TRIM(FileString),MPI_SingleFile=.TRUE.)
CASE(OUTPUTFORMAT_PARAVIEW_2D_MULTI)
  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileTypeStrIn)//'2D',OutputTime)) !without .vtu
  CALL WriteDataToVTK(2,3,(/NVisu,Nvisu/),nVisuElems,nOutVars,VarNamesIn,Coords_NVisu(1:3,:,:,0,:), &
                              U_Nvisu(:,:,:,0,:),TRIM(FileString),MPI_SingleFile=.FALSE.)
END SELECT

DEALLOCATE(VisuElem)

DEALLOCATE(U_Nvisu)
DEALLOCATE(Coords_NVisu)

END SUBROUTINE VisualizeAny


!==================================================================================================================================
!> Creates or opens file for structured output of data time series in comma seperated value or Tecplot ASCII format.
!> Default is comma seperated value.
!> Searches the file for a dataset at restart time and tries to resume at this point.
!> Otherwise a new file is created.
!==================================================================================================================================
SUBROUTINE InitOutputToFile(Filename,ZoneName,nVar,VarNames,lastLine)
! MODULES
USE MOD_Globals
USE MOD_Restart_Vars, ONLY: RestartTime
USE MOD_Output_Vars,  ONLY: ProjectName,ASCIIOutputFormat
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileName                 !< file to be written, without data type extension
CHARACTER(LEN=*),INTENT(IN)   :: ZoneName                 !< name of zone (e.g. names of boundary conditions), used for tecplot
INTEGER,INTENT(IN)            :: nVar                     !< number of variables
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVar)           !< variable names to be written
REAL,INTENT(OUT),OPTIONAL     :: lastLine(nVar+1)         !< last written line to search for, when appending to the file 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: stat                         !< File IO status
INTEGER                        :: ioUnit,i
REAL                           :: dummytime                    !< Simulation time read from file
LOGICAL                        :: fileExists                   !< marker if file exists and is valid
CHARACTER(LEN=255)             :: FileName_loc                 ! FileName with data type extension
!==================================================================================================================================
IF(.NOT.MPIRoot) RETURN
IF(PRESENT(lastLine)) lastLine=-HUGE(1.)

! Append data type extension to FileName
IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
  FileName_loc = TRIM(FileName)//'.csv'
ELSE
  FileName_loc = TRIM(FileName)//'.dat'
END IF

! Check for file
INQUIRE(FILE = TRIM(Filename_loc), EXIST = fileExists)
IF(RestartTime.LT.0.0) fileExists=.FALSE.
!! File processing starts here open old and extratct information or create new file.
ioUnit=GETFREEUNIT()

IF(fileExists)THEN ! File exists and append data
  OPEN(UNIT     = ioUnit             , &
       FILE     = TRIM(Filename_loc) , &
       FORM     = 'FORMATTED'        , &
       STATUS   = 'OLD'              , &
       POSITION = 'APPEND'           , &
       RECL     = 50000              , &
       IOSTAT = stat                 )
  IF(stat.NE.0)THEN
    WRITE(UNIT_stdOut,*)' File '//TRIM(FileName_loc)// ' is invalid. Rewriting file...'
    fileExists=.FALSE.
  END IF
END IF

IF(fileExists)THEN
  ! If we have a restart we need to find the position from where to move on.
  ! Read the values from the previous analyse interval, get the CPUtime
  WRITE(UNIT_stdOut,*)' Opening file '//TRIM(FileName_loc)
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'Searching for time stamp...'

  REWIND(ioUnit)
  DO i=1,4
    READ(ioUnit,*,IOSTAT=stat)
    IF(stat.NE.0)THEN
      ! file is broken, rewrite
      fileExists=.FALSE.
      WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' failed. Writing new file.'
      EXIT
    END IF
  END DO
END IF

IF(fileExists)THEN
  ! Loop until we have found the position
  Dummytime = 0.0
  stat=0
  DO WHILE ((Dummytime.LT.RestartTime) .AND. (stat.EQ.0))
    READ(ioUnit,*,IOSTAT=stat) Dummytime
  END DO
  IF(stat.EQ.0)THEN
    ! read final dataset
    IF(PRESENT(lastLine))THEN
      BACKSPACE(ioUnit)
      READ(ioUnit,*,IOSTAT=stat) lastLine
    END IF
    BACKSPACE(ioUnit)
    ENDFILE(ioUnit) ! delete from here to end of file
    WRITE(UNIT_stdOut,'(A,ES15.5)',ADVANCE='YES')' successfull. Resuming file at time ',Dummytime
  ELSE
    WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' failed. Appending data to end of file.'
  END IF
END IF
CLOSE(ioUnit) ! outputfile

IF(.NOT.fileExists)THEN ! No restart create new file
  ioUnit=GETFREEUNIT()
  OPEN(UNIT   = ioUnit             ,&
       FILE   = TRIM(Filename_loc) ,&
       STATUS = 'UNKNOWN'          ,&
       ACCESS = 'SEQUENTIAL'       ,&
       IOSTAT = stat               )
  IF (stat.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR: cannot open '//TRIM(Filename_loc))
  END IF
  ! Create a new file with the CSV or Tecplot header
  IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
    WRITE(ioUnit,'(A)',ADVANCE='NO') 'Time'
    DO i=1,nVar
      WRITE(ioUnit,'(A)',ADVANCE='NO') ','//TRIM(VarNames(i))
    END DO
  ELSE
    WRITE(ioUnit,*)'TITLE="'//TRIM(ZoneName)//','//TRIM(ProjectName)//'"'
    WRITE(ioUnit,'(A)',ADVANCE='NO')'VARIABLES = "Time"'
    DO i=1,nVar
      WRITE(ioUnit,'(A)',ADVANCE='NO') ',"'//TRIM(VarNames(i))//'"'
    END DO
    WRITE(ioUnit,'(A)',ADVANCE='YES')
    WRITE(ioUnit,'(A)') 'ZONE T="'//TRIM(ZoneName)//','//TRIM(ProjectName)//'"'
  END IF
  CLOSE(ioUnit) ! outputfile
END IF
END SUBROUTINE InitOutputToFile


!==================================================================================================================================
!> Outputs formatted data into a text file.
!> Format is either comma seperated value or tecplot.
!==================================================================================================================================
SUBROUTINE OutputToFile(FileName,time,nVar,output)
! MODULES
USE MOD_Globals
USE MOD_Output_Vars,  ONLY: ASCIIOutputFormat
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileName                 !< name of file to be written, without data type extension
INTEGER,INTENT(IN)            :: nVar(2)                  !< 1: number of variables, 2: number of time samples
REAL,INTENT(IN)               :: time(nVar(2))            !< array of output times
REAL,INTENT(IN)               :: output(nVar(1)*nVar(2))  !< array containing one dataset vector per output time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: openStat                !< File IO status
CHARACTER(LEN=50)              :: formatStr               !< format string for the output and Tecplot header
INTEGER                        :: ioUnit,i,j
CHARACTER(LEN=255)             :: FileName_loc            ! FileName with data type extension
!==================================================================================================================================
! Append data type extension to FileName
IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
  FileName_loc = TRIM(FileName)//'.csv'
ELSE
  FileName_loc = TRIM(FileName)//'.dat'
END IF

ioUnit=GETFREEUNIT()
OPEN(UNIT     = ioUnit             , &
     FILE     = TRIM(Filename_loc) , &
     FORM     = 'FORMATTED'        , &
     STATUS   = 'OLD'              , &
     POSITION = 'APPEND'           , &
     RECL     = 50000              , &
     IOSTAT = openStat             )
IF(openStat.NE.0) THEN
  CALL abort(__STAMP__, &
    'ERROR: cannot open '//TRIM(Filename_loc))
END IF
! Choose between CSV and tecplot output format
IF (ASCIIOutputFormat.EQ.ASCIIOUTPUTFORMAT_CSV) THEN
  DO i=1,nVar(2) ! Loop over all output times (equals lines to write)
    ! Write time stamp
    WRITE(ioUnit,'(E23.14E5)',ADVANCE='NO') time(i)
    DO j=1,nVar(1)-1 ! Loop over variables
      WRITE(ioUnit,'(1A,E23.14E5)',ADVANCE='NO') ",",output(nVar(1)*(i-1)+j)
    END DO ! 
    ! New line at last variable
    WRITE(ioUnit,'(1A,E23.14E5)',ADVANCE='YES') ",",output(nVar(1)*(i-1)+nVar(2))
  END DO
ELSE
  ! Create format string for the variable output
  WRITE(formatStr,'(A10,I2,A14)')'(E23.14E5,',nVar(1),'(1X,E23.14E5))'
  DO i=1,nVar(2)
    WRITE(ioUnit,formatstr) time(i),output(nVar(1)*(i-1)+1:nVar(1)*i)
  END DO
END IF
CLOSE(ioUnit) ! outputfile
END SUBROUTINE OutputToFile


!==================================================================================================================================
!> Deallocate arrays of output routines
!==================================================================================================================================
SUBROUTINE FinalizeOutput()
! MODULES
USE MOD_Output_Vars,ONLY:Vdm_GaussN_NVisu,OutputInitIsDone
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(Vdm_GaussN_NVisu)
OutputInitIsDone = .FALSE.
END SUBROUTINE FinalizeOutput

END MODULE MOD_Output
