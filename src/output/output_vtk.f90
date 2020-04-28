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
!> Module for generic data output in vtk xml fromat
!> WARNING: WriteDataToVTK works only for POSTPROCESSING
!==================================================================================================================================
MODULE MOD_VTK
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDataToVTK3D
  MODULE PROCEDURE WriteDataToVTK3D
END INTERFACE

PUBLIC::WriteDataToVTK3D
PUBLIC::WriteDataToVTK
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Subroutine to write 3D point data to VTK format
!==================================================================================================================================
SUBROUTINE WriteDataToVTK3D(NPlot,nElems,nVal,VarNames,Coord,Value,FileString)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nVal                                        !< Number of nodal output variables
INTEGER,INTENT(IN)          :: NPlot                                       !< Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)          :: nElems                                      !< Number of output elements
REAL,INTENT(IN)             :: Coord(3,0:NPlot,0:NPlot,0:NPlot,nElems)     !< CoordsVector
CHARACTER(LEN=*),INTENT(IN) :: VarNames(nVal)                              !< Names of all variables that will be written out
REAL,INTENT(IN)             :: Value(nVal,0:NPlot,0:NPlot,0:NPlot,nElems)  !< Statevector
CHARACTER(LEN=*),INTENT(IN) :: FileString                                  !< Output file name
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iVal,iElem,Offset,nBytes,nVTKPoints,nVTKCells,ivtk
INTEGER            :: nGlobalElems_loc
INTEGER            :: INTdummy
INTEGER            :: NPlot_p1_3,NPlot_p1_2,PointID,CellID,ElemType
INTEGER,ALLOCATABLE:: Vertex(:,:)
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
CHARACTER(LEN=255) :: tmpVarName,tmpVarNameY,tmpVarNameZ
REAL(KIND=8)       :: FLOAT64dummy
INTEGER            :: StrLen,iValVec,nValVec,nVal_loc,VecOffset(0:nVal)
#if MPI
INTEGER            :: iProc,nElems_proc,nElemsMax
REAL,ALLOCATABLE   ::  buf2(:,:,:,:,:)
#endif /*MPI*/
INTEGER            :: nElems_glob(0:nProcessors-1)
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY FILE: " //TRIM(FileString)//" ..."
ivtk=GETFREEUNIT()
NPlot_p1_3=(NPlot+1)**3
NPlot_p1_2=(NPlot+1)**2

#if MPI
CALL MPI_GATHER(nElems,1,MPI_INTEGER,nElems_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
#else
nElems_glob(0) = nElems
#endif

IF(MPIROOT)THEN
  ! here comes the MPI stuff
  nGlobalElems_loc=nElems
#if MPI
  !ALLOCATE receive buffer for Root
  nElemsMax=MAXVAL(nElems_glob)
  ALLOCATE(buf2(3,0:Nplot,0:Nplot,0:Nplot,nElemsMax))
  nGlobalElems_loc=SUM(nElems_glob)
#endif /*MPI*/

  ! Line feed character
  lf = char(10)

  ! Write file
  OPEN(UNIT=ivtk,FILE=TRIM(FileString),STATUS='REPLACE',ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  nVTKPoints=NPlot_p1_3*nGlobalElems_loc
  nVTKCells =NPlot**3  *nGlobalElems_loc
  Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  WRITE(TempStr1,'(I16)')nVTKPoints
  WRITE(TempStr2,'(I16)')nVTKCells
  Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
         &NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset

  !accout for vectors: 
  ! if Variable Name ends with an X and the following have the same name with Y and Z 
  ! then it forms a vector variable (X is omitted for the name) 
  
  iVal=0 !scalars
  iValVec=0 !scalars & vectors
  VecOffset(0)=0
  DO WHILE(iVal.LT.nVal)
    iVal=iVal+1
    iValVec=iValVec+1
    tmpVarName=TRIM(VarNames(iVal)) 
    StrLen=LEN(TRIM(tmpVarName))
    IF(iVal+2.LE.nVal)THEN !variable could be a vector
      tmpVarNameY=TRIM(VarNames(iVal+1)) 
      tmpVarNameZ=TRIM(VarNames(iVal+2)) 
    END IF
    IF((iVal+2.LE.nVal).AND.(INDEX(tmpVarName( StrLen:StrLen),"X").NE.0).AND. &
       (INDEX(tmpVarNameY(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Y").NE.0).AND. &
       (INDEX(tmpVarNameZ(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Z").NE.0))THEN !variable is a vector!
      tmpVarName=tmpVarName(:StrLen-1)
      Buffer='        <DataArray type="Float64" Name="'//TRIM(tmpVarName)//'" NumberOfComponents="3" &
             &format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+SIZEOF_F(INTdummy)+3*nVTKPoints*SIZEOF_F(FLOAT64dummy)
      WRITE(StrOffset,'(I16)')Offset
      VecOffset(iValVec)=VecOffset(iValVec-1)+3
      iVal=iVal+2 !skip the Y & Z components
    ELSE
      Buffer='        <DataArray type="Float64" Name="'//TRIM(tmpVarName)//'" &
             &format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+SIZEOF_F(INTdummy)+nVTKPoints*SIZEOF_F(FLOAT64dummy)
      WRITE(StrOffset,'(I16)')Offset
      VecOffset(iValVec)=VecOffset(iValVec-1)+1
    END IF
  END DO !iVal <=nVal
  nValVec=iValVec



  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+3*nVTKPoints*SIZEOF_F(FLOAT64dummy)
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify necessary cell data
  Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Connectivity
  Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+8*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Offsets
  Buffer='        <DataArray type="Int32" Name="offsets" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF_F(INTdummy)+nVTKCells*SIZEOF_F(INTdummy)
  WRITE(StrOffset,'(I16)')Offset
  ! Elem types
  Buffer='        <DataArray type="Int32" Name="types" format="appended" &
         &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)

END IF !MPIROOT

#if MPI
CALL MPI_BCAST(nValVec,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
CALL MPI_BCAST(vecOffset(0:nValVec),nValVec+1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
#endif /*MPI*/
! Write binary raw data into append section
! Solution data
DO iValVec=1,nValVec
  iVal    = vecOffset(iValVec-1)
  nVal_loc= vecOffset(iValVec)-vecOffset(iValVec-1)
  IF(MPIroot)THEN    
    nBytes = nVTKPoints*SIZEOF_F(FLOAT64dummy)
    WRITE(ivtk) nVal_loc*nBytes
    WRITE(ivtk)REAL(Value(iVal+1:iVal+nVal_loc,0:NPlot,0:NPlot,0:Nplot,1:nElems),8)
#if MPI
    DO iProc=1,nProcessors-1
      IF(nElems_glob(iProc).EQ.0) CYCLE      
      nElems_proc=nElems_glob(iProc)
      CALL MPI_RECV(   buf2(1:nVal_loc,0:NPlot,0:NPlot,0:Nplot,1:nElems_proc), &
                                 nVal_loc*NPlot_p1_3*nElems_proc,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
      WRITE(ivtk) REAL(buf2(1:nVal_loc,0:NPlot,0:NPlot,0:Nplot,1:nElems_proc),8)
    END DO !iProc
  ELSE
    IF(nElems.GT.0)THEN
      CALL MPI_SEND(  Value(iVal+1:iVal+nVal_loc,0:NPlot,0:NPlot,0:Nplot,1:nElems), &
                           nVal_loc*NPlot_p1_3*nElems,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_WORLD,iError)
    END IF !nElems>0
#endif /*MPI*/
  END IF !MPIroot
END DO       ! iValVec=1,nValVec

! Coordinates
IF(MPIRoot)THEN
  nBytes = nVTKPoints*SIZEOF_F(FLOAT64dummy) * 3
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL( Coord(1:3,0:NPlot,0:NPlot,0:Nplot,1:nElems),8)
#if MPI
  DO iProc=1,nProcessors-1
    IF(nElems_glob(iProc).EQ.0) CYCLE      
    nElems_proc=nElems_glob(iProc)
    CALL MPI_RECV(   buf2(1:3,0:NPlot,0:NPlot,0:Nplot,1:nElems_proc), &
                            3*NPlot_p1_3*nElems_proc,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
    WRITE(ivtk) REAL(buf2(1:3,0:NPlot,0:NPlot,0:Nplot,1:nElems_proc),8)
  END DO !iProc
ELSE
  IF(nElems.GT.0)THEN
    CALL MPI_SEND(    Coord(1:3,0:NPlot,0:NPlot,0:Nplot,1:nElems), &
                      3*NPlot_p1_3*nElems,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_WORLD,iError)
  END IF !MPIroot
#endif /*MPI*/
END IF !MPIroot


! Connectivity
IF(MPIROOT)THEN
  PointID=0
  CellID=0
  ALLOCATE(Vertex(8,nVTKCells))
  DO iElem=1,nGlobalElems_loc
    DO k=1,NPlot; DO j=1,NPlot; DO i=1,NPlot
      CellID=CellID+1
      !
      Vertex(:,CellID)=(/                                       &
        PointID+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P4(CGNS=tecplot standard)
        PointID+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P1
        PointID+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P2
        PointID+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P3
        PointID+i+   j   *(NPlot+1)+ k   *NPlot_p1_2-1,      & !P8
        PointID+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P5
        PointID+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P6
        PointID+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2-1      /) !P7
    END DO; END DO; END DO
    PointID=PointID+NPlot_p1_3
  END DO
  nBytes = 8*nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) Vertex(:,:)
  ! Offset
  nBytes = nVTKCells*SIZEOF_F(INTdummy)
  WRITE(ivtk) nBytes
  WRITE(ivtk) (Offset,Offset=8,8*nVTKCells,8)
  ! Elem type
  ElemType = 12 ! VTK_HEXAHEDRON
  WRITE(ivtk) nBytes
  WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
  ! Write footer
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
  SDEALLOCATE(Vertex)
#if MPI
  SDEALLOCATE(buf2)
#endif /*MPI*/
ENDIF
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK3D


!===================================================================================================================================
!> Subroutine to write 2D or 3D unstructured element data to VTK format, with different Nplot(1:dim1) in the directions of the element. 
!> If in the variable name list, 3 consecutive entries have the same name with only the ending being X,Y (vecdim=2) or X,Y,Z(vecdim=3), they are described in paraview as a vector.
!!
!===================================================================================================================================
SUBROUTINE WriteDataToVTK(dim1,vecDim,NPlot,nElems,nVal,VarNames,Coord,Values,FileString_in,MPI_SingleFile)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim1                    !! dimension of the data (either 2=quads or 3=hexas)
INTEGER,INTENT(IN)            :: vecdim                  !! dimension of coordinates 
INTEGER,INTENT(IN)            :: NPlot(dim1)             !! Number of output points per element : (nPlot+1)**dim1
INTEGER,INTENT(IN)            :: nElems                  !! Number of output elements
INTEGER,INTENT(IN)            :: nVal                    !! Number of nodal output variables
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          !! Names of all variables that will be written out
REAL   ,INTENT(IN)            :: Coord(vecdim,1:PRODUCT(Nplot+1),nElems)      ! CoordinatesVector 
REAL   ,INTENT(IN)            :: Values(nVal,1:PRODUCT(Nplot+1),nElems)   !! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString_in           !! Output file name (without .vtu!!)
LOGICAL,INTENT(IN),OPTIONAL   :: MPI_singleFile          !! for MPI, TRUE: gather data of all procs and write single file, FALSE: write one file per proc and an additional link file. default=TRUE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER     :: kindFloat=8  !set floating point accuracy single (4) double (8), should be equal or lower than input data!
REAL(KIND=kindFloat)  :: FLOATdummy
CHARACTER(LEN=7)      :: strfloat
INTEGER               :: INTdummy
INTEGER               :: sizefloat,sizeInt
INTEGER               :: i,j,k,iVal,iElem,Offset,nBytes,nVTKPoints,nVTKCells,ivtk
INTEGER,ALLOCATABLE   :: Vertex(:,:)
INTEGER               :: ProdNplot,ProdNplot_p1,NPlot_p1(dim1),CellID,PointID,ElemType  ! ?
CHARACTER(LEN=35)     :: StrOffset,TempStr1,TempStr2  ! ?
CHARACTER(LEN=300)    :: Buffer,FileString
CHARACTER(LEN=255)    :: tmpVarName,tmpVarNameY,tmpVarNameZ,strabort
INTEGER               :: StrLen,iValVec,nValVec,nVal_loc,VecOffset(0:nVal)
LOGICAL               :: isVector,maybeVector
CHARACTER(LEN=1)      :: strvecdim
CHARACTER(LEN=1)      :: lf
INTEGER               :: nGlobalElems_loc
LOGICAL               :: WriteToFile
#if MPI               
LOGICAL               :: SingleFile
INTEGER               :: iProc,nElems_proc,nElemsMax
REAL,ALLOCATABLE      :: buf2(:,:,:)
#endif /*MPI*/        
INTEGER               :: nElems_glob(0:nProcessors-1)
!===================================================================================================================================
#if MPI
IF(PRESENT(MPI_SingleFile))THEN
  singleFile=MPI_singleFile
ELSE 
  singleFile=.TRUE.
END IF
IF(nProcessors.EQ.1) singleFile=.TRUE. !overwrite input if only 1 proc is used!
IF(SingleFile)THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'  WRITE DATA TO SINGLE XML BINARY (VTU) FILE "'//TRIM(FileString_in)//'.vtu" ...'
ELSE
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'  WRITE DATA TO MULTIPLE XML BINARY (VTU) FILES "'//TRIM(FileString_in)//'_Proc*_.vtu" ...'
END IF
#else
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')'  WRITE DATA TO SINGLE XML BINARY (VTU) FILE "'//TRIM(FileString_in)//'.vtu" ...'
#endif

NPlot_p1  =(Nplot(:)+1)
ProdNPlot  =PRODUCT(Nplot(:))
ProdNPlot_p1  =PRODUCT(Nplot_p1(:))

IF(kindFloat.EQ.4) THEN 
  strfloat='Float32'
ELSEIF(kindFloat.EQ.8)THEN
  strfloat='Float64'
ELSE 
  CALL abort(__STAMP__, &
       'kindFloat not implemented in output vtk')
END IF
sizefloat=SIZEOF_F(FLOATdummy)
sizeInt  =SIZEOF_F(INTdummy)
! Line feed character
lf = char(10)
WRITE(strvecdim,'(I1)') vecdim

#if MPI
CALL MPI_GATHER(nElems,1,MPI_INTEGER,nElems_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
IF(singleFile)THEN
  writeToFile=MPIroot
  IF(MPIROOT)THEN
    !ALLOCATE receive buffer for Root
    nElemsMax=MAXVAL(nElems_glob)
    ALLOCATE(buf2(vecdim,1:ProdNplot_p1,nElemsMax))
    nGlobalElems_loc=SUM(nElems_glob)
    FileString=TRIM(FileString_in)//'.vtu'
  END IF ! mpiroot
ELSE
  !every processor writes his own file
  writeToFile=(nElems.GT.0)
  FileString=TRIM(INTSTAMP(TRIM(FileString_in),myRank))//'_.vtu'
  nElems_glob(0) = nElems
  nGlobalElems_loc=nElems
END IF !singlefile
#else 
writeToFile=(nElems.GT.0)
FileString=TRIM(FileString_in)//'.vtu'
nElems_glob(0) = nElems
nGlobalElems_loc=nElems
#endif

IF(vecdim.LT.dim1) THEN
  WRITE(strabort,*)'WARNING: data dimension should be <= vecdim! dim1= ',dim1,' vecdim= ',vecdim
  CALL abort(__STAMP__, strabort)
END IF
IF(ANY(Nplot.EQ.0)) THEN
  WRITE(strabort,*)'WARNING: nplot= ',nplot,', all dimensions must be >0!'
  CALL abort(__STAMP__, strabort)
END IF
!accout for vectors: 
! if Variable Name ends with an X and the following have the same name with Y and Z 
! then it forms a vector variable (X is omitted for the name) 

iVal=0 !scalars
iValVec=0 !scalars & vectors
VecOffset(0)=0
DO WHILE(iVal.LT.nVal)
  iVal=iVal+1
  iValVec=iValVec+1
  tmpVarName=TRIM(VarNames(iVal)) 
  StrLen=LEN(TRIM(tmpVarName))
  maybeVector=(iVal+vecdim-1.LE.nVal)
  isVector=.FALSE.
  IF(maybeVector)THEN
    SELECT CASE(vecdim)
    CASE(2)
      tmpVarNameY=TRIM(VarNames(iVal+1))
      isVector=((iVal+2.LE.nVal).AND.(INDEX(tmpVarName( StrLen:StrLen),"X").NE.0) &
                                .AND.(INDEX(tmpVarNameY(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Y").NE.0))
    CASE(3)
      tmpVarNameY=TRIM(VarNames(iVal+1))
      tmpVarNameZ=TRIM(VarNames(iVal+2)) 
      isVector=((iVal+2.LE.nVal).AND.(INDEX(tmpVarName( StrLen:StrLen),"X").NE.0) &
                                .AND.(INDEX(tmpVarNameY(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Y").NE.0) &
                                .AND.(INDEX(tmpVarNameZ(:StrLen),TRIM(tmpVarName(:StrLen-1))//"Z").NE.0))
    END SELECT
  END IF !maybevector

  IF(isvector)THEN !variable is a vector!
    VecOffset(iValVec)=VecOffset(iValVec-1)+vecdim
    iVal=iVal+vecdim-1 !skip the Y (& Z) components
  ELSE
    VecOffset(iValVec)=VecOffset(iValVec-1)+1
  END IF !isvector
END DO !iVal <=nVal
nValVec=iValVec


IF(writeToFile)THEN  
  ivtk=GETFREEUNIT()
  ! Write file
  OPEN(UNIT=ivtk,FILE=TRIM(FileString),STATUS='REPLACE',ACCESS='STREAM')

  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  nVTKPoints= ProdNPlot_p1*nGlobalElems_loc
  nVTKCells = ProdNPlot   *nGlobalElems_loc
  Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  WRITE(TempStr1,'(I16)')nVTKPoints
  WRITE(TempStr2,'(I16)')nVTKCells
  Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
         &NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify point data
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset

  DO iValVec=1,nValVec
    iVal    = vecOffset(iValVec-1)
    nVal_loc= vecOffset(iValVec)-vecOffset(iValVec-1)
    tmpVarName=TRIM(VarNames(iVal+1)) 
    IF(nVal_loc.EQ.vecDim)THEN !vector
      StrLen=LEN(TRIM(tmpVarName))
      tmpVarName=tmpVarName(:StrLen-1)
      Buffer='        <DataArray type="'//strfloat//'" Name="'//TRIM(tmpVarName)//'" NumberOfComponents="'//strvecdim// &
             &'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+sizeInt+vecdim*nVTKPoints*sizefloat
      WRITE(StrOffset,'(I16)')Offset
    ELSE !scalar
      Buffer='        <DataArray type="'//strfloat//'" Name="'//TRIM(tmpVarName)// &
             &'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+sizeInt+nVTKPoints*sizeFloat
      WRITE(StrOffset,'(I16)')Offset
    END IF
  END DO

  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify coordinate data
  Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="'//strfloat//'" Name="Coordinates" NumberOfComponents="'//strvecdim// &
  '" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+sizeInt+vecdim*nVTKPoints*sizeFloat
  WRITE(StrOffset,'(I16)')Offset
  Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify necessary cell data
  Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Connectivity
  Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" &
           &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+sizeInt+2**dim1*nVTKCells*sizeInt
  WRITE(StrOffset,'(I16)')Offset
  ! Offset in connectivity data
  Buffer='        <DataArray type="Int32" Name="offsets" format="appended" &
           &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+sizeInt+nVTKCells*sizeInt
  WRITE(StrOffset,'(I16)')Offset
  ! Elem types
  Buffer='        <DataArray type="Int32" Name="types" format="appended" &
           &offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Prepare append section
  Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Write leading data underscore
  Buffer='_';WRITE(ivtk) TRIM(Buffer)

END IF !writeToFile

! Write binary raw data into append section
! Point data
DO iValVec=1,nValVec
  iVal    = vecOffset(iValVec-1)
  nVal_loc= vecOffset(iValVec)-vecOffset(iValVec-1)
  IF(writeToFile)THEN
    nBytes = nVTKPoints*sizeFloat
    WRITE(ivtk) nVal_loc*nBytes
    WRITE(ivtk) REAL(Values(iVal+1:iVal+nVal_loc,:,1:nElems),kindFloat)
  END IF !writeToFile
#if MPI
  IF(SingleFile)THEN
    IF(MPIroot)THEN
      DO iProc=1,nProcessors-1
        nElems_proc=nElems_glob(iProc)
        IF(nElems_proc.EQ.0) CYCLE      
        CALL MPI_RECV(   buf2( 1:nVal_loc,1:ProdNplot_p1,1:nElems_proc), &
                                 nVal_loc  *ProdNPlot_p1  *nElems_proc,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
        WRITE(ivtk) REAL(buf2( 1:nVal_loc,1:ProdNplot_p1,1:nElems_proc),kindFloat)
      END DO !iProc
    ELSE
      IF(nElems.GT.0)THEN
        CALL MPI_SEND(Values(iVal+1:iVal+nVal_loc,1:ProdNplot_p1,1:nElems), &
                                         nVal_loc  *ProdNPlot_p1  *nElems,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_WORLD,iError)
      END IF !nElems>0
    END IF !MPIroot
  END IF!SingleFile
#endif /*MPI*/
END DO !iValVec
! Point coordinates
IF(writeToFile)THEN
  nBytes = nVTKPoints * vecdim*sizeFloat
  WRITE(ivtk) nBytes
  WRITE(ivtk) REAL(Coord(1:vecdim,1:ProdNplot_p1,1:nElems),kindFloat)
END IF !writeToFile
#if MPI
IF(SingleFile)THEN
  IF(MPIroot)THEN
    DO iProc=1,nProcessors-1
      nElems_proc=nElems_glob(iProc)
      IF(nElems_proc.EQ.0) CYCLE      
      CALL MPI_RECV(   buf2(1:vecDim,1:ProdNPlot_p1,1:nElems_proc), &
                              vecDim*  ProdNPlot_p1  *nElems_proc,MPI_DOUBLE_PRECISION,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
      WRITE(ivtk) REAL(buf2(1:vecDim,1:ProdNPlot_p1,1:nElems_proc),kindFloat)
    END DO !iProc
  ELSE
    IF(nElems.GT.0)THEN
      CALL MPI_SEND(Coord(1:vecDim,1:ProdNPlot_p1,1:nElems), &
                            vecDim  *ProdNPlot_p1  *nElems,MPI_DOUBLE_PRECISION, 0,0,MPI_COMM_WORLD,iError)
    END IF !nElems>0
  END IF !MPIroot
END IF !SingleFile
#endif /*MPI*/
! Connectivity
IF(writeToFile)THEN
  ALLOCATE(Vertex(2**dim1,nVTKCells))
  SELECT CASE(dim1)
  CASE(2)
    CellID = 0
    PointID= 0
    DO iElem=1,nGlobalElems_loc
      DO j=1,NPlot(2)
        DO i=1,NPlot(1)
          CellID = CellID+1
          !visuQuadElem
          Vertex(:,CellID) = (/                  &
            PointID+(i-1)+  j   * NPlot_p1(1) ,    & !P4
            PointID+(i-1)+ (j-1)* NPlot_p1(1) ,    & !P1(CGNS=tecplot standard)
            PointID+ i   + (j-1)* NPlot_p1(1) ,    & !P2
            PointID+ i   +  j   * NPlot_p1(1)     /) !P3
        END DO
      END DO
      PointID=PointID+ProdNPlot_p1
    END DO
  CASE(3)
    CellID=0
    PointID=0
    DO iElem=1,nGlobalElems_loc
      DO k=1,NPlot(3)
        DO j=1,NPlot(2)
          DO i=1,NPlot(1)
            CellID=CellID+1
            !
            Vertex(:,CellID)=(/                                       &
              PointID+(i-1)+( j   +(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P4(CGNS=tecplot standard)
              PointID+(i-1)+((j-1)+(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P1
              PointID+ i   +((j-1)+(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P2
              PointID+ i   +( j   +(k-1)*NPlot_p1(2))*NPlot_p1(1),      & !P3
              PointID+(i-1)+( j   + k   *NPlot_p1(2))*NPlot_p1(1),      & !P8
              PointID+(i-1)+((j-1)+ k   *NPlot_p1(2))*NPlot_p1(1),      & !P5
              PointID+ i   +((j-1)+ k   *NPlot_p1(2))*NPlot_p1(1),      & !P6
              PointID+ i   +( j   + k   *NPlot_p1(2))*NPlot_p1(1)      /) !P7
          END DO
        END DO
      END DO
      !
      PointID=PointID+ProdNPlot_p1
    END DO
  END SELECT
  nBytes = 2**dim1*nVTKCells*sizeInt
  WRITE(ivtk) nBytes
  WRITE(ivtk) Vertex(:,:)
  DEALLOCATE(Vertex)
  ! Offset in connectivity
  nBytes = nVTKCells*sizeInt
  WRITE(ivtk) nBytes
  WRITE(ivtk) (Offset,Offset=2**dim1,2**dim1*nVTKCells,2**dim1)
  ! Elem type
  ElemType =3+3*dim1 !9 VTK_QUAD 12  VTK_HEXAHEDRON
  WRITE(ivtk) nBytes
  WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
  ! Write footer
  Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
END IF!writeToFile
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"

#if MPI
IF((.NOT.SingleFile).AND.(MPIroot))THEN
  FileString=TRIM(FileString_in)//'.pvtu'
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   LINKING VTK FILES INTO "//TRIM(FileString)//" ..."
  ivtk=GETFREEUNIT()    
  OPEN(UNIT=ivtk,FILE=TRIM(FileString),STATUS='REPLACE',ACCESS='STREAM')
  ! Write header
  Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify file type
  Buffer='  <PUnstructuredGrid GhostLevel="0">'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    <PPointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  DO iValVec=1,nValVec
    iVal    = vecOffset(iValVec-1)
    nVal_loc= vecOffset(iValVec)-vecOffset(iValVec-1)
    tmpVarName=TRIM(VarNames(iVal+1)) 
    IF(nVal_loc.EQ.vecDim)THEN !vector
      StrLen=LEN(TRIM(tmpVarName))
      tmpVarName=tmpVarName(:StrLen-1)
      Buffer='        <DataArray type="'//strfloat//'" Name="'//TRIM(tmpVarName)//'" NumberOfComponents="'//strvecdim//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    ELSE !scalar
      Buffer='        <DataArray type="'//strfloat//'" Name="'//TRIM(tmpVarName)//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    END IF
  END DO
  !DO iVal=1,nVal
  !  Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  !END DO
  Buffer='    </PPointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    <PCellData> </PCellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    <PPoints>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='        <DataArray type="'//strfloat//'" Name="Coordinates" NumberOfComponents="'//strvecdim//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='    </PPoints>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Link files
  DO iProc=0,nProcessors-1
    IF(nElems_glob(iProc).EQ.0) CYCLE      
    FileString=TRIM(INTSTAMP(TRIM(FileString_in),iProc))//'_.vtu'
    Buffer='    <Piece Source="'//TRIM(FileString)//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  END DO
  ! Write footer
  Buffer='  </PUnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
  Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
  CLOSE(ivtk)
  WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END IF !not singleFile & MPIroot
#endif 
END SUBROUTINE WriteDataToVTK

END MODULE MOD_VTK
