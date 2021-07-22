!==================================================================================================================================
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
!> IO routines to gather informations from HDF5 files
!==================================================================================================================================
MODULE MOD_HDF5_Input
! MODULES
USE MOD_IO_HDF5
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ISVALIDHDF5FILE
  MODULE PROCEDURE ISVALIDHDF5FILE
END INTERFACE

INTERFACE ISVALIDMESHFILE
  MODULE PROCEDURE ISVALIDMESHFILE
END INTERFACE

INTERFACE GetNextFileName
  MODULE PROCEDURE GetNextFileName
END INTERFACE

INTERFACE DatasetExists
  MODULE PROCEDURE DatasetExists
END INTERFACE

INTERFACE GetDataSize
  MODULE PROCEDURE GetDataSize
END INTERFACE

INTERFACE GetDataProps
  MODULE PROCEDURE GetDataProps
END INTERFACE

INTERFACE ReadAttribute
  MODULE PROCEDURE ReadAttribute
END INTERFACE


PUBLIC :: File_ID,HSize,nDims        ! Variables from MOD_IO_HDF5 that need to be public
PUBLIC :: OpenDataFile,CloseDataFile ! Subroutines from MOD_IO_HDF5 that need to be public
PUBLIC :: ISVALIDHDF5FILE,ISVALIDMESHFILE,GetDataSize,GetDataProps,GetNextFileName
PUBLIC :: ReadArray,ReadAttribute
PUBLIC :: DatasetExists
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Subroutine to check if a file is a valid FLUXO HDF5 file
!==================================================================================================================================
FUNCTION ISVALIDHDF5FILE(FileName,FileVersionOpt)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName        !< name of file to be checked
REAL,INTENT(IN),OPTIONAL       :: FileVersionOpt  !< desired version
LOGICAL                        :: isValidHDF5File !< result: file is valid HDF5 file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: FileVersion,FileVersionRef
INTEGER(HID_T)                 :: Plist_ID
CHARACTER(LEN=255)             :: ProgramName
LOGICAL                        :: fileExists
!==================================================================================================================================
isValidHDF5File=.TRUE.
iError=0
FileVersionRef=0.1
IF(PRESENT(FileVersionOpt)) FileVersionRef=FileVersionOpt

! Disable error messages
CALL H5ESET_AUTO_F(0, iError)
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Create property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if MPI
! Setup file access property list with parallel I/O access (MPI)
CALL H5PSET_FAPL_MPIO_F(Plist_ID,MPI_COMM_WORLD, MPIInfo, iError)
#endif /* MPI */

! Check if file exists
INQUIRE(FILE=TRIM(FileName),EXIST=fileExists)
IF(.NOT.fileExists) THEN
  CALL abort(__STAMP__,'ERROR: HDF5 file '//TRIM(FileName)//' does not exist.')
  RETURN
END IF


! Open HDF5 file
CALL H5FOPEN_F(TRIM(FileName), H5F_ACC_RDONLY_F, File_ID, iError,access_prp = Plist_ID)
CALL H5PCLOSE_F(Plist_ID, iError)
IF(iError.EQ.0) THEN
  isValidHDF5File=.TRUE.
  ! Check program name -------------------------------------------------------------------------------------------------------------
  ! Open the attribute "Program" of root group
  CALL ReadAttribute(File_ID,'Program',1,StrScalar=ProgramName)
  IF(TRIM(ProgramName) .NE. 'Fluxo') isValidHDF5File=.FALSE.
  IF (isValidHDF5File) THEN
    ! Check file version -------------------------------------------------------------------------------------------------------------
    ! Open the attribute "File_Version" of root group
    CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersion)
    IF(FileVersion .LT. FileVersionRef)THEN
      isValidHDF5File=.FALSE.
      SWRITE(UNIT_stdOut,'(A)')' ERROR: FILE VERSION TOO OLD! FileName: '//TRIM(FileName)
    END IF
  END IF
  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close the file.
  CALL H5FCLOSE_F(File_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  isValidHDF5File=.FALSE.
  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
END IF
END FUNCTION ISVALIDHDF5FILE

!==================================================================================================================================
!> Subroutine to check if a file is a valid mesh file
!==================================================================================================================================
FUNCTION ISVALIDMESHFILE(MeshFileName)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName    !< name of mesh file to be checked
LOGICAL                        :: isValidMeshFile !< result: file is valid mesh file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: NGeoExists,fileExists
INTEGER(HID_T)                 :: Plist_ID
!==================================================================================================================================
! Disable error messages
CALL H5ESET_AUTO_F(0, iError)

! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Create property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if MPI
! Setup file access property list with parallel I/O access (MPI)
CALL H5PSET_FAPL_MPIO_F(Plist_ID,MPI_COMM_WORLD, MPIInfo, iError)
#endif /* MPI */

! Check if file exists
INQUIRE(FILE=TRIM(MeshFileName),EXIST=fileExists)
IF(.NOT.fileExists) THEN
  CALL abort(__STAMP__,'ERROR: Mesh file '//TRIM(MeshFileName)//' does not exist.')
  isValidMeshFile = .FALSE.
  RETURN
END IF

! Open HDF5 file
CALL H5FOPEN_F(TRIM(MeshFileName), H5F_ACC_RDONLY_F, File_ID, iError,access_prp = Plist_ID)
IF(iError.NE.0) THEN
  CALL abort(__STAMP__,'ERROR: Mesh file '//TRIM(MeshFileName)//' cannot be opened, problem with HDF5.')
  isValidMeshFile = .FALSE.
  RETURN
END IF

IF(iError.EQ.0) THEN
  isValidMeshFile=.TRUE.

  ! Check NGeo attribute --------------------------------------------------------------------------------------------------------
  CALL DatasetExists(File_ID,'Ngeo',NGeoExists,attrib=.TRUE.)
  IF (.NOT.NGeoExists) isValidMeshFile = .FALSE.

  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close the file.
  CALL H5FCLOSE_F(File_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  isValidMeshFile=.FALSE.
  ! Close property list
  CALL H5PCLOSE_F(Plist_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
END IF
END FUNCTION ISVALIDMESHFILE

!==================================================================================================================================
!> Subroutine to determine HDF5 datasize
!==================================================================================================================================
SUBROUTINE GetDataSize(Loc_ID,DSetName,nDims,Size)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName !< name if dataset to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID   !< ID of dataset
INTEGER,INTENT(OUT)                  :: nDims    !< found data size dimensions
INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: Size(:)  !< found data size
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: DSet_ID,FileSpace
INTEGER(HSIZE_T), POINTER            :: SizeMax(:)
!==================================================================================================================================
! Open the dataset with default properties.
CALL H5DOPEN_F(Loc_ID, TRIM(DSetName) , DSet_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
! Get size and max size of data space
ALLOCATE(Size(nDims),SizeMax(nDims))
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Size, SizeMax, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(DSet_ID, iError)
DEALLOCATE(SizeMax)
END SUBROUTINE GetDataSize


!==================================================================================================================================
!> @brief Subroutine to check wheter a dataset on the hdf5 file exists
!>
!> We have no "h5dexists_f", so we use the error given by h5dopen_f.
!> this produces hdf5 error messages even if everything is ok, so we turn the error msgs off
!> during this operation.
!> auto error messages off
!==================================================================================================================================
SUBROUTINE DatasetExists(Loc_ID,DSetName,Exists,attrib)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName !< name if dataset to be checked
INTEGER(HID_T),INTENT(IN)            :: Loc_ID   !< ID of dataset
LOGICAL,INTENT(IN),OPTIONAL          :: attrib   !< check dataset or attribute 
LOGICAL,INTENT(OUT)                  :: Exists   !< result: dataset exists
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: DSet_ID
INTEGER                              :: hdferr
LOGICAL                              :: isAttrib
!==================================================================================================================================
CALL h5eset_auto_f(0, hdferr)
! Open the dataset with default properties.
isAttrib = .FALSE.
IF(PRESENT(attrib)) isAttrib=attrib
IF(isAttrib)THEN
  CALL H5AOPEN_F(Loc_ID, TRIM(DSetName), DSet_ID, iError)
  CALL H5ACLOSE_F(DSet_ID, iError)
ELSE
  CALL H5DOPEN_F(Loc_ID, TRIM(DSetName), DSet_ID, iError)
  CALL H5DCLOSE_F(DSet_ID, iError)
END IF
Exists=.TRUE.
IF(iError.LT.0) Exists=.FALSE.
! auto error messages on
CALL h5eset_auto_f(1, hdferr)
END SUBROUTINE DatasetExists



!==================================================================================================================================
!> Subroutine to determine HDF5 dataset properties in FLUXO terminology
!==================================================================================================================================
SUBROUTINE GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(OUT)                     :: nVar_HDF5     !< number of variables
INTEGER,INTENT(OUT)                     :: N_HDF5        !< polynomial degree
INTEGER,INTENT(OUT)                     :: nElems_HDF5   !< inumber of elements
CHARACTER(LEN=255),OPTIONAL,INTENT(OUT) :: NodeType_HDF5 !< nodetype string
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                 :: Rank
INTEGER(HID_T)                          :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,A)')' GET SIZE OF DATA IN HDF5 FILE... '

! Read in attributes
! Open the dataset with default properties.
CALL H5DOPEN_F(File_ID, 'DG_Solution', Dset_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(Dset_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, Rank, iError)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Rank of database',' | ',Rank,' | HDF5    | '
! Get size and max size of data space
Dims   =0
DimsMax=0
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Dims(1:Rank), DimsMax(1:Rank), iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(Dset_ID, iError)
IF(PRESENT(NodeType_HDF5)) THEN
  ! Read in NodeType
  CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType_HDF5)
END IF

! Display data
! nVar = first array index
!nVar_HDF5 = INT(Dims(1),4)
CHECKSAFEINT(Dims(1),4)
nVar_HDF5 = INT(Dims(1),4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Number of variables nVar',' | ',nVar_HDF5,' | HDF5    | '
! N = index 2-4 of array, is expected to have the same value for each direction
CHECKSAFEINT(Dims(2)-1,4)
N_HDF5 = INT(Dims(2)-1,4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','Polynomial degree N',' | ',N_HDF5,' | HDF5    | '
IF(PRESENT(NodeType_HDF5)) THEN
  SWRITE(UNIT_stdOut,'(A3,A30,A3,A33,A13)')' | ','          Node type',' | ',TRIM(NodeType_HDF5),' | HDF5    | '
END IF
! nElems = index 5 of array
CHECKSAFEINT(Dims(5),4)
nElems_HDF5 = INT(Dims(5),4)
SWRITE(UNIT_stdOut,'(A3,A30,A3,I33,A13)')' | ','GeometricnElems',' | ',nElems_HDF5,' | HDF5    | '

SWRITE(UNIT_stdOut,'(A)')' DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE GetDataProps


!==================================================================================================================================
!> Subroutine to read arrays of rank "Rank" with dimensions "Dimsf(1:Rank)".
!==================================================================================================================================
SUBROUTINE ReadArray(ArrayName,Rank,nVal,Offset_in,Offset_dim,RealArray,IntegerArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER                        :: Rank                  !< number of dimensions of the array
INTEGER                        :: offset_in             !< offset =0, start at beginning of the array
INTEGER                        :: offset_dim            !< which dimension is the offset (only one dimension possible here)
INTEGER                        :: nVal(Rank)            !< size of complete (local) array to write
CHARACTER(LEN=*),INTENT(IN)    :: ArrayName             !< name of array to be read
REAL              ,DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET :: RealArray    !< only if real array shall be read
INTEGER           ,DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET :: IntegerArray !< only if integer array shall be read
CHARACTER(LEN=255),DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT),TARGET :: StrArray     !< only if real string shall be read
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                                                   :: DSet_ID,Type_ID,MemSpace,FileSpace,PList_ID
INTEGER(HSIZE_T)                                                 :: Offset(Rank),Dimsf(Rank)
TYPE(C_PTR)                                                      :: buf
!==================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')'    READ ',Rank,'D ARRAY "',TRIM(ArrayName),'"'
Dimsf=nVal
LOGWRITE(*,*)'Dimsf,Offset=',Dimsf,Offset_in
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
CALL H5DOPEN_F(File_ID, TRIM(ArrayName) , DSet_ID, iError)
! Define and select the hyperslab to use for reading.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
Offset(:)=0
Offset(offset_dim)=Offset_in
CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
! Create property list
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#if MPI
! Set property list to collective dataset read
CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F, iError)
#endif
CALL H5DGET_TYPE_F(DSet_ID, Type_ID, iError)

! Read the data
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(RealArray))THEN
  CALL H5DREAD_F(DSet_ID,Type_ID,RealArray,Dimsf,&
                 iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
END IF
IF(PRESENT(IntegerArray))THEN
  CALL H5DREAD_F(DSet_ID,Type_ID,IntegerArray,Dimsf,&
                 iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DREAD_F(DSet_ID,Type_ID,StrArray,Dimsf,&
                 iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
END IF
#else /*HDF5_F90*/
IF(PRESENT(RealArray))    buf=C_LOC(RealArray)
IF(PRESENT(IntegerArray)) buf=C_LOC(IntegerArray)
IF(PRESENT(StrArray))     buf=C_LOC(StrArray(1))
CALL H5DREAD_F(DSet_ID,Type_ID,buf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
#endif /*HDF5_F90*/

! Close the datatype, property list, dataspaces and dataset.
CALL H5TCLOSE_F(Type_ID, iError)
CALL H5PCLOSE_F(PList_ID,iError)
CALL H5SCLOSE_F(FileSpace,iError)
CALL H5DCLOSE_F(DSet_ID, iError)
CALL H5SCLOSE_F(MemSpace,iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadArray



!==================================================================================================================================
!> Subroutine to read attributes from HDF5 file.
!==================================================================================================================================
SUBROUTINE ReadAttribute(Loc_ID_in,AttribName,nVal,DatasetName,RealScalar,IntegerScalar,&
                                 StrScalar,LogicalScalar,RealArray,IntegerArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)                  :: Loc_ID_in         !< HDF5 file id of opened file
INTEGER           ,INTENT(IN)                  :: nVal              !< number of attributes in case an array is expected
CHARACTER(LEN=*)  ,INTENT(IN)                  :: AttribName        !< name of attribute to be read
CHARACTER(LEN=*)  ,INTENT(IN) ,OPTIONAL        :: DatasetName       !< dataset name in case attribute is located in a dataset
REAL              ,INTENT(OUT),OPTIONAL,TARGET :: RealArray(nVal)   !< Array of real array attributes
INTEGER           ,INTENT(OUT),OPTIONAL,TARGET :: IntegerArray(nVal)!< Array for integer array for attributes
REAL              ,INTENT(OUT),OPTIONAL,TARGET :: RealScalar        !< Scalar real attribute
INTEGER           ,INTENT(OUT),OPTIONAL,TARGET :: IntegerScalar     !< Scalar integer attribute
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL,TARGET :: StrScalar         !< Scalar string attribute
CHARACTER(LEN=255),INTENT(OUT),OPTIONAL,TARGET :: StrArray(nVal)    !< Array for character array attributes
LOGICAL           ,INTENT(OUT),OPTIONAL        :: LogicalScalar     !< Scalar logical attribute
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Attr_ID,Type_ID,Loc_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER                        :: i
INTEGER,TARGET                 :: IntToLog
CHARACTER(LEN=255),TARGET      :: StrTmp(1)
TYPE(C_PTR)                    :: buf
!==================================================================================================================================
LOGWRITE(*,*)' READ ATTRIBUTE "',TRIM(AttribName),'" FROM HDF5 FILE...'
StrTmp=''
Dimsf(1)=nVal
Loc_ID=Loc_ID_in
IF(PRESENT(DatasetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
END IF
! Create scalar data space for the attribute.
! Create the attribute for group Loc_ID.
CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)
CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)

! Nullify
IF(PRESENT(RealArray))     RealArray=0.
IF(PRESENT(RealScalar))    RealScalar=0.
IF(PRESENT(IntegerArray))  IntegerArray=0
IF(PRESENT(IntegerScalar)) IntegerScalar=0
IF(PRESENT(LogicalScalar)) LogicalScalar=.FALSE.
IF(PRESENT(StrScalar))     StrScalar=''
IF(PRESENT(StrArray))THEN
  DO i=1,nVal
    StrArray(i)=''
  END DO
END IF

! Read the attribute data.
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(RealArray))      CALL H5AREAD_F(Attr_ID, Type_ID, RealArray,     Dimsf, iError)
IF(PRESENT(RealScalar))     CALL H5AREAD_F(Attr_ID, Type_ID, RealScalar,    Dimsf, iError)
IF(PRESENT(IntegerArray))   CALL H5AREAD_F(Attr_ID, Type_ID, IntegerArray,  Dimsf, iError)
IF(PRESENT(IntegerScalar))  CALL H5AREAD_F(Attr_ID, Type_ID, IntegerScalar, Dimsf, iError)
IF(PRESENT(LogicalScalar))  CALL H5AREAD_F(Attr_ID, Type_ID, IntToLog,      Dimsf, iError)
IF(PRESENT(StrScalar))      CALL H5AREAD_F(Attr_ID, Type_ID, StrScalar,     Dimsf, iError)
IF(PRESENT(StrArray))       CALL H5AREAD_F(Attr_ID, Type_ID, StrArray,      Dimsf, iError)
#else /* HDF5_F90 */
StrTmp(1) = ''
IF(PRESENT(RealArray))      buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))     buf=C_LOC(RealScalar)
IF(PRESENT(IntegerArray))   buf=C_LOC(IntegerArray)
IF(PRESENT(IntegerScalar))  buf=C_LOC(IntegerScalar)
IF(PRESENT(LogicalScalar))  buf=C_LOC(IntToLog)
IF(PRESENT(StrScalar))      buf=C_LOC(StrTmp(1))
IF(PRESENT(StrArray))       buf=C_LOC(StrArray(1))
CALL H5AREAD_F(Attr_ID, Type_ID, buf, iError)
IF(PRESENT(StrScalar))      StrScalar=StrTmp(1)
#endif /* HDF5_F90 */
IF(PRESENT(LogicalScalar)) LogicalScalar=(IntToLog.EQ.1)

CALL H5TCLOSE_F(Type_ID, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadAttribute



!==================================================================================================================================
!> Subroutine to determine filename of next HDF5 file for FlushFiles
!==================================================================================================================================
SUBROUTINE GetNextFileName(FileName,NextFileName_HDF5,single)
! MODULES
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName                !< filename to check
LOGICAL,INTENT(IN)             :: single                  !< switch whether file is being accessed in parallel my MPI_COMM_WORLD
CHARACTER(LEN=255),INTENT(OUT) :: NextFileName_HDF5       !< output: follow up file according to checked file opened
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: ReadError
INTEGER(HID_T)                 :: File_ID_loc,Plist_ID
!==================================================================================================================================
LOGWRITE(*,*)' GET NEXT FILE NAME FROM HDF5 FILE ', TRIM (FileName),' ...'
ReadError=0
NextFileName_HDF5=''
! Disable error messages
CALL H5ESET_AUTO_F(0, iError)
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Setup file access property list
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#if MPI
IF(.NOT.single)THEN
  ! Set property list to MPI IO
  CALL H5PSET_FAPL_MPIO_F(Plist_ID, MPI_COMM_WORLD, MPI_INFO_NULL, iError)
END IF
#endif /* MPI */
! Open file
CALL H5FOPEN_F(TRIM(FileName), H5F_ACC_RDONLY_F, File_ID_loc, iError,access_prp = Plist_ID)
ReadError=iError
CALL H5PCLOSE_F(Plist_ID, iError)
iError=ReadError
IF (iError .EQ. 0) THEN
  ! Get Name of the mesh file, stored as third atrribute with name "NextFile"
  ! Open the attribute "NextFile" of opened file
  CALL ReadAttribute(File_ID_loc,'NextFile',1,StrScalar=NextFileName_HDF5)
  ! Close the file.
  CALL H5FCLOSE_F(File_ID_loc, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
  iError=-1
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE GetNextFileName

END MODULE MOD_HDF5_Input
