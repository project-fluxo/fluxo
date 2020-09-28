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
!> Contains control routines to read high-order meshes,
!> provide mesh data to the solver, 
!> build the metrics, partition the domain.
!==================================================================================================================================
MODULE MOD_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::FinalizeMesh
!==================================================================================================================================

PUBLIC::DefineParametersMesh
CONTAINS


!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateStringOption(  'MeshFile',            "(relative) path to meshfile (mandatory).")
CALL prms%CreateLogicalOption( 'useCurveds',          "Controls usage of high-order information in mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes." //&
                                                      "With mortar meshes, it should always be true!(watertight constraint)" &
                                                      , '.TRUE.')
CALL prms%CreateRealOption(    'meshScale',           "Scale the mesh by this factor (shrink/enlarge).",&
                                                      '1.0')
CALL prms%CreateLogicalOption( 'meshdeform',          "Apply simple sine-shaped deformation on cartesion mesh (for testing).",&
                                                      '.FALSE.')
CALL prms%CreateLogicalOption( 'crossProductMetrics', "Compute mesh metrics using cross product form. Caution: in this case "//&
                                                      "free-stream preservation is only guaranteed for N=3*NGeo.",&
                                                      '.FALSE.')
CALL prms%CreateIntOption(     'debugmesh',           "Output file with visualization and debug information for the mesh. "//&
                                                      "0: no visualization, 1: Paraview vtk",'0')
CALL prms%CreateStringOption(  'BoundaryName',        "Names of boundary conditions to be set (must be present in the mesh!)."//&
                                                      "For each BoundaryName a BoundaryType needs to be specified.",&
                                                      multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType',        "Type of boundary conditions to be set. Format: (BC_TYPE,BC_STATE)",&
                                                      multiple=.TRUE.)
CALL prms%CreateLogicalOption( 'writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
END SUBROUTINE DefineParametersMesh


!==================================================================================================================================
!> Routine controlling the initialization of the mesh.
!>
!> - parameter and mesh reading
!> - domain partitioning
!> - allocation of mesh arrays
!> - build mesh mappings to handle volume/surface operations
!> - compute the mesh metrics
!> - provide mesh metrics for overintegration
!==================================================================================================================================
SUBROUTINE InitMesh()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_HDF5_Input
USE MOD_Interpolation_Vars, ONLY:NodeType,NodeTypeGL,InterpolationInitIsDone
USE MOD_Interpolation,      ONLY:GetVandermonde
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR,GETREAL,GETINT
USE MOD_Metrics,            ONLY:CalcMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Mappings,           ONLY:buildMappings
#if USE_AMR
USE MOD_AMR_Vars,           ONLY: p4estFileExist, UseAMR
#endif /*USE_AMR*/
#if MPI
USE MOD_Prepare_Mesh,       ONLY:exchangeFlip
#endif
USE MOD_IO_HDF5,            ONLY:AddToElemData
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: x(3),meshScale
INTEGER           :: iElem,i,j,k
#if USE_AMR
INTEGER           :: iMortar
#endif /*USE_AMR*/
LOGICAL           :: validMesh, dsExists
INTEGER           :: firstMasterSide     !< lower side ID of array U_master/gradUx_master...
INTEGER           :: lastMasterSide      !< upper side ID of array U_master/gradUx_master...
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'InitMesh not ready to be called or already called.')
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

! prepare pointer structure (get nElems, etc.)
! MeshFile = GETSTR('MeshFile')
IF (ICHAR(Meshfile(1:1))==0) THEN !MeshFile not defined in AMR
  MeshFile = GETSTR('MeshFile')
ENDIF
validMesh = ISVALIDMESHFILE(MeshFile)
IF(.NOT.validMesh) &
    CALL CollectiveStop(__STAMP__,'ERROR - Mesh file not a valid HDF5 mesh.')

useCurveds=GETLOGICAL('useCurveds','.TRUE.')
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
#if USE_AMR

IF (UseAMR) THEN
  iMortar = 0
  dsExists = .FALSE.
  CALL DatasetExists(File_ID,'isMortarMesh',dsExists,.TRUE.)
  IF(dsExists)&
    CALL ReadAttribute(File_ID,'isMortarMesh',1,IntegerScalar=iMortar)
  isMortarMesh=(iMortar.EQ.1)

  IF ((iMortar .EQ. 1) .AND. (.NOT. p4estFileExist)) THEN
    ! Error, we can't Run AMR on Mortar Mesh without p4est file
    CALL CollectiveStop(__STAMP__,&
    "Error, we cannot use AMR on Mortar Mesh without p4est file.")
  ENDIF
ENDIF
#endif /*USE_AMR*/
CALL CloseDataFile()

CALL readMesh(MeshFile) !set nElems

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
CALL setLocalSideIDs()

#if MPI
! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
CALL exchangeFlip()
#endif

! Side Datastructure: RANGES
!-----------------|-----------------|-------------------|
!    U_master     | U_slave         |    Flux           |
!-----------------|-----------------|-------------------|
!  BCsides        |                 |    BCSides        |
!  InnerMortars   |                 |    InnerMortars   | < only the big mortar sides
!  InnerSides     | InnerSides      |    InnerSides     |
!  MPI_MINE sides | MPI_MINE sides  |    MPI_MINE sides |
!                 | MPI_YOUR sides  |    MPI_YOUR sides |
!  MPIMortars     |                 |    MPIMortars     | < only the big mortar sides 
!-----------------|-----------------|-------------------|
!   ...small mortar sides are treated like normal inner or MPI sides

firstBCSide          = 1
firstMortarInnerSide = firstBCSide         +nBCSides
firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
firstMPISide_MINE    = firstInnerSide      +nInnerSides
firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR

lastBCSide           = firstMortarInnerSide-1
lastMortarInnerSide  = firstInnerSide    -1
lastInnerSide        = firstMPISide_MINE -1
lastMPISide_MINE     = firstMPISide_YOUR -1
lastMPISide_YOUR     = firstMortarMPISide-1
lastMortarMPISide    = nSides


firstMasterSide = 1
lastMasterSide  = nSides
firstSlaveSide  = firstInnerSide
lastSlaveSide   = lastMPISide_YOUR
nSidesMaster    = lastMasterSide-firstMasterSide+1
nSidesSlave     = lastSlaveSide -firstSlaveSide+1

LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A25,I8)')   'first/lastMasterSide     ', firstMasterSide,lastMasterSide
LOGWRITE(*,'(A25,I8)')   'first/lastSlaveSide      ', firstSlaveSide, lastSlaveSide
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A25,I8,I8)')'first/lastBCSide         ', firstBCSide         ,lastBCSide
LOGWRITE(*,'(A25,I8,I8)')'first/lastMortarInnerSide', firstMortarInnerSide,lastMortarInnerSide
LOGWRITE(*,'(A25,I8,I8)')'first/lastInnerSide      ', firstInnerSide      ,lastInnerSide
LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_MINE   ', firstMPISide_MINE   ,lastMPISide_MINE
LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_YOUR   ', firstMPISide_YOUR   ,lastMPISide_YOUR
LOGWRITE(*,'(A30,I8,I8)')'first/lastMortarMPISide  ', firstMortarMPISide  ,lastMortarMPISide
LOGWRITE(*,*)'-------------------------------------------------------'

! fill ElemToSide, SideToElem,BC
ALLOCATE(ElemToSide(2,6,nElems))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(BC(1:nBCSides))
ALLOCATE(AnalyzeSide(1:nSides))
ElemToSide  = 0
SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
BC          = 0
AnalyzeSide = 0

!NOTE: nMortarSides=nMortarInnerSides+nMortarMPISides
ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo
ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides)) ! [1]: 1: Neighbour sides, 2: Flip, [2]: small sides
MortarType=-1
MortarInfo=-1

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
CALL fillMeshInfo()

! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

! Build necessary mappings 
CALL buildMappings(PP_N,V2S=V2S,V2S2=V2S2,S2V=S2V,S2V2=S2V2,S2V3=S2V3,CS2V2=CS2V2,FS2M=FS2M)

! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----
! scale and deform mesh if desired (warning: no mesh output!)
meshScale=GETREAL('meshScale','1.0')
IF(ABS(meshScale-1.).GT.1e-14)&
  NodeCoords = NodeCoords*meshScale

IF(GETLOGICAL('meshdeform','.FALSE.'))THEN
  DO iElem=1,nElems
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      x(:)=NodeCoords(:,i,j,k,iElem)
      NodeCoords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))*SIN(PP_Pi*x(3))
    END DO; END DO; END DO;
  END DO
END IF

! volume data
ALLOCATE(Vdm_GLN_N(0:PP_N,0:PP_N))
ALLOCATE(Vdm_N_GLN(0:PP_N,0:PP_N))
CALL GetVandermonde(    PP_N, NodeTypeGL , PP_N, NodeType, Vdm_GLN_N , Vdm_N_GLN , modal=.FALSE. )
ALLOCATE(      Elem_xGP(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(dXGL_N      (3,3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(            sJ(  0:PP_N,0:PP_N,0:PP_N,nElems))
NGeoRef=3*NGeo ! build jacobian at higher degree
ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

ALLOCATE(    Elem_inCyl(nElems))

! surface data
ALLOCATE(      Face_xGP(3,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(       NormVec(3,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(      TangVec1(3,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(      TangVec2(3,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(      SurfElem(  0:PP_N,0:PP_N,1:nSides))


! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
! assign 1/detJ (sJ)
! assign normal and tangential vectors and surfElems on faces

! compute metrics using cross product instead of curl form (warning: no free stream preservation!)
crossProductMetrics=GETLOGICAL('crossProductMetrics','.FALSE.')
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."
CALL CalcMetrics()     ! DG metrics

DEALLOCATE(NodeCoords)

! debugmesh: param specifies format to output, 0: no output, 1: tecplot ascii, 2: tecplot binary, 3: paraview binary
CALL WriteDebugMesh(GETINT('debugmesh','0'))

CALL AddToElemData('myRank',IntScalar=myRank)

Elem_inCyl=.FALSE.

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE FinalizeMesh()
! MODULES
USE MOD_Mesh_Vars
IMPLICIT NONE
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)

! Volume
SDEALLOCATE(Elem_xGP)
SDEALLOCATE(Metrics_fTilde)
SDEALLOCATE(Metrics_gTilde)
SDEALLOCATE(Metrics_hTilde)
SDEALLOCATE(dXGL_N)
SDEALLOCATE(Vdm_GLN_N)
SDEALLOCATE(Vdm_N_GLN)
SDEALLOCATE(sJ)
SDEALLOCATE(detJac_Ref)

! surface
SDEALLOCATE(Face_xGP)
SDEALLOCATE(NormVec)
SDEALLOCATE(TangVec1)
SDEALLOCATE(TangVec2)
SDEALLOCATE(SurfElem)


! mappings
SDEALLOCATE(FS2M)
SDEALLOCATE(V2S)
SDEALLOCATE(V2S2)
SDEALLOCATE(S2V3)
SDEALLOCATE(CS2V2)

MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
