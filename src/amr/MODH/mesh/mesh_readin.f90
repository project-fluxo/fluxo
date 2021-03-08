#include "../hopest_f.h"

MODULE MODH_Mesh_ReadIn
    !===================================================================================================================================
    ! Add comments please!
    !===================================================================================================================================
    ! MODULES
    USE MOD_HDF5_Input
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! Private Part ---------------------------------------------------------------------------------------------------------------------

    ! Public Part ----------------------------------------------------------------------------------------------------------------------
    INTERFACE ReadMeshFromHDF5
        MODULE PROCEDURE ReadMeshFromHDF5
    END INTERFACE

    INTERFACE ReadMeshHeader
        MODULE PROCEDURE ReadMeshHeader
    END INTERFACE

    PUBLIC :: ReadMeshHeader
    !===================================================================================================================================

CONTAINS

    SUBROUTINE ReadBCs()
        !===================================================================================================================================
        ! Read boundary conditions from data file
        !===================================================================================================================================
        ! MODULES
        USE MOD_Globals
        USE MOD_HDF5_Input
        USE MODH_Mesh_Vars, ONLY : BoundaryName, BoundaryType, nBCs
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: iBC, Offset
        !===================================================================================================================================
        ! Get number of boundary condtions
        CALL GetDataSize(File_ID, 'BCNames', nDims, HSize)
        nBCs = INT(HSize(1), 4)
        DEALLOCATE(HSize)
        CALL GetDataSize(File_ID, 'BCType', nDims, HSize)
        IF(HSize(2).NE.nBCs) STOP 'Problem in readBC'
        DEALLOCATE(HSize)

        ALLOCATE(BoundaryName(nBCs))
        ALLOCATE(BoundaryType(BC_SIZE, nBCs))
        offset = 0
        CALL ReadArray('BCNames', 1, (/nBCs/), Offset, 1, StrArray = BoundaryName)

        CALL ReadArray('BCType', 2, (/BC_SIZE, nBCs/), Offset, 1, IntegerArray = BoundaryType)

        SWRITE(UNIT_StdOut, '(132("."))')
        SWRITE(UNIT_StdOut, '(A,A16,A20,A9,A9,A9,A9)')'BOUNDARY CONDITIONS', '|', 'Name', 'Type', 'Curved', 'State', 'Alpha'
        DO iBC = 1, nBCs
            SWRITE(UNIT_StdOut, '(A,A33,A20,I9,I9,I9,I9)')' |', '|', TRIM(BoundaryName(iBC)), BoundaryType(:, iBC)
        END DO
        SWRITE(UNIT_StdOut, '(132("."))')
    END SUBROUTINE ReadBCs

    ! BC must be set by FLUXO
    SUBROUTINE SetUserBCs()
        !===================================================================================================================================
        ! The user can redefine boundaries in the ini file. We create the mappings for the boundaries.
        !===================================================================================================================================
        ! MODULES
        USE MOD_Globals
        ! USE MODH_Mesh_Vars,  ONLY: BoundaryName,BoundaryType,nBCs,nUserBCs
        USE MOD_Mesh_Vars, ONLY : BoundaryName, BoundaryType, nBCs, nUserBCs

        USE MOD_ReadinTools, ONLY : GETSTR, GETINTARRAY
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: BCMapping(nBCs)
        CHARACTER(LEN = 255), ALLOCATABLE :: BoundaryNameUser(:)
        INTEGER, ALLOCATABLE :: BoundaryTypeUser(:, :)
        INTEGER :: iBC, iUserBC
        !===================================================================================================================================
        ! read in boundary conditions, will overwrite BCs from meshfile!
        ! nUserBCs = CNTSTR('BoundaryName',0)
        ! print *, "=============> HERE! <===================", nUserBCs
        ! IF(nUserBCs.EQ.0)
        RETURN

        ! Read user BC
        ALLOCATE(BoundaryNameUser(nUserBCs))
        ALLOCATE(BoundaryTypeUser(nUserBCs, 2))
        DO iBC = 1, nUserBCs
            BoundaryNameUser(iBC) = GETSTR('BoundaryName')
            BoundaryTypeUser(iBC, :) = GETINTARRAY('BoundaryType', 2) !(/Type,State/)
        END DO

        ! Override BCs
        BCMapping = 0
        DO iBC = 1, nBCs
            DO iUserBC = 1, nUserBCs
                IF(INDEX(TRIM(BoundaryNameUser(iUserBC)), TRIM(BoundaryName(iBC))) .NE.0)THEN
                    SWRITE(Unit_StdOut, '(A,A)')    ' |     Boundary in HDF file found | ', TRIM(BoundaryName(iBC))
                    SWRITE(Unit_StdOut, '(A,I2,I2)')' |                            was | ', BoundaryType(1, iBC), BoundaryType(3, iBC)
                    SWRITE(Unit_StdOut, '(A,I2,I2)')' |                      is set to | ', BoundaryTypeUser(iUserBC, 1:2)
                    BoundaryType(BC_TYPE, iBC) = BoundaryTypeUser(iUserBC, 1)
                    BoundaryType(BC_STATE, iBC) = BoundaryTypeUser(iUserBC, 2)
                END IF
            END DO
        END DO

        SWRITE(UNIT_StdOut, '(132("."))')
        DEALLOCATE(BoundaryNameUser, BoundaryTypeUser)
    END SUBROUTINE SetUserBCs


    SUBROUTINE ReadMeshHeader(FileString)
        !===================================================================================================================================
        ! Subroutine to read the mesh from a mesh data file
        !===================================================================================================================================
        ! MODULES
        USE MOD_Globals
        USE MOD_HDF5_Input
#if USE_AMR
        USE MOD_AMR_Vars,           ONLY: p4estFileExist
#endif /*USE_AMR*/
        USE MOD_Globals, ONLY : myrank, MPIRoot
        USE MODH_Mesh_Vars, ONLY : nGlobalTrees
        USE MOD_Mesh_Vars, ONLY : nElems, Ngeo, isMortarMesh
        ! USE MODH_Mesh,     ONLY: SetCurvedInfo
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        CHARACTER(LEN = *), INTENT(IN) :: FileString
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        LOGICAL :: isMesh, dsExists
        INTEGER :: NGEO1, nGlobalTrees1, iMortar
        !===================================================================================================================================
        CALL CheckIfMesh(FileString, isMesh)
        IF(.NOT.isMesh) CALL abort(__STAMP__, &
                'ERROR: Given file ' // TRIM(FileString) // ' is no valid mesh file.')

        ! print *, "READ IN PROGRESS"
        !Create a Communicator with MPI root processor
            

       
       
        CALL OpenDataFile(FileString, create = .FALSE., single = .FALSE., readOnly = .TRUE.)
        CALL ReadAttribute(File_ID, 'nElems', 1, IntegerScalar = nGlobalTrees1) !global number of elements
        CALL ReadAttribute(File_ID, 'Ngeo', 1, IntegerScalar = NGeo1)

#if USE_AMR
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
#endif /*USE_AMR*/

        nGlobalTrees = nGlobalTrees1
        nElems = nGlobalTrees1
        nGeo = Ngeo1
        ! CALL SetCurvedInfo()
        ! useCurveds=.TRUE. ! TODO: maybe implement as optional parameter

        CALL GetDataSize(File_ID, 'NodeCoords', nDims, HSize)
        ! nDims=3
        IF(HSize(2).NE.(NGeo + 1)**3 * nGlobalTrees) CALL abort(__STAMP__, &
                'ERROR: Number of nodes in NodeCoords is not consistent with nTrees and NGeo.')
        DEALLOCATE(HSize)

        ! RETURN

        CALL readBCs()

        CALL setUserBCs()

        CALL CloseDataFile()

    END SUBROUTINE ReadMeshHeader


    SUBROUTINE ReadMeshFromHDF5(FileString)
        !===================================================================================================================================
        ! Subroutine to read the mesh from a mesh data file and build p4est_connectivity
        !===================================================================================================================================
        ! MODULES
        USE MOD_Globals, ONLY : MPIRoot
        USE MOD_HDF5_Input
        USE MOD_Globals
        USE MODH_Mesh_Vars
        USE MOD_Mesh_Vars, ONLY : NGeo
        USE MOD_AMR_Vars, ONLY : H2P_VertexMap, H2P_FaceMap
        USE MOD_AMR_Vars, ONLY : connectivity_ptr
        ! USE MODH_P4EST_Binding, ONLY: p4_connectivity_treevertex
        USE MOD_P4EST, ONLY : p4_connectivity_treevertex, p4_build_bcs
        ! USE MODH_P4EST,           ONLY: getHFlip
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        CHARACTER(LEN = *), INTENT(IN) :: FileString
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: i, j, k, l
        INTEGER :: FirstElemInd, LastElemInd
        INTEGER :: FirstSideInd, LastSideInd
        INTEGER :: nbLocSide, offsetSideID
        INTEGER :: BCindex
        INTEGER :: iElem, iTree, ElemID, iNode
        INTEGER :: iLocSide, iSide
        INTEGER :: nSideIDs
        INTEGER :: C2V(3, 8) ! Volume to 1D corner node mapping
        INTEGER, ALLOCATABLE :: ElemInfo(:, :), SideInfo(:, :), GlobalNodeIDs(:, :, :, :)
        LOGICAL :: oriented, fileExists, doConnection
        TYPE(tElem), POINTER :: aElem
        TYPE(tSide), POINTER :: aSide, bSide
        TYPE(tNode), POINTER :: aNode
        ! p4est interface
        INTEGER :: num_vertices, num_trees, num_periodics, iPeriodic
        INTEGER :: PFlip, HFlip, HFlip_test
        INTEGER, ALLOCATABLE :: tree_to_vertex(:, :), JoinFaces(:, :)
        REAL, ALLOCATABLE :: vertices(:, :)
        INTEGER(KIND = C_INT32_T), ALLOCATABLE :: TreeToBC(:, :)
        !===================================================================================================================================
        INQUIRE (FILE = TRIM(FileString), EXIST = fileExists)
        IF(.NOT.FileExists)  &
                CALL abort(__STAMP__, &
                        'readMesh from data file "' // TRIM(FileString) // '" does not exist')

        SWRITE(UNIT_stdOut, '(A)')'READ MESH FROM DATA FILE FOR AMR"'//TRIM(FileString)//'" ...'
        SWRITE(UNIT_StdOut, '(132("-"))')

        ! map variables to use as much of Flexi codebase as possible
        offsetElem = 0
        nElems = nGlobalTrees   !local number of Elements
        nTrees = nGlobalTrees   !local number of Elements
        ! Open data file
        CALL OpenDataFile(FileString, create = .FALSE., single = .FALSE., readOnly = .TRUE.)
        ! IF (MPIRoot) THEN
        !----------------------------------------------------------------------------------------------------------------------------
        !                              NODES
        !----------------------------------------------------------------------------------------------------------------------------
        ! Define corner to volume mapping: C2V in CGNS ordering
        C2V(:, 1) = (/0, 0, 0   /)
        C2V(:, 2) = (/NGeo, 0, 0   /)
        C2V(:, 3) = (/NGeo, NGeo, 0   /)
        C2V(:, 4) = (/0, NGeo, 0   /)
        C2V(:, 5) = (/0, 0, NGeo/)
        C2V(:, 6) = (/NGeo, 0, NGeo/)
        C2V(:, 7) = (/NGeo, NGeo, NGeo/)
        C2V(:, 8) = (/0, NGeo, NGeo/)

        ! get unique node mapping
        nNodes = nElems * (NGeo + 1)**3
        ALLOCATE(GlobalNodeIDs(0:NGeo, 0:NGeo, 0:NGeo, nElems))
        ALLOCATE(XGeo(3, 0:NGeo, 0:NGeo, 0:NGeo, nElems))

        CALL ReadArray('GlobalNodeIDs', 1, (/nNodes/), offsetElem * (NGeo + 1)**3, 1, IntegerArray = GlobalNodeIDs)

        CALL ReadArray('NodeCoords', 2, (/3, nNodes/), offsetElem * (NGeo + 1)**3, 2, RealArray = XGeo)

        ! Unique nodes for elements
        CALL ReadAttribute(File_ID, 'nUniqueNodes', 1, IntegerScalar = nUniqueNodes)
        ALLOCATE(UniqueNodes(1:nUniqueNodes))
        DO iNode = 1, nUniqueNodes
            NULLIFY(UniqueNodes(iNode)%np)
        END DO

        num_vertices = 0
        DO iElem = 1, nElems
            DO iNode = 1, 8
                i = C2V(1, iNode);j = C2V(2, iNode);k = C2V(3, iNode);
                l = GlobalNodeIDs(i, j, k, iElem)
                IF(.NOT.ASSOCIATED(UniqueNodes(l)%np))THEN
                    ALLOCATE(UniqueNodes(l)%np)
                    num_vertices = num_vertices + 1
                END IF
                UniqueNodes(l)%np%x = XGeo(:, i, j, k, iElem)
                UniqueNodes(l)%np%ind = l
                UniqueNodes(l)%np%tmp = -1
            END DO
        END DO

        !----------------------------------------------------------------------------------------------------------------------------
        !                              ELEMENTS
        !----------------------------------------------------------------------------------------------------------------------------

        !read local ElemInfo from data file
        FirstElemInd = offsetElem + 1
        LastElemInd = offsetElem + nElems
        ALLOCATE(Trees(FirstElemInd:LastElemInd))
        ALLOCATE(ElemInfo(ElemInfoSize, FirstElemInd:LastElemInd))
        CALL ReadArray('ElemInfo', 2, (/ElemInfoSize, nElems/), offsetElem, 2, IntegerArray = ElemInfo)

        DO iElem = FirstElemInd, LastElemInd
            iSide = ElemInfo(ELEM_FirstSideInd, iElem) !first index -1 in Sideinfo
            Trees(iElem)%ep => GETNEWELEM()
            aElem => Trees(iElem)%ep
            aElem%Ind = iElem
            aElem%Type = ElemInfo(ELEM_Type, iElem)
            aElem%Zone = ElemInfo(ELEM_Zone, iElem)
            ! WARNING: THIS IS CARTESIAN ORDERING NOT CGNS ORDERING AS IN FLEXI !!!
            DO iNode = 1, 8
                i = C2V(1, iNode);j = C2V(2, iNode);k = C2V(3, iNode);
                aElem%Node(iNode)%np => UniqueNodes(GlobalNodeIDs(i, j, k, iElem))%np
            END DO
            CALL createSides(aElem)
        END DO

        DEALLOCATE(GlobalNodeIDs)

        !----------------------------------------------------------------------------------------------------------------------------
        !                              SIDES
        !----------------------------------------------------------------------------------------------------------------------------

        offsetSideID = ElemInfo(ELEM_FirstSideInd, FirstElemInd) ! hdf5 array starts at 0-> -1
        nSideIDs = ElemInfo(ELEM_LastSideInd, LastElemInd) - ElemInfo(ELEM_FirstSideInd, FirstElemInd)
        !read local SideInfo from data file
        FirstSideInd = offsetSideID + 1
        LastSideInd = offsetSideID + nSideIDs
        ALLOCATE(SideInfo(SideInfoSize, FirstSideInd:LastSideInd))
        CALL ReadArray('SideInfo', 2, (/SideInfoSize, nSideIDs/), offsetSideID, 2, IntegerArray = SideInfo)

        DO iElem = FirstElemInd, LastElemInd
            aElem => Trees(iElem)%ep
            iSide = ElemInfo(ELEM_FirstSideInd, iElem) !first index -1 in Sideinfo
            !build up sides of the element according to CGNS standard
            ! assign flip
            DO iLocSide = 1, 6
                aSide => aElem%Side(iLocSide)%sp
                iSide = iSide + 1
                aSide%Elem => aElem
                oriented = (Sideinfo(SIDE_ID, iSide).GT.0)
                aSide%Ind = ABS(SideInfo(SIDE_ID, iSide))
                IF(oriented)THEN !oriented side
                    aSide%flip = 0
                ELSE !not oriented
                    aSide%flip = MOD(Sideinfo(SIDE_Flip, iSide), 10)
                    IF((aSide%flip.LT.0).OR.(aSide%flip.GT.4)) STOP 'NodeID doesnt belong to side'
                END IF

                ! Check if mortar element
                ElemID = SideInfo(SIDE_nbElemID, iSide) !IF nbElemID <0, this marks a mortar master side.
                ! The number (-1,-2,-3) is the Type of mortar
                IF(ElemID.LT.0)THEN ! mortar Sides attached!
                    CALL abort(__STAMP__, &
                            'Only conforming meshes in readin.')
                END IF
            END DO !i=1,locnSides
        END DO !iElem

        ! build up side connection
        DO iElem = FirstElemInd, LastElemInd
            aElem => Trees(iElem)%ep
            iSide = ElemInfo(ELEM_FirstSideInd, iElem) !first index -1 in Sideinfo
            DO iLocSide = 1, 6
                aSide => aElem%Side(iLocSide)%sp
                iSide = iSide + 1
                elemID = SideInfo(SIDE_nbElemID, iSide)
                BCindex = SideInfo(SIDE_BCID, iSide)
                doConnection = .TRUE. ! for periodic sides if BC is reassigned as non periodic
                IF(BCindex.NE.0)THEN !BC
                    aSide%BCindex = BCindex
                    IF((BoundaryType(BC_TYPE, aSide%BCindex).NE.1).AND.&
                            (BoundaryType(BC_TYPE, aSide%BCindex).NE.100))THEN ! Reassignement from periodic to non-periodic
                        doConnection = .FALSE.
                        aSide%flip = 0
                    END IF
                ELSE
                    aSide%BCindex = 0
                END IF

                IF(.NOT.doConnection) CYCLE
                IF(ASSOCIATED(aSide%connection)) CYCLE

                ! check if neighbor on local proc or MPI connection
                IF((elemID.LE.LastElemInd).AND.(elemID.GE.FirstElemInd))THEN !local
                    nbLocSide = Sideinfo(SIDE_Flip, iSide) / 10
                    IF((nbLocSide.LT.1).OR.(nbLocSide.GT.6))&
                            CALL abort(__STAMP__, 'SideInfo: Index of local side must be between 1 and 6!')
                    bSide => Trees(elemID)%ep%Side(nbLocSide)%sp
                    aSide%connection => bSide
                    bSide%connection => aSide
                    IF(bSide%ind.NE.aSide%ind)&
                            CALL abort(__STAMP__, &
                                    'SideInfo: Index of side and neighbor side have to be identical!')
                ELSE !MPI
                    ! #ifdef MPI
                    !       aSide%connection=>GETNEWSIDE()
                    !       aSide%connection%flip=aSide%flip
                    !       aSide%connection%Elem=>GETNEWELEM()
                    !       aSide%NbProc = ELEMIPROC(elemID)
                    ! #else
                    !       CALL abort(__STAMP__, &
                    !         ' elemID of neighbor not in global Elem list ')
                    ! #endif
                END IF
            END DO !iLocSide
        END DO !iElem

        DEALLOCATE(ElemInfo, SideInfo)
        ! ENDIF ! MPIRoot
        CALL CloseDataFile()
        ! IF (MPIRoot) THEN


        !----------------------------------------------------------------------------------------------------------------------------
        !                  P4EST MESH CONNECTIVITY
        !----------------------------------------------------------------------------------------------------------------------------
        ! should be replaced by connectivity information ?
        ! needs unique corner nodes for mesh connectivity

        num_vertices = 0
        DO iTree = 1, nTrees
            aElem => Trees(iTree)%ep
            DO iNode = 1, 8
                aNode => aElem%Node(iNode)%np
                IF(aNode%tmp.EQ.-1)THEN
                    num_vertices = num_vertices + 1
                    aElem%Node(iNode)%np%tmp = num_vertices
                END IF
            END DO
        END DO !iTree

        ALLOCATE(Vertices(3, num_vertices))
        DO iNode = 1, nUniqueNodes
            aNode => UniqueNodes(iNode)%np
            IF(ASSOCIATED(aNode))THEN ! only corner nodes are associated
                IF(aNode%tmp.GT.0) Vertices(:, aNode%tmp) = aNode%x
            END IF
        END DO

        num_trees = nTrees
        ALLOCATE(tree_to_vertex(8, num_trees))
        DO iTree = 1, nTrees
            aElem => Trees(iTree)%ep
            DO iNode = 1, 8
                tree_to_vertex(iNode, iTree) = aElem%Node(H2P_VertexMap(iNode) + 1)%np%tmp - 1
            END DO
        END DO

        !periodic Boundaries
        num_periodics = 0
        DO iTree = 1, nTrees
            aElem => Trees(iTree)%ep
            DO iLocSide = 1, 6
                aSide => aElem%Side(iLocSide)%sp
                IF(aSide%BCIndex.EQ.0) CYCLE ! NO Boundary Condition
                IF((BoundaryType(BC_TYPE, aSide%BCIndex).EQ.1).AND.(aSide%flip.EQ.0))THEN !periodic side
                    num_periodics = num_periodics + 1
                END IF
            END DO !iLocSide
        END DO !iTree

        ALLOCATE(JoinFaces(5, num_periodics))
        IF(num_periodics.GT.0) THEN
            DO iTree = 1, nTrees
                aElem => Trees(iTree)%ep
                DO iLocSide = 1, 6
                    aElem%Side(iLocSide)%sp%tmp = H2P_FaceMap(iLocSide)  !local Face ID in p4est
                END DO
            END DO

            iperiodic = 0
            DO iTree = 1, nTrees
                aElem => Trees(iTree)%ep
                DO iLocSide = 1, 6
                    aSide => aElem%Side(iLocSide)%sp
                    IF(aSide%BCIndex.EQ.0) CYCLE ! NO Boundary Condition
                    IF((BoundaryType(BC_TYPE, aSide%BCIndex).EQ.1).AND.(aSide%flip.EQ.0))THEN !periodic side
                        HFlip = aSide%connection%flip
                        iperiodic = iperiodic + 1
                        bSide => aSide%connection
                        IF(aSide%tmp.GT.bSide%tmp)THEN
                            aSide => aSide%connection
                        END IF
                        bSide => aSide%connection
                        JoinFaces(1, iPeriodic) = aSide%Elem%ind - 1  !treeID of face with smaller p4est locfaceID
                        JoinFaces(2, iPeriodic) = bSide%Elem%ind - 1  ! neighbor tree id
                        JoinFaces(3, iPeriodic) = aSide%tmp         ! p4est locSideID
                        JoinFaces(4, iPeriodic) = bSide%tmp         ! p4est neighbor locSideID
                        DO PFlip = 0, 3
                            Hflip_test = getHflip(aSide%tmp, bSide%tmp, PFlip)
                            IF(HFlip_test.EQ.HFlip) EXIT
                        END DO
                        JoinFaces(5, iPeriodic) = PFlip
                    END IF
                END DO !iLocSide
            END DO !iTree
        END IF !num_periodics>0

        ! IF (MPIRoot) THEN
        CALL p4_connectivity_treevertex(num_vertices, num_trees, vertices, tree_to_vertex, &
                num_periodics, JoinFaces, connectivity_ptr)
        ! ENDIF

        DEALLOCATE(Vertices, tree_to_vertex)
        IF(num_periodics.GT.0) DEALLOCATE(JoinFaces)

        ! INTEGER(KIND=C_INT32_T) :: TreeToBC(0:5,nTrees)
        ALLOCATE(TreeToBC(0:5, nTrees))
        ! Now pack BC to the tree_to_attr in Connectivity
        TreeToBC = -1
        DO iTree = 1, nTrees
            aElem => Trees(iTree)%ep
            DO iSide = 1, 6
                aSide => aElem%side(iSide)%sp
                TreeToBC(H2P_FaceMap(iSide), iTree) = aSide%BCIndex
            END DO
        END DO
        !Print *, TreeToBC(:,1)
        CALL p4_build_bcs(connectivity_ptr, nTrees, TreeToBC)

        DEALLOCATE(TreeToBC)
        ! COUNT SIDES

        ! nBCSides=0
        ! nSides=0
        ! nPeriodicSides=0
        ! DO iTree=1,nTrees
        !   aElem=>Trees(iTree)%ep
        !   DO iLocSide=1,6
        !     aSide=>aElem%Side(iLocSide)%sp
        !     aSide%tmp=0
        !   END DO !iLocSide
        ! END DO !iTree
        ! DO iTree=1,nTrees
        !   aElem=>Trees(iTree)%ep
        !   DO iLocSide=1,6
        !     aSide=>aElem%Side(iLocSide)%sp

        !     IF(aSide%tmp.EQ.0)THEN
        !       nSides=nSides+1
        !       aSide%tmp=-1 !used as marker
        !       IF(ASSOCIATED(aSide%connection)) aSide%connection%tmp=-1
        !       IF(aSide%BCindex.NE.0)THEN !side is BC or periodic side
        !         IF(ASSOCIATED(aSide%connection))THEN
        !           nPeriodicSides=nPeriodicSides+1
        !         ELSE
        !           nBCSides=nBCSides+1
        !         END IF
        !       END IF
        !     END IF
        !   END DO !iLocSide
        ! END DO !iTree


        ! WRITE(*,*)'-------------------------------------------------------'
        ! WRITE(*,'(A22,I8)' )'NGeo:',NGeo
        ! WRITE(*,'(A22,X7L)')'useCurveds:',useCurveds
        ! WRITE(*,'(A22,I8)' )'nTrees:',nTrees
        ! WRITE(*,'(A22,I8)' )'nNodes:',nNodes
        ! WRITE(*,'(A22,I8)' )'nSides:',nSides
        ! WRITE(*,'(A22,I8)' )'nBCSides:',nBCSides
        ! WRITE(*,'(A22,I8)' )'nPeriodicSides:',nPeriodicSides
        ! WRITE(*,*)'-------------------------------------------------------'
        ! ENDIF !MPIRoot
    END SUBROUTINE ReadMeshFromHDF5


    SUBROUTINE CheckIfMesh(MeshFile_in, isMeshFile)
        !===================================================================================================================================
        ! Check if the file is a mesh file
        !===================================================================================================================================
        ! MODULES
        !-----------------------------------------------------------------------------------------------------------------------------------
        USE MOD_HDF5_Input
        IMPLICIT NONE
        ! INPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        CHARACTER(LEN = 255), INTENT(IN) :: MeshFile_in
        ! INPUT/OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        LOGICAL, INTENT(OUT) :: isMeshFile
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        REAL :: FileVersion
        INTEGER :: i, nDataSets
        CHARACTER(LEN = 255) :: DataSets(9)
        !===================================================================================================================================
        CALL OpenDataFile(MeshFile_in, create = .FALSE., single = .FALSE., readOnly = .TRUE.)
        
        ! It is used short version check. P4fluxo created the short version of mesh file
        nDataSets = 5
        DataSets(1) = 'BCNames'
        DataSets(2) = 'BCType'
        DataSets(3) = 'ElemInfo'
        DataSets(4) = 'NodeCoords'
        DataSets(5) = 'SideInfo'
        CALL DatasetExists(File_ID, 'Version', isMeshFile, attrib = .TRUE.)
        IF(isMeshFile)THEN
            CALL ReadAttribute(File_ID, 'Version', 1, RealScalar = FileVersion)
            IF(FileVersion.LT.0.999) isMeshFile = .FALSE.
        END IF
        DO i = 1, nDataSets
            IF(.NOT.isMeshFile) EXIT
            CALL DatasetExists(File_ID, DataSets(i), isMeshFile)
        END DO! i=1,nDataSets
        CALL CloseDataFile()
    END SUBROUTINE CheckIfMesh

    FUNCTION GetHFlip(PSide0, PSide1, PFlip)
        !===================================================================================================================================
        ! transform an p4est orientation (r = PFlip)  in HOPR Flip, using local Side and neighbor side
        !===================================================================================================================================
        ! MODULES
        USE MOD_AMR_Vars, ONLY : P2H_FaceMap, H2P_FaceNodeMap, P2H_FaceNodeMap, P4R, P4Q, P4P
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        INTEGER, INTENT(IN) :: PSide0  ! local side ID of HOPEST
        INTEGER, INTENT(IN) :: PSide1  ! P4EST neighbour local side id
        INTEGER, INTENT(IN) :: PFlip   ! Neighbour side in p4est convention: 0..5
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        INTEGER :: GetHFlip   ! Neighbour side in p4est convention: 0..5
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: HSide0, PNode0, PNode1 ! p4est Side,Flip,Mortar encoding
        !-----------------------------------------------------------------------------------------------------------------------------------
        !1. Get CGNS side from P4 side
        HSide0 = P2H_FaceMap(PSide0)
        !2. First node CGNS -> p4est
        PNode0 = H2P_FaceNodeMap(1, HSide0)
        !3. Get oriented node on neighbour side, formula and matrices see paper Burstedde p4est, 2011
        PNode1 = P4P(P4Q(P4R(PSide0, PSide1), PFlip), PNode0)
        !4. P4EST node -> CGNS
        GetHFlip = P2H_FaceNodeMap(PNode1, PSide1)

    END FUNCTION GetHFlip


END MODULE MODH_Mesh_ReadIn
