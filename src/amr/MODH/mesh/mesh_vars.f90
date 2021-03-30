!===================================================================================================================================
! Copyright (c) 2018 - 2020 Alexander Astanin
!
! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
!===================================================================================================================================
#include "../hopest_f.h"
MODULE MODH_Mesh_Vars
    !===================================================================================================================================
    ! Contains global variables provided by the mesh routines
    !===================================================================================================================================
    ! MODULES
    USE, INTRINSIC :: ISO_C_BINDING
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    PUBLIC
    SAVE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! basis
    !-----------------------------------------------------------------------------------------------------------------------------------

    INTEGER :: NGeo_out                    ! polynomial degree of geometric transformation for output
    REAL, ALLOCATABLE :: XGeo(:, :, :, :, :)             ! High order geometry nodes, per element (1:3,0:Ngeo,0:Ngeo,0:Ngeo,nTrees)

    !-----------------------------------------------------------------------------------------------------------------------------------
    ! GLOBAL VARIABLES
    !-----------------------------------------------------------------------------------------------------------------------------------
    INTEGER, ALLOCATABLE :: BC(:)
    INTEGER, ALLOCATABLE :: BoundaryType(:, :)
    CHARACTER(LEN = 255), ALLOCATABLE :: BoundaryName(:)
    !-----------------------------------------------------------------------------------------------------------------------------------
    INTEGER :: nGlobalTrees = 0      ! number of elements in mesh
    INTEGER :: nTrees = 0            ! number of local elements
    INTEGER :: offsetElem = 0
    INTEGER :: nGlobalElems = 0      ! number of elements / quadrants in mesh
    INTEGER :: nElems = 0            ! local number of elements/ quadrants
    INTEGER :: nBCs = 0              ! number of BCs in mesh
    INTEGER :: nUserBCs = 0          ! number of BC in inifile

    INTEGER :: nNodes = 0            ! number of nodes in mesh (=nElems*(NGeo+1)**3)
    INTEGER :: nUniqueNodes = 0      ! number of unique nodes in mesh
    INTEGER :: nCurvedNodes = 0      ! number of curved nodes per element = (Ngeo+1)^3
    !-----------------------------------------------------------------------------------------------------------------------------------
    
    ! USER DEFINED TYPES
    TYPE tNodePtr
        TYPE(tNode), POINTER :: np                     ! node pointer
    END TYPE tNodePtr

    TYPE tSidePtr
        TYPE(tSide), POINTER :: sp                     ! side pointer
    END TYPE tSidePtr

    TYPE tElemPtr
        TYPE(tElem), POINTER :: ep                     ! Local element pointer
    END TYPE tElemPtr

    TYPE tElem
        INTEGER :: ind             ! global element index
        INTEGER :: Type            ! element type (linear/bilinear/curved)
        INTEGER :: Zone
        INTEGER :: treeID
        TYPE(tNodePtr) :: Node(8)
        TYPE(tSidePtr) :: Side(6)
    END TYPE tElem

    TYPE tSide
        INTEGER :: ind             ! global side ID
        INTEGER :: sideID          ! local side ID on Proc
        INTEGER :: locSide         ! local side in element [1..6]
        INTEGER :: tmp
        INTEGER :: NbProc
        INTEGER :: BCindex         ! index in BoundaryType array!
        INTEGER :: flip
        INTEGER :: nMortars        ! number of slave mortar sides associated with master mortar
        INTEGER :: MortarType      ! type of mortar: Type1 : 1-4 , Type 2: 1-2 in eta, Type 3: 1-2 in xi
        !TODO: if tSide is small side of a mortar group, mortar type is -1
        TYPE(tNodePtr) :: Node(4)
        TYPE(tElem), POINTER :: Elem
        TYPE(tSide), POINTER :: connection
        TYPE(tSidePtr), POINTER :: MortarSide(:)   ! array of side pointers to slave mortar sides
    END TYPE tSide

    TYPE tNode
        INTEGER :: ind = 0           ! global unique node index
        INTEGER :: tmp = 0
        REAL :: x(3) = 0.
    END TYPE tNode
    !-----------------------------------------------------------------------------------------------------------------------------------
    TYPE(tElemPtr), POINTER :: Trees(:)        ! list of tree elements (coarsest level)
    TYPE(tNodePtr), POINTER :: UniqueNodes(:)  ! pointer array of unique nodes in mesh
    INTEGER, ALLOCATABLE :: HexMap(:, :, :)   ! for input: 0:Ngeo,0:Ngeo,0:Ngeo -> i [0;(Ngeo+1)^3]
    INTEGER, ALLOCATABLE :: HexMapInv(:, :)
    INTEGER, ALLOCATABLE :: HexMap_out(:, :, :) ! for output 0:Ngeo_out,0:Ngeo_out,0:Ngeo_out -> i [0;(Ngeo_out+1)^3]
    ! DATA STRUCTURES BUILT USING P4EST CONNECTIVITY
    TYPE(tElemPtr), POINTER :: Elems(:)        ! new element list elements are "quadrants/octants"
    !===================================================================================================================================
    INTERFACE GETNEWSIDE
        MODULE PROCEDURE GETNEWSIDE
    END INTERFACE

    INTERFACE GETNEWELEM
        MODULE PROCEDURE GETNEWELEM
    END INTERFACE

    INTERFACE deleteMeshPointer
        MODULE PROCEDURE deleteMeshPointer
    END INTERFACE

CONTAINS


    FUNCTION GETNEWSIDE()
        !===================================================================================================================================
        !
        !===================================================================================================================================
        ! MODULES
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        TYPE(tSide), POINTER :: getNewSide
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: iNode
        !===================================================================================================================================
        ALLOCATE(getNewSide)
        DO iNode = 1, 4
            NULLIFY(getNewSide%Node(iNode)%np)
        END DO
        NULLIFY(getNewSide%Elem)
        NULLIFY(getNewSide%MortarSide)
        NULLIFY(getNewSide%connection)
        getNewSide%sideID = 0
        getNewSide%ind = 0
        getNewSide%tmp = 0
        getNewSide%NbProc = -1
        getNewSide%BCindex = 0
        getNewSide%flip = 0
        getNewSide%nMortars = 0
        getNewSide%MortarType = 0
    END FUNCTION GETNEWSIDE

    FUNCTION GETNEWELEM()
        !===================================================================================================================================
        !
        !===================================================================================================================================
        ! MODULES
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        TYPE(tElem), POINTER :: getNewElem
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: iNode, iLocSide
        !===================================================================================================================================
        ALLOCATE(getNewElem)
        DO iNode = 1, 8
            NULLIFY(getNewElem%Node(iNode)%np)
        END DO
        DO iLocSide = 1, 6
            getNewElem%Side(iLocSide)%sp => getNewSide()
        END DO
        getNewElem%ind = 0
        getNewElem%Zone = 0
        getNewElem%Type = 0
    END FUNCTION GETNEWELEM


    SUBROUTINE createSides(Elem)
        !===================================================================================================================================
        ! if element nodes already assigned, create Sides using CGNS standard
        !===================================================================================================================================
        ! MODULES
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        TYPE(tElem), POINTER :: Elem
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        !===================================================================================================================================
        !side 1
        Elem%Side(1)%sp%Node(1)%np => Elem%Node(1)%np
        Elem%Side(1)%sp%Node(2)%np => Elem%Node(4)%np
        Elem%Side(1)%sp%Node(3)%np => Elem%Node(3)%np
        Elem%Side(1)%sp%Node(4)%np => Elem%Node(2)%np
        Elem%Side(1)%sp%elem => Elem
        !side 2
        Elem%Side(2)%sp%Node(1)%np => Elem%Node(1)%np
        Elem%Side(2)%sp%Node(2)%np => Elem%Node(2)%np
        Elem%Side(2)%sp%Node(3)%np => Elem%Node(6)%np
        Elem%Side(2)%sp%Node(4)%np => Elem%Node(5)%np
        Elem%Side(2)%sp%elem => Elem
        !side 3
        Elem%Side(3)%sp%Node(1)%np => Elem%Node(2)%np
        Elem%Side(3)%sp%Node(2)%np => Elem%Node(3)%np
        Elem%Side(3)%sp%Node(3)%np => Elem%Node(7)%np
        Elem%Side(3)%sp%Node(4)%np => Elem%Node(6)%np
        Elem%Side(3)%sp%elem => Elem
        !side 4
        Elem%Side(4)%sp%Node(1)%np => Elem%Node(3)%np
        Elem%Side(4)%sp%Node(2)%np => Elem%Node(4)%np
        Elem%Side(4)%sp%Node(3)%np => Elem%Node(8)%np
        Elem%Side(4)%sp%Node(4)%np => Elem%Node(7)%np
        Elem%Side(4)%sp%elem => Elem
        !side 5
        Elem%Side(5)%sp%Node(1)%np => Elem%Node(1)%np
        Elem%Side(5)%sp%Node(2)%np => Elem%Node(5)%np
        Elem%Side(5)%sp%Node(3)%np => Elem%Node(8)%np
        Elem%Side(5)%sp%Node(4)%np => Elem%Node(4)%np
        Elem%Side(5)%sp%elem => Elem
        !side 6
        Elem%Side(6)%sp%Node(1)%np => Elem%Node(5)%np
        Elem%Side(6)%sp%Node(2)%np => Elem%Node(6)%np
        Elem%Side(6)%sp%Node(3)%np => Elem%Node(7)%np
        Elem%Side(6)%sp%Node(4)%np => Elem%Node(8)%np
        Elem%Side(6)%sp%elem => Elem
    END SUBROUTINE createSides


    SUBROUTINE deleteMeshPointer()
        !===================================================================================================================================
        ! Deallocates all pointers used for the mesh readin
        !===================================================================================================================================
        ! MODULES
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: iTree, iElem, iLocSide, iNode
        TYPE(tElem), POINTER :: aTree, aElem
        TYPE(tSide), POINTER :: aSide
        !===================================================================================================================================
        IF(ASSOCIATED(Trees))THEN
            DO iTree = 1, nTrees
                aTree => Trees(iTree)%ep
                DO iLocSide = 1, 6
                    aSide => aTree%Side(iLocSide)%sp
                    DEALLOCATE(aSide)
                END DO
                DEALLOCATE(aTree)
            END DO
            DEALLOCATE(Trees)
            DO iNode = 1, nUniqueNodes
                IF(ASSOCIATED(UniqueNodes(iNode)%np))THEN
                    DEALLOCATE(UniqueNodes(iNode)%np)
                END IF
            END DO
            DEALLOCATE(UniqueNodes)
        END IF
        IF(ASSOCIATED(Elems))THEN
            DO iElem = 1, nElems
                aElem => Elems(iElem)%ep
                DO iLocSide = 1, 6
                    aSide => aElem%Side(iLocSide)%sp
                    DEALLOCATE(aSide)
                END DO
                DEALLOCATE(aElem)
            END DO
            DEALLOCATE(Elems)
        END IF
    END SUBROUTINE deleteMeshPointer


END MODULE MODH_Mesh_Vars
