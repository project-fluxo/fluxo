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

MODULE MODH_Mesh
    !===================================================================================================================================
    ! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
    !===================================================================================================================================
    ! MODULES
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    PRIVATE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! GLOBAL VARIABLES (PUBLIC)
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! Public Part ----------------------------------------------------------------------------------------------------------------------

    INTERFACE InitMesh
        MODULE PROCEDURE InitMesh
    END INTERFACE

    INTERFACE FinalizeMesh
        MODULE PROCEDURE FinalizeMesh
    END INTERFACE

    PUBLIC :: InitMesh
    PUBLIC :: FinalizeMesh
    !===================================================================================================================================

CONTAINS


    ! Actually we don't need it, ProjectName and MeshFile have Already been defined
    SUBROUTINE InitMesh()
        !===================================================================================================================================
        ! Read Parameter from inputfile
        !===================================================================================================================================
        ! MODULES
        USE MOD_Globals
        USE MOD_ReadInTools, ONLY : GETSTR
        USE MOD_Mesh_Vars, ONLY : MeshFile
        ! USE MOD_Output_vars,       ONLY: ProjectName

        IMPLICIT NONE
        ! INPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT/OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !-----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        !===================================================================================================================================
        SWRITE(UNIT_StdOut, '(132("-"))')
        SWRITE(UNIT_stdOut, '(A)') ' INIT AMR MESH ...'
        ! prepare pointer structure (get nTrees, etc.)

        MeshFile = GETSTR('MeshFile')

        ! Deform = GETINT('Deform','0')


        ! doSplineInterpolation = GETLOGICAL('doSplineInterpolation','.FALSE.')

        SWRITE(UNIT_stdOut, '(A)')' INIT AMR MESH DONE!'
        SWRITE(UNIT_StdOut, '(132("-"))')
    END SUBROUTINE InitMesh

    SUBROUTINE FinalizeMesh()
        !============================================================================================================================
        ! Deallocate all global interpolation variables.
        !============================================================================================================================
        ! MODULES
        ! USE MOD_Globals
        USE MODH_Mesh_Vars
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !----------------------------------------------------------------------------------------------------------------------------
        !input parameters
        !----------------------------------------------------------------------------------------------------------------------------
        !output parameters
        !----------------------------------------------------------------------------------------------------------------------------
        !local variables
        INTEGER :: iTree, iLocSide, iNode
        !============================================================================================================================
        ! Deallocate global variables, needs to go somewhere else later
        DO iTree = 1, nTrees
            DO iLocSide = 1, 6
                DEALLOCATE(Trees(iTree)%ep%Side(iLocSide)%sp)
            END DO
            DEALLOCATE(Trees(iTree)%ep)
        END DO
        DEALLOCATE(Trees)
        DO iNode = 1, nUniqueNodes
            ADEALLOCATE(UniqueNodes(iNode)%np)
        END DO
        DEALLOCATE(UniqueNodes)
        SDEALLOCATE(XGeo)
        SDEALLOCATE(HexMap)
        SDEALLOCATE(HexMap_out)
        SDEALLOCATE(HexMapInv)
        SDEALLOCATE(BoundaryName)
        SDEALLOCATE(BoundaryType)
    END SUBROUTINE FinalizeMesh

END MODULE MODH_Mesh
