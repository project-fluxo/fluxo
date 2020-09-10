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
#include "amr_f.h"

#if USE_AMR
!==================================================================================================================================
!> Contains control routines for AMR
!==================================================================================================================================
MODULE MOD_P4EST

#if MPI
USE mpi
#endif

USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR, C_F_POINTER, C_LOC
! MODULES
IMPLICIT NONE
PUBLIC
!  PRIVATE
INTERFACE
    FUNCTION p4est_init(COMM) BIND(C, NAME = 'p4est_init_f')
        IMPORT :: C_INT
        INTEGER(C_INT), VALUE :: COMM
        INTEGER(C_INT) :: P4EST_INIT
    END FUNCTION P4EST_INIT

    FUNCTION P8EST_CONNECTIVITY_NEW_PERIODIC() BIND(C, NAME = 'p8est_connectivity_new_periodic_f')
        IMPORT :: C_PTR
        TYPE(C_PTR) :: P8EST_CONNECTIVITY_NEW_PERIODIC
    END FUNCTION P8EST_CONNECTIVITY_NEW_PERIODIC


    FUNCTION P4EST_NEW(CONN) BIND(C, NAME = 'p4est_new_f')
        IMPORT :: C_PTR
        TYPE(C_PTR) :: CONN
        TYPE(C_PTR) :: P4EST_NEW
    END FUNCTION P4EST_NEW


    FUNCTION p4estGetMPIData(p4est) BIND(C, NAME = 'p4estGetMPIData_f')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: p4est
        TYPE(C_PTR) :: p4estGetMPIData
    END FUNCTION p4estGetMPIData

    SUBROUTINE p4estDelMPIData(data_ptr) BIND(C, NAME = 'p4estDelMPIData_f')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: data_ptr
        !TYPE(C_PTR) ::p4est_new
    END SUBROUTINE p4estDelMPIData

    SUBROUTINE p4est_destroy(P4) BIND(C, NAME = 'p4est_destroy_f')
        IMPORT :: C_PTR
        TYPE(C_PTR) :: P4
        !TYPE(C_PTR) ::p4est_new
    END SUBROUTINE p4est_destroy


    SUBROUTINE FillElemsChanges(P4EST, P4fortrandata) BIND(C, NAME = 'FillElemsChange')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4EST
        TYPE(C_PTR), VALUE :: P4fortrandata
    END SUBROUTINE FillElemsChanges


    SUBROUTINE p4est_ResetElementNumber(P4) BIND(C, NAME = 'ResetElementNumber')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4
        !TYPE(C_PTR) ::p4est_new
    END SUBROUTINE p4est_ResetElementNumber

    SUBROUTINE p4est_finalize() BIND(C, NAME = 'p4est_finalize')
        IMPORT :: C_INT
        !TYPE(C_PTR) :: P4
        !TYPE(C_PTR) ::p4est_new
        !INTEGER(C_INT) ::p4est_finalize
    END SUBROUTINE p4est_finalize

    SUBROUTINE p4est_loadbalancing(P4, user_pointer)  BIND(C, NAME = 'p4est_loadbalancing')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4
        TYPE(C_PTR), VALUE :: user_pointer
    END SUBROUTINE p4est_loadbalancing


    SUBROUTINE p4est_loadbalancing_init(P4, user_pointer)  BIND(C, NAME = 'p4est_loadbalancing_init')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4
        TYPE(C_PTR), VALUE :: user_pointer
    END SUBROUTINE p4est_loadbalancing_init

    SUBROUTINE p4est_loadbalancing_go(P4, user_pointer)  BIND(C, NAME = 'p4est_loadbalancing_go')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4
        TYPE(C_PTR), VALUE :: user_pointer
    END SUBROUTINE p4est_loadbalancing_go

    SUBROUTINE p4est_connectivity_destroy(CONN) BIND(C, NAME = 'p4est_connectivity_destroy_f')
        IMPORT :: C_PTR
        TYPE(C_PTR) :: CONN
        !TYPE(C_PTR) ::p4est_new
    END SUBROUTINE p4est_connectivity_destroy

    subroutine free_data_memory(Memory) BIND(C, NAME = 'free_data_memory')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: Memory
    END subroutine free_data_memory

    subroutine free_balance_memory(Memory) BIND(C, NAME = 'free_balance_memory')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: Memory
    END subroutine free_balance_memory


    subroutine free_savemesh_memory(Memory) BIND(C, NAME = 'savef_data_destroy')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: Memory
    END subroutine free_savemesh_memory

    FUNCTION P4EST_CONN_BCAST(CONN, ROOT, COMM) BIND(C, NAME = 'p8est_conn_bcast')
        IMPORT :: C_PTR, C_INT
        TYPE(C_PTR), VALUE :: CONN
        INTEGER(C_INT), VALUE :: ROOT
        INTEGER(C_INT), VALUE :: COMM
        TYPE(C_PTR) :: P4EST_CONN_BCAST
    END FUNCTION P4EST_CONN_BCAST


    ! This function is called in RefineAndCourse
    SUBROUTINE GetData(P4EST, P4fortrandata) BIND(C, NAME = 'GetData')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4EST
        ! TYPE(C_PTR) ::GetData
        TYPE(C_PTR), VALUE :: P4fortrandata
    END SUBROUTINE GetData


    SUBROUTINE SetEtSandStE(P4EST, P4fortrandata) BIND(C, NAME = 'SetEtSandStE')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4EST
        ! TYPE(C_PTR) ::GetData
        TYPE(C_PTR), VALUE :: P4fortrandata
    END SUBROUTINE SetEtSandStE

    SUBROUTINE GetnNBProcs(P4EST, P4fortrandata) BIND(C, NAME = 'GetnNBProcs')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4EST
        TYPE(C_PTR), VALUE :: P4fortrandata
    END SUBROUTINE GetnNBProcs

    FUNCTION GetNElems(P4EST) BIND(C, NAME = 'GetNElems')
        IMPORT :: C_PTR, C_INT
        TYPE(C_PTR), VALUE :: P4EST
        INTEGER :: GetNElems
    END FUNCTION GetNElems

    SUBROUTINE RefineCoarse(P4EST, ElemToRC) BIND(C, NAME = 'RefineCoarse')
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4EST, ElemToRC
        ! TYPE(C_PTR) ::RefineCoarse
    END SUBROUTINE RefineCoarse

    FUNCTION SaveMeshP4(P4EST) BIND(C, NAME = 'save_mesh') !
        IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4EST
        TYPE(C_PTR) :: SaveMeshP4
    END FUNCTION SaveMeshP4

    subroutine SaveP4(P4EST, FNAME) BIND(C, NAME = 'save_p4est') !
        USE ISO_C_BINDING
        ! IMPORT :: C_PTR
        TYPE(C_PTR), VALUE :: P4EST
        CHARACTER(KIND = C_CHAR), DIMENSION(*) :: FNAME

        ! TYPE(C_PTR) ::SaveMeshP4
    END subroutine SaveP4

    FUNCTION LoadP4(COMM, FNAME) BIND(C, NAME = 'load_p4est') !
        USE ISO_C_BINDING
        INTEGER(C_INT), VALUE :: COMM
        CHARACTER(KIND = C_CHAR), DIMENSION(*) :: FNAME
        ! TYPE(C_PTR),VALUE ::
        TYPE(C_PTR) :: LoadP4
    END FUNCTION LoadP4


    FUNCTION GetConnectivity(P4EST) BIND(C, NAME = 'GetConnectivity')
        IMPORT :: C_PTR
        TYPE(C_PTR) :: P4EST
        TYPE(C_PTR) :: GetConnectivity
    END FUNCTION GetConnectivity
    !This is from Hopest!
    SUBROUTINE p4_connectivity_treevertex(num_vertices, num_trees, vertices, tree_to_vertex, &
            num_periodics, JoinFaces, connectivity) BIND(C, NAME = "p4_connectivity_treevertex")
        !=================================================================================================================================
        ! builds up p4est connectivit, using only element connectivity and vertex positions
        !=================================================================================================================================
        ! MODULES
        USE, INTRINSIC :: ISO_C_BINDING
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !---------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        INTEGER(KIND = C_INT32_T), VALUE :: num_vertices
        INTEGER(KIND = C_INT32_T), VALUE :: num_trees
        REAL(KIND = C_DOUBLE) :: Vertices(3, num_vertices)
        INTEGER(KIND = C_INT32_T) :: tree_to_vertex(8 * num_trees)
        INTEGER(KIND = C_INT32_T), VALUE :: num_periodics
        INTEGER(KIND = C_INT32_T) :: JoinFaces(5 * num_periodics)
        !---------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        TYPE(C_PTR) :: connectivity
        !=================================================================================================================================
    END SUBROUTINE p4_connectivity_treevertex

    SUBROUTINE p4_build_bcs(conn, num_trees, bcelemmap) BIND(C, NAME = "p4_build_bcs")
        !=================================================================================================================================
        ! function to store BCs in Connectivity of p4est.
        !=================================================================================================================================
        ! MODULES
        USE, INTRINSIC :: ISO_C_BINDING
        ! IMPLICIT VARIABLE HANDLING
        IMPLICIT NONE
        !---------------------------------------------------------------------------------------------------------------------------------
        ! INPUT VARIABLES
        TYPE(C_PTR), VALUE, INTENT(IN) :: conn !Connectivity
        INTEGER(KIND = C_INT), VALUE :: num_trees
        INTEGER(KIND = C_INT32_T), INTENT(IN) :: bcelemmap(0:5, num_trees)
        !---------------------------------------------------------------------------------------------------------------------------------
        ! OUTPUT VARIABLES
        !=================================================================================================================================
    END SUBROUTINE p4_build_bcs

END INTERFACE


!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================


CONTAINS

    !==================================================================================================================================
    SUBROUTINE InitAMR_P4est()
        USE MOD_AMR_vars, ONLY : P4EST_PTR, CONNECTIVITY_PTR

        ! CONNECTIVITY_PTR = P8EST_CONNECTIVITY_NEW_PERIODIC();
        P4EST_PTR = P4EST_NEW(CONNECTIVITY_PTR);

    END SUBROUTINE InitAMR_P4est


    SUBROUTINE SaveP4est(FileString)
        USE MOD_AMR_vars, ONLY : P4EST_PTR, CONNECTIVITY_PTR
        USE, INTRINSIC :: ISO_C_BINDING
        IMPLICIT NONE
        !----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT/OUTPUT VARIABLES
        CHARACTER(LEN = *), INTENT(IN) :: FileString !< (IN) mesh filename
        !----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        CHARACTER(LEN = 10, KIND = C_CHAR) :: DIGIT_STRING = '123456789' // C_NULL_CHAR
        ! CONNECTIVITY_PTR = P8EST_CONNECTIVITY_NEW_PERIODIC();
        DIGIT_STRING = FileString // C_NULL_CHAR
        CALL SaveP4(P4EST_PTR, DIGIT_STRING);
    END SUBROUTINE SaveP4est


    SUBROUTINE LoadP4est(FileString)
        USE MOD_AMR_vars, ONLY : P4EST_PTR, CONNECTIVITY_PTR
        USE, INTRINSIC :: ISO_C_BINDING
        IMPLICIT NONE
        !----------------------------------------------------------------------------------------------------------------------------------
        ! INPUT/OUTPUT VARIABLES
        CHARACTER(LEN = *), INTENT(IN) :: FileString !< (IN) mesh filename
        !----------------------------------------------------------------------------------------------------------------------------------
        ! LOCAL VARIABLES
        CHARACTER(LEN = 10, KIND = C_CHAR) :: DIGIT_STRING = '123456789' // C_NULL_CHAR
        ! CONNECTIVITY_PTR = P8EST_CONNECTIVITY_NEW_PERIODIC();
        DIGIT_STRING = FileString // C_NULL_CHAR
        P4EST_PTR = LoadP4(MPI_COMM_WORLD, DIGIT_STRING);

        CONNECTIVITY_PTR = GetConnectivity(P4EST_PTR)

        !write(*,"(Z32)")  CONNECTIVITY_PTR
        CALL EXIT()
    END SUBROUTINE LoadP4est


    SUBROUTINE p4estSetMPIData()
        USE MOD_AMR_vars, ONLY : P4EST_PTR, p4est_mpi_data
        ! USE MOD_MPI_Vars,            ONLY:
        USE MOD_Mesh_Vars, ONLY : offsetElem, nElems, nGlobalElems
        USE MOD_MPI_Vars, ONLY : offsetElemMPI
        USE, INTRINSIC :: ISO_C_BINDING
        IMPLICIT NONE
        TYPE(p4est_mpi_data), POINTER :: DATAF
        INTEGER(KIND = C_INT64_T), POINTER :: Offset_f(:)
        TYPE(C_PTR) :: DataPtr;
        DataPtr = p4estGetMPIData(p4est_ptr)

        CALL C_F_POINTER(DataPtr, DATAF)
        nElems = DATAF%local_num_quad
        nGlobalElems = DATAF%global_num_quad

        ! DATAF%mpisize
        ! DATAF%mpirank
        ! DATAF%local_num_quad
        ! DATAF%global_num_quad
        ! TYPE(C_PTR) ::  offsetMPI;
        IF (DataF%mpisize.EQ.1) THEN
            nElems = nGlobalElems   !local number of Elements
            offsetElem = 0
        ELSE
            CALL C_F_POINTER(DataF%offsetMPI, Offset_f, [DataF%mpisize + 1])
            SDEALLOCATE(offsetElemMPI)
            ALLOCATE(offsetElemMPI(0:DATAF%mpisize))
            offsetElemMPI(0:DataF%mpisize) = offset_f(1:DataF%mpisize + 1)
            offsetElemMPI(DataF%mpisize) = nGlobalElems
            offsetElem = offsetElemMPI(DATAF%mpirank)
        ENDIF
        NULLIFY(Offset_f)
        NULLIFY(DATAF)
        CALL p4estDelMPIData(DataPtr)
        ! NULLIFY(DataPtr)
    END SUBROUTINE p4estSetMPIData

END MODULE MOD_P4EST
#endif
