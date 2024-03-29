!===================================================================================================================================
! Copyright (c) 2018 - 2020 Alexander Astanin
! Copyright (c) 2020 - 2021 Florian Hindenlang
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
#include "amr_f.h"
MODULE MOD_AMR_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
use MOD_Indicators, only: Indicator_PerssonPeraire
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

! !-----------------------------------------------------------------------------------------------------------------------------------
! ! P4EST related data structures 
! !-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)          :: p4estFile          ! name of hdf5 meshfile (write with ending .h5!)
LOGICAL                     :: UseAMR
LOGICAL                     :: AMRInitIsDone
LOGICAL                     :: p4estFileExist
INTEGER                     :: MaxLevel             ! Loaded from .ini file
INTEGER                     :: MinLevel             ! Loaded from .ini file
INTEGER                     :: nWriteDataAMR             ! Loaded from .ini file
INTEGER                     :: nDoAMR               ! Time-step interval to do AMR
INTEGER                     :: nDoAMRShift          ! Initial shift fot time-step interval to do AMR
integer                     :: InitialRefinement    ! Initial refinement switch
integer                     :: N_2                  ! Half of polynomial degree (interpolation operators)
real                        :: IniHalfwidthAMR
REAL                        :: RefineVal            ! Loaded from .ini file
REAL                        :: CoarseVal            ! Loaded from .ini file 
real, allocatable           :: Vdm_Interp_1_0(:,:)    ! Interpolation Vandermonde matrix from [-1,1] to first half of [-1,3]
real, allocatable           :: Vdm_Interp_2_0(:,:)    ! Interpolation Vandermonde matrix from [-1,1] to second half of [-3,1]
real, allocatable           :: M_0_1(:,:),M_0_2(:,:)  ! 1D Interpolation from big to small [-1,1] to [-1,0] / [0,1] 
real, allocatable           :: M_1_0(:,:),M_2_0(:,:)  ! 1D L2 Projection matrix for small to big [-1,0]/[0,1] to [-1,1]  
TYPE(C_PTR)                 :: P4EST_PTR              ! c pointers to p4est structures
TYPE(C_PTR)                 :: connectivity_ptr       ! c pointer to p4est connectivity 
type(Indicator_PerssonPeraire) :: AMR_Indicator

TYPE p4est_save_data
    INTEGER ::  nGlobalSides
    INTEGER ::  nLocalSides
    INTEGER ::  nSidesArrIndex;
    TYPE(C_PTR) ::  OffsetSideMPI;
    TYPE(C_PTR) ::  OffsetSideArrIndexMPI;
    TYPE(C_PTR) ::  ElemInfo;
    TYPE(C_PTR) ::  SideInfo;
END TYPE p4est_save_data

! Data strcuture to pass data from p4est to FLUXO
TYPE P4EST_FORTRAN_DATA
    INTEGER     ::  nGlobalElems;
    INTEGER     ::  nElems
    INTEGER     ::  nSides
    INTEGER     ::  nBCSides
    INTEGER     ::  nMortarInnerSides
    INTEGER     ::  nInnerSides
    INTEGER     ::  nMPISides
    INTEGER     ::  nMortarMPISides
    TYPE(C_PTR) ::  EtSPtr; !
    TYPE(C_PTR) ::  StEPtr; !
    TYPE(C_PTR) ::  MTPtr; !
    TYPE(C_PTR) ::  MIPtr; !
    TYPE(C_PTR) ::  ChngElmPtr; !
    TYPE(C_PTR) ::  ChngSidePtr; !NOT
    TYPE(C_PTR) ::  nNbProc; !
    TYPE(C_PTR) ::  nMPISides_Proc; !
    INTEGER     ::   nNBProcs;
    TYPE(C_PTR) ::  nMPISides_MINE_Proc; !
    TYPE(C_PTR) ::  nMPISides_YOUR_Proc; !
    INTEGER     ::   nMPISides_YOUR;
    INTEGER     ::   nMPISides_MINE; 
    TYPE(C_PTR) ::  offsetMPISides_MINE; !
    TYPE(C_PTR) ::  offsetMPISides_YOUR; ! 
    TYPE(C_PTR) ::  BCs;
    TYPE(C_PTR) ::  GhostPtr; ! Used in C only
    TYPE(C_PTR) ::  GhostDataPtr; !Used in C only
    TYPE(C_PTR) ::  ghost_to_proc_Ptr; !Used in C only
END TYPE P4EST_FORTRAN_DATA

TYPE p4est_balance_data
    INTEGER     ::  nVar
    INTEGER     ::  PP
    INTEGER     ::  nElems;
    INTEGER     ::  DataSize
    TYPE(C_PTR) ::  DataSetU;
    TYPE(C_PTR) ::  DataSetElem_xGP;
END TYPE p4est_balance_data


TYPE p4est_balance_datav2
    INTEGER     ::  DataSize !DataSize of 1 Elem U = number of REALs
    INTEGER     ::  GPSize !DataSize of 1 ElemxGP = number of REALs
    INTEGER     ::  CoarseSize !DataSize of 1 ElemxGP = number of REALs
    INTEGER     ::  nElems ! Number of new elements
    TYPE(C_PTR) ::  GlbIdxData; ! For transfer data. for p4est only
    TYPE(C_PTR) ::  ElemxGPnew_Ptr; ! new Array ElemxGP
    TYPE(C_PTR) ::  ElemxGPold_Ptr; ! old Array ElemxGP
    TYPE(C_PTR) ::  Unew_Ptr; !new U
    TYPE(C_PTR) ::  Uold_Ptr; !old U
    TYPE(C_PTR) ::  ElemWasCoarsened_new; 
    TYPE(C_PTR) ::  ElemWasCoarsened_old; 
#if SHOCK_NFVSE
    TYPE(C_PTR) ::  AlphaNew_Ptr; !new alpha
    TYPE(C_PTR) ::  AlphaOld_Ptr; !old alpha
#endif /*SHOCK_NFVSE*/
END TYPE p4est_balance_datav2

TYPE p4est_mpi_data
    INTEGER ::  mpisize
    INTEGER ::  mpirank
    P4EST_F90_GLOIDX ::  local_num_quad
    P4EST_F90_GLOIDX ::  global_num_quad
    TYPE(C_PTR) ::  offsetMPI;
END TYPE p4est_mpi_data

INTEGER,PARAMETER   :: H2P_FaceMap(1:6)     =  (/4,2,1,3,0,5/)     !mapping from local face order (CGNS) to p4est face
INTEGER,PARAMETER   :: P2H_FaceMap(0:5)     =  (/5,3,2,4,1,6/)     !mapping from p4est face to local face order (CGNS) 
INTEGER,PARAMETER   :: H2P_VertexMap(1:8)   =  (/0,1,3,2,4,5,7,6/) !mapping from local node order (CGNS) to p4est node order 

! ! mapping from CGNS node of local sides to P4EST nodes of local sides
INTEGER,PARAMETER   :: H2P_FaceNodeMap(1:4,1:6) = &
                                      RESHAPE((/ 0,2,3,1,&
                                                 0,1,3,2,&
                                                 0,1,3,2,&
                                                 1,0,2,3,&
                                                 0,2,3,1,&
                                                 0,1,3,2 /),(/4,6/))

! ! mapping from P4EST node of local sides to CGNS node of local sides
INTEGER,PARAMETER   :: P2H_FaceNodeMap(0:3,0:5) = &
                                      RESHAPE((/ 1,4,2,3,&
                                                 1,2,4,3,&
                                                 1,2,4,3,&
                                                 2,1,3,4,&
                                                 1,4,2,3,&
                                                 1,2,4,3 /),(/4,6/))

! ! Mapping matrices for computation of same node on adjacent face, see paper Burstedde p4est, 2011
! ! Node1= P4P(P4Q(P4R(Face0,Face1),orientation),Node0)
INTEGER,PARAMETER   :: P4R(0:5,0:5) = TRANSPOSE(RESHAPE((/ 0,1,1,0,0,1,&
                                                           2,0,0,1,1,0,&
                                                           2,0,0,1,1,0,&
                                                           0,2,2,0,0,1,&
                                                           0,2,2,0,0,1,&
                                                           2,0,0,2,2,0 /),(/6,6/)))

INTEGER,PARAMETER   :: P4Q(0:2,0:3) = TRANSPOSE(RESHAPE((/ 1,2,5,6,&
                                                           0,3,4,7,&
                                                           0,4,3,7 /),(/4,3/)))

INTEGER,PARAMETER   :: P4P(0:7,0:3) = TRANSPOSE(RESHAPE((/ 0,1,2,3,&
                                                           0,2,1,3,&
                                                           1,0,3,2,&
                                                           1,3,0,2,&
                                                           2,0,3,1,&
                                                           2,3,0,1,&
                                                           3,1,2,0,&
                                                           3,2,1,0 /),(/4,8/)))
! !===================================================================================================================================
    TYPE(P4EST_FORTRAN_DATA), TARGET    :: FortranData

END MODULE MOD_AMR_Vars
