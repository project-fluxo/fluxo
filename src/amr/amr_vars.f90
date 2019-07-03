#include "amr_f.h"
MODULE MOD_AMR_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
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

TYPE(C_PTR)                 :: P4EST_PTR              ! c pointers to p4est structures
TYPE(C_PTR)                 :: connectivity_ptr       !

! typedef struct p4est_savef_data
! {   
!     int nGlobalSides;
!     int nLocalSides;
!     int nSidesArrIndex;
!     int *OffsetSideMPI;
!     int *OffsetSideArrIndexMPI;
!     int *ElemInfo; // ElemInfo in HDF5
!     int *SideInfo; // SideInfo in HDF5

! } p4est_savef_data_t;

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
    INTEGER ::  nGlobalElems;
    INTEGER ::  nElems
    INTEGER ::  nSides
    INTEGER ::  nBCSides
    INTEGER ::  nMortarInnerSides
    INTEGER ::  nInnerSides
    INTEGER ::  nMPISides
    INTEGER ::  nMortarMPISides
    TYPE(C_PTR) ::  EtSPtr; !
    TYPE(C_PTR) ::  StEPtr; !
    TYPE(C_PTR) ::  MTPtr; !
    TYPE(C_PTR) ::  MIPtr; !
    TYPE(C_PTR) ::  ChngElmPtr; !
    TYPE(C_PTR) ::  ChngSidePtr; !NOT
    TYPE(C_PTR) ::  nNbProc; !
    TYPE(C_PTR) ::  nMPISides_Proc; !
    INTEGER ::   nNBProcs;
    TYPE(C_PTR) ::  nMPISides_MINE_Proc; !
    TYPE(C_PTR) ::  nMPISides_YOUR_Proc; !
    INTEGER ::   nMPISides_YOUR;
    INTEGER ::   nMPISides_MINE; 
    TYPE(C_PTR) ::  offsetMPISides_MINE; !
    TYPE(C_PTR) ::  offsetMPISides_YOUR; ! 
    TYPE(C_PTR) ::  BCs;
END TYPE P4EST_FORTRAN_DATA

TYPE p4est_balance_data
    INTEGER ::  nVar
    INTEGER ::  PP_N
    INTEGER ::  nElems;
    INTEGER ::  DataSize
    TYPE(C_PTR) ::  DataSetU;
    TYPE(C_PTR) ::  DataSetElem_xGP;
END TYPE p4est_balance_data

TYPE p4est_mpi_data
    INTEGER ::  mpisize
    INTEGER ::  mpirank
    P4EST_F90_GLOIDX ::  local_num_quad
    P4EST_F90_GLOIDX ::  global_num_quad
    TYPE(C_PTR) ::  offsetMPI;
END TYPE p4est_mpi_data
! !-----------------------------------------------------------------------------------------------------------------------------------
! INTEGER,PARAMETER   :: EdgeToElemNode(1:2,1:12) = RESHAPE((/ 1, 2,&  ! CGNS corner nodes mapped 
!                                                              4, 3,&  ! to p4est edges
!                                                              5, 6,&
!                                                              8, 7,&
!                                                              1, 4,&
!                                                              2, 3,&
!                                                              5, 8,&
!                                                              6, 7,&
!                                                              1, 5,&
!                                                              2, 6,&
!                                                              4, 8,&
!                                                              3, 7 /),(/2,12/))
INTEGER,PARAMETER   :: H2P_FaceMap(1:6)     =  (/4,2,1,3,0,5/)     !mapping from local face order (CGNS) to p4est face
INTEGER,PARAMETER   :: P2H_FaceMap(0:5)     =  (/5,3,2,4,1,6/)     !mapping from p4est face to local face order (CGNS) 
INTEGER,PARAMETER   :: H2P_VertexMap(1:8)   =  (/0,1,3,2,4,5,7,6/) !mapping from local node order (CGNS) to p4est node order 
! INTEGER,PARAMETER   :: P2H_VertexMap(0:7)   =  (/1,2,4,3,5,6,8,7/) !mapping from local node order (CGNS) to p4est node order 

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
! INTEGER,PARAMETER   :: H_MortarCase(1:4,1:4) = &                              !  first CGNS node and second CGNS node->Mortar Case [1:8]
!                                       TRANSPOSE(RESHAPE((/  0,1,0,2,&                         ! (1,2)->1, (1,4)->2
!                                                             3,0,4,0,&                         ! (2,1)->3, (2,3)->4
!                                                             0,5,0,6,&                         ! (3,2)->5, (3,4)->6
!                                                             7,0,8,0 /),(/4,4/)))              ! (4,1)->7, (4,3)->8
! INTEGER,PARAMETER   :: P2H_MortarMap(0:3,1:8) = &                 !p4est mortar ID, MortarCase -> iMortar CGNS
!                                       RESHAPE((/ 1,2,3,4,&        ! iMortar = P2H_MortarMap(iPMortar, H_MortarCase( node1, node2) ) 
!                                                  1,3,2,4,&
!                                                  2,1,4,3,&
!                                                  2,4,1,3,&
!                                                  4,2,3,1,&
!                                                  4,3,2,1,&
!                                                  3,1,4,2,&
!                                                  3,4,1,2 /),(/4,8/))
! INTEGER,PARAMETER   :: P_FaceToEdge(0:3,0:5) = &  !mapping from face edges 0...3 (zordered) for each face 0..5 -> element edges  0..11
!                                       RESHAPE((/  4, 6, 8,10,&     
!                                                   5, 7, 9,11,&
!                                                   0, 2, 8, 9,&
!                                                   1, 3,10,11,&
!                                                   0, 1, 4, 5,&
!                                                   2, 3, 6, 7 /),(/4,6/))
! INTEGER,PARAMETER   :: P_EdgeToFaces(1:6,0:11) = & !mapping from element edges  0..11 -> first and second adjacent face 0...6 , i/j direction(0/1) and lower/upper bound(0/1)
!                                       RESHAPE((/  2,1,0,4,1,0,&    !edge 0: Face2,j=0-Face4,j=0
!                                                   3,1,0,4,1,1,&    !edge 1: Face3,j=0-Face4,j=N
!                                                   2,1,1,5,1,0,&    !edge 2: Face2,j=N-Face5,j=0
!                                                   3,1,1,5,1,1,&    !edge 3: Face3,j=N-Face5,j=N
!                                                   0,1,0,4,0,0,&    !edge 4: Face0,j=0-Face4,i=0
!                                                   1,1,0,4,0,1,&    !edge 5: Face1,j=0-Face4,i=N
!                                                   0,1,1,5,0,0,&    !edge 6: Face0,j=N-Face5,i=0
!                                                   1,1,1,5,0,1,&    !edge 7: Face1,j=N-Face5,i=N
!                                                   0,0,0,2,0,0,&    !edge 8: Face0,i=0-Face2,i=0
!                                                   1,0,0,2,0,1,&    !edge 9: Face1,i=0-Face2,i=N
!                                                   0,0,1,3,0,0,&    !edge10: Face0,i=N-Face3,i=0
!                                                   1,0,1,3,0,1 &    !edge11: Face1,i=N-Face3,i=N
!                                                   /),(/6,12/))

! !===================================================================================================================================


END MODULE MOD_AMR_Vars
