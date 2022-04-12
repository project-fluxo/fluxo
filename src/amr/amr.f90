!==================================================================================================================================
! Copyright (c) 2018 - 2020 Alexander Astanin
! Copyright (c) 2020 - 2020 Andrés Rueda
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
!==================================================================================================================================
#include "defines.h"
#if USE_AMR
!==================================================================================================================================
!> Contains control routines for AMR
!==================================================================================================================================
MODULE MOD_AMR
! MODULES
#if MPI
USE mpi
#endif
IMPLICIT NONE
! PUBLIC
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitAMR
  MODULE PROCEDURE InitAMR
END INTERFACE

INTERFACE RunAMR
  MODULE PROCEDURE RunAMR
END INTERFACE


INTERFACE InitAMR_Connectivity
  MODULE PROCEDURE InitAMR_Connectivity
END INTERFACE

INTERFACE WriteStateAMR
  MODULE PROCEDURE WriteStateAMR
END INTERFACE

INTERFACE FinalizeAMR
  MODULE PROCEDURE FinalizeAMR
END INTERFACE

INTERFACE SaveMesh
  MODULE PROCEDURE SaveMesh
END INTERFACE

PUBLIC::InitAMR, SaveMesh
PUBLIC::FinalizeAMR
PUBLIC :: WriteStateAMR
!==================================================================================================================================
PUBLIC::RunAMR
PUBLIC::InitAMR_Connectivity
PUBLIC::DefineParametersAMR
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAMR()
  ! MODULES
  USE MOD_Globals
  USE MOD_ReadInTools   , ONLY: prms
  use MOD_Equation_Vars , only: IndicatorQuantityNames, nIndVar
  IMPLICIT NONE
  !-local-variables----------------------------
  integer            :: i
  character(len=255) :: IndicatorQuantities,fmt
  !==================================================================================================================================
  CALL prms%SetSection("AMR")
  CALL prms%CreateLogicalOption( 'UseAMR',              " Use AMR?",&
                                                      '.FALSE.')
  CALL prms%CreateStringOption('p4estFile',             " Path to p4est's connectivity file (mandatory).")
  
  write(fmt,'(A,I0,A)') '(',nIndVar,'A)'
  write(IndicatorQuantities,fmt) ('  * '//trim(IndicatorQuantityNames(i))//'\n', i=1, nIndVar)

  CALL prms%CreateStringOption("AMRIndicatorQuantity",  " Specifies the quantity to be used for the AMR indicator. One of the following:\n"//&
                                                          trim(IndicatorQuantities)&
                                             ,                "DensityTimesPressure")
 CALL prms%CreateRealOption(  'RefineVal',              " The value of the AMR indicator above which an element has to be subdivided (refinement)", "-1.0")
 CALL prms%CreateRealOption(  'CoarseVal',              " The value of the AMR indicator below which elements can be merged (coarsening)", "-2.0")
 
 CALL prms%CreateIntOption(  'MinLevel',                " Minimum refinment level of the mesh with respect to base mesh", "0")
 CALL prms%CreateIntOption(  'MaxLevel',                " Maximum refinemen level of the mesh with respect to base mesh", "0")
 CALL prms%CreateIntOption(  'nWriteDataAMR',           " Interval as multiple of nWriteData at which Mesh and p4est files"//&
                                                        " (_mesh.h5 and .p4est) are written.",&
                                                              '1')
 CALL prms%CreateIntOption(  'nDoAMR'           ,       " Time-step interval to call the AMR routines.",&
                                                  '1')
 CALL prms%CreateIntOption(  'nDoAMRShift'      ,       " Initial shift for time-step interval to call the AMR routines.",&
                                                  '0')
 CALL prms%CreateIntOption(  'InitialRefinement',       " Initial refinement to be used\n"//&
                                                        "  0: Use the custom indicator\n"//&
                                                        "  1: Refine the elements in the sphere with r=IniHalfwidthAMR\n"//&
                                                        "  2: Do nothing (needed to restart AMR simulation properly)",&
                                                  '0')
 CALL prms%CreateRealOption(   'IniHalfwidthAMR',       " Parameter for InitialRefinement.","0.1")
END SUBROUTINE DefineParametersAMR




  

!==================================================================================================================================
!> Routine controlling the initialization of the AMR.
!==================================================================================================================================
SUBROUTINE InitAMR()
    ! MODULES
    USE MOD_Globals
    USE MOD_PreProc
    USE MOD_AMR_Vars
    USE MOD_P4EST
    USE MOD_ReadInTools         ,ONLY:GETLOGICAL,GETSTR, GETINT, GETREAL
    use MOD_Interpolation_Vars  ,only: NodeType
    use MOD_Basis               ,only: InitializeVandermonde
    use MOD_Interpolation       ,only: GetNodesAndWeights
    use MOD_Mortar_vars         ,only: MortarBasis_BigToSmall,MortarBasis_SmallToBig_projection
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    real :: xi_In(0:PP_N),w_in(0:PP_N),wBary_In(0:PP_N)
    character(len=255) :: AMRIndicatorQuantity
    !==================================================================================================================================
    INTEGER RET

    IF(AMRInitIsDone) THEN
      CALL CollectiveStop(__STAMP__,&
        'InitAMR not ready to be called or already called.')
    END IF

    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT AMR...'

    UseAMR = GETLOGICAL('UseAMR','.FALSE.')
    IF (UseAMR) THEN
      p4estFile = GETSTR('p4estFile')
    ELSE
      SWRITE(UNIT_stdOut,'(A)') ' AMR is not used'
      RETURN;
    ENDIF
    p4estFileExist = CheckP4estFileExist(trim(p4estFile))
    MinLevel = GetINT('MinLevel',"0")
    MaxLevel = GetINT('MaxLevel',"0")
    RefineVal = GetReal('RefineVal',"0.")
    CoarseVal = GetREal('CoarseVal',"0.")
    nWriteDataAMR = GetINT('nWriteDataAMR',"1")
    nDoAMR = GetINT('nDoAMR',"1")
    nDoAMRShift = GetINT('nDoAMRShift',"0")
    InitialRefinement = GETINT('InitialRefinement','0')
    IniHalfwidthAMR   = GetREal('IniHalfwidthAMR',"0.1")
    AMRIndicatorQuantity = GETSTR('AMRIndicatorQuantity','DensityTimesPressure')
    
    ! Do some checks
    if (MinLevel > MaxLevel) then
      CALL CollectiveStop(__STAMP__,&
          'AMR: MinLevel must be less than or equal to MaxLevel.')
    end if
    if (MinLevel < 0) then
      CALL CollectiveStop(__STAMP__,&
          'AMR: MinLevel must be >=0.')
    end if
    if (RefineVal < CoarseVal) then
      CALL CollectiveStop(__STAMP__,&
          'AMR: RefineVal must be greater than or equal to MaxLevel.')
    end if
    
    ! Construct AMR indicator
    call AMR_Indicator % construct(AMRIndicatorQuantity,.TRUE.)
    
#if MPI 
    RET=P4EST_INIT(MPI_COMM_WORLD); 
#else  /*MPI*/
    RET=P4EST_INIT(INT(Z'44000000',KIND=4));      
#endif  /*MPI*/
    
    IF  (p4estFileExist) THEN
      
      IF(MPIRoot)THEN
        WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')'THE P4EST FILE HAS BEEN FOUND: '
        WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')TRIM(p4estFile)
      END IF
      CALL LoadP4est(trim(p4estFile))
    ELSE
      IF(MPIRoot)THEN
        WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')'THE P4EST FILE HAS NOT BEEN FOUND : '
        WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')TRIM(p4estFile)
        WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')'CREATE THE FOREST'
      END IF
      CALL InitAMR_Connectivity()
      CALL InitAMR_P4est()
    ENDIF
    
    ! Refinement and coarsening interpolation operators
    ! -------------------------------------------------
    CALL GetNodesAndWeights(PP_N,NodeType,xi_In,wIP=w_in,wIPBary=wBary_In)
    
    ! Refinement operators, simply interpolation  (use 1D mortar routine)
    allocate (M_0_1(0:PP_N,0:PP_N), M_0_2(0:PP_N,0:PP_N))
    CALL MortarBasis_BigToSmall(PP_N,NodeType,   M_0_1,M_0_2)

    
    ! Coarsening operators, interpolation for the mapping
    N_2 = PP_N/2
    allocate( Vdm_Interp_1_0(0:N_2,0:PP_N),Vdm_Interp_2_0(N_2+1:PP_N,0:PP_N) ) ! Temporary
    
    CALL InitializeVandermonde(PP_N,N_2              ,wBary_In,xi_In,2.0*xi_In(0    :N_2 )+1.,Vdm_Interp_1_0)
    CALL InitializeVandermonde(PP_N,N_2-mod(PP_N+1,2),wBary_In,xi_In,2.0*xi_In(N_2+1:PP_N)-1.,Vdm_Interp_2_0)

    ! USE MORTAR L2 PROJECTION 1D MATRIX for projection when coarsening of solution
    ALLOCATE(M_1_0(0:PP_N,0:PP_N), M_2_0(0:PP_N,0:PP_N))
    CALL MortarBasis_SmallToBig_Projection(PP_N,NodeType,   M_1_0,M_2_0)
        
    AMRInitIsDone=.TRUE.
    SWRITE(UNIT_stdOut,'(A)')' INIT AMR DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')
    
  
END SUBROUTINE InitAMR

!==================================================================================================================================
!> Routine creating and 
!==================================================================================================================================
SUBROUTINE InitAMR_Connectivity()
   ! MODULES
   USE MOD_Globals
   USE MOD_Mesh_Vars,             ONLY: MeshFile
   USE MOD_P4EST
   USE MOD_AMR_Vars,              ONLY: connectivity_ptr
   USE MODH_Mesh,                 ONLY: InitMesh, FinalizeMesh
   USE MODH_Mesh_ReadIn,        ONLY: ReadMeshHeader,ReadMeshFromHDF5
   USE, INTRINSIC :: ISO_C_BINDING
   IMPLICIT NONE
   INTEGER CONN_OWNER
   !----------------------------------------------------------------------------------------------------------------------------------
   ! INPUT/OUTPUT VARIABLES
   !----------------------------------------------------------------------------------------------------------------------------------
   ! LOCAL VARIABLES
   !==================================================================================================================================
   
   CONN_OWNER=0;
    CALL InitMesh()
   ! Create Connectivity
   CALL ReadMeshHeader(MeshFile)   ! read mesh header file including BCs
   CALL ReadMeshFromHDF5(MeshFile)
  
   ! The result is connectivity PTR
#if MPI 
    connectivity_ptr=P4EST_CONN_BCAST(connectivity_ptr, CONN_OWNER, MPI_COMM_WORLD)
#else  /*MPI*/
    connectivity_ptr=P4EST_CONN_BCAST(connectivity_ptr, CONN_OWNER,  INT(Z'44000000',KIND=4))
#endif  /*MPI*/
  
   CALL FinalizeMesh()
  
END SUBROUTINE InitAMR_Connectivity

!==================================================================================================================================
!> write mesh and p4est files to the 
!==================================================================================================================================
SUBROUTINE WriteStateAMR(OutputTime,isErrorFile)
  ! MODULES
  USE MOD_PreProc
  USE MOD_Globals
  USE MOD_Output_Vars  ,ONLY: ProjectName
  USE MOD_P4EST         ,ONLY: SaveP4est
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  REAL,INTENT(IN)                :: OutputTime   !< simulation time when output is performed
  LOGICAL,INTENT(IN)             :: isErrorFile  !< indicate whether an error file is written in case of a crashed simulation
  !----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  !----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  CHARACTER(LEN=255)             :: FileType
  CHARACTER(LEN=255)             :: FileNameAMRMesh ! Mesh File
  CHARACTER(LEN=255)             :: FileNameP !P4est File
  
  !==================================================================================================================================
  IF(isErrorFile) THEN
    FileType='ERROR_State'
  ELSE
    FileType='State'
  END IF
  FileNameAMRMesh=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileType),OutputTime))//'_mesh.h5'
  FileNameP=  TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(FileType),OutputTime))//'.p4est'
  CALL SaveMesh(FileNameAMRMesh)
  CALL SaveP4est(FileNameP)
  END SUBROUTINE WriteStateAMR
  
LOGICAL FUNCTION CheckP4estFileExist(FileString)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString !< (IN) mesh filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
inquire( file=trim(FileString), exist=CheckP4estFileExist )
END FUNCTION CheckP4estFileExist

!==================================================================================================================================
!>  The main SUBROUTINE used to adapt the mesh (Coarsening and refining).
!==================================================================================================================================
SUBROUTINE RunAMR(ElemToRefineAndCoarse)
  USE MOD_Globals
  USE MOD_PreProc,            ONLY: PP_N
  USE MOD_Analyze_Vars,       ONLY: ElemVol
  USE MOD_AMR_Vars,           ONLY: P4EST_FORTRAN_DATA, P4est_ptr, UseAMR, FortranData
  USE MOD_Mesh_Vars,          ONLY: Elem_xGP, ElemToSide, SideToElem, Face_xGP, NormVec, TangVec1, TangVec2
  USE MOD_Mesh_Vars,          ONLY: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde,dXGL_N, sJ, SurfElem
  USE MOD_P4est,              ONLY: free_data_memory, RefineCoarse, GetData, p4estSetMPIData, GetnNBProcs, SetEtSandStE
  USE MOD_P4est,              ONLY: FillElemsChanges, GetNElems
  USE MOD_Metrics,            ONLY: CalcMetrics
  USE MOD_DG_Vars,            ONLY: U,Ut,Source,nTotalU, nTotal_vol, nTotal_IP, nTotal_face, nDOFElem, U_master, U_SLAVE, Flux_master, Flux_slave
  USE MOD_Mesh_Vars,          ONLY: AnalyzeSide, MortarInfo, MortarType, NGeo, DetJac_Ref, BC
  USE MOD_TimeDisc_Vars,      ONLY:   dtElem
  USE MOD_Mesh_Vars,          ONLY: LastSlaveSide, firstSlaveSide, nSides, nElems, firstMortarInnerSide !, lastMortarInnerSide
#if MPI 
  USE MOD_MPI_Vars,           ONLY: NbProc , nMPISides_MINE_Proc, nMPISides_YOUR_Proc, offsetMPISides_YOUR, offsetMPISides_MINE 
  USE MOD_MPI_Vars,           ONLY: nMPISides_Proc, nMPISides_send, nMPISides_rec, OffsetMPISides_send, OffsetMPISides_rec
  USE MOD_MPI_Vars,           ONLY: MPIRequest_U, MPIRequest_Flux, nNbProcs
#endif  /*MPI*/
  USE MOD_Globals ,           ONLY: nProcessors, MPIroot, myrank
  use MOD_GetBoundaryFlux,    only: InitBC,FinalizeBC
  use MOD_Mortar  ,           only: FinalizeMortarArrays,InitMortarArrays
#if POSITIVITYPRES
  use MOD_PositivityPreservation, only: FinalizePositivityPreservation, InitPositivityPreservation
#endif /*POSITIVITYPRES*/
#if PARABOLIC
  USE  MOD_Lifting_Vars
#if MPI
  USE  MOD_MPI_Vars,          ONLY: MPIRequest_Lifting
#endif /*MPI*/
#endif /* PARABOLIC */
#if SHOCKCAPTURE
  use MOD_ShockCapturing,     only: InitShockCapturingAfterAdapt1, InitShockCapturingAfterAdapt2
#endif /*SHOCKCAPTURE*/
  use MOD_Equation,           only: InitEquationAfterAdapt
  USE, INTRINSIC :: ISO_C_BINDING
! ----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENTS
  INTEGER, ALLOCATABLE, TARGET :: ElemToRefineAndCoarse(:) ! positive Number - refine, negative - coarse, 0 - do nothing
! ----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL,ALLOCATABLE :: Elem_xGP_New(:,:,:,:,:), U_New(:,:,:,:,:)
  INTEGER :: Ie
  TYPE(C_PTR) :: DataPtr;
  INTEGER, POINTER :: MInfo(:,:,:), ChangeElem(:,:)
  INTEGER, POINTER :: nBCsF(:)
  INTEGER :: i,j,k,iElem, nMortarSides, NGeoRef
  INTEGER :: nElemsOld, nSidesOld, LastSlaveSideOld, firstSlaveSideOld, firstMortarInnerSideOld
  integer, allocatable :: ElemWasCoarsened(:)
  integer :: max_nElems, min_nElems, sum_nElems
  integer :: new_nElems
!==================================================================================================================================
  IF (.NOT. UseAMR) THEN
    RETURN;
  ENDIF

! ==============
! Store old data  
! ==============
  nElemsOld = nElems;
  nSidesOld = nSides
  LastSlaveSideOld = LastSlaveSide;
  firstSlaveSideOld = firstSlaveSide;
  firstMortarInnerSideOld = firstMortarInnerSide
  
! ================================================
! Compute which elements will be refined/coarsened
! ================================================
  CALL RefineCoarse(p4est_ptr,C_LOC(ElemToRefineAndCoarse))
  DATAPtr = C_LOC(FortranData)
  FortranData%nElems = GetNElems(p4est_ptr)

! =========================================================================================================================
! Get the element changes into the ChangeElem variable
!       ChangeElem(1:8,1:nElems)
!                    |        |_ New element index (in the adapted mesh)
!                    |__________ Stores old element index(es)
!                                * Refinement: First position = minus the old element index, second to eighth positions = 0
!                                * Coarsening: All poisitive old element indices
!                                * Nothing:    First position = the old element index (positive), second to eighth = 0
! =========================================================================================================================
  ALLOCATE(ChangeElem(8,FortranData%nElems))
  FortranData%ChngElmPtr = C_LOC(ChangeElem)
  
  CALL FillElemsChanges(p4est_ptr, DATAPtr)
  
! =====================================================================================================
! Transfer the solution (U) and the node coordinates (Elem_xGP) to the new mesh
! ATTENTION!!: 1) For the elements that are coarsened, we transfer J*U. This is needed to obtain an 
!                 accurate L2 projection in physical space!!
!              2) U is recovered in a later stage after the metric terms are computed for the new mesh.
! =====================================================================================================
  ALLOCATE(Elem_xGP_New(3      ,0:PP_N,0:PP_N,0:PP_N,FortranData%nElems))
  ALLOCATE(U_New       (PP_nVar,0:PP_N,0:PP_N,0:PP_N,FortranData%nElems))
  allocate(ElemWasCoarsened(FortranData%nElems) )
  ElemWasCoarsened = 0
  iElem=0;
  DO 
    iElem=iElem+1;
    IF (iElem .GT. FortranData%nElems) EXIT
    
    Ie= ChangeElem(1,iElem);    
    IF (Ie .LT. 0) THEN
      !This is refine and this and next 7 elements [iElem: iElem+7] number negative and 
      ! contains the number of child element 
      call InterpolateSolution_Refinement(3      ,Elem_xGP_New(:,:,:,:,iElem:iElem+7), Elem_xGP(:,:,:,:,-Ie))
      call InterpolateSolution_Refinement(PP_nVar,U_New       (:,:,:,:,iElem:iElem+7), U(:,:,:,:,-Ie))
      iElem=iElem+7;
    ELSE IF (ChangeElem(2,iElem) .GT. 0) THEN
      !  This is COARSE. Array ChangeElem(:,iElem) Contains 
      !  8 Element which must be COARSED to the new number iElem
      ElemWasCoarsened(iElem) = 1
      call InterpolateCoords_Coarsening( Elem_xGP_New(:,:,:,:,iElem),&
                                         Elem_xGP    (:,:,:,:,ChangeElem(:,iElem)) )
      call ProjectSolution_Coarsening(U_New(:,:,:,:,iElem), U(:,:,:,:,ChangeElem(:,iElem)),sJ(:,:,:,ChangeElem(:,iElem)))
    ELSE
      IF (iE .LE. 0) THEN
        print *, "Error, iE = 0!, iElem = ", ielem
        CALL EXIT()
      ENDIF
      ! This is simple case of renumeration of Elements
      Elem_xGP_New(:,:,:,:,iElem)= Elem_xGP(:,:,:,:,Ie)
      U_New       (:,:,:,:,iElem)= U       (:,:,:,:,Ie)
    ENDIF
  END DO
  CALL MOVE_ALLOC(Elem_xGP_New, Elem_xGP)
  CALL MOVE_ALLOC(U_New, U)
! =====================================================================================================
! Transfer the additional variables to the new mesh
! ATTENTION!!: All variables that are load-balanced must be transferred!
! =====================================================================================================
#if SHOCKCAPTURE
  call InitShockCapturingAfterAdapt1(FortranData%nElems,ChangeElem)
#endif /*SHOCKCAPTURE*/

  
! =============================================
! Perform load balancing and display statistics
! =============================================
  
  CALL LoadBalancingAMR(ElemWasCoarsened,new_nElems)

#if MPI 
  
          
  call MPI_Reduce(new_nElems, max_nElems, 1 , MPI_INT, MPI_MAX, 0, MPI_Comm_WORLD, i)
  call MPI_Reduce(new_nElems, min_nElems, 1 , MPI_INT, MPI_MIN, 0, MPI_Comm_WORLD, i)
  call MPI_Reduce(new_nElems, sum_nElems, 1 , MPI_INT, MPI_SUM, 0, MPI_Comm_WORLD, i)

#endif  /*MPI*/        
  IF (MPIRoot) THEN
    WRITE(*,'(A,I0,A,I0,A,F0.2,A,I0)') "LoadBalance: Done! nGlobalElems=", sum_nElems, ", min_nElems=", min_nElems, ", avg_nElems=", sum_nElems/real(nProcessors), ", max_nElems=", max_nElems
  ENDIF

! ==============================
! Reallocate MPI data structures
! ==============================
 
 CALL GetnNBProcs(p4est_ptr, DATAPtr)

#if MPI 
  IF (nProcessors .GT. 1) THEN
    nNbProcs=FortranData%nNBProcs
    IF (ALLOCATED(NbProc)) THEN 
      SDEALLOCATE(NbProc); ALLOCATE(NbProc(1:nNbProcs))
    ENDIF
    FortranData%nNbProc = C_LOC(NbProc)
    
    SDEALLOCATE(nMPISides_Proc)
    ALLOCATE(nMPISides_Proc(1:nNbProcs))
    FortranData%nMPISides_Proc = C_LOC(nMPISides_Proc)
    
    
    SDEALLOCATE(nMPISides_MINE_Proc)
    ALLOCATE(nMPISides_MINE_Proc(1:nNbProcs))
    FortranData%nMPISides_MINE_Proc = C_LOC(nMPISides_MINE_Proc)
    
    
    SDEALLOCATE(nMPISides_YOUR_Proc)
    ALLOCATE(nMPISides_YOUR_Proc(1:nNbProcs))
    FortranData%nMPISides_YOUR_Proc = C_LOC(nMPISides_YOUR_Proc)
    
    
    SDEALLOCATE(offsetMPISides_YOUR)
    ALLOCATE(offsetMPISides_YOUR(0:nNbProcs))
    FortranData%offsetMPISides_YOUR = C_LOC(offsetMPISides_YOUR)
    
    SDEALLOCATE(offsetMPISides_MINE)
    ALLOCATE(offsetMPISides_MINE(0:nNbProcs))
    FortranData%offsetMPISides_MINE = C_LOC(offsetMPISides_MINE)
    
  ENDIF
          
#endif  /*MPI*/
  
! ==========================================================================
! Count and organize sides (nSides, nMPISides_MINE, nMPISides_YOUR, etc.)
! ==========================================================================
  CALL GetData(p4est_ptr,DATAPtr)

#if MPI
  IF (nProcessors .GT. 1) THEN
    nNbProcs=FortranData%nNBProcs
    ! Here reallocate all arrays and redefine  all parameters
    
    SDEALLOCATE(nMPISides_send)
    ALLOCATE(nMPISides_send(       nNbProcs,2))
    
    nMPISides_send(:,1)     =nMPISides_MINE_Proc
    nMPISides_send(:,2)     =nMPISides_YOUR_Proc
    
    SDEALLOCATE(nMPISides_rec)
    ALLOCATE(nMPISides_rec(        nNbProcs,2))
    nMPISides_rec(:,1)      =nMPISides_YOUR_Proc
    nMPISides_rec(:,2)      =nMPISides_MINE_Proc
    
    SDEALLOCATE(OffsetMPISides_send)
    ALLOCATE(OffsetMPISides_send(0:nNbProcs,2))
    OffsetMPISides_send(:,1)=OffsetMPISides_MINE
    OffsetMPISides_send(:,2)=OffsetMPISides_YOUR
    
    SDEALLOCATE(OffsetMPISides_rec)
    ALLOCATE(OffsetMPISides_rec( 0:nNbProcs,2))
    OffsetMPISides_rec(:,1) =OffsetMPISides_YOUR
    OffsetMPISides_rec(:,2) =OffsetMPISides_MINE
    
    SDEALLOCATE(MPIRequest_U)
    SDEALLOCATE(MPIRequest_Flux)
    ALLOCATE(MPIRequest_U(nNbProcs,2)    )
    ALLOCATE(MPIRequest_Flux(nNbProcs,2) )
    MPIRequest_U      = MPI_REQUEST_NULL
    MPIRequest_Flux   = MPI_REQUEST_NULL
#if PARABOLIC && MPI
    SDEALLOCATE(MPIRequest_Lifting)
    ALLOCATE(MPIRequest_Lifting(nNbProcs,3,2))
    MPIRequest_Lifting = MPI_REQUEST_NULL
#endif /*PARABOLIC*/
  ENDIF
#endif  /*MPI*/

! ==========================================================================
! Reallocate the connectivity arrays (ElemToSide, SideToElem and MortarType)
! ... And fill them with updated information
! ==========================================================================

  CALL p4estSetMPIData()
  nElems=FortranData%nElems
  nSides=FortranData%nSides

  IF (nElemsOld .NE. nElems) THEN
    deallocate(ElemToSide)
    ALLOCATE(ElemToSide(2,6,FortranData%nElems))
  ENDIF

  IF (nSidesOld .NE. nSides) THEN
    deallocate(SideToElem)
    ALLOCATE(SideToElem(5,nSides))

    deallocate(MortarType)
    ALLOCATE(MortarType(2,nSides))
  ENDIF
  FortranData%EtSPtr = C_LOC(ElemToSide)
  FortranData%StEPtr = C_LOC(SideToElem)
  FortranData%MTPtr = C_LOC(MortarType)

  CALL SetEtSandStE(p4est_ptr,DATAPtr)

! ==========================================
! Reallocate the BC and MortarInfo arrays
! ... And fill them with updated information
! ==========================================
 
  CALL C_F_POINTER(FortranData%BCs, nBCsF,[FortranData%nBCSides])
  SDEALLOCATE(BC)
  ALLOCATE(BC(FortranData%nBCSides),SOURCE = nBCsF)

  nMortarSides    = FortranData%nMortarInnerSides +  FortranData%nMortarMPISides

  CALL C_F_POINTER(FortranData%MIPtr, MInfo,[2,5,nMortarSides])
  deallocate(MortarInfo)
  ALLOCATE(MortarInfo(MI_FLIP,0:4,nMortarSides),SOURCE = MInfo)

! ==========================================
! Reallocate remaining arrays of size nElems
! ==========================================

  ! Reallocate Arrays if the nElems was changed
  IF (nElemsOld .NE. nElems) THEN
    SDEALLOCATE(Ut); ALLOCATE(Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
    Ut = 0.0
    SDEALLOCATE(Source); ALLOCATE(Source(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
    Source = 0.0
    SDEALLOCATE(Metrics_fTilde); ALLOCATE(Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(Metrics_gTilde); ALLOCATE(Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(Metrics_hTilde); ALLOCATE(Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(dXGL_N); ALLOCATE(dXGL_N(3,3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(sJ); ALLOCATE(            sJ(  0:PP_N,0:PP_N,0:PP_N,nElems))
    
    SDEALLOCATE(DetJac_Ref); 
    NGeoRef=3*NGeo ! build jacobian at higher degree
    ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

    IF (ALLOCATED(dtElem))  THEN
      DEALLOCATE(dtElem); ALLOCATE(dtElem(nElems)); 
    ENDIF
    
    IF (ALLOCATED(ElemVol))  THEN 
      DEALLOCATE(ElemVol); ALLOCATE(ElemVol(nElems)); 
    ENDIF 
  
#if PARABOLIC
    IF (ALLOCATED(gradPx))  THEN 
      DEALLOCATE(gradPx); ALLOCATE(gradPx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
    ENDIF

    IF (ALLOCATED(gradPy))  THEN 
      DEALLOCATE(gradPy); ALLOCATE(gradPy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
    ENDIF

    IF (ALLOCATED(gradPz))  THEN 
      DEALLOCATE(gradPz); ALLOCATE(gradPz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
    ENDIF
#endif /* PARABOLIC */
  ENDIF  !IF (nElemsOld .NE. nElems)

! ==========================================
! Reallocate remaining arrays of size nSides
! ==========================================
  
  IF (nSidesOld .NE. nSides) THEN
! surface data
    SDEALLOCATE(Face_xGP); ALLOCATE(      Face_xGP(3,0:PP_N,0:PP_N,1:nSides))
    ! Face_xGP=0;
    SDEALLOCATE(NormVec); ALLOCATE(       NormVec(3,0:PP_N,0:PP_N,1:nSides))
    
    SDEALLOCATE(TangVec1); ALLOCATE(      TangVec1(3,0:PP_N,0:PP_N,1:nSides))
    
    SDEALLOCATE(TangVec2); ALLOCATE(      TangVec2(3,0:PP_N,0:PP_N,1:nSides))
    SDEALLOCATE(SurfElem); ALLOCATE(      SurfElem(  0:PP_N,0:PP_N,1:nSides))

#if PARABOLIC
   
    IF (ALLOCATED(FluxX))  THEN 
      DEALLOCATE(FluxX); ALLOCATE(FluxX        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
    ENDIF

    IF (ALLOCATED(FluxY))  THEN 
      DEALLOCATE(FluxY); ALLOCATE(FluxY        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
    ENDIF

    IF (ALLOCATED(FluxZ))  THEN 
      DEALLOCATE(FluxZ); ALLOCATE(FluxZ        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
    ENDIF


    IF (ALLOCATED(gradPx_master))  THEN 
      DEALLOCATE(gradPx_master); ALLOCATE(gradPx_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
    ENDIF

    IF (ALLOCATED(gradPy_master))  THEN 
      DEALLOCATE(gradPy_master); ALLOCATE(gradPy_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
    ENDIF

    IF (ALLOCATED(gradPz_master))  THEN 
      DEALLOCATE(gradPz_master); ALLOCATE(gradPz_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
    ENDIF
#endif /* PARABOLIC */
    IF (ALLOCATED(AnalyzeSide))  THEN 
        DEALLOCATE(AnalyzeSide); ALLOCATE(AnalyzeSide(1:nSides))
        AnalyzeSide=0;
    ENDIF
    DEALLOCATE(U_master)
    ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
  
    DEALLOCATE(Flux_master)
    ALLOCATE(Flux_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
  ENDIF !  IF (nSidesOld .NE. nSides) THEN

! ======================
! Recalculate parameters
! ======================
  ! Mesh parameters
  CALL RecalculateParameters(FortranData)
  
  ! Parameters of MOD_DG_Vars
  nDOFElem=(PP_N+1)**3
  nTotalU=PP_nVar*nDOFElem*nElems
  nTotal_face=(PP_N+1)*(PP_N+1)
  nTotal_vol=nTotal_face*(PP_N+1)
  nTotal_IP=nTotal_vol*nElems
  nTotalU=PP_nVar*nTotal_vol*nElems
  
! ======================================================
! Reallocate arrays of size firstSlaveSide:LastSlaveSide
! ======================================================
  
  IF ((LastSlaveSideOld .NE. LastSlaveSide) .OR. (firstSlaveSideOld .NE. firstSlaveSide)) THEN
#if PARABOLIC
    IF (ALLOCATED(gradPx_slave))  THEN 
      DEALLOCATE(gradPx_slave); 
      ALLOCATE(gradPx_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
    ENDIF
   
    IF (ALLOCATED(gradPy_slave))  THEN 
      DEALLOCATE(gradPy_slave); 
      ALLOCATE(gradPy_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
    ENDIF
   
    IF (ALLOCATED(gradPz_slave))  THEN 
      DEALLOCATE(gradPz_slave); 
      ALLOCATE(gradPz_slave (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
    ENDIF
#endif /* PARABOLIC */
    DEALLOCATE(U_SLAVE)
    ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
    DEALLOCATE(Flux_SLAVE)
    
    ALLOCATE(Flux_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
  ENDIF ! IF ((LastSlaveSideOld .NE. LastSlaveSide) .OR. firstSlaveSideOld .NE. firstSlaveSide)) 
  
! ======================
! Recompute metric terms
! ======================
  CALL CalcMetrics((/0/))
  
! ================================================================================================
! Recover U in the elements that were coarsened (their U is currently scaled by the Jacobian: J*U)
! ================================================================================================
  ! Scale coarsened elements
  do iElem=1, nElems
    if ( ElemWasCoarsened(iElem)>0 ) then
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        U(:,i,j,k,iElem) = U(:,i,j,k,iElem) * sJ(i,j,k,iElem)
      end do       ; end do       ; end do
    end if
  end do
  
! ===================================================
! Re-initialize other modules that depend on the mesh
! ===================================================
  call FinalizeMortarArrays()
  call InitMortarArrays()
#if SHOCKCAPTURE
  call InitShockCapturingAfterAdapt2(nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld)
#endif /*SHOCKCAPTURE*/
#if POSITIVITYPRES
  call FinalizePositivityPreservation()
  call InitPositivityPreservation()
#endif /*POSITIVITYPRES*/
  call InitEquationAfterAdapt()
! =========
! Finish up
! =========  
  CALL FinalizeBC()
  CALL InitBC()
  call free_data_memory(DataPtr)
  DEALLOCATE(ChangeElem)
  SDEALLOCATE(ElemWasCoarsened)
  NULLIFY(MInfo)
  NULLIFY(nBCsF)
END SUBROUTINE RunAMR


  SUBROUTINE RecalculateParameters(FortranData)
        USE MOD_AMR_Vars,            ONLY: P4EST_FORTRAN_DATA
        USE MOD_Mesh_Vars
#if MPI 
        
        USE MOD_MPI_Vars,             ONLY:  offsetElemMPI
        
#endif  /*MPI*/        
        USE MOD_Globals,             ONLY:  myrank
        USE MOD_P4est,               ONLY: p4estGetMPIData

        IMPLICIT NONE
        TYPE(p4est_fortran_data) :: FortranData
        INTEGER     ::firstMasterSide, lastMasterSide

        nBCSides            =   FortranData%nBCSides
        nElems              =   FortranData%nElems
        nGlobalElems        =   FortranData%nGlobalElems
        nSides              =   FortranData%nSides
        nMortarInnerSides   =   FortranData%nMortarInnerSides
        nInnerSides         =   FortranData%nInnerSides
        ! Must be set later, for parallel version
        nMortarMPISides        = FortranData%nMortarMPISides
        nMPISides_MINE         = FortranData%nMPISides_MINE
        nMPISides_YOUR         = FortranData%nMPISides_YOUR


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
        nMortarSides    = nMortarInnerSides +  nMortarMPISides
#if MPI 
        offsetElem  = offsetElemMPI(myrank)
#endif  /*MPI*/        
        
    END SUBROUTINE RecalculateParameters
!===================================================================================================================================
!> Does an L2 projection of the solution from 8 elements to 1
!> ATTENTION: The projected solution must be scaled with the Jacobians because the projection integrals are evaluated on the fine and coarse elements
!===================================================================================================================================
  subroutine ProjectSolution_Coarsening(Unew, Uold, sJold)
    use MOD_PreProc     , only: PP_N
    use MOD_AMR_Vars    , only: M_1_0,M_2_0
    implicit none
    !-arguments-----------------------------------------------
    real, intent(inout) :: Unew(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    real, intent(in)    :: Uold(PP_nVar,0:PP_N,0:PP_N,0:PP_N,8)
    real, intent(in)    :: sJold       (0:PP_N,0:PP_N,0:PP_N,8)
    !-local-variables-----------------------------------------
    integer :: l,m,s,p,q,r
    real :: JUold(PP_nVar,0:PP_N,0:PP_N,0:PP_N,8)
    !---------------------------------------------------------
    
    ! Scale U with the Jacobian
    do r=1, 8 ; do s=0, PP_N ; do m=0, PP_N ; do l=0, PP_N
      JUold(:,l,m,s,r) = Uold(:,l,m,s,r) / sJold(l,m,s,r)
    end do       ; end do       ; end do       ; end do
    
    ! Project to new (big) element (The scaling with the new Jacobian is missing)
    Unew = 0.0
    do r=0, PP_N ; do q=0, PP_N ; do p=0, PP_N
      do s=0, PP_N ; do m=0, PP_N ; do l=0, PP_N
        Unew(:,p,q,r) = Unew(:,p,q,r) + M_1_0(l,p)*M_1_0(m,q)*M_1_0(s,r) * JUold(:,l,m,s,1) &
                                      + M_2_0(l,p)*M_1_0(m,q)*M_1_0(s,r) * JUold(:,l,m,s,2) &
                                      + M_1_0(l,p)*M_2_0(m,q)*M_1_0(s,r) * JUold(:,l,m,s,3) &
                                      + M_2_0(l,p)*M_2_0(m,q)*M_1_0(s,r) * JUold(:,l,m,s,4) &
                                      + M_1_0(l,p)*M_1_0(m,q)*M_2_0(s,r) * JUold(:,l,m,s,5) &
                                      + M_2_0(l,p)*M_1_0(m,q)*M_2_0(s,r) * JUold(:,l,m,s,6) &
                                      + M_1_0(l,p)*M_2_0(m,q)*M_2_0(s,r) * JUold(:,l,m,s,7) &
                                      + M_2_0(l,p)*M_2_0(m,q)*M_2_0(s,r) * JUold(:,l,m,s,8) 
      end do       ; end do       ; end do
    end do       ; end do       ; end do
    
  end subroutine ProjectSolution_Coarsening
!===================================================================================================================================
!> Projects a field from 1 element to 8. Can be used for
!>  1. Solution (L2 projection): The projected solution DOES NOT HAVE to be scaled with the Jacobians because the projection 
!>                               integrals are evaluated only on the fine elements
!>  2. Coordinates: Using the interpolation matrices
!===================================================================================================================================
  subroutine InterpolateSolution_Refinement(nVar,Unew, Uold)
    use MOD_PreProc     , only: PP_N
    use MOD_AMR_vars    , only: M_0_1,M_0_2
    implicit none
    !-arguments-----------------------------------------------
    integer, intent(in)    :: nVar
    real   , intent(inout) :: Unew(nVar,0:PP_N,0:PP_N,0:PP_N,8)
    real   , intent(in)    :: Uold(nVar,0:PP_N,0:PP_N,0:PP_N)
    !-local-variables-----------------------------------------
    integer :: l,m,s,p,q,r
    !---------------------------------------------------------
    
    Unew = 0.0
    do r=0, PP_N ; do q=0, PP_N ; do p=0, PP_N
      do s=0, PP_N ; do m=0, PP_N ; do l=0, PP_N
        Unew(:,p,q,r,1) = Unew(:,p,q,r,1) + M_0_1(l,p)*M_0_1(m,q)*M_0_1(s,r) * Uold(:,l,m,s)
        Unew(:,p,q,r,2) = Unew(:,p,q,r,2) + M_0_2(l,p)*M_0_1(m,q)*M_0_1(s,r) * Uold(:,l,m,s)
        Unew(:,p,q,r,3) = Unew(:,p,q,r,3) + M_0_1(l,p)*M_0_2(m,q)*M_0_1(s,r) * Uold(:,l,m,s)
        Unew(:,p,q,r,4) = Unew(:,p,q,r,4) + M_0_2(l,p)*M_0_2(m,q)*M_0_1(s,r) * Uold(:,l,m,s)
        Unew(:,p,q,r,5) = Unew(:,p,q,r,5) + M_0_1(l,p)*M_0_1(m,q)*M_0_2(s,r) * Uold(:,l,m,s)
        Unew(:,p,q,r,6) = Unew(:,p,q,r,6) + M_0_2(l,p)*M_0_1(m,q)*M_0_2(s,r) * Uold(:,l,m,s)
        Unew(:,p,q,r,7) = Unew(:,p,q,r,7) + M_0_1(l,p)*M_0_2(m,q)*M_0_2(s,r) * Uold(:,l,m,s)
        Unew(:,p,q,r,8) = Unew(:,p,q,r,8) + M_0_2(l,p)*M_0_2(m,q)*M_0_2(s,r) * Uold(:,l,m,s)
      end do       ; end do       ; end do
    end do       ; end do       ; end do
    
  end subroutine InterpolateSolution_Refinement
!===================================================================================================================================
!> Does a simple interpolation (extrapolation) from element 1 to the big element. This works if Ngeo<=N
!===================================================================================================================================
  subroutine InterpolateCoords_Coarsening(Xnew, Xold)
    use MOD_PreProc     , only: PP_N
    use MOD_AMR_Vars    , only: Vdm_Interp_1_0, Vdm_Interp_2_0, N_2
    implicit none
    !-arguments-----------------------------------------------
    real, intent(inout) :: Xnew (1:3,0:PP_N,0:PP_N,0:PP_N)
    real, intent(in)    :: Xold (1:3,0:PP_N,0:PP_N,0:PP_N,1:8)
    !-local-variables-----------------------------------------
    integer :: l,m,s,p,q,r
    !---------------------------------------------------------
    
    ! Fill Xnew (subcell by subcell)
    Xnew = 0.0
    do s=0, PP_N ; do m=0, PP_N ; do l=0, PP_N
      ! Subcell 1
      do r=0, N_2 ; do q=0, N_2 ; do p=0, N_2
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_1_0(p,l)*Vdm_Interp_1_0(q,m)*Vdm_Interp_1_0(r,s) * Xold(:,l,m,s,1) 
      end do       ; end do       ; end do
      ! Subcell 2
      do r=0, N_2 ; do q=0, N_2 ; do p=N_2+1, PP_N
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_2_0(p,l)*Vdm_Interp_1_0(q,m)*Vdm_Interp_1_0(r,s) * Xold(:,l,m,s,2)
      end do       ; end do       ; end do
      ! Subcell 3
      do r=0, N_2 ; do q=N_2+1, PP_N ; do p=0, N_2
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_1_0(p,l)*Vdm_Interp_2_0(q,m)*Vdm_Interp_1_0(r,s) * Xold(:,l,m,s,3)
      end do       ; end do       ; end do
      ! Subcell 4
      do r=0, N_2 ; do q=N_2+1,PP_N ; do p=N_2+1,PP_N
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_2_0(p,l)*Vdm_Interp_2_0(q,m)*Vdm_Interp_1_0(r,s) * Xold(:,l,m,s,4)
      end do       ; end do       ; end do
      ! Subcell 5
      do r=N_2+1, PP_N ; do q=0, N_2 ; do p=0, N_2
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_1_0(p,l)*Vdm_Interp_1_0(q,m)*Vdm_Interp_2_0(r,s) * Xold(:,l,m,s,5)
      end do       ; end do       ; end do
      ! Subcell 6
      do r=N_2+1, PP_N ; do q=0, N_2 ; do p=N_2+1, PP_N
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_2_0(p,l)*Vdm_Interp_1_0(q,m)*Vdm_Interp_2_0(r,s) * Xold(:,l,m,s,6)
      end do       ; end do       ; end do
      ! Subcell 7
      do r=N_2+1, PP_N ; do q=N_2+1, PP_N ; do p=0, N_2
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_1_0(p,l)*Vdm_Interp_2_0(q,m)*Vdm_Interp_2_0(r,s) * Xold(:,l,m,s,7)
      end do       ; end do       ; end do
      ! Subcell 8
      do r=N_2+1, PP_N ; do q=N_2+1, PP_N ; do p=N_2+1, PP_N
        Xnew(:,p,q,r) = Xnew(:,p,q,r) + Vdm_Interp_2_0(p,l)*Vdm_Interp_2_0(q,m)*Vdm_Interp_2_0(r,s) * Xold(:,l,m,s,8)
      end do       ; end do       ; end do
    end do       ; end do       ; end do
    
  end subroutine InterpolateCoords_Coarsening

!============================================================================================================================
!> Balances the loads using p4est
!============================================================================================================================
SUBROUTINE LoadBalancingAMR(ElemWasCoarsened,new_nElems)
  ! MODULES
  USE MOD_Globals
  USE MOD_AMR_Vars
  USE MOD_P4EST
  USE MOD_DG_Vars,            ONLY: U
  USE MOD_Mesh_Vars,          ONLY: Elem_xGP
#if SHOCK_NFVSE
  use MOD_NFVSE_Vars        , only: alpha
#endif /*SHOCK_NFVSE*/
  USE MOD_PreProc,            ONLY: PP_N
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  !----------------------------------------------------------------------------------------------------------------------------------
  ! ARGUMENTS
  integer, intent(inout), allocatable, target :: ElemWasCoarsened(:)
  integer, intent(out) :: new_nElems
  !----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  REAL,ALLOCATABLE, TARGET :: Elem_xGP_New(:,:,:,:,:), U_New(:,:,:,:,:), Alpha_New(:)
  integer,allocatable,target :: ElemWasCoarsened_new(:)
  !============================================================================================================================
  TYPE(p4est_balance_datav2), TARGET :: BalanceData;
  
  IF (.NOT. UseAMR) THEN
    RETURN;
  ENDIF
  
  
  BalanceData%DataSize = sizeof(U(:,:,:,:,1))
  BalanceData%GPSize = sizeof(Elem_xGP(:,:,:,:,1))
  BalanceData%CoarseSize = sizeof(ElemWasCoarsened(1))

  
  BalanceData%Uold_Ptr = C_LOC(U)
  BalanceData%ElemxGPold_Ptr = C_LOC(Elem_xGP)
#if SHOCK_NFVSE
  BalanceData%AlphaOld_Ptr = C_LOC(alpha)
#endif /*SHOCK_NFVSE*/
  BalanceData%ElemWasCoarsened_old = C_LOC(ElemWasCoarsened)
  
  CALL p4est_loadbalancing_init(P4EST_PTR, C_LOC(BalanceData))
  
  
  new_nElems = BalanceData%nElems
  ALLOCATE(U_New(PP_nVar,0:PP_N,0:PP_N,0:PP_N,BalanceData%nElems))
  BalanceData%Unew_Ptr = C_LOC(U_New)
  
  ALLOCATE(Elem_xGP_New(3,0:PP_N,0:PP_N,0:PP_N,BalanceData%nElems))
  BalanceData%ElemxGPnew_Ptr = C_LOC(Elem_xGP_New)
  
  allocate(ElemWasCoarsened_new(BalanceData%nElems))
  BalanceData%ElemWasCoarsened_new = C_LOC(ElemWasCoarsened_new)
  
#if SHOCK_NFVSE
  ALLOCATE(Alpha_New(BalanceData%nElems))
  BalanceData%AlphaNew_Ptr = C_LOC(Alpha_New)
#endif /*SHOCK_NFVSE*/
  
  CALL p4est_loadbalancing_go(P4EST_PTR, C_LOC(BalanceData))

  CALL MOVE_ALLOC(Elem_xGP_New, Elem_xGP)
  CALL MOVE_ALLOC(U_New, U)
  CALL MOVE_ALLOC(ElemWasCoarsened_new, ElemWasCoarsened)
  
#if SHOCK_NFVSE
  CALL MOVE_ALLOC(Alpha_New, alpha)
#endif /*SHOCK_NFVSE*/
  
  CALL p4est_ResetElementNumber(P4EST_PTR)
  
END SUBROUTINE LoadBalancingAMR
!============================================================================================================================
!> Save mesh to HDF5 file
!============================================================================================================================
SUBROUTINE SaveMesh(FileString)
  ! MODULES
  USE MOD_AMR_Vars
  USE MOD_P4EST
  USE MOD_PreProc,            ONLY: PP_N
  USE MOD_Mesh_vars,           ONLY: nElems, nGlobalElems, offsetElem, BoundaryName, BoundaryType, nBCs, Elem_xGP, NGeo
  USE MOD_Globals,             only: myrank,nProcessors, MPIRoot, UNIT_stdOut
  USE MOD_IO_HDF5
  USE MOD_HDF5_Output,            only: WriteHeader, WriteAttribute, WriteArray, GatheredWriteArray
  USE MOD_Interpolation_Vars, ONLY: NodeType,NodeTypeVISU
  USE MOD_Interpolation,      ONLY: GetVandermonde
  USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN)  :: FileString !< (IN) mesh filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  TYPE(C_PTR) :: DataPtr;
  TYPE(p4est_save_data), POINTER :: FortranDataSave
  INTEGER, POINTER :: ElInfoF(:,:), SiInfoF(:,:)  ! ElemInfo and SideInfo to pass the data
  INTEGER, POINTER :: OffSideMPIF(:), OffSideArrIndMPIF(:)  ! OffsetSideMPI and OffsetSideArrIndexMPI to pass the data
  INTEGER, ALLOCATABLE :: OffsetSideMPI(:), OffsetSideArrIndexMPI(:), ElemInfoW(:,:)
  INTEGER          :: iElem, nIndexSide, mpisize, FirstElemInd, LastElemInd, nLocalIndSides, nGlobalIndSides, OffsetIndSides
  INTEGER(HID_T)                 :: DSet_ID,FileSpace,HDF5DataType
  INTEGER(HSIZE_T)               :: DimsM(2)
  CHARACTER(LEN=255), ALLOCATABLE:: BCNames(:)
  INTEGER, ALLOCATABLE           :: BCType(:,:)
  INTEGER                        :: i,j,k, index
  REAL,ALLOCATABLE               :: NodeCoords(:,:)
  REAL,ALLOCATABLE               :: NodeCoordsTMP(:,:,:,:)
  REAL,ALLOCATABLE               :: Vdm_FromNodeType_toNVisu(:,:)
  integer                        :: NGeo_new
  !============================================================================================================================
  
  ! Set NGeo of the mesh to save
  NGeo_new = min(PP_N,NGeo)
  
  IF (.NOT. UseAMR) RETURN

  IF(MPIRoot)THEN
    WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE MESH TO _MESH.H5 FILE: '
    WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')TRIM(FileString)
  END IF

  mpisize = nProcessors
  DATAPtr=SaveMeshP4(p4est_ptr)
  CALL C_F_POINTER(DataPtr, FortranDataSave)

  CALL C_F_POINTER(FortranDataSave%ElemInfo, ElInfoF,[6,nElems])
  !read local ElemInfo from data file
  FirstElemInd=offsetElem+1
  LastElemInd=offsetElem+nElems
  ALLOCATE(ElemInfoW(6,nElems), SOURCE  = ElInfoF)
  ElemInfoW = ElInfoF
  ! ALLOCATE(ElemInfo1(6,nElems))
  
  if (nElems > 0) then
    nIndexSide = ElInfoF(4,nElems)
    CALL C_F_POINTER(FortranDataSave%SideInfo, SiInfoF,[5,nIndexSide])
  else
    allocate ( SiInfoF(5,0) )
  end if
  
  CALL C_F_POINTER(FortranDataSave%OffsetSideMPI, OffSideMPIF,[mpisize+1])
  
  ALLOCATE(OffsetSideMPI(0:mpisize), SOURCE  = OffSideMPIF)

  CALL C_F_POINTER(FortranDataSave%OffsetSideArrIndexMPI, OffSideArrIndMPIF,[mpisize+1])
  ALLOCATE(OffsetSideArrIndexMPI(0:mpisize), SOURCE  = OffSideArrIndMPIF)
  

  ElInfoF(3,:) = ElInfoF(3,:) + OffsetSideArrIndexMPI(Myrank)
  ElInfoF(4,:) = ElInfoF(4,:) + OffsetSideArrIndexMPI(Myrank)
  
  ! 1. Create Mesh file
  !Generate skeleton for the file with all relevant data on a single proc (MPIRoot)

  IF(MPIRoot)THEN
   ! Create file
    CALL OpenDataFile(TRIM(FileString),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'Version',1,RealScalar=1.0)
    CALL WriteAttribute(File_ID,'Ngeo',1,IntScalar=NGeo_new)
    CALL WriteAttribute(File_ID,'nElems',1,IntScalar=nGlobalElems)
    CALL WriteAttribute(File_ID,'nSides',1,IntScalar=OffsetSideArrIndexMPI(mpisize))
    CALL WriteAttribute(File_ID,'nNodes',1,IntScalar=nGlobalElems*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1))
    CALL WriteAttribute(File_ID,'nUniqueSides',1,IntScalar=OffsetSideMPI(mpisize))
    CALL WriteAttribute(File_ID,'nUniqueNodes',1,IntScalar=343)
    CALL WriteAttribute(File_ID,'nBCs',1,IntScalar=nBCs)
    
    CALL WriteAttribute(File_ID,'hasAMRmortars',1,IntScalar=1)

     DimsM=(/6, nGlobalElems/)
    CALL H5SCREATE_SIMPLE_F(2, DimsM, FileSpace, iError)
    HDF5DataType=H5T_NATIVE_INTEGER
    CALL H5DCREATE_F(File_ID,'ElemInfo', HDF5DataType, FileSpace, DSet_ID, iError)
    ! Close the filespace and the dataset
    CALL H5DCLOSE_F(Dset_id, iError)
    CALL H5SCLOSE_F(FileSpace, iError)  

    DimsM=(/5, OffsetSideArrIndexMPI(mpisize)/)
    CALL H5SCREATE_SIMPLE_F(2, DimsM, FileSpace, iError)
    HDF5DataType=H5T_NATIVE_INTEGER
    CALL H5DCREATE_F(File_ID,'SideInfo', HDF5DataType, FileSpace, DSet_ID, iError)
    ! Close the filespace and the dataset
    CALL H5DCLOSE_F(Dset_id, iError)
    CALL H5SCLOSE_F(FileSpace, iError)  


    ! ALLOCATE(NodeCoords)
    DimsM=(/3, nGlobalElems*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1)/)
    CALL H5SCREATE_SIMPLE_F(2, DimsM, FileSpace, iError)
    HDF5DataType=H5T_NATIVE_DOUBLE
    CALL H5DCREATE_F(File_ID,'NodeCoords', HDF5DataType, FileSpace, DSet_ID, iError)
    ! Close the filespace and the dataset
    CALL H5DCLOSE_F(Dset_id, iError)
    CALL H5SCLOSE_F(FileSpace, iError)  
  
  ALLOCATE(BCType(4,nBCs))

  ALLOCATE(BCNames(nBCs))
  BCType=0
  BCNames = BoundaryName
  BCType(1,:) = BoundaryType(:,BC_TYPE)  
  BCType(3,:) = BoundaryType(:,BC_STATE)
  BCType(4,:) = BoundaryType(:,BC_ALPHA)
  DimsM=(/4, nBCs/)

    CALL H5SCREATE_SIMPLE_F(2, DimsM, FileSpace, iError)
    HDF5DataType=H5T_NATIVE_INTEGER
    CALL H5DCREATE_F(File_ID,'BCType', HDF5DataType, FileSpace, DSet_ID, iError)
    ! Close the filespace and the dataset
    CALL H5DCLOSE_F(Dset_id, iError)
    CALL H5SCLOSE_F(FileSpace, iError)  

  CALL WriteArray('BCType',2,(/4, nBCs/),(/4, nBCs/),&
                                                (/0, 0/),collective = .FALSE.,IntArray =BCType)

  CALL WriteArray('BCNames',1,(/nBCs/),(/nBCs/),&
                                                (/0/),collective = .FALSE.,StrArray =BCNames)
                                                
  DEALLOCATE(BCNames,BCType)

    CALL CloseDataFile()
  END IF !  IF(MPIRoot)THEN
#if MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

DO iElem = 1, nElems
  ElInfoF(5,iElem) = (iElem + offsetElem - 1)*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1)
  ElInfoF(6,iElem) = (iElem + offsetElem)*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1)
ENDDO

IF (NGeo_new .EQ. 1) THEN
  ElInfoF(1,:)=108
ENDIF
CALL GatheredWriteArray(FileString,create=.FALSE.,&
                        DataSetName='ElemInfo', rank=2,  &
                        nValGlobal=(/6,nGlobalElems/),&
                        nVal=      (/6,nElems      /),&
                        offset=    (/0   ,offSetElem  /),&
                        collective=.TRUE.,IntArray=ElInfoF)
#if MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

nLocalIndSides = FortranDataSave%nSidesArrIndex
nGlobalIndSides = OffsetSideArrIndexMPI(mpisize)
OffsetIndSides = OffsetSideArrIndexMPI(myrank)

CALL GatheredWriteArray(FileString,create=.FALSE.,&
                        DataSetName='SideInfo', rank=2,  &
                        nValGlobal=(/5,nGlobalIndSides/),&
                        nVal=      (/5,nLocalIndSides      /),&
                        offset=    (/0   ,OffsetIndSides  /),&
                        collective=.TRUE.,IntArray=SiInfoF)

ALLOCATE(NodeCoords(3,nElems*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1)))
ALLOCATE(NodeCoordsTMP(3,0:NGeo_new,0:NGeo_new,0:NGeo_new))
ALLOCATE(Vdm_FromNodeType_toNVisu(0:NGeo_new,0:PP_N))

CALL GetVandermonde(PP_N, NodeType, NGeo_new,      NodeTypeVISU, Vdm_FromNodeType_toNVisu)
    

index = 1
DO iElem=1,nElems
 
  CALL ChangeBasis3D(3,PP_N,NGeo_new,Vdm_FromNodeType_toNVisu,Elem_xGP(:,:,:,:,iElem),NodeCoordsTMP(:,:,:,:))
  DO k=0,NGeo_new
    DO j=0,NGeo_new
      DO i=0,NGeo_new

        NodeCoords(:,index) = NodeCoordsTMP(:,i,j,k)
        index = index + 1
      ENDDO
    ENDDO
  ENDDO
ENDDO
#if MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
CALL GatheredWriteArray(FileString,create=.FALSE.,&
                        DataSetName='NodeCoords', rank=2,  &
                        nValGlobal=(/3,nGlobalElems*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1)/),&
                        nVal=      (/3,nElems*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1)/),&
                        offset=    (/0   ,offSetElem*(NGeo_new+1)*(NGeo_new+1)*(NGeo_new+1)/),&
                        collective=.TRUE.,RealArray=NodeCoords)

SDEALLOCATE(NodeCoords)
SDEALLOCATE(NodeCoordsTMP)
SDEALLOCATE(Vdm_FromNodeType_toNVisu)

#if MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

  DEALLOCATE(OffsetSideMPI)
  DEALLOCATE(OffsetSideArrIndexMPI)
  DEALLOCATE(ElemInfoW)
  CALL free_savemesh_memory(DATAPtr)
  NULLIFY(ElInfoF)
  NULLIFY(SiInfoF)
  NULLIFY(OffSideMPIF)
  NULLIFY(OffSideArrIndMPIF)
END SUBROUTINE SaveMesh

!============================================================================================================================
!> Deallocate AMR data.
!============================================================================================================================
SUBROUTINE FinalizeAMR()
! MODULES
USE MOD_AMR_Vars
USE MOD_P4EST
use MOD_Indicators
IMPLICIT NONE
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
IF (.NOT. UseAMR) THEN
  RETURN;
ENDIF
SDEALLOCATE(Vdm_Interp_1_0)
SDEALLOCATE(Vdm_Interp_2_0)
SDEALLOCATE(M_0_1)
SDEALLOCATE(M_0_2)
SDEALLOCATE(M_1_0)
SDEALLOCATE(M_2_0)
CALL p4est_destroy(P4EST_PTR);
CALL p4est_connectivity_destroy(CONNECTIVITY_PTR)
CALL p4est_finalize()
call AMR_Indicator % destruct
AMRInitIsDone = .FALSE.
END SUBROUTINE FinalizeAMR

END MODULE MOD_AMR
#endif
