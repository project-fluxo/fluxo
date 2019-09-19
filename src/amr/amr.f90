#include "defines.h"
! #include "amr_f.h"
! #include "mpif.h"
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


INTERFACE FinalizeAMR
  MODULE PROCEDURE FinalizeAMR
END INTERFACE

INTERFACE LoadBalancingAMR
  MODULE PROCEDURE LoadBalancingAMR
END INTERFACE

INTERFACE LoadBalancingAMRold
  MODULE PROCEDURE LoadBalancingAMRold
END INTERFACE

INTERFACE SaveMesh
  MODULE PROCEDURE SaveMesh
END INTERFACE

PUBLIC::InitAMR, SaveMesh
PUBLIC::FinalizeAMR

!==================================================================================================================================
PUBLIC::RunAMR
PUBLIC::InitAMR_Connectivity
PUBLIC::DefineParametersAMR
PUBLIC::LoadBalancingAMR
PUBLIC::LoadBalancingAMRold
! INTEGER :: COUNT =0 
CONTAINS


!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAMR()
  ! MODULES
  USE MOD_Globals
  USE MOD_ReadInTools ,ONLY: prms
  IMPLICIT NONE
  !==================================================================================================================================
  CALL prms%SetSection("AMR")
  CALL prms%CreateLogicalOption( 'UseAMR',          "Use AMR for solution.",&
                                                      '.FALSE.')
  CALL prms%CreateStringOption('p4estFile', "(relative) path to p4ests connectivity file (mandatory).")

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
    USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR
    ! USE MOD_P4EST_Binding, ONLY: p4_initvars
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !==================================================================================================================================
    INTEGER RET

    IF(AMRInitIsDone) THEN
      CALL CollectiveStop(__STAMP__,&
        'InitAMR not ready to be called or already called.')
    END IF

    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT AMR...'

    ! CALL p4_initvars(IntSize)

    UseAMR = GETLOGICAL('UseAMR','.FALSE.')
    IF (UseAMR) THEN
      p4estFile = GETSTR('p4estFile')
      ! PRINT *, "!!!!!!!!!!!!!!!!!!!!!!!!USE  AMR = ", UseAMR
      ! PRINT *, "!!!!!!!!!!!!!!!!!!!!!!!!p4estFile  AMR = ", p4estFile
    ELSE
      SWRITE(UNIT_stdOut,'(A)') ' AMR will not be used! UseAMR is FALSE!'
      RETURN;
    ENDIF
    RET=P4EST_INIT(MPI_COMM_WORLD);

    AMRInitIsDone=.TRUE.
    SWRITE(UNIT_stdOut,'(A)')' INIT AMR DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')
    CALL InitAMR_Connectivity()
    CALL InitAMR_P4est()
    ! CALL EXIT()
    ! CALL EXIT()
END SUBROUTINE InitAMR

!==================================================================================================================================
!> Routine creating and 
!==================================================================================================================================
SUBROUTINE InitAMR_Connectivity()
    ! MODULES
    USE MOD_Globals
    USE MOD_Mesh_Vars,             ONLY: MeshFile
    ! USE MOD_P4Mesh_ReadIn,         ONLY: P4ReadMesh,ReadMeshHeader
    USE MOD_P4EST
    USE MOD_AMR_Vars,              ONLY: connectivity_ptr
    USE MODH_Mesh,                 ONLY: InitMesh, FinalizeMesh
    USE MODH_Mesh_ReadIn,        ONLY: ReadMeshHeader,ReadMeshFromHDF5
    ! USE MOD_IO_HDF5,              ONLY: nDims
    ! USE MODH_Mesh_ReadIn,        ONLY: ReadMeshFromHDF5,ReadMeshHeader
    !  USE MOD_P4EST_Binding, ONLY: p4_initvars
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    INTEGER CONN_OWNER
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !==================================================================================================================================

    ! MeshFile = GETSTR('MeshFile')
        ! print *, "~~~~~~~~~~~~~~~~~~~~~~ n Dims ~~~~~~~~~~~~~~~~~~:",nDims
        CONN_OWNER=0;
    ! IF(MPIroot)THEN
    !   CONN_OWNER=myrank
    ! ENDIF
        CALL InitMesh()
    ! ! Create Connectivity
        CALL ReadMeshHeader(MeshFile)   ! read mesh header file including BCs
    !  print *, " myrank = ", myrank
        CALL ReadMeshFromHDF5(MeshFile)
        ! The result is connectivity PTR


    connectivity_ptr=P4EST_CONN_BCAST(connectivity_ptr, CONN_OWNER, MPI_COMM_WORLD)

    CALL FinalizeMesh()
   
END SUBROUTINE InitAMR_Connectivity




SUBROUTINE RunAMR(ElemToRefineAndCoarse)
  USE MOD_Globals
  USE MOD_Analyze_Vars,        ONLY: ElemVol
  USE MOD_AMR_Vars,           ONLY: P4EST_FORTRAN_DATA, P4est_ptr, UseAMR
  USE MOD_Mesh_Vars,          ONLY: Elem_xGP, ElemToSide, SideToElem, Face_xGP, NormVec, TangVec1, TangVec2
  USE MOD_Mesh_Vars,          ONLY: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde,dXGL_N, sJ, SurfElem, nBCs
  USE MOD_P4est,              ONLY: free_data_memory, RefineCoarse, GetData, p4estSetMPIData
  USE MOD_Metrics,            ONLY: CalcMetrics
  USE MOD_DG_Vars,            ONLY: U,Ut,nTotalU, nTotal_vol, nTotal_IP, nTotal_face, nDOFElem, U_master, U_SLAVE, Flux_master, Flux_slave
  USE MOD_Mesh_Vars,          ONLY: AnalyzeSide, MortarInfo, MortarType, NGeo, DetJac_Ref, BC
  USE MOD_TimeDisc_Vars,      ONLY:   dtElem
  USE MOD_Mesh_Vars,          ONLY: LastSlaveSide, firstSlaveSide, nSides, nElems!, firstMortarInnerSide, lastMortarInnerSide
  USE MOD_MPI_Vars,           ONLY: NbProc , nMPISides_MINE_Proc, nMPISides_YOUR_Proc, offsetMPISides_YOUR, offsetMPISides_MINE 
  USE MOD_MPI_Vars,           ONLY: nMPISides_Proc, nMPISides_send, nMPISides_rec, OffsetMPISides_send, OffsetMPISides_rec
  USE MOD_MPI_Vars,           ONLY: MPIRequest_U, MPIRequest_Flux, nNbProcs, offsetElemMPI
  USE MOD_Globals ,           ONLY: nProcessors, MPIroot, myrank
  USE MOD_Mesh_Vars,          ONLY: nMPISides_MINE, nMPISides_YOUR,nGlobalElems, firstMPISide_YOUR, firstMPISide_MINE,firstMortarMPISide
#if PARABOLIC
  USE  MOD_Lifting_Vars
  USE  MOD_MPI_Vars,          ONLY: MPIRequest_Lifting
#endif /* PARABOLIC */
  USE, INTRINSIC :: ISO_C_BINDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  REAL,ALLOCATABLE :: Elem_xGP_New(:,:,:,:,:), U_New(:,:,:,:,:)
  REAL,ALLOCATABLE :: Ut_New(:,:,:,:,:), Metrics_fTilde_New(:,:,:,:,:)
  REAL,ALLOCATABLE :: Metrics_gTilde_New(:,:,:,:,:), Metrics_hTilde_New(:,:,:,:,:)

  REAL,ALLOCATABLE :: dXGL_N_New(:,:,:,:,:,:), sJ_new(:,:,:,:)

  REAL,ALLOCATABLE :: DetJac_Ref_New(:,:,:,:,:)
  REAL,ALLOCATABLE :: Face_xGP_New(:,:,:,:)
  
  REAL,ALLOCATABLE :: NormVec_New(:,:,:,:), TangVec1_New(:,:,:,:)
  REAL,ALLOCATABLE :: TangVec2_New(:,:,:,:), SurfElem_New(:,:,:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, ALLOCATABLE, TARGET, Optional  :: ElemToRefineAndCoarse(:) ! positive Number - refine, negative - coarse, 0 - do nothing
  INTEGER :: PP, NELM, Ie, nVar;
  INTEGER, ALLOCATABLE, TARGET ::  ElementToCalc(:)
  TYPE(C_PTR) :: DataPtr;
  TYPE(p4est_fortran_data), POINTER :: DataF
  ! TYPE (P4EST_FORTRAN_DATA) :: 
  INTEGER, POINTER :: EtSF(:,:,:), MType(:,:), MInfo(:,:,:), StEF(:,:), ChangeElem(:,:), ChangeSide(:,:)
  INTEGER, POINTER :: nNBProcF(:), nMPISides_ProcF(:), nMPISides_MINE_ProcF(:), nMPISides_YOUR_ProcF(:)
  INTEGER, POINTER :: offsetMPISides_MINEF(:), offsetMPISides_YOURF(:), nBCsF(:)
  INTEGER :: i,j,iElem, PP_N, nMortarSides, NGeoRef


  IF (.NOT. UseAMR) THEN
    RETURN;
  ENDIF
  ! COUNT = COUNT + 1
  nVar=size(U(:,1,1,1,1))
  PP=size(U(1,:,1,1,1))-1
  NELM=size(U(1,1,1,1,:))
  PP_N=PP
  
  IF (PRESENT(ElemToRefineAndCoarse)) THEN
    DATAPtr=RefineCoarse(p4est_ptr,C_LOC(ElemToRefineAndCoarse))
  ELSE
    DATAPtr = GetData(p4est_ptr)
  ENDIF 

  
  CALL p4estSetMPIData()
  
  CALL C_F_POINTER(DataPtr, DataF)
  
  ! PRINT *, "DataF%nElems = ", DataF%nElems, myrank


  CALL C_F_POINTER(DataF%EtSPtr, EtSF,[2,6,DataF%nElems])
  deallocate(ElemToSide)
  ! ALLOCATE(ElemToSide(2,6,DataF%nElems),SOURCE=EtSF)
  ALLOCATE(ElemToSide(2,6,DataF%nElems), SOURCE  = EtSF)


  ! i=DataF%nBCSides;
 
  CALL C_F_POINTER(DataF%BCs, nBCsF,[DataF%nBCSides])
  SDEALLOCATE(BC)
  ALLOCATE(BC(DataF%nBCSides),SOURCE = nBCsF)


 


  CALL C_F_POINTER(DataF%StEPtr, StEF,[5,DataF%nSides])
  deallocate(SideToElem)
  ALLOCATE(SideToElem(5,DataF%nSides),SOURCE = StEF)


  CALL C_F_POINTER(DataF%MTPtr, MType,[2,DataF%nSides])
  deallocate(MortarType)
  ALLOCATE(MortarType(2,DataF%nSides),SOURCE = MType)
!   ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo


    nMortarSides    = DataF%nMortarInnerSides +  DataF%nMortarMPISides

  CALL C_F_POINTER(DataF%MIPtr, MInfo,[2,4,nMortarSides])
  deallocate(MortarInfo)
  ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides),SOURCE = MInfo)

  



  CALL C_F_POINTER(DataF%ChngElmPtr, ChangeElem,[8,DataF%nElems])



  IF (nProcessors .GT. 1) THEN
        nNbProcs=DataF%nNBProcs

     

      ! Here reallocate all arrays and redefine  all parameters
      CALL C_F_POINTER(DataF%nNbProc, nNbProcF,[nNbProcs])
      SDEALLOCATE(NbProc)
      ALLOCATE(NbProc(1:nNbProcs), SOURCE = nNbProcF)
   
      CALL C_F_POINTER(DataF%nMPISides_Proc, nMPISides_ProcF,[nNbProcs])
      SDEALLOCATE(nMPISides_Proc)
      ALLOCATE(nMPISides_Proc(1:nNbProcs),SOURCE = nMPISides_ProcF)

      CALL C_F_POINTER(DataF%nMPISides_MINE_Proc, nMPISides_MINE_ProcF,[nNbProcs])
      SDEALLOCATE(nMPISides_MINE_Proc)
      ALLOCATE(nMPISides_MINE_Proc(1:nNbProcs),SOURCE = nMPISides_MINE_ProcF)
   

      CALL C_F_POINTER(DataF%nMPISides_YOUR_Proc, nMPISides_YOUR_ProcF,[nNbProcs])
      SDEALLOCATE(nMPISides_YOUR_Proc)
      ALLOCATE(nMPISides_YOUR_Proc(1:nNbProcs),SOURCE = nMPISides_YOUR_ProcF)
      
      ! ALLOCATE(nMPISides_MINE_Proc(1:nNbProcs),nMPISides_YOUR_Proc(1:nNbProcs))

    ! ALLOCATE(offsetMPISides_YOUR(0:nNbProcs),offsetMPISides_MINE(0:nNbProcs))
      CALL C_F_POINTER(DataF%offsetMPISides_YOUR, offsetMPISides_YOURF,[nNbProcs+1])
      SDEALLOCATE(offsetMPISides_YOUR)
      ALLOCATE(offsetMPISides_YOUR(0:nNbProcs), SOURCE = offsetMPISides_YOURF)

    ! print *,  "shape(offsetMPISides_MINE) =",shape(offsetMPISides_MINE)
      CALL C_F_POINTER(DataF%offsetMPISides_MINE, offsetMPISides_MINEF,[nNbProcs+1])
      SDEALLOCATE(offsetMPISides_MINE)
      ALLOCATE(offsetMPISides_MINE(0:nNbProcs),SOURCE = offsetMPISides_MINEF)


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
#if PARABOLIC
      SDEALLOCATE(MPIRequest_Lifting)
      ALLOCATE(MPIRequest_Lifting(nNbProcs,3,2))
      MPIRequest_Lifting = MPI_REQUEST_NULL
#endif /*PARABOLIC*/
      IF (nNbProcs .EQ. 0) nNbProcs =1;
    ENDIF
    
    nElems=DataF%nElems
    nSides=DataF%nSides

    PP_N=PP





    SDEALLOCATE(Ut_new)
    ALLOCATE(Ut_new(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(Metrics_fTilde_new)
    ALLOCATE(Metrics_fTilde_new(3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(Metrics_gTilde_new)
    ALLOCATE(Metrics_gTilde_new(3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(Metrics_hTilde_new)
    ALLOCATE(Metrics_hTilde_new(3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(dXGL_N_new)
    ALLOCATE(dXGL_N_new      (3,3,0:PP_N,0:PP_N,0:PP_N,nElems))

    SDEALLOCATE(sJ_new)
    ALLOCATE(            sJ_new(  0:PP_N,0:PP_N,0:PP_N,nElems))
    
    SDEALLOCATE(DetJac_Ref_new)
    NGeoRef=3*NGeo ! build jacobian at higher degree
    ALLOCATE(    DetJac_Ref_new(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

! surface data
    SDEALLOCATE(Face_xGP_new)
    ALLOCATE(      Face_xGP_new(3,0:PP_N,0:PP_N,1:nSides))
    ! Face_xGP=0;
    SDEALLOCATE(NormVec_new)
    ALLOCATE(       NormVec_new(3,0:PP_N,0:PP_N,1:nSides))
    
    SDEALLOCATE(TangVec1_new)
    ALLOCATE(      TangVec1_new(3,0:PP_N,0:PP_N,1:nSides))
    
    SDEALLOCATE(TangVec2_new)
    
    ALLOCATE(      TangVec2_new(3,0:PP_N,0:PP_N,1:nSides))
    SDEALLOCATE(SurfElem_new)
    ALLOCATE(      SurfElem_new(  0:PP_N,0:PP_N,1:nSides))


    CALL RecalculateParameters(DataF)
        !From DG Vars ????
        nDOFElem=(PP_N+1)**3
        nTotalU=PP_nVar*nDOFElem*nElems
        nTotal_face=(PP_N+1)*(PP_N+1)
        nTotal_vol=nTotal_face*(PP_N+1)
        nTotal_IP=nTotal_vol*nElems
        nTotalU=PP_nVar*nTotal_vol*nElems

        DEALLOCATE(U_master)
        DEALLOCATE(U_SLAVE)
        DEALLOCATE(Flux_master)
        DEALLOCATE(Flux_SLAVE)
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

        IF (ALLOCATED(gradPx_master))  THEN 
          DEALLOCATE(gradPx_master); 
          ALLOCATE(gradPx_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ENDIF

        IF (ALLOCATED(gradPy_master))  THEN 
          DEALLOCATE(gradPy_master); 
          ALLOCATE(gradPy_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ENDIF

        IF (ALLOCATED(gradPz_master))  THEN 
          DEALLOCATE(gradPz_master); 
          ALLOCATE(gradPz_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ENDIF
        
        IF (ALLOCATED(gradPx))  THEN 
          DEALLOCATE(gradPx); 
          ALLOCATE(gradPx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
        ENDIF
        
        IF (ALLOCATED(gradPy))  THEN 
          DEALLOCATE(gradPy); 
          ALLOCATE(gradPy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
        ENDIF
        
        IF (ALLOCATED(gradPz))  THEN 
          DEALLOCATE(gradPz); 
          ALLOCATE(gradPz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
        ENDIF

        IF (ALLOCATED(FluxX))  THEN 
          DEALLOCATE(FluxX); 
          ALLOCATE(FluxX        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ENDIF
        
        IF (ALLOCATED(FluxY))  THEN 
          DEALLOCATE(FluxY); 
          ALLOCATE(FluxY        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ENDIF
        
        IF (ALLOCATED(FluxZ))  THEN 
          DEALLOCATE(FluxZ); 
          ALLOCATE(FluxZ        (PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ENDIF

#endif /* PARABOLIC */

        IF (ALLOCATED(dtElem)) THEN  
          DEALLOCATE(dtElem); 
          ALLOCATE(dtElem(nElems)); 
        ENDIF
        
        
        IF (ALLOCATED(ElemVol)) THEN  
          DEALLOCATE(ElemVol); 
          ALLOCATE(ElemVol(nElems)); 
        ENDIF
        

        DEALLOCATE(AnalyzeSide)
        ALLOCATE(AnalyzeSide(1:nSides))
        AnalyzeSide=0;
        ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))

        ! ALLOCATE(Flux(PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ALLOCATE(Flux_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
        ALLOCATE(Flux_slave(PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))

!   ALLOCATE(ElementToCalc(DataF%nElems))
  ALLOCATE(ElementToCalc(DataF%nElems))
  do j=1,DataF%nElems
    ElementToCalc(j)=j;
  enddo

  ALLOCATE(Elem_xGP_New(3,0:PP,0:PP,0:PP,DataF%nElems))
  ALLOCATE(U_New(nVar,0:PP,0:PP,0:PP,DataF%nElems))
  iElem=0;
  !  DO iElem=1,DataF%nElems
   DO 
   iElem=iElem+1;
   if (iElem .GT. DataF%nElems) EXIT
          i=1;
           Ie= ChangeElem(i,iElem);
        !    print *, "iElem = ", iElem
           IF (Ie .LT. 0) THEN
            !This is refine and this and next 7 elements [iElem: iElem+7] number negative and 
            ! contains the number of child element 
            ! print *, "Ie = ", Ie
        !   IF ((Elem_xGP(1,PP,0,0,-iE) - Elem_xGP(1,0,0,0,-iE)) .LT. 0.015*2.) THEN

         !  Print *, "Error!!!", "COUNT =",COUNT
         !  Print *, "ElemToRefineAndCoarse =" ,ElemToRefineAndCoarse(-iE)
         !  Print *, "ChangeElem(:,iElem);=", ChangeElem(:,iElem);
         !  CALL EXIT()
         !  ENDIF
            CALL InterpolateCoarseRefine(U_new(:,:,:,:,iElem:iElem+7), &
                                        U(:,:,:,:,-Ie:-Ie),&
                                        Elem_xGP_New(:,:,:,:,iElem:iElem+7),&
                                        Elem_xGP(:,:,:,:,-Ie:-Ie))
            !It is also possible to use (/1,5,7,8/) instead of iElem:iElem+7
            iElem=iElem+7;
           ELSE IF (ChangeElem(2,iElem) .GT. 0) THEN
           !This is COARSE. Array ChangeElem(:,iElem) Contains 
           ! 8 Element which must be COARSED to the new number iElem
            CALL InterpolateCoarseRefine(U_new(:,:,:,:,iElem:iElem), &
                                        U(:,:,:,:,ChangeElem(:,iElem)),&
                                        Elem_xGP_New(:,:,:,:,iElem:iElem),&
                                        Elem_xGP(:,:,:,:,ChangeElem(:,iElem)))

         ! else IF (ChangeElem(2,iElem) .EQ. -1) THEN
          !     ! Reserve for creating the array for calculating Sides 
          !     ! Neighbours of refined Elements
          !     Elem_xGP_New(:,:,:,:,iElem)= Elem_xGP(:,:,:,:,Ie)
          !     U_New(:,:,:,:,iElem)= U(:,:,:,:,Ie)

          !     Ut_new(:,:,:,:,iElem) = Ut(:,:,:,:,Ie)
          !     Metrics_fTilde_new(:,:,:,:,iElem) = Metrics_fTilde(:,:,:,:,Ie)
          !     Metrics_gTilde_new(:,:,:,:,iElem) = Metrics_gTilde(:,:,:,:,Ie)
          !     Metrics_hTilde_new(:,:,:,:,iElem) = Metrics_hTilde(:,:,:,:,Ie)
          !     dXGL_N_new(:,:,:,:,:,iElem)=dXGL_N(:,:,:,:,:,Ie)

          !     sJ_new(:,:,:,iElem) = sJ(:,:,:,Ie)
          !     DetJac_Ref_new(:,:,:,:,iElem) = DetJac_Ref(:,:,:,:,Ie)
          !     ! ElemToSide(1,:, iElem) - Sides
    
          !     DO i=1,6
          !       IF (ChangeElem(i+2,iElem) .GT. nSides) THEN 
          !         PRINT *, "ERROR: Side number > nSides"
          !         CALL EXIT()
          !       ENDIF
          !       Face_xGP_new(:,:,:,ElemToSide(1,i,iElem)) = Face_xGP(:,:,:,ChangeElem(i+2,iElem))
          !       NormVec_new(:,:,:,ElemToSide(1,i,iElem)) = NormVec(:,:,:,ChangeElem(i+2,iElem))
          !       TangVec1_new(:,:,:,ElemToSide(1,i,iElem)) = TangVec1(:,:,:,ChangeElem(i+2,iElem))
          !       TangVec2_new(:,:,:,ElemToSide(1,i,iElem)) = TangVec2(:,:,:,ChangeElem(i+2,iElem))
          !       SurfElem_new(:,:,ElemToSide(1,i,iElem)) = SurfElem(:,:,ChangeElem(i+2,iElem))

          !       ! Print *, ",ElemToSide(1,i,iElem) = ", ElemToSide(1,i,iElem), "ChangeElem(i+2,iElem) = ",ChangeElem(i+2,iElem)
          !     ENDDO
          !     ! ElementToCalc(iElem)=0;
            ELSE
              IF (iE .LE. 0) THEN
                print *, "Error, iE = 0!, iElem = ", ielem
                ! print *, Count
                CALL EXIT()
              ENDIF
            ! This is simple case of renumeration of Elements
              ! IF (myrank .EQ. 1) THEN
              !   PRINT *, "Renumeration new iElem = ", iElem, "Old iElem = ", iE
              ! ENDIF
               Elem_xGP_New(:,:,:,:,iElem)= Elem_xGP(:,:,:,:,Ie)
               U_New(:,:,:,:,iElem)= U(:,:,:,:,Ie)
            ENDIF
  
   END DO
  CALL MOVE_ALLOC(Elem_xGP_New, Elem_xGP)
  CALL MOVE_ALLOC(U_New, U)




    CALL MOVE_ALLOC(Ut_new, Ut)

  CALL MOVE_ALLOC(Metrics_fTilde_new,Metrics_fTilde)
  CALL MOVE_ALLOC(Metrics_gTilde_new,Metrics_gTilde)
  CALL MOVE_ALLOC(Metrics_hTilde_new,Metrics_hTilde)

  CALL MOVE_ALLOC(dXGL_N_new,dXGL_N)
  CALL MOVE_ALLOC(sJ_new,sJ)
  CALL MOVE_ALLOC(DetJac_Ref_new,DetJac_Ref)

  CALL MOVE_ALLOC(Face_xGP_new,Face_xGP)
  
  CALL MOVE_ALLOC(NormVec_new,NormVec)
  CALL MOVE_ALLOC(TangVec1_new,TangVec1)
    CALL MOVE_ALLOC(TangVec2_new,TangVec2)
      CALL MOVE_ALLOC(SurfElem_new,SurfElem)

            ! PRint *,  ElementToCalc(:);
            ! CALL EXIT()


  ! CALL EXIT()

    !IF (PRESENT(ElemToRefineAndCoarse) ) THEN
      ! IF (myrank .EQ. 0) THEN 
        ! CALL CalcMetrics(ElementToCalc)
      ! ENDIF
        ! CALL EXIT()
    !ELSE 
      CALL CalcMetrics(ElementToCalc)

  Deallocate(ElementToCalc)



  call free_data_memory(DataPtr)
 
  NULLIFY(nMPISides_YOUR_ProcF)
  NULLIFY(nMPISides_MINE_ProcF)
  NULLIFY(offsetMPISides_MINEF)
  NULLIFY(offsetMPISides_YOURF)

  NULLIFY(nNbProcF)
  NULLIFY(nMPISides_ProcF)
  NULLIFY(ChangeSide)
  NULLIFY(ChangeElem)
  NULLIFY(EtSF) 
  NULLIFY(StEF)
  NULLIFY(MInfo)
  NULLIFY(MType)
  NULLIFY(nBCsF)
  
  NULLIFY(DataF)
    !  print *," ERRRORRRR !!!!!!!!!!!!"
!   IF (MPIroot) THEN
! ! 
!   ! ELSE 


END SUBROUTINE RunAMR


  SUBROUTINE RecalculateParameters(DataF)
        USE MOD_AMR_Vars,            ONLY: P4EST_FORTRAN_DATA
        USE MOD_Mesh_Vars
        USE MOD_MPI_Vars,             ONLY:  offsetElemMPI
        ! USE MOD_AMR_vars,            ONLY: P4EST_PTR, p4est_mpi_data
        USE MOD_Globals,             ONLY:  myrank
        USE MOD_P4est,               ONLY: p4estGetMPIData

        IMPLICIT NONE
        TYPE(p4est_fortran_data), POINTER :: DataF
        ! TYPE(C_PTR) :: DataPtr;
        INTEGER     ::firstMasterSide, lastMasterSide, nMortarMPISide

        ! DataPtr = p4estGetMPIData(p4est_ptr)
        ! CALL C_F_POINTER(DataPtr, DATAF)
  

        nBCSides            =   DataF%nBCSides
        nElems              =   DataF%nElems
        nGlobalElems        =   DataF%nGlobalElems
        nSides              =   DataF%nSides
        nMortarInnerSides   =   DataF%nMortarInnerSides
        nInnerSides         =   DataF%nInnerSides
        ! Must be set later, for parallel version
        nMortarMPISide         = DataF%nMortarMPISides
        nMPISides_MINE         = DataF%nMPISides_MINE
        nMPISides_YOUR         = DataF%nMPISides_YOUR


        firstBCSide          = 1
        firstMortarInnerSide = firstBCSide         +nBCSides
        firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
        firstMPISide_MINE    = firstInnerSide      +nInnerSides
        firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
        firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR
        ! print *, "nMPISides_MINE = ", nMPISides_MINE
        ! print *, "nMPISides_YOUR = ", nMPISides_YOUR
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
        
        offsetElem  = offsetElemMPI(myrank)
        
    END SUBROUTINE RecalculateParameters


SUBROUTINE InterpolateCoarseRefine(Unew, Uold,Elem_xGPnew,Elem_xGPold)
    ! Enumerate the cells number as in P4est
    USE MOD_Interpolation_Vars, ONLY: NodeType,NodeTypeVISU
    USE MOD_Interpolation,      ONLY: GetVandermonde
    USE MOD_ChangeBasis,        ONLY: ChangeBasis3D

    IMPLICIT NONE
    REAL, INTENT(INOUT)          :: Unew(:,:,:,:,:),Elem_xGPnew(:,:,:,:,:)
    REAL, INTENT(IN)             :: Uold(:,:,:,:,:),Elem_xGPold(:,:,:,:,:)
    INTEGER                      :: SizeNew, SizeOld
    INTEGER                      :: i,PP_N, nVar_local
    REAL,ALLOCATABLE                :: Unew_big(:,:,:,:,:),Elem_xGP_big(:,:,:,:,:)
    REAL,ALLOCATABLE                :: Vdm_fromSmalltoBig(:,:),Vdm_FromNVisu_toNodeType(:,:)
    ! REAL,ALLOCATABLE                :: Unew_big(:,:,:,:,:) , Elem_xGP_big(:,:,:,:,:)
    REAL,ALLOCATABLE                :: Vdm_fromBigtoSmall(:,:),Vdm_FromNodeType_toNVisu(:,:)

    sizenew=size(Unew(1,1,1,1,:))
    sizeold=size(Uold(1,1,1,1,:))


    nVar_local=size(Unew(:,1,1,1,1))
    PP_N=size(Unew(1,:,1,1,1))-1
    
    !     ALLOCATE(Vdm_fromSmalltoBig(0:2*PP_N,0:PP_N))
    !     ALLOCATE(Vdm_FromNVisu_toNodeType(0:PP_N,0:PP_N))

    IF ((sizenew .EQ. 1 .AND. sizeold .EQ. 8) .OR. (sizenew .EQ.  8 .AND. sizeold .EQ. 1)) THEN
      IF (sizenew > sizeold) THEN
      !This is refine


        ALLOCATE(Vdm_fromSmalltoBig(0:2*PP_N,0:PP_N))
        ALLOCATE(Vdm_FromNVisu_toNodeType(0:PP_N,0:PP_N))

    !   print *, "Size!!!!!!!11 U old ", nVar_local
        CALL GetVandermonde(PP_N, NodeType, 2*PP_N,      NodeTypeVISU, Vdm_fromSmallToBig)
        CALL GetVandermonde(PP_N, NodeTypeVISU, PP_N,      NodeType, Vdm_FromNVisu_toNodeType)

        ALLOCATE(Unew_big(PP_nVar,0:2*PP_N,0:2*PP_N,0:2*PP_N,1))
        Unew_big=0
        CALL ChangeBasis3D(PP_nVar,PP_N,2*PP_N,Vdm_fromSmallToBig,Uold(:,:,:,:,1),Unew_big(:,:,:,:,1))

        !Map the data from Big Cell into the 8 small cells
        i=2*PP_N
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,0:PP_N,0:PP_N,0:PP_N,1),Unew(:,:,:,:,1))
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,PP_N:i,0:PP_N,0:PP_N,1),Unew(:,:,:,:,2))
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,PP_N:i,PP_N:i,0:PP_N,1),Unew(:,:,:,:,4))
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,0:PP_N,PP_N:i,0:PP_N,1),Unew(:,:,:,:,3))
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,0:PP_N,0:PP_N,PP_N:i,1),Unew(:,:,:,:,5))
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,PP_N:i,0:PP_N,PP_N:i,1),Unew(:,:,:,:,6))
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,PP_N:i,PP_N:i,PP_N:i,1),Unew(:,:,:,:,8))
        CALL ChangeBasis3D(PP_nVar,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Unew_big(:,0:PP_N,PP_N:i,PP_N:i,1),Unew(:,:,:,:,7))


        DEALLOCATE(Unew_big)
        ALLOCATE(Elem_xGP_big(size(Elem_xGPnew(:,1,1,1,1)),0:2*PP_N,0:2*PP_N,0:2*PP_N,1))
        !Add interpolation of this elements


        !CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
        CALL ChangeBasis3D(3,PP_N,2*PP_N,Vdm_fromSmallToBig,Elem_xGPold(:,:,:,:,1),Elem_xGP_big(:,:,:,:,1))



        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,0:PP_N,0:PP_N,0:PP_N,1),Elem_xGPnew(1:3,:,:,:,1))! Octant%Elems(1)%ElemID))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,PP_N:i,0:PP_N,0:PP_N,1),Elem_xGPnew(1:3,:,:,:,2))! Octant%Elems(2)%ElemID))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,PP_N:i,PP_N:i,0:PP_N,1),Elem_xGPnew(1:3,:,:,:,4))! Octant%Elems(3)%ElemID))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,0:PP_N,PP_N:i,0:PP_N,1),Elem_xGPnew(1:3,:,:,:,3))! Octant%Elems(4)%ElemID))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,0:PP_N,0:PP_N,PP_N:i,1),Elem_xGPnew(1:3,:,:,:,5))! Octant%Elems(5)%ElemID))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,PP_N:i,0:PP_N,PP_N:i,1),Elem_xGPnew(1:3,:,:,:,6))! Octant%Elems(6)%ElemID))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,PP_N:i,PP_N:i,PP_N:i,1),Elem_xGPnew(1:3,:,:,:,8))! Octant%Elems(7)%ElemID))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNVisu_toNodeType,Elem_xGP_big(1:3,0:PP_N,PP_N:i,PP_N:i,1),Elem_xGPnew(1:3,:,:,:,7))! Octant%Elems(8)%ElemID))

        !~
        DEALLOCATE(Elem_xGP_big)
        DEALLOCATE(Vdm_fromSmalltoBig)
        DEALLOCATE(Vdm_FromNVisu_toNodeType)


      ELSE
      !This is Coarse
  
        ALLOCATE(Vdm_fromBigtoSmall(0:PP_N,0:2*PP_N))
        ALLOCATE(Vdm_FromNodeType_toNVisu(0:PP_N,0:PP_N))
         
        CALL GetVandermonde(2*PP_N, NodeTypeVISU, PP_N, NodeType,      Vdm_fromBigToSmall)
        CALL GetVandermonde(PP_N,      NodeType,PP_N, NodeTypeVISU,  Vdm_FromNodeType_toNVisu)


        ALLOCATE(Unew_big(PP_nVar,0:2*PP_N,0:2*PP_N,0:2*PP_N,1))
        Unew_big=0
        i=2*PP_N  

      

        CALL ChangeBasis3D(PP_nVar ,PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,1), Unew_big(:,0:PP_N,0:PP_N,0:PP_N,1))
        CALL ChangeBasis3D(PP_nVar, PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,2), Unew_big(:,PP_N:i,0:PP_N,0:PP_N,1))
        CALL ChangeBasis3D(PP_nVar, PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,4), Unew_big(:,PP_N:i,PP_N:i,0:PP_N,1))
        CALL ChangeBasis3D(PP_nVar, PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,3), Unew_big(:,0:PP_N,PP_N:i,0:PP_N,1))
        CALL ChangeBasis3D(PP_nVar, PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,5), Unew_big(:,0:PP_N,0:PP_N,PP_N:i,1))
        CALL ChangeBasis3D(PP_nVar, PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,6), Unew_big(:,PP_N:i,0:PP_N,PP_N:i,1))
        CALL ChangeBasis3D(PP_nVar, PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,8), Unew_big(:,PP_N:i,PP_N:i,PP_N:i,1))
        CALL ChangeBasis3D(PP_nVar ,PP_N, PP_N, Vdm_FromNodeType_toNVisu, Uold(:,:,:,:,7), Unew_big(:,0:PP_N,PP_N:i,PP_N:i,1))
        
        CALL ChangeBasis3D(PP_nVar,2*PP_N,PP_N,Vdm_fromBigToSmall,Unew_big(:,:,:,:,1),Unew(:,:,:,:,1))
        
        DEALLOCATE(Unew_big)
        ALLOCATE(Elem_xGP_big(size(Elem_xGPnew(:,1,1,1,1)),0:2*PP_N,0:2*PP_N,0:2*PP_N,1))
        !Add interpolation of this elements
        
     
        
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,1),Elem_xGP_big(1:3,0:PP_N,0:PP_N,0:PP_N,1))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,2),Elem_xGP_big(1:3,PP_N:i,0:PP_N,0:PP_N,1))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,4),Elem_xGP_big(1:3,PP_N:i,PP_N:i,0:PP_N,1))  
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,3),Elem_xGP_big(1:3,0:PP_N,PP_N:i,0:PP_N,1))  
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,5),Elem_xGP_big(1:3,0:PP_N,0:PP_N,PP_N:i,1))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,6),Elem_xGP_big(1:3,PP_N:i,0:PP_N,PP_N:i,1))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,8),Elem_xGP_big(1:3,PP_N:i,PP_N:i,PP_N:i,1))
        CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_FromNodeType_toNVisu,Elem_xGPold(1:3,:,:,:,7),Elem_xGP_big(1:3,0:PP_N,PP_N:i,PP_N:i,1))
        
        CALL ChangeBasis3D(3,2*PP_N,PP_N,Vdm_fromBigtoSmall,Elem_xGP_big(:,:,:,:,1),Elem_xGPnew(:,:,:,:,1))
        !~     
        DEALLOCATE(Elem_xGP_big)
        DEALLOCATE(Vdm_fromBigtoSmall)
        DEALLOCATE(Vdm_FromNodeType_toNVisu)

      ENDIF  !if (sizenew > sizeold) THEN
    ELSE 
      !There is an Error !!!
      print*, "Error in =InterpolateCoarseRefine= the number of cells not 1 to 8 or 8 to 1"
    ENDIF

  

    END SUBROUTINE InterpolateCoarseRefine



!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE LoadBalancingAMR()
  ! MODULES
  USE MOD_Globals
  USE MOD_AMR_Vars
  USE MOD_P4EST
  USE MOD_DG_Vars,            ONLY: U
  USE MOD_Mesh_Vars,          ONLY: Elem_xGP, nElems
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  
  REAL,ALLOCATABLE, TARGET :: Elem_xGP_New(:,:,:,:,:), U_New(:,:,:,:,:)
  ! REAL,ALLOCATABLE, TARGET :: ExGP_New(:,:), ExGP_old(:,:)
  INTEGER :: PP, nVar
  ! REAL,ALLOCATABLE :: Elem_xGP(:,:,:,:,:)
  ! REAL, POINTER :: Elem_xGPP(:,:,:,:,:)
  ! REAL, POINTER :: UP(:,:,:,:,:)
  !============================================================================================================================
  ! Deallocate global variables, needs to go somewhere else later
  TYPE(p4est_balance_datav2), TARGET :: BalanceData;
  
  IF (.NOT. UseAMR) THEN
    RETURN;
  ENDIF
  
  
  BalanceData%DataSize = sizeof(U(:,:,:,:,1))
  BalanceData%GPSize = sizeof(Elem_xGP(:,:,:,:,1))
  PP = size(U(1,:,0,0,1))-1
  nVar = size(U(:,0,0,0,1))
  ! PRINT *, "BalanceData%DataSize =", BalanceData%DataSize
  ! PRINT *, "BalanceData%GPSize =", BalanceData%GPSize
  ! PRINT *, "PP =", PP
  ! PRINT *, "nVar =", nVar
  BalanceData%Uold_Ptr = C_LOC(U)
  BalanceData%ElemxGPold_Ptr = C_LOC(Elem_xGP)
  
  ! BalanceData%ElemxGPold_Ptr = C_LOC(ExGP_old)
  
  
  
  CALL p4est_loadbalancing_init(P4EST_PTR, C_LOC(BalanceData))
  ! PRINT *, "BalanceData%nElems =", BalanceData%nElems

  ! CALL p4est_loadbalancing(P4EST_PTR, C_LOC(BalanceData))
  ALLOCATE(U_New(PP_nVar,0:PP,0:PP,0:PP,BalanceData%nElems))
  BalanceData%Unew_Ptr = C_LOC(U_New)
  ALLOCATE(Elem_xGP_New(3,0:PP,0:PP,0:PP,BalanceData%nElems))
  ! Elem_xGP_New = 0.
  ! U_new = 0.
  BalanceData%ElemxGPnew_Ptr = C_LOC(Elem_xGP_New)
  ! ALLOCATE(ExGP_New(1000,1:BalanceData%nElems))
  ! BalanceData%ElemxGPnew_Ptr = C_LOC(ExGP_New)
  CALL p4est_loadbalancing_go(P4EST_PTR, C_LOC(BalanceData))

  IF (Myrank .EQ. 1) THEN
    ! PRINT *, "1-177", Elem_xGP_New(:,:,1,1,177)
  ENDIF
  CALL MOVE_ALLOC(Elem_xGP_New, Elem_xGP)
  CALL MOVE_ALLOC(U_New, U)
  CALL p4est_ResetElementNumber(P4EST_PTR)
  !-- 
  ! ! print *, "BalanceData%nElemsNew = ",BalanceData%nElems
  ! nElemsNew=BalanceData%nElems;
  ! CALL C_F_POINTER(BalanceData%DataSetU, U_New,[nVar,PP+1,PP+1,PP+1,nElemsNew])
  ! CALL C_F_POINTER(BalanceData%DataSetElem_xGP, Elem_xGP_New,[3,PP+1,PP+1,PP+1,nElemsNew])
  ! SDEALLOCATE(Elem_xGP)
  ! SDEALLOCATE(U)
  ! ALLOCATE(Elem_xGP(1:3,0:PP,0:PP,0:PP,1:nElemsNew))
  ! Elem_xGP(1:3,0:PP,0:PP,0:PP,1:nElemsNew) = Elem_xGP_New(1:3,1:PP+1,1:PP+1,1:PP+1,1:nElemsNew)
  
  ! ALLOCATE(U(1:nVar,0:PP,0:PP,0:PP,1:nElemsNew))
  ! U(1:nVar,0:PP,0:PP,0:PP,1:nElemsNew) = U_new(1:nVar,1:PP+1,1:PP+1,1:PP+1,1:nElemsNew)
  
  ! CALL free_balance_memory(C_LOC(BalanceData))
  
  ! DEALLOCATE(U_new)
  ! NULLIFY(U_New)
  ! NULLIFY(Elem_xGP_New)
  CALL RunAMR()
  ! PRINT *,"YAHOOOO!!!! "
  ! CALL EXIT()
 
END SUBROUTINE LoadBalancingAMR
      

!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE LoadBalancingAMRold()
! MODULES
USE MOD_AMR_Vars
USE MOD_P4EST
USE MOD_DG_Vars,            ONLY: U
USE MOD_Mesh_Vars,          ONLY: Elem_xGP
USE, INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE

REAL,POINTER :: Elem_xGP_New(:,:,:,:,:), U_New(:,:,:,:,:)
INTEGER :: DataSize, nElemsNew, PP, nVar
! REAL,ALLOCATABLE :: Elem_xGP(:,:,:,:,:)
! REAL, POINTER :: Elem_xGPP(:,:,:,:,:)
! REAL, POINTER :: UP(:,:,:,:,:)
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
TYPE(p4est_balance_data), TARGET :: BalanceData;

IF (.NOT. UseAMR) THEN
  RETURN;
ENDIF
! UP=>U
! Elem_xGPP=>Elem_xGP
BalanceData%nVar = size(U(:,0,0,0,1))
BalanceData%PP = size(U(1,:,0,0,1))-1
nVar = BalanceData%nVar
PP = BalanceData%PP
BalanceData%nElems = size(U(1,0,0,0,:))
BalanceData%DataSize=INT(sizeof(U(:,:,:,:,1)) + sizeof(Elem_xGP(:,:,:,:,1)))


BalanceData%DataSetU = C_LOC(U)
BalanceData%DataSetElem_xGP = C_LOC(Elem_xGP)

CALL p4est_loadbalancing(P4EST_PTR, C_LOC(BalanceData))

! print *, "BalanceData%nElemsNew = ",BalanceData%nElems
nElemsNew=BalanceData%nElems;
CALL C_F_POINTER(BalanceData%DataSetU, U_New,[nVar,PP+1,PP+1,PP+1,nElemsNew])
CALL C_F_POINTER(BalanceData%DataSetElem_xGP, Elem_xGP_New,[3,PP+1,PP+1,PP+1,nElemsNew])
SDEALLOCATE(Elem_xGP)
SDEALLOCATE(U)
ALLOCATE(Elem_xGP(1:3,0:PP,0:PP,0:PP,1:nElemsNew))
Elem_xGP(1:3,0:PP,0:PP,0:PP,1:nElemsNew) = Elem_xGP_New(1:3,1:PP+1,1:PP+1,1:PP+1,1:nElemsNew)

ALLOCATE(U(1:nVar,0:PP,0:PP,0:PP,1:nElemsNew))
U(1:nVar,0:PP,0:PP,0:PP,1:nElemsNew) = U_new(1:nVar,1:PP+1,1:PP+1,1:PP+1,1:nElemsNew)

CALL free_balance_memory(C_LOC(BalanceData))

NULLIFY(U_New)
NULLIFY(Elem_xGP_New)
CALL RunAMR()
! PRINT *,"YAHOOOO!!!! "
! CALL EXIT()


END SUBROUTINE LoadBalancingAMRold



!============================================================================================================================
!> Save mesh to HDF5 file
!============================================================================================================================
SUBROUTINE SaveMesh(FileString)
  ! MODULES
  USE MOD_AMR_Vars
  USE MOD_P4EST
  USE MOD_Mesh_vars,           ONLY: nElems, nGlobalElems, offsetElem, BoundaryName, BoundaryType, nBCs, Elem_xGP
  ! DEbug
  USE MOD_Mesh_vars,           ONLY: ElemToSide, SideToElem, nSides,MortarType, MortarInfo
  !EndDebug
  USE MOD_Globals,                only: myrank,nProcessors, MPIRoot
  USE MOD_IO_HDF5
  USE MOD_HDF5_Output,            only: WriteHeader, WriteAttribute, WriteArray, GatheredWriteArray
  USE MOD_DG_Vars,            ONLY: U
  USE MOD_Interpolation_Vars, ONLY: NodeType,NodeTypeVISU, NodeTypeGL
  USE MOD_Interpolation,      ONLY: GetVandermonde
  USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
  ! USE MOD_Mesh_Vars,           ONLY:Vdm_GLN_N
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
    CHARACTER(LEN=*),INTENT(IN)  :: FileString !< (IN) mesh filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  TYPE(C_PTR) :: DataPtr;
  TYPE(p4est_save_data), POINTER :: DataF
  INTEGER, POINTER :: ElInfoF(:,:), SiInfoF(:,:)  ! ElemInfo and SideInfo to pass the data
  INTEGER, POINTER :: OffSideMPIF(:), OffSideArrIndMPIF(:)  ! OffsetSideMPI and OffsetSideArrIndexMPI to pass the data
  INTEGER, ALLOCATABLE :: OffsetSideMPI(:), OffsetSideArrIndexMPI(:), ElemInfoW(:,:)
  INTEGER          :: iElem, nIndexSide, mpisize, FirstElemInd, LastElemInd, nLocalIndSides, nGlobalIndSides, OffsetIndSides
  CHARACTER(LEN=255)             :: FileName,TypeString
  INTEGER(HID_T)                 :: DSet_ID,FileSpace,HDF5DataType
  INTEGER(HSIZE_T)               :: Dimsf(5)
  INTEGER(HSIZE_T)               :: DimsM(2)
  CHARACTER(LEN=255), ALLOCATABLE:: BCNames(:)
  INTEGER, ALLOCATABLE           :: BCMapping(:),BCType(:,:)
  INTEGER                        :: offset = 0, PP, i,j,k, index
  REAL,ALLOCATABLE               :: NodeCoords(:,:)
  REAL,ALLOCATABLE               :: NodeCoordsTMP(:,:,:,:)
  REAL,ALLOCATABLE               :: Vdm_FromNodeType_toNVisu(:,:), Vdm_GLN_EQN(:,:)
  REAL,ALLOCATABLE               :: Vdm_N_GLN(    :   ,:)
  REAL,ALLOCATABLE                :: XGL_N(      :,  :,:,:)          ! mapping X(xi) P\in N
  ! REAL,ALLOCATABLE,TARGET :: NodeCoords(:,:,:,:,:) !< XYZ positions (equidistant,NGeo) of element interpolation points from meshfile
  ! CHARACTER(LEN=255)             :: MeshFile255
  !============================================================================================================================
  ! 


  
  IF (.NOT. UseAMR) RETURN;

  mpisize = nProcessors
  DATAPtr=SaveMeshP4(p4est_ptr)
  CALL C_F_POINTER(DataPtr, DataF)

  CALL C_F_POINTER(DataF%ElemInfo, ElInfoF,[6,nElems])
  !read local ElemInfo from data file
FirstElemInd=offsetElem+1
LastElemInd=offsetElem+nElems
ALLOCATE(ElemInfoW(6,nElems), SOURCE  = ElInfoF)
ElemInfoW = ElInfoF
  ! ALLOCATE(ElemInfo1(6,nElems))
  nIndexSide = ElInfoF(4,nElems)
  
  CALL C_F_POINTER(DataF%SideInfo, SiInfoF,[5,nIndexSide])
  
  CALL C_F_POINTER(DataF%OffsetSideMPI, OffSideMPIF,[mpisize+1])
  ALLOCATE(OffsetSideMPI(0:mpisize), SOURCE  = OffSideMPIF)

  CALL C_F_POINTER(DataF%OffsetSideArrIndexMPI, OffSideArrIndMPIF,[mpisize+1])
  ALLOCATE(OffsetSideArrIndexMPI(0:mpisize), SOURCE  = OffSideArrIndMPIF)
  

  ElInfoF(3,:) = ElInfoF(3,:) + OffsetSideArrIndexMPI(Myrank)
  ElInfoF(4,:) = ElInfoF(4,:) + OffsetSideArrIndexMPI(Myrank)
  ! IF(nVar_Avg.GT.0)THEN
  ! FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_TimeAvg',OutputTime))//'.h5'
! 1. Create Mesh file
  !Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
  PP = size(U(1,:,1,1,1))-1
  IF(MPIRoot)THEN
   ! Create file
   FileName=TRIM(FileString)//'.h5'
  !  print *, "FileName = ", TRIM(FileName), FileName, TRIM(FileName)
   TypeString = FileName
    CALL OpenDataFile(TRIM(FileString),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    ! CALL WriteHeader(TRIM(TypeString),File_ID)
    CALL WriteAttribute(File_ID,'Version',1,RealScalar=1.0)
    CALL WriteAttribute(File_ID,'Ngeo',1,IntScalar=PP)
    CALL WriteAttribute(File_ID,'nElems',1,IntScalar=nGlobalElems)
    CALL WriteAttribute(File_ID,'nSides',1,IntScalar=OffsetSideArrIndexMPI(mpisize))
    CALL WriteAttribute(File_ID,'nNodes',1,IntScalar=nGlobalElems*(PP+1)*(PP+1)*(PP+1))
    CALL WriteAttribute(File_ID,'nUniqueSides',1,IntScalar=OffsetSideMPI(mpisize))
    CALL WriteAttribute(File_ID,'nUniqueNodes',1,IntScalar=343)
    CALL WriteAttribute(File_ID,'nBCs',1,IntScalar=nBCs)
    
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
    DimsM=(/3, nGlobalElems*(PP+1)*(PP+1)*(PP+1)/)
    CALL H5SCREATE_SIMPLE_F(2, DimsM, FileSpace, iError)
    HDF5DataType=H5T_NATIVE_DOUBLE
    CALL H5DCREATE_F(File_ID,'NodeCoords', HDF5DataType, FileSpace, DSet_ID, iError)
    ! Close the filespace and the dataset
    CALL H5DCLOSE_F(Dset_id, iError)
    CALL H5SCLOSE_F(FileSpace, iError)  

    ! CALL WriteAttribute(File_ID,'nNodes',1,IntScalar=5)
  !   CALL GenerateFileSkeleton(TRIM(FileName),'TimeAvg',nVar_Avg,PP_N,VarNamesAvg,MeshFileName,OutputTime,FutureTime)
  !   CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  !   CALL WriteAttribute(File_ID,'AvgTime',1,RealScalar=dtAvg)


  ! ALLOCATE(BoundaryName(nBCs))
! ALLOCATE(BoundaryType(nBCs,3))
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

! SWRITE(UNIT_StdOut,'(132("."))')
! SWRITE(Unit_StdOut,'(A,A16,A20,A10,A10,A10)')'BOUNDARY CONDITIONS','|','Name','Type','State','Alpha'
! DO iBC=1,nBCs
!   SWRITE(*,'(A,A33,A20,I10,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:)
! END DO
! SWRITE(UNIT_StdOut,'(132("."))')
  CALL WriteArray('BCType',2,(/4, nBCs/),(/4, nBCs/),&
                                                (/0, 0/),collective = .FALSE.,IntArray =BCType)

  ! CALL ReadArray('BCNames',1,(/nBCs/),Offset,1,StrArray=BCNames)  ! Type is a dummy type only
  CALL WriteArray('BCNames',1,(/nBCs/),(/nBCs/),&
                                                (/0/),collective = .FALSE.,StrArray =BCNames)
                                                
  DEALLOCATE(BCNames,BCType)

    CALL CloseDataFile()
  END IF !  IF(MPIRoot)THEN
#if MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
! 2. Save File
  ! Reopen file and write ElemInfo
  ! CALL OpenDataFile(TRIM(FileString),create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
  ! ! IF(PRESENT(RealArray)) CALL WriteArray(DataSetName,rank,nValGlobal,nVal,&
  !                                             !  offset,collective,RealArray=RealArray)
  ! ! IF(PRESENT(IntArray)) 
  !  CALL WriteArray('ElemInfo',  2,                (/6,nGlobalElems/),&
  !                                               (/6,nElems/),&
  !                                               (/0, offsetElem/),&
  !                                              .TRUE.,IntArray =ElInfoF)
  ! CALL CloseDataFile()

DO iElem = 1, nElems
  ElInfoF(5,iElem) = (iElem + offsetElem - 1)*(PP+1)*(PP+1)*(PP+1)
  ElInfoF(6,iElem) = (iElem + offsetElem)*(PP+1)*(PP+1)*(PP+1)
ENDDO

IF (PP .EQ. 1) THEN
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

nLocalIndSides = DataF%nSidesArrIndex
nGlobalIndSides = OffsetSideArrIndexMPI(mpisize)
OffsetIndSides = OffsetSideArrIndexMPI(myrank)

CALL GatheredWriteArray(FileString,create=.FALSE.,&
                        DataSetName='SideInfo', rank=2,  &
                        nValGlobal=(/5,nGlobalIndSides/),&
                        nVal=      (/5,nLocalIndSides      /),&
                        offset=    (/0   ,OffsetIndSides  /),&
                        collective=.TRUE.,IntArray=SiInfoF)

ALLOCATE(NodeCoords(3,nElems*(PP+1)*(PP+1)*(PP+1)))
ALLOCATE(NodeCoordsTMP(3,0:PP,0:PP,0:PP))
ALLOCATE(Vdm_FromNodeType_toNVisu(0:PP,0:PP))
! ALLOCATE(Vdm_GLN_EQN(0:PP,0:PP))
! ALLOCATE(Vdm_N_GLN(0:PP,0:PP))
! ALLOCATE(XGL_N(3,  0:PP,0:PP,0:PP))          
    
    ! 1.a) NodeCoords: EQUI Ngeo to GLNgeo and GLN
! CALL GetVandermonde(    PP   ,NodeTypeGL , PP    ,NodeTypeVISU, Vdm_GLN_EQN , modal=.FALSE.)

! CALL GetVandermonde(    PP   ,NodeType , PP    ,NodeTypeGL, Vdm_N_GLN , modal=.FALSE.)

CALL GetVandermonde(PP, NodeType, PP,      NodeTypeVISU, Vdm_FromNodeType_toNVisu)
    

index = 1
! print *, "index = ", index
DO iElem=1,nElems
 
  CALL ChangeBasis3D(3,PP,PP,Vdm_FromNodeType_toNVisu,Elem_xGP(:,:,:,:,iElem),NodeCoordsTMP(:,:,:,:))! Octant%Elems(1)%ElemID))
  ! CALL ChangeBasis3D(3,PP,PP,Vdm_N_GLN,Elem_xGP(:,:,:,:,iElem),XGL_N)! Octant%Elems(1)%ElemID))
  ! CALL ChangeBasis3D(3,PP,PP,Vdm_GLN_EQN,XGL_N, NodeCoordsTMP(:,:,:,:))
  DO k=0,PP
    DO j=0,PP
      DO i=0,PP

        NodeCoords(:,index) = NodeCoordsTMP(:,i,j,k)
        index = index + 1
      ENDDO
    ENDDO
  ENDDO
ENDDO
CALL GatheredWriteArray(FileString,create=.FALSE.,&
                        DataSetName='NodeCoords', rank=2,  &
                        nValGlobal=(/3,nGlobalElems*(PP+1)*(PP+1)*(PP+1)/),&
                        nVal=      (/3,nElems*(PP+1)*(PP+1)*(PP+1)/),&
                        offset=    (/0   ,offSetElem*(PP+1)*(PP+1)*(PP+1)/),&
                        collective=.TRUE.,RealArray=NodeCoords)

SDEALLOCATE(NodeCoords)
SDEALLOCATE(NodeCoordsTMP)
SDEALLOCATE(Vdm_FromNodeType_toNVisu)

#if MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

  DEALLOCATE(OffsetSideMPI)
  DEALLOCATE(OffsetSideArrIndexMPI)
  CALL free_savemesh_memory(DATAPtr)
  NULLIFY(ElInfoF)
  NULLIFY(SiInfoF)
  NULLIFY(OffSideMPIF)
  NULLIFY(OffSideArrIndMPIF)
  ! CALL EXIT()
END SUBROUTINE SaveMesh

!============================================================================================================================
!> Deallocate AMR data.
!============================================================================================================================
SUBROUTINE FinalizeAMR()
! MODULES
USE MOD_AMR_Vars
USE MOD_P4EST
IMPLICIT NONE
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
  !INTEGER A=0
! returngit ..
IF (.NOT. UseAMR) THEN
  RETURN;
ENDIF
CALL p4est_destroy(P4EST_PTR);
CALL p4est_connectivity_destroy(CONNECTIVITY_PTR)
CALL p4est_finalize()
AMRInitIsDone = .FALSE.
END SUBROUTINE FinalizeAMR





END MODULE MOD_AMR
#endif
