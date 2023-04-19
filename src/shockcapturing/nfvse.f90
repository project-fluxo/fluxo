!==================================================================================================================================
! Copyright (c) 2020 - 2020 Andrés Rueda
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
! Use bar states by default, except if the equation has non-conservative terms
!TODO: define  barstates in another manner!!!
#if NONCONS
#define barStates 0
#else
#define barStates 1
#endif /*NONCONS*/
module MOD_NFVSE
#if SHOCK_NFVSE
  use MOD_NFVSE_Vars
  use MOD_PreProc
  implicit none
!==================================================================================================================================
!> Computes the spatial operator using Native Finite Volume Sub-Elements (NFVSE)
!> Based on: * Hennemann et al. (2020). "A provably entropy stable subcell shock capturing approach for high order split form 
!>                                       DG for the compressible Euler Equations"
!>           * Rueda-Ramírez et al. (2021) "An Entropy Stable Nodal Discontinuous Galerkin Method for the resistive MHD Equations. 
!>                                          Part II: Subcell Finite Volume Shock Capturing"
!> Attention 1: The current strategy is to perform the operations on "inner faces". This might not be the most optimal choice, 
!>              as the inner element data has to be copied each time to improve performance. This can be improved, but the Riemann  
!>              solver routines have to be adapted
!> Attention 2: The inner faces do not have the same indexing as the inter-element faces (sides). The reshape routines, and the 
!>              definition of TanDirs1 and TanDirs2, account for that...
!==================================================================================================================================
!
! Example for N=4 in 1D
! ---------------------
!
! O Gauss-Lobatto point
! · FV boundary
! 0 Both a Gauss-Lobatto point and a FV boundary
!
! Three different kinds of indexes
!
!(1) DG DOF idx:     0    1      2      3    4
!                    |    |      |      |    |
!                    0-·--O--·---O---·--O--·-0
!                    | |     |       |     | |
!(2) FV BC idx:     -1 0     1       2     3 4 (inner faces)
!                    └─┴─────┴───────┴─────┴─┘
!(3) FV subcell idx:  0   1      2      3   4
!
!   * (1) matches (3) for consistency
!===================================================================================================================================
  private
  public :: VolInt_NFVSE, InitNFVSE, FinalizeNFVSE
  public :: DefineParametersNFVSE
  public :: InitNFVSEAfterAdaptation1, InitNFVSEAfterAdaptation2
  public :: CalcBlendingCoefficient
  
contains
!===================================================================================================================================
!> Defines parameters for NFVSE (to be called in shockcapturing.f90)
!===================================================================================================================================
  subroutine DefineParametersNFVSE()
    use MOD_ReadInTools,  only: prms
#if NFVSE_CORR
    use MOD_IDP        ,  only: DefineParameters_IDP
#endif /*NFVSE_CORR*/
    implicit none
    
    call prms%CreateIntOption     (   "SubFVMethod",  " Specifies subcell Finite-Volume method to be used:\n"//&
                                              "   1: 1st order FV\n"//&
                                              "   2: TVD method (not ES)\n"//&
                                              "   3: TVD-ES method (entropy fix)\n"//&
                                              "   4: TVD-ES method (à la Fjordhom)\n"&
                                             ,"1")
                                             
    call prms%CreateIntOption     (  "ComputeAlpha",  " Specifies how to compute the blending coefficient:\n"//&
                                              "   1: Use the shock indicator with heuristic mapping of Hennemann et al. \n"//&
                                              "  11: Use the shock indicator capped in interval [0,1] \n"//&
                                              "   2: Randomly assign the blending coef.,changes over time\n"//&
                                              "   20: Randomly assign the blending coef.0<alpha<alpha_max, fixed over time\n"//&
                                              "   3: Fixed blending coef. (alpha=ShockBlendCoef)"&
                                             ,"1")
    call prms%CreateRealOption(   "ShockBlendCoef",  " Fixed blending coefficient to be used with ComputeAlpha=3", "0.0")

    call prms%CreateIntOption(     "ModalThreshold",  " Threshold function to be used for the indicator (for ComputeAlpha=1) "//&
                                                  "  1: 0.5 * 10.0 ** (-1.8 * (PP_N+1)**0.25)"//&
                                                  "  2: 0.5 * 10.0 ** (-1.8 * PP_N**0.25)"&
                                                 ,"1")
    call prms%CreateRealOption    (      "alpha_max", "Maximum value for the blending coefficient", "0.5")
    call prms%CreateRealOption    (      "alpha_min", "Minimum value for the blending coefficient (below this, alpha=0)", "0.01")
    call prms%CreateRealOption    ("SpacePropFactor", "Space propagation factor", "0.5")
    call prms%CreateIntOption     ("SpacePropSweeps", "Number of space propagation sweeps (MPI-optimized only for 0 or 1)", "1")
    call prms%CreateRealOption    (  "TimeRelFactor", "Time relaxation factor", "0.0")
    call prms%CreateIntOption    ("ReconsBoundaries", "Reconstruction procedure on boundary subcells (only for SubFVMethod=2,3,4):\n"//&
                                                       "   1: No reconstruction\n"//&
                                                       "   2: Central reconstruction\n"//&
                                                       "   3: Neighbor element reconstruction"&
                                                       ,"1")

#if NFVSE_CORR
    call DefineParameters_IDP()
#endif /*NFVSE_CORR*/
   
  end subroutine DefineParametersNFVSE
!===================================================================================================================================
!> Initializes the NFVSE module
!===================================================================================================================================
  subroutine InitNFVSE()
    USE MOD_Globals
    use MOD_NFVSE_Vars
    use MOD_ReadInTools        , only: GETINT, GETREAL, GETLOGICAL
    USE MOD_Mesh_Vars          , only: nElems,nSides,firstSlaveSide,LastSlaveSide, MeshIsNonConforming, firstMortarInnerSide
    use MOD_Interpolation_Vars , only: wGP, xGP
    use MOD_Equation_Vars      , only: RiemannVolFluxAndDissipMatrices
#if USE_AMR
    use MOD_AMR_Vars           , only: UseAMR
#endif /*USE_AMR*/
#if MPI
    use MOD_MPI_Vars           , only: nNbProcs
    use MOD_IDP_Vars           , only: IDPneedsUprev_ext, IDPneedsUsafe_ext
#endif /*MPI*/
#if NFVSE_CORR
    use MOD_IDP, only: Init_IDP
#endif /*NFVSE_CORR*/
    implicit none
    !-local-variables---------------------------------------------------------------------------------------------------------------
    integer :: i
    real    :: sumWm1
    logical :: MeshNonConforming
    !--------------------------------------------------------------------------------------------------------------------------------
    
    ! Safety deallocations
    ! --------------------
    SDEALLOCATE(SubCellMetrics)
    SDEALLOCATE(alpha)
    SDEALLOCATE(alpha_Master)
    SDEALLOCATE(alpha_Slave)
    
    ! Get parameters
    ! --------------
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' NFVSE specific parameters: '
    
    SubFVMethod      = GETINT    ('SubFVMethod','1')
    ComputeAlpha     = GETINT    ('ComputeAlpha','1')
    ShockBlendCoef   = GETREAL   ('ShockBlendCoef','0.0')
    ModalThreshold   = GETINT    ('ModalThreshold','1')
    alpha_max        = GETREAL   ('alpha_max','0.5')
    alpha_min        = GETREAL   ('alpha_min','0.01')
    SpacePropFactor  = GETREAL   ('SpacePropFactor','0.5')
    SpacePropSweeps  = GETINT    ('SpacePropSweeps','1')
    TimeRelFactor    = GETREAL   ('TimeRelFactor'  ,'0.0')
    ! ReconsBoundaries is read afterwards only if needed
    
    ! Initialize everything
    ! ---------------------
    
    ! Check if the mesh can be non-conforming
    MeshNonConforming = MeshIsNonConforming
#if USE_AMR
    MeshNonConforming = UseAMR .or. MeshNonConforming
#endif /*USE_AMR*/
    
    ! Select the sub-FV method
    select case(SubFVMethod)
      case(1) ; Compute_FVFluxes => Compute_FVFluxes_1st_Order
      case(2,3,4) 
        ReconsBoundaries = GETINT('ReconsBoundaries','1')
        
        if (ReconsBoundaries >= RECONS_NEIGHBOR .and. MeshNonConforming) then
          SWRITE(*,*) 'WARNING :: ReconsBoundaries>=3 is only possible for conforming meshes. Setting to 1.'
          ReconsBoundaries = RECONS_NONE
        end if
        
        select case(SubFVMethod)
          case(2) ; Compute_FVFluxes => Compute_FVFluxes_TVD
          case(3) 
            Compute_FVFluxes => Compute_FVFluxes_TVD2ES
            if ( .not. associated(RiemannVolFluxAndDissipMatrices)) then
              SWRITE(*,*) 'ERROR :: SubFVMethod=4 needs a Riemann solver with RiemannVolFluxAndDissipMatrices.'
              stop
            end if
          case(4) 
            Compute_FVFluxes => Compute_FVFluxes_TVD_Fjordholm
            if ( .not. associated(RiemannVolFluxAndDissipMatrices)) then
              SWRITE(*,*) 'ERROR :: SubFVMethod=4 needs a Riemann solver with RiemannVolFluxAndDissipMatrices.'
              stop
            end if
          case default
            SWRITE(*,*) 'ERROR :: SubFVMethod not implemented.'
            stop
        end select
      case default
        SWRITE(*,*) 'ERROR :: SubFVMethod not implemented.'
        stop
    end select
    
    ! Allocate storage
    allocate ( alpha(nElems) )
    allocate ( alpha_Master(firstMortarInnerSide:nSides ) )
    allocate ( alpha_Slave (firstSlaveSide:LastSlaveSide) )
    allocate ( sWGP(0:PP_N) )
    allocate ( SubCellMetrics(nElems) )
    
    ! Some initializations
    alpha        = 0.0
    alpha_Master = 0.0
    alpha_Slave  = 0.0
    select case (ModalThreshold)
      case(1) ; threshold = 0.5 * 10.0 ** (-1.8 * (PP_N+1.)**0.25) ! Sebastian's thresold (Euler and MHD paper)
      case(2) ; threshold = 0.5 * 10.0 ** (-1.8 * PP_N**0.25)      ! New threshold
    end select
    call SubCellMetrics % construct(PP_N)
    call ComputeSubcellMetrics()
    sWGP = 1.d0 / wGP ! Inverse of the quadrature weights (sub-cell dimensions)
    
    ! Compute variables for 2nd order FV
    ! **********************************
    if (SubFVMethod>=2) then
      if (ReconsBoundaries >= RECONS_NEIGHBOR) allocate ( U_ext(1:PP_nVar,0:PP_N,0:PP_N,6,nElems) )
      allocate ( sdxR(1:PP_N-1) )
      allocate ( sdxL(1:PP_N-1) )
      allocate ( rR  (1:PP_N-1) )
      allocate ( rL  (1:PP_N-1) )
      
      sumWm1 = wGP(0) - 1.
      do i=1, PP_N-1
        ! dx calculation (with nodes position)
        ! ------------------------------------
        
        sdxR(i) = 1.0/(xGP(i+1)-xGP(i))
        sdxL(i) = 1.0/(xGP(i)-xGP(i-1))
        
        ! r calculation
        ! -------------
        rL(i) = sumWm1 - xGP(i)
        sumWm1 = sumWm1 + wGP(i)
        rR(i) = sumWm1 - xGP(i)
      end do
    end if
    
    ! For Gauss shock capturing, we initialize an array to store the surface contribution of the FV part
#if (PP_NodeType==1)
    allocate( Ut_DGGauss(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) )
#if FV_BLENDSURFACE
    if (MeshNonConforming) then
      CALL abort(__STAMP__,'Gauss shock capturing with FV_BLENDSURFACE cannot be used with nonconforming meshes!',999,999.)
      RETURN
    end if
    allocate( Ut_FVGauss(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) )
    allocate ( U_master_FV(PP_nVar,0:PP_N,0:PP_N,1:nSides) )
    allocate ( U_slave_FV (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide) )
    allocate (Flux_master_FV(PP_nVar,0:PP_N,0:PP_N,1:nSides))
    allocate (Flux_slave_FV (PP_nVar,0:PP_N,0:PP_N,firstSlaveSide:LastSlaveSide))
#if MPI
    ALLOCATE(MPIRequest_U_FV(nNbProcs,2)    )
    ALLOCATE(MPIRequest_Flux_FV(nNbProcs,2)    )
    MPIRequest_U_FV      = MPI_REQUEST_NULL
    MPIRequest_Flux_FV   = MPI_REQUEST_NULL
#endif /*MPI*/
#endif /*FV_BLENDSURFACE*/
#endif /*(PP_NodeType==1)*/
    
#if NFVSE_CORR
    call Init_IDP()
#endif /*NFVSE_CORR*/
    
#if MPI
    allocate(MPIRequest_alpha(nNbProcs,4)    ) ! 1: send slave, 2: send master, 3: receive slave, 4, receive master
    
    if (ReconsBoundaries >= RECONS_NEIGHBOR .or. IDPneedsUprev_ext .or. IDPneedsUsafe_ext) then
      allocate(MPIRequest_Umaster(nNbProcs,2)) ! 1: send master, 2: receive master
    end if
#endif
  end subroutine InitNFVSE
!===================================================================================================================================
!> Computes the subcell metrics for one element
!===================================================================================================================================
  subroutine ComputeSubcellMetrics()
    USE MOD_Globals
    use MOD_NFVSE_Vars         , only: SubCellMetrics
    use MOD_Mesh_Vars          , only: nElems, Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
    use MOD_Interpolation_Vars , only: wGP
    use MOD_DG_Vars            , only: D
    implicit none
    !-local-variables-----------------------------
    integer :: iElem, i, j, k,l,m         !DOF counters
    real :: Metrics_fCont(3,0:PP_N,0:PP_N,0:PP_N) ! Container for the (reshaped) xi metrics
    real :: Metrics_gCont(3,0:PP_N,0:PP_N,0:PP_N) ! Container for the (reshaped) eta metrics
    real :: Metrics_hCont(3,0:PP_N,0:PP_N,0:PP_N) ! Container for the (reshaped) zeta metrics
    !---------------------------------------------
    
    ! Compute inner normal and tangent vectors
    do iElem=1, nElems
      
!     Xi planes
!     ---------
      Metrics_fCont = reshape(Metrics_fTilde(:,:,:,:,iElem) , shape(Metrics_fCont), order = [1,4,2,3])
      Metrics_gCont = reshape(Metrics_gTilde(:,:,:,:,iElem) , shape(Metrics_gCont), order = [1,4,2,3])
      
      do i=-1, PP_N
        ! Compute vectors
        SubCellMetrics(iElem) % xi   % nv(:,:,:,i) = Metrics_fCont(:,:,:,0)
        SubCellMetrics(iElem) % xi   % t1(:,:,:,i) = Metrics_gCont(:,:,:,0)
        
        do m=0, PP_N  ; do l=0, i
          SubCellMetrics(iElem) % xi   % nv(:,:,:,i) = SubCellMetrics(iElem) % xi   % nv(:,:,:,i) + wGP(l)*D(l,m) * Metrics_fCont(:,:,:,m)
          SubCellMetrics(iElem) % xi   % t1(:,:,:,i) = SubCellMetrics(iElem) % xi   % t1(:,:,:,i) + wGP(l)*D(l,m) * Metrics_gCont(:,:,:,m)
        end do        ; end do
        
        ! Normalize each
        do k=0, PP_N ; do j=0, PP_N
          associate (norm => SubCellMetrics(iElem) % xi   % norm(j,k,i), &
                     nv   => SubCellMetrics(iElem) % xi   % nv(:,j,k,i), &
                     t1   => SubCellMetrics(iElem) % xi   % t1(:,j,k,i), &
                     t2   => SubCellMetrics(iElem) % xi   % t2(:,j,k,i) )
          norm = norm2 (nv)
          nv   = nv / norm
          
          t1   = t1 - dot_product(t1, nv) * nv
          t1   = t1 / norm2(t1)
          
          t2   = CROSS(nv, t1)
          end associate
        end do       ; end do
      end do
      
!     Eta planes
!     ----------
      Metrics_gCont = reshape(Metrics_gTilde(:,:,:,:,iElem) , shape(Metrics_gCont), order = [1,2,4,3])
      Metrics_hCont = reshape(Metrics_hTilde(:,:,:,:,iElem) , shape(Metrics_hCont), order = [1,2,4,3])
      
      do i=-1, PP_N
        ! Compute vectors
        SubCellMetrics(iElem) % eta  % nv(:,:,:,i) = Metrics_gCont(:,:,:,0)
        SubCellMetrics(iElem) % eta  % t1(:,:,:,i) = Metrics_hCont(:,:,:,0)
        
        do m=0, PP_N  ; do l=0, i
          SubCellMetrics(iElem) % eta  % nv(:,:,:,i) = SubCellMetrics(iElem) % eta  % nv(:,:,:,i) + wGP(l)*D(l,m) * Metrics_gCont(:,:,:,m)
          SubCellMetrics(iElem) % eta  % t1(:,:,:,i) = SubCellMetrics(iElem) % eta  % t1(:,:,:,i) + wGP(l)*D(l,m) * Metrics_hCont(:,:,:,m)
        end do        ; end do
        
        ! Normalize each
        do k=0, PP_N ; do j=0, PP_N
          associate (norm => SubCellMetrics(iElem) % eta  % norm(j,k,i), &
                     nv   => SubCellMetrics(iElem) % eta  % nv(:,j,k,i), &
                     t1   => SubCellMetrics(iElem) % eta  % t1(:,j,k,i), &
                     t2   => SubCellMetrics(iElem) % eta  % t2(:,j,k,i) )
          norm = norm2 (nv)
          nv   = nv / norm
          
          t1   = t1 - dot_product(t1, nv) * nv
          t1   = t1 / norm2(t1)
          
          t2   = CROSS(nv, t1)
          end associate
        end do       ; end do
      end do
      
!     Zeta planes
!     -----------
      ! (here we don't have to reshape)
      do i=-1, PP_N
        ! Compute vectors
        SubCellMetrics(iElem) % zeta % nv(:,:,:,i) = Metrics_hTilde(:,:,:,0,iElem)
        SubCellMetrics(iElem) % zeta % t1(:,:,:,i) = Metrics_fTilde(:,:,:,0,iElem)
        
        do m=0, PP_N  ; do l=0, i
          SubCellMetrics(iElem) % zeta % nv(:,:,:,i) = SubCellMetrics(iElem) % zeta % nv(:,:,:,i) + wGP(l)*D(l,m) * Metrics_hTilde(:,:,:,m,iElem)
          SubCellMetrics(iElem) % zeta % t1(:,:,:,i) = SubCellMetrics(iElem) % zeta % t1(:,:,:,i) + wGP(l)*D(l,m) * Metrics_fTilde(:,:,:,m,iElem)
        end do        ; end do
        
        ! Normalize each
        do k=0, PP_N ; do j=0, PP_N
          associate (norm => SubCellMetrics(iElem) % zeta % norm(j,k,i), &
                     nv   => SubCellMetrics(iElem) % zeta % nv(:,j,k,i), &
                     t1   => SubCellMetrics(iElem) % zeta % t1(:,j,k,i), &
                     t2   => SubCellMetrics(iElem) % zeta % t2(:,j,k,i) )
          norm = norm2 (nv)
          nv   = nv / norm
          
          t1   = t1 - dot_product(t1, nv) * nv
          t1   = t1 / norm2(t1)
          
          t2   = CROSS(nv, t1)
          end associate
        end do       ; end do
      end do
      
    end do !iElem
    
  end subroutine ComputeSubcellMetrics
!===================================================================================================================================
!> Reinitializes all variables that need reinitialization after the h-adaptation
!> Transfers the blending coefficient to the new mesh after adaptation (but before load balancing!)
!===================================================================================================================================  
  subroutine InitNFVSEAfterAdaptation1(nElems,ChangeElem)
    use MOD_NFVSE_Vars         , only: alpha
    implicit none
    !-arguments-----------------------------------
    integer, intent(in) :: nElems
    integer, intent(in) :: ChangeElem(8,nElems)
    !-local-variables-----------------------------
    integer          :: eID
    real,allocatable,target :: alphaNew(:)
    !---------------------------------------------
    allocate ( alphaNew(nElems) )
    ! Set with old values
    do eID=1, nElems
      if (ChangeElem(1,eID) < 0) then
        ! refinement
        alphaNew(eID) = alpha(-ChangeElem(1,eID))
      elseif (ChangeElem(2,eID) > 0) then
        ! coarsening
        alphaNew(eID) = maxval(alpha(ChangeElem(1:8,eID)))
      else
        ! simple reasignment
        alphaNew(eID) = alpha(ChangeElem(1,eID))
      endif
    end do
    call move_alloc(alphaNew,alpha)
  end subroutine InitNFVSEAfterAdaptation1
!===================================================================================================================================
!> Reinitializes all variables that need reinitialization after the h-adaptation
!> ATTENTION: The subcell metrics are always recomputed, as the metrics of the high-order DG elements
!===================================================================================================================================
  subroutine InitNFVSEAfterAdaptation2(nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld)
    USE MOD_Globals
    use MOD_NFVSE_Vars         , only: SubCellMetrics, alpha, alpha_Master, alpha_Slave, TimeRelFactor, alpha_max, alpha_min, ComputeAlpha, ShockBlendCoef
    use MOD_Mesh_Vars          , only: nElems,nSides,firstSlaveSide,LastSlaveSide,firstMortarInnerSide
    use MOD_IDP_Vars
#if NFVSE_CORR
    use MOD_NFVSE_Vars         , only: alpha_old
#endif /*NFVSE_CORR*/
#if MPI
    use MOD_MPI_Vars           , only: nNbProcs
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha, MPIRequest_Umaster, ReconsBoundaries, RECONS_NEIGHBOR
#endif /*MPI*/
    implicit none
    !-arguments-----------------------------------
    integer, intent(in) :: nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld
    !-local-variables-----------------------------
    integer          :: eID
    real,allocatable,target :: alphaNew(:)
    !---------------------------------------------
    
!   Reallocate storage
!   ------------------
  
    if ( (firstMortarInnerSideOld .ne. firstMortarInnerSide) .or. (nSides .ne. nSidesOld) ) then
      SDEALLOCATE(alpha_Master)
      allocate ( alpha_Master(firstMortarInnerSide:nSides ) )
    end if
    if ( (firstSlaveSide .ne. firstSlaveSideOld) .or. (LastSlaveSide .ne. LastSlaveSideOld) ) then
      SDEALLOCATE(alpha_Slave)
      allocate ( alpha_Slave (firstSlaveSide:LastSlaveSide) )
    end if
    if (nElems /= nElemsOld) then
      SDEALLOCATE(SubCellMetrics)
      allocate ( SubCellMetrics(nElems) )
      call SubCellMetrics % construct(PP_N)
#if NFVSE_CORR
      ! TODO: Call IDP method
      SDEALLOCATE(FFV_m_FDG)
      SDEALLOCATE(alpha_old)
      allocate ( FFV_m_FDG(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) )
#if LOCAL_ALPHA
      SDEALLOCATE(alpha_loc)
      allocate ( alpha_loc(0:PP_N,0:PP_N,0:PP_N,nElems) )
      do eID=1, nElems
        alpha_loc(:,:,:,eID) = alpha(eID)
      end do
      
      
      allocate ( alpha_old(0:PP_N,0:PP_N,0:PP_N,nElems) )
#else
      allocate ( alpha_old(nElems) )
#endif /*LOCAL_ALPHA*/
      FFV_m_FDG = 0.
      alpha_old = 0.
      
      ! Usafe = U_FV (with external DOFs)
    if (IDPneedsUsafe) then
      SDEALLOCATE(Usafe)
      SDEALLOCATE(p_safe)
      allocate ( Usafe      (PP_nVar,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,nElems) )
      allocate ( p_safe             ( 0:PP_N  , 0:PP_N  , 0:PP_N  ,nElems) )
    end if
    
    ! Solution in the previous step (with external DOFs)
    if (IDPneedsUprev) then
      SDEALLOCATE(Uprev)
      allocate ( Uprev      (PP_nVar,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,nElems) )
    end if
    
    ! Container for external Uprev
    if (IDPneedsUprev_ext) then
      SDEALLOCATE(Uprev_ext)
      allocate( Uprev_ext    (PP_nVar, 0:PP_N, 0:PP_N, 6,nElems) )
    end if
    
    ! Allocate bar states if needed
#if barStates
    if (IDPneedsUbar) then
      SDEALLOCATE(Ubar_xi  )
      SDEALLOCATE(Ubar_eta )
      SDEALLOCATE(Ubar_zeta)
      allocate( Ubar_xi      (PP_nVar,-1:PP_N  , 0:PP_N  , 0:PP_N  ,nElems) )
      allocate( Ubar_eta     (PP_nVar, 0:PP_N  ,-1:PP_N  , 0:PP_N  ,nElems) )
      allocate( Ubar_zeta    (PP_nVar, 0:PP_N  , 0:PP_N  ,-1:PP_N  ,nElems) )
    end if
#else
    if (IDPneedsUsafe_ext) then
      SDEALLOCATE(Usafe_ext)
      allocate( Usafe_ext    (PP_nVar, 0:PP_N  , 0:PP_N          ,6,nElems) )
    end if
#endif /*barStates*/
    
    ! Variables for local alpha
#if LOCAL_ALPHA
    SDEALLOCATE(f_antidiff)
    SDEALLOCATE(g_antidiff)
    SDEALLOCATE(h_antidiff)
    allocate ( f_antidiff(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N,nElems) )
    allocate ( g_antidiff(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N,nElems) )
    allocate ( h_antidiff(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N,nElems) )
    f_antidiff = 0.0
    g_antidiff = 0.0
    h_antidiff = 0.0
#if NONCONS
    SDEALLOCATE(f_antidiffR)
    SDEALLOCATE(g_antidiffR)
    SDEALLOCATE(h_antidiffR)
    allocate ( f_antidiffR(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N,nElems) )
    allocate ( g_antidiffR(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N,nElems) )
    allocate ( h_antidiffR(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N,nElems) )
    f_antidiffR = 0.0
    g_antidiffR = 0.0
    h_antidiffR = 0.0
#endif /*NONCONS*/
    SDEALLOCATE(dalpha_loc)
    allocate ( dalpha_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,nElems) )
#endif /*LOCAL_ALPHA*/
    SDEALLOCATE(dalpha)
    allocate ( dalpha(nElems) )
    
    ! Bounds containers
    if (IDPStateTVD .or. IDPPositivity) then
      SDEALLOCATE(state_min)
      allocate ( state_min   (PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) )
    end if
    if (IDPStateTVD) then
      SDEALLOCATE(state_max)
      allocate ( state_max   (PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) )
    end if
    if (IDPSpecEntropy) then
      SDEALLOCATE(s_min)
      allocate ( s_min       (0:PP_N,0:PP_N,0:PP_N,nElems) )
    end if
    if (IDPMathEntropy) then
      SDEALLOCATE(s_max)
      allocate ( s_max       (0:PP_N,0:PP_N,0:PP_N,nElems) )
    end if
    if (IDPPositivity) then
      SDEALLOCATE(p_min)
      allocate ( p_min       (0:PP_N,0:PP_N,0:PP_N,nElems) )
    end if
#endif /*NFVSE_CORR*/
    end if
    
!   Initialize values
!   -----------------
    alpha_Master = 0.0
    alpha_Slave  = 0.0
    
!   Compute Subcell Metrics
!   -----------------------
    call ComputeSubcellMetrics()
    
!   Reallocate MPI variables
!   ------------------------
#if MPI
    SDEALLOCATE(MPIRequest_alpha)
    allocate(MPIRequest_alpha(nNbProcs,4)    ) ! 1: send slave, 2: send master, 3: receive slave, 4, receive master
    
    if (ReconsBoundaries >= RECONS_NEIGHBOR .or. IDPneedsUprev_ext .or. IDPneedsUsafe_ext) then
      SDEALLOCATE(MPIRequest_Umaster)
      allocate(MPIRequest_Umaster(nNbProcs,2)) ! 1: send master, 2: receive master
    end if
#endif
  end subroutine InitNFVSEAfterAdaptation2
!===================================================================================================================================
!> Computes the "volume integral": Spatial contribution to Ut by the subcell finite volumes
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: This routine has to be called after VolInt_adv_SplitForm, since here Ut is updated with the finite volume contribution
!===================================================================================================================================
  subroutine VolInt_NFVSE(Ut,tIn)
    use MOD_PreProc
    use MOD_DG_Vars            , only: U, U_master, U_slave
    use MOD_Mesh_Vars          , only: nElems, sJ
    use MOD_NFVSE_Vars         , only: SubCellMetrics, sWGP, Compute_FVFluxes, ReconsBoundaries, SpacePropSweeps, RECONS_NEIGHBOR, alpha, U_ext
    use MOD_IDP_Vars           , only: IDPneedsUprev, IDPneedsUprev_ext, Uprev
#if NFVSE_CORR
    use MOD_IDP_Vars           , only: FFV_m_FDG, IDPafterIndicator, IDPPositivity
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars         , only: alpha_loc, f_antidiff, g_antidiff, h_antidiff
#if NONCONS
    use MOD_NFVSE_Vars         , only: f_antidiffR, g_antidiffR, h_antidiffR
#endif /*NONCONS*/
#endif /*LOCAL_ALPHA*/
#endif /*NFVSE_CORR*/
    use MOD_Basis              , only: ALMOSTEQUAL
    use MOD_NFVSE_MPI          , only: PropagateBlendingCoeff, ProlongBlendingCoeffToFaces
#if MPI
    use MOD_NFVSE_Vars         , only: MPIRequest_Umaster
    use MOD_MPI                , only: FinishExchangeMPIData
    use MOD_MPI_Vars           , only: nNbProcs
#endif /*MPI*/
    use MOD_NFVSE_MPI          , only: Get_externalU
#if (PP_NodeType==1)
    use MOD_NFVSE_Vars         , only: Ut_DGGauss
#if FV_BLENDSURFACE
    use MOD_NFVSE_Vars         , only: Ut_FVGauss
#endif /*FV_BLENDSURFACE*/
#endif /*(PP_NodeType==1)*/
    implicit none
    !-arguments---------------------------------------------------------------------------------------------------------------------
    real,intent(inout)                              :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
    real,intent(in)                                 :: tIn
    !-local-variables---------------------------------------------------------------------------------------------------------------
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N) :: ftilde   ! transformed inter-subcell flux in xi (with ghost cells)
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N) :: gtilde   ! transformed inter-subcell flux in eta (with ghost cells)
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N) :: htilde   ! transformed inter-subcell flux in zeta (with ghost cells)
#if NONCONS
    ! For the right interfaces:
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N) :: ftildeR  ! transformed inter-subcell flux in xi (with ghost cells)
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N) :: gtildeR  ! transformed inter-subcell flux in eta (with ghost cells)
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N) :: htildeR  ! transformed inter-subcell flux in zeta (with ghost cells)
#endif /*NONCONS*/
    real,dimension(PP_nVar)                         :: F_FV
    integer                                         :: i,j,k,iElem, sweep
    !===============================================================================================================================
    
    !if reconstruction:
    if (ReconsBoundaries >= RECONS_NEIGHBOR .or. IDPneedsUprev_ext) then
#if MPI
      call FinishExchangeMPIData(2*nNbProcs,MPIRequest_Umaster) 
#endif /*MPI*/
      if (ReconsBoundaries >= RECONS_NEIGHBOR) call Get_externalU(PP_nVar,U_ext,U,U_master,U_slave,tIn)
    end if
    
    if (SpacePropSweeps > 0) then
      ! Receive alpha (MPI) and propagate
      call PropagateBlendingCoeff()
      ! Do furher sweeps if needed (not MPI-optimized)
      do sweep=2, SpacePropSweeps
        call ProlongBlendingCoeffToFaces()
        call PropagateBlendingCoeff()
      end do
    end if
    
#if NFVSE_CORR
    if (IDPneedsUprev) then ! IDP correction needs the solution at the previous step:
      Uprev(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) = U
    end if
#if LOCAL_ALPHA
    do iElem=1, nElems ! Store alpha in local array
      alpha_loc(:,:,:,iElem) = alpha(iElem)
    end do
#endif /*LOCAL_ALPHA*/
#endif /*NFVSE_CORR*/
    
    do iElem=1,nElems
      
#if NFVSE_CORR
      if ( IDPafterIndicator .and. (.not. IDPPositivity) .and. ALMOSTEQUAL(alpha(iElem),0.0) ) cycle
#else
      if ( ALMOSTEQUAL(alpha(iElem),0.0) ) cycle
#endif /*NFVSE_CORR*/
      
!     Compute the finite volume fluxes
!     --------------------------------
      call Compute_FVFluxes( U(:,:,:,:,iElem), ftilde , gtilde , htilde , &
#if NONCONS
                                               ftildeR, gtildeR, htildeR, &
#endif /*NONCONS*/
                                               SubCellMetrics(iElem), iElem )
!     Compute antidiffusive fluxes in the case of local alpha
#if LOCAL_ALPHA
      f_antidiff(:,:,:,:,iElem) = f_antidiff(:,:,:,:,iElem) - ftilde
      g_antidiff(:,:,:,:,iElem) = g_antidiff(:,:,:,:,iElem) - gtilde
      h_antidiff(:,:,:,:,iElem) = h_antidiff(:,:,:,:,iElem) - htilde
#if NONCONS
      f_antidiffR(:,:,:,:,iElem) = f_antidiffR(:,:,:,:,iElem) - ftildeR
      g_antidiffR(:,:,:,:,iElem) = g_antidiffR(:,:,:,:,iElem) - gtildeR
      h_antidiffR(:,:,:,:,iElem) = h_antidiffR(:,:,:,:,iElem) - htildeR
#endif /*NONCONS*/
#endif /*LOCAL_ALPHA*/
!     Update Ut
!     ---------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        ! Get Finite Volume Ut
#if NONCONS
        F_FV =   sWGP(i) * ( ftildeR(:,i,j,k) - ftilde(:,i-1,j  ,k  ) ) &
               + sWGP(j) * ( gtildeR(:,i,j,k) - gtilde(:,i  ,j-1,k  ) ) &
               + sWGP(k) * ( htildeR(:,i,j,k) - htilde(:,i  ,j  ,k-1) )
#else
        F_FV =   sWGP(i) * ( ftilde(:,i,j,k) - ftilde(:,i-1,j  ,k  ) ) &
               + sWGP(j) * ( gtilde(:,i,j,k) - gtilde(:,i  ,j-1,k  ) ) &
               + sWGP(k) * ( htilde(:,i,j,k) - htilde(:,i  ,j  ,k-1) )
#endif /*NONCONS*/
        
#if (PP_NodeType==1)
        ! For Gauss shock capturing, add the surface contribution
        ! -------------------------------------------------------
#if FV_BLENDSURFACE
        ! If FV_BLENDSURFACE is active, use the FV surface contribution (temporary)
        F_FV = F_FV + Ut_FVGauss(:,i,j,k,iElem)
#if !(LOCAL_ALPHA)
        ! Impose FV flux at the border if not using LOCAL_ALPHA
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) - Ut_DGGauss(:,i,j,k,iElem) + Ut_FVGauss(:,i,j,k,iElem)
#endif /*!(LOCAL_ALPHA)*/
#else
        ! use the DG contribution for the FV method as well
        F_FV = F_FV + Ut_DGGauss(:,i,j,k,iElem)
#endif
#endif /*(PP_NodeType==1)*/
        
#if NFVSE_CORR
        FFV_m_FDG(:,i,j,k,iElem) = (F_FV - Ut(:,i,j,k,iElem)) * (-sJ(i,j,k,iElem)) ! Account for the sign change and the Jacobian division
#endif /*NFVSE_CORR*/
        
        Ut   (:,i,j,k,iElem) = (1. - alpha(iElem)) * Ut(:,i,j,k,iElem) +  alpha(iElem) * F_FV

      end do       ; end do       ; end do ! i,j,k
      
    end do ! iElem
  end subroutine VolInt_NFVSE
  
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a first-order FV consistent flux
!===================================================================================================================================
  subroutine Compute_FVFluxes_1st_Order (U, F , G , H , &
#if NONCONS
                                            FR, GR, HR, &
#endif /*NONCONS*/
                                            sCM, iElem )
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: SubCellMetrics_t
#if NONCONS
    USE MOD_Riemann   , only: AddNonConsFlux
#endif /*NONCONS*/
    implicit none
    !-arguments---------------------------------------------------------------
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), intent(in)    :: U   !< The element solution
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: F   !< Left flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: G   !< Left flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: H   !< Left flux in zeta
#if NONCONS
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: FR  !< Right flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: GR  !< Right flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: HR  !< Right flux in zeta
#endif /*NONCONS*/
    type(SubCellMetrics_t)                         , intent(in)    :: sCM       !< Sub-cell metric terms
    integer                                        , intent(in)    :: iElem
    !-local-variables---------------------------------------------------------
    real :: U_ (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)
    real :: F_ (PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
#if NONCONS
    real :: FR_(PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
#endif /*NONCONS*/
    !-------------------------------------------------------------------------
    
!   *****************
!   Initialize values
!   *****************
    
    F   = 0.0
    G   = 0.0
    H   = 0.0
#if NONCONS
    FR  = 0.0
    GR  = 0.0
    HR  = 0.0
#endif /*NONCONS*/
    
!   *********
!   Xi-planes
!   *********
    F_  = 0.0
    U_ = reshape(U , shape(U_), order = [1,4,2,3])
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    call Compute_FVFluxes_1D(U_,U_,F_ , &
#if NONCONS
                                U_,FR_, &
#endif /*NONCONS*/
                                   sCM % xi)
    
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    F  = reshape(F_ , shape(F ), order = [1,3,4,2])
#if NONCONS
    FR = reshape(FR_, shape(FR), order = [1,3,4,2])
#endif /*NONCONS*/

!   **********    
!   Eta-planes
!   **********
    F_ = 0.0
    U_ = reshape(U , shape(U_), order = [1,2,4,3])
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    call Compute_FVFluxes_1D(U_,U_,F_ , &
#if NONCONS
                                U_,FR_, &
#endif /*NONCONS*/
                                   sCM % eta)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    G  = reshape(F_ , shape(G ), order = [1,2,4,3])
#if NONCONS
    GR = reshape(FR_, shape(GR), order = [1,2,4,3])
#endif /*NONCONS*/

!   ***********    
!   Zeta-planes
!   ***********    
!   Fill inner interfaces
!   ---------------------
    call Compute_FVFluxes_1D(U , U , H , &
#if NONCONS
                                 U , HR, &
#endif /*NONCONS*/
                                   sCM % zeta)
  end subroutine Compute_FVFluxes_1st_Order
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a FV consistent flux using a second-order reconstruction of the primitive variables
!===================================================================================================================================
  subroutine Compute_FVFluxes_TVD(U, F , G , H , &
#if NONCONS
                                 FR, GR, HR, &
#endif /*NONCONS*/
                                              sCM, iElem )
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: SubCellMetrics_t
    use MOD_DG_Vars       , only: nTotal_vol
    use MOD_Equation_Vars , only: ConsToPrimVec,PrimToConsVec
    implicit none
    !-arguments---------------------------------------------------------------
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), intent(in)    :: U   !< The element solution
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: F   !< Left flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: G   !< Left flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: H   !< Left flux in zeta
#if NONCONS
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: FR  !< Right flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: GR  !< Right flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: HR  !< Right flux in zeta
#endif /*NONCONS*/
    type(SubCellMetrics_t)                         , intent(in)    :: sCM       !< Sub-cell metric terms
    integer                                        , intent(in)    :: iElem
    !-local-variables---------------------------------------------------------
    real :: UL   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the left
    real :: UR   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the right
    real :: prim (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables
    real :: prim_(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables after reshape
    real :: primL(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables left
    real :: primR(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables right
    real :: F_ (PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
#if NONCONS
    real :: FR_(PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
    real :: U_ (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)
#endif /*NONCONS*/
    !-------------------------------------------------------------------------
    
!   *****************
!   Initialize values
!   *****************
    
    F   = 0.0
    G   = 0.0
    H   = 0.0
#if NONCONS
    FR  = 0.0
    GR  = 0.0
    HR  = 0.0
#endif /*NONCONS*/
    
    call ConsToPrimVec(nTotal_vol,prim,U)
    
!   *********
!   Xi-planes
!   *********
    F_  = 0.0
    prim_ = reshape(prim , shape(prim_), order = [1,4,2,3])
#if NONCONS
    U_    = reshape(U    , shape(U_)   , order = [1,4,2,3])
#endif /*NONCONS*/
    
!   Do the solution reconstruction
!   ------------------------------
    call ReconstructSolution_2nd_Order(prim_,primL,primR,5,3,iElem)
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill fluxes
!   -----------
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
    call Compute_FVFluxes_1D(UL,UR,F_ , &
#if NONCONS
                                U_,FR_, &
#endif /*NONCONS*/
                                sCM % xi)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    F  = reshape(F_ , shape(F ), order = [1,3,4,2])
#if NONCONS
    FR = reshape(FR_, shape(FR), order = [1,3,4,2])
#endif /*NONCONS*/

!   **********    
!   Eta-planes
!   **********
    F_ = 0.0
    prim_ = reshape(prim , shape(prim_), order = [1,2,4,3])
#if NONCONS
    U_    = reshape(U    , shape(U_)   , order = [1,2,4,3])
#endif /*NONCONS*/
    
!   Do the solution reconstruction
!   ------------------------------
    call ReconstructSolution_2nd_Order(prim_,primL,primR,2,4,iElem)
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill fluxes
!   -----------
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
    call Compute_FVFluxes_1D(UL,UR,F_ , &
#if NONCONS
                                U_,FR_, &
#endif /*NONCONS*/
                                         sCM % eta)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    G  = reshape(F_ , shape(G ), order = [1,2,4,3])
#if NONCONS
    GR = reshape(FR_, shape(GR), order = [1,2,4,3])
#endif /*NONCONS*/


!   ***********    
!   Zeta-planes
!   ***********
    
!   Do the solution reconstruction
!   ------------------------------
    call ReconstructSolution_2nd_Order(prim,primL,primR,1,6,iElem)
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill fluxes
!   -----------
    call Compute_FVFluxes_1D(UL,UR,H , &
#if NONCONS
                                U ,HR, &
#endif /*NONCONS*/
                                sCM % zeta)
    
  end subroutine Compute_FVFluxes_TVD
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs an entropy-stable FV consistent flux that is obtained with a central (not dissipative) 
!> two-point flux evaluated at the nodal values (i.e. a first-order reconstruction) plus a dissipation term that uses a second-order 
!> reconstruction of the primitive variables.
!> ATTENTION: 1) These FV fluxes are only ES if the central two-point flux is EC.
!>            2) Since the second-order reconstruction of the primitive variables may break the entropy stability (ES), a check is 
!>               performed in AdvRiemannRecons and the scheme falls to first-order dissipation locally if needed to ensure ES.
!>            3) A Riemann solver that implements RiemannVolFluxAndDissipMatrices is needed.
!>            4) This routine has its own Compute_FVFluxes_1D
!===================================================================================================================================
  subroutine Compute_FVFluxes_TVD2ES(U, F , G , H , &
#if NONCONS
                                        FR, GR, HR, &
#endif /*NONCONS*/
                                              sCM, iElem )
    use MOD_PreProc
    use MOD_NFVSE_Vars    , only: SubCellMetrics_t
    use MOD_Equation_Vars , only: ConsToPrimVec, PrimToConsVec
    use MOD_DG_Vars       , only: nTotal_vol
    implicit none
    !-arguments---------------------------------------------------------------
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), intent(in)    :: U   !< The element solution
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: F   !< Left flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: G   !< Left flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: H   !< Left flux in zeta
#if NONCONS
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: FR  !< Right flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: GR  !< Right flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: HR  !< Right flux in zeta
#endif /*NONCONS*/
    type(SubCellMetrics_t)                         , intent(in)    :: sCM       !< Sub-cell metric terms
    integer                                        , intent(in)    :: iElem
    !-local-variables---------------------------------------------------------
    real :: U_   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)
    real :: UL   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the left
    real :: UR   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the right
    real :: prim (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables
    real :: prim_(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables after reshape
    real :: primL(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables left
    real :: primR(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables right
    real :: F_   (PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
#if NONCONS
    real :: FR_(PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
#endif /*NONCONS*/
    !-------------------------------------------------------------------------
    
!   *****************
!   Initialize values
!   *****************
    
    F   = 0.0
    G   = 0.0
    H   = 0.0
#if NONCONS
    FR  = 0.0
    GR  = 0.0
    HR  = 0.0
#endif /*NONCONS*/
    
    call ConsToPrimVec(nTotal_vol,prim,U)
    
!   *********
!   Xi-planes
!   *********
    F_  = 0.0
    prim_ = reshape(prim , shape(prim_), order = [1,4,2,3])
    U_    = reshape(U    , shape(U_)   , order = [1,4,2,3])
    
!   Do the solution reconstruction
!   ------------------------------
    call ReconstructSolution_2nd_Order(prim_,primL,primR,5,3,iElem)
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
        
!   Fill inner interfaces
!   ---------------------
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
    call Compute_FVFluxes_1D_TVD2ES(U_,UL,UR,F_ , &
#if NONCONS
                                             FR_, &
#endif /*NONCONS*/
                                             sCM % xi)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    F  = reshape(F_ , shape(F ), order = [1,3,4,2])
#if NONCONS
    FR = reshape(FR_, shape(FR), order = [1,3,4,2])
#endif /*NONCONS*/

!   **********    
!   Eta-planes
!   **********
    F_ = 0.0
    prim_ = reshape(prim , shape(prim_), order = [1,2,4,3])
    U_    = reshape(U    , shape(U_)   , order = [1,2,4,3])
    
!   Do the solution reconstruction
!   ------------------------------
    call ReconstructSolution_2nd_Order(prim_,primL,primR,2,4,iElem)
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
        
!   Fill inner interfaces
!   ---------------------
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
    call Compute_FVFluxes_1D_TVD2ES(U_,UL,UR,F_ , &
#if NONCONS
                                             FR_, &
#endif /*NONCONS*/
                                             sCM % eta)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    G  = reshape(F_ , shape(G ), order = [1,2,4,3])
#if NONCONS
    GR = reshape(FR_, shape(GR), order = [1,2,4,3])
#endif /*NONCONS*/

!   ***********    
!   Zeta-planes
!   ***********
    primL = 0.0
    primR = 0.0
    
!   Do the solution reconstruction
!   ------------------------------
    call ReconstructSolution_2nd_Order(prim,primL,primR,1,6,iElem)
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill inner interfaces
!   ---------------------
    call Compute_FVFluxes_1D_TVD2ES(U,UL,UR,H , &
#if NONCONS
                                            HR, &
#endif /*NONCONS*/
                                            sCM % zeta)
! ========
  contains
! ========
!===================================================================================================================================
!> Computes the numerical fluxes and the numerical non-conservative terms for all the DOFs of an inner interface
!> ATTENTION: 1) This works only for Compute_FVFluxes_TVD2ES
!             1) The U_, UL, UR and F_ arrays have to be reshaped before calling this routine
!===================================================================================================================================
    pure subroutine Compute_FVFluxes_1D_TVD2ES(U_,UL,UR,F_ , &
#if NONCONS
                                                 FR_, &
#endif /*NONCONS*/
                                                 metrics)
      use MOD_NFVSE_Vars, only: InnerFaceMetrics_t
#if NONCONS
      use MOD_Riemann   , only: AddNonConsFlux
#endif /*NONCONS*/
      implicit none
      !-arguments---------------------------------------------------------------
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: U_       !< The nodal solution
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: UL       !< The solution on the left of the subcells
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: UR       !< The solution on the right of the subcells
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N ), intent(inout) :: F_       !< Left flux
#if NONCONS
      real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N ), intent(inout) :: FR_      !< Right flux
#endif /*NONCONS*/
      type(InnerFaceMetrics_t)                        , intent(in)    :: metrics  !< Sub-cell metric terms
      !-local-variables---------------------------------------------------------
      integer  :: i,j,k
      !-------------------------------------------------------------------------
      
!     Loop over inner interfaces
!     --------------------------
      do i=0, PP_N-1
        
        call AdvRiemannRecons(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1),UR(:,:,:,i),UL(:,:,:,i+1), &
                              metrics%nv(:,:,:,i),metrics%t1(:,:,:,i),metrics%t2(:,:,:,i))
        
#if NONCONS
        ! Copy conservative part
        FR_(:,:,:,i) = F_(:,:,:,i)
        
        ! Add nonconservative fluxes
        CALL AddNonConsFlux(F_(:,:,:,i), &
                            U_(:,:,:,i+1), U_(:,:,:,i),&
                            metrics%nv(:,:,:,i),metrics%t1(:,:,:,i),metrics%t2(:,:,:,i))
        CALL AddNonConsFlux(FR_(:,:,:,i), &
                            U_(:,:,:,i  ), U_(:,:,:,i+1),&
                            metrics%nv(:,:,:,i),metrics%t1(:,:,:,i),metrics%t2(:,:,:,i))
        
        ! Scale right flux
        do k=0, PP_N ; do j=0, PP_N
          FR_(:,j,k,i) = FR_(:,j,k,i) * metrics%norm(j,k,i)
        end do       ; end do
#endif /*NONCONS*/
        
        ! Scale flux
        do k=0, PP_N ; do j=0, PP_N
          F_ (:,j,k,i) = F_ (:,j,k,i) * metrics%norm(j,k,i)
        end do       ; end do
      end do
    end subroutine Compute_FVFluxes_1D_TVD2ES
!===================================================================================================================================
!> Computes the advective Riemann solver for an inner interface with reconstructed states. If the flux is not EC, it returns the 
!> Riemann solution computed with the nodal states
!===================================================================================================================================
    pure subroutine AdvRiemannRecons(F,UL,UR,UL_r,UR_r,nv,t1,t2)
      use MOD_Riemann       , only: RotateState, RotateFluxBack
      use MOD_Equation_Vars , only: RiemannVolFluxAndDissipMatrices, ConsToEntropy
      implicit none
      !-arguments-----------------------------------------------------------
      real, intent(out):: F(       PP_nVar,0:PP_N,0:PP_N) !< numerical flux on face
      real, intent(in) :: UL(      PP_nVar,0:PP_N,0:PP_N) !<  left state on face
      real, intent(in) :: UR(      PP_nVar,0:PP_N,0:PP_N) !< right state on face
      real, intent(in) :: UL_r(    PP_nVar,0:PP_N,0:PP_N) !<  left state on face
      real, intent(in) :: UR_r(    PP_nVar,0:PP_N,0:PP_N) !< right state on face
      real, intent(in) :: nv(            3,0:PP_N,0:PP_N) !< normal vector of face
      real, intent(in) :: t1(            3,0:PP_N,0:PP_N) !< 1st tangential vector of face
      real, intent(in) :: t2(            3,0:PP_N,0:PP_N) !< 2nd tangential vector of face
      !-local-variables-----------------------------------------------------
      integer :: j,k
      real :: U_L  (PP_nVar), U_L_r  (PP_nVar)       ! Rotated nodal and reconstructed state on the left
      real :: U_R  (PP_nVar), U_R_r  (PP_nVar)       ! Rotated nodal and reconstructed state on the right
      real :: RT_Dv(PP_nVar), RT_Dv_r(PP_nVar)       ! Jump of entropy variables on the left (also multiplied by the R^T matrix) of rotated states
      real :: Dmat (PP_nVar)
      real :: Rmat (PP_nVar,PP_nVar)
      real :: RmatT(PP_nVar,PP_nVar)
      !---------------------------------------------------------------------
      
      ! Loop over the face and get the ES flux
      do k=0, PP_N ; do j=0, PP_N
        ! Rotate states to 1D frame
        ! -------------------------
        U_L   = RotateState(UL  (:,j,k), nv(:,j,k), t1(:,j,k), t2(:,j,k))
        U_R   = RotateState(UR  (:,j,k), nv(:,j,k), t1(:,j,k), t2(:,j,k))
        U_L_r = RotateState(UL_r(:,j,k), nv(:,j,k), t1(:,j,k), t2(:,j,k))
        U_R_r = RotateState(UR_r(:,j,k), nv(:,j,k), t1(:,j,k), t2(:,j,k))
        ! Compute EC flux and dissipation matrices with nodal states
        ! ----------------------------------------------------------
        call RiemannVolFluxAndDissipMatrices(U_L,U_R,F(:,j,k),Dmat,Rmat)
        RmatT = transpose(Rmat)
        ! Get jump of entropy vars of rotated states
        ! ------------------------------------------
        RT_Dv   = ConsToEntropy(U_R)   - ConsToEntropy(U_L)
        RT_Dv_r = ConsToEntropy(U_R_r) - ConsToEntropy(U_L_r)
        
        ! Check sign condition and compute flux
        ! First variant: Use the dissipation matrices computed with the nodal states
        ! **************************************************************************
        
        ! Scale entropy jump with the right-eigV mat
        RT_Dv   = matmul(RmatT,RT_Dv)
        RT_Dv_r = matmul(RmatT,RT_Dv_r)
        
        ! Discrete switch    
        if ( any(RT_Dv*RT_Dv_r < -1.e-10) ) then ! Reconstructed is not ES
          F(:,j,k) = F(:,j,k) - 0.5*matmul(Rmat,Dmat*RT_Dv)
        else  ! Reconstructed is ES
          F(:,j,k) = F(:,j,k) - 0.5*matmul(Rmat,Dmat*RT_Dv_r)
        end if
        
        ! Rotate flux back to 3D frame
        ! ----------------------------
        call RotateFluxBack(F(:,j,k), nv(:,j,k), t1(:,j,k), t2(:,j,k))
      end do       ; end do !j,k

    end subroutine AdvRiemannRecons
    
  end subroutine Compute_FVFluxes_TVD2ES
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs an entropy-stable FV consistent flux à la Fjordholm that is obtained with a central  
!> (not dissipative) two-point flux evaluated at the nodal values (i.e. a first-order reconstruction) plus a dissipation term that 
!> uses a second-order reconstruction of the "scaled entropy variables".
!>       See: 1) Fjordholm, U. S., Mishra, S., & Tadmor, E. (2012). Arbitrarily high-order accurate entropy stable essentially nonoscillatory schemes for systems of conservation laws. SIAM Journal on Numerical Analysis, 50(2), 544-573. 
!>            2) Rueda-Ramírez, A. M., Hennemann, S., Hindenlang, F. J., Winters, A. R., & Gassner, G. (2020). An Entropy Stable Nodal Discontinuous Galerkin Method for the resistive MHD Equations. Part II: Subcell Finite Volume Shock Capturing. arXiv preprint arXiv:2012.12040.
!> ATTENTION: 1) A Riemann solver that implements RiemannVolFluxAndDissipMatrices is needed.
!>            2) The resulting FV fluxes are only ES if VolumeFluxAverage is an EC flux.
!===================================================================================================================================
  subroutine Compute_FVFluxes_TVD_Fjordholm(U, F , G , H , &
#if NONCONS
                                               FR, GR, HR, &
#endif /*NONCONS*/
                                              sCM, iElem )
    use MOD_PreProc
    use MOD_NFVSE_Vars        , only: SubCellMetrics_t
    implicit none
    !-arguments---------------------------------------------------------------
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), intent(in)    :: U   !< The element solution
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: F   !< Left flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: G   !< Left flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: H   !< Left flux in zeta
#if NONCONS
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: FR  !< Right flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: GR  !< Right flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: HR  !< Right flux in zeta
#endif /*NONCONS*/
    type(SubCellMetrics_t)                         , intent(in)    :: sCM       !< Sub-cell metric terms
    integer                                        , intent(in)    :: iElem
    !-local-variables---------------------------------------------------------
    real :: U_   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)
    real :: F_   (PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
#if NONCONS
    real :: FR_          (PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
#endif /*NONCONS*/
    !-------------------------------------------------------------------------
    
!   *****************
!   Initialize values
!   *****************
    
    F   = 0.0
    G   = 0.0
    H   = 0.0
#if NONCONS
    FR  = 0.0
    GR  = 0.0
    HR  = 0.0
#endif /*NONCONS*/
    
!   *********
!   Xi-planes
!   *********
    F_  = 0.0
    U_    = reshape(U    , shape(U_)   , order = [1,4,2,3])
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
    call Compute_FVFluxes_1D_TVD_Fjordholm(U_,F_ , &
#if NONCONS
                                              FR_, &
#endif /*NONCONS*/
                                              sCM % xi  ,5,3,iElem)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    F  = reshape(F_ , shape(F ), order = [1,3,4,2])
#if NONCONS
    FR = reshape(FR_, shape(FR), order = [1,3,4,2])
#endif /*NONCONS*/

!   **********    
!   Eta-planes
!   **********
    F_ = 0.0
    U_    = reshape(U    , shape(U_)   , order = [1,2,4,3])
#if NONCONS
    FR_ = 0.0
#endif /*NONCONS*/
    
    call Compute_FVFluxes_1D_TVD_Fjordholm(U_,F_ , &
#if NONCONS
                                              FR_, &
#endif /*NONCONS*/
                                              sCM % eta ,2,4,iElem)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    G  = reshape(F_ , shape(G ), order = [1,2,4,3])
#if NONCONS
    GR = reshape(FR_, shape(GR), order = [1,2,4,3])
#endif /*NONCONS*/

!   ***********    
!   Zeta-planes
!   ***********
    
    call Compute_FVFluxes_1D_TVD_Fjordholm(U ,H , &
#if NONCONS
                                              HR, &
#endif /*NONCONS*/
                                              sCM % zeta,1,6,iElem)
! ========
  contains
! ========
!===================================================================================================================================
!> Computes the numerical fluxes and the numerical non-conservative terms for all the DOFs of an inner interface
!> ATTENTION: 1) This works only for Compute_FVFluxes_TVD_Fjordholm
!             2) The reconstruction in scaled entropy variables is done here!
!             3) The U_ array has to be reshaped before calling this routine
!===================================================================================================================================
    pure subroutine Compute_FVFluxes_1D_TVD_Fjordholm(U_,F_ , &
#if NONCONS
                                                         FR_, &
#endif /*NONCONS*/
                                                         metrics,sideL,sideR,iElem)
      use MOD_Riemann           , only: RotateState, RotateFluxBack
      use MOD_NFVSE_Vars        , only: InnerFaceMetrics_t, sdxR, sdxL, rL, rR, U_ext
      use MOD_NFVSE_Vars        , only: ReconsBoundaries, RECONS_CENTRAL, RECONS_NEIGHBOR, RECONS_NONE
      use MOD_Interpolation_Vars, only: wGP
      use MOD_Equation_Vars     , only: RiemannVolFluxAndDissipMatrices, ConsToEntropy
#if NONCONS
      use MOD_Riemann           , only: AddNonConsFlux
#endif /*NONCONS*/
      implicit none
      !-arguments---------------------------------------------------------------      
      real,dimension        (PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: U_       !< The nodal solution
      real,dimension        (PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N ), intent(inout) :: F_       !< Left flux
#if NONCONS      
      real,dimension        (PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N ), intent(inout) :: FR_      !< Right flux
#endif /*NONCONS*/
      type(InnerFaceMetrics_t)                                , intent(in)    :: metrics  !< Sub-cell metric terms
      integer                                                 , intent(in)    :: sideL    !< Side on the left
      integer                                                 , intent(in)    :: sideR    !< Side on the right
      integer                                                 , intent(in)    :: iElem    !< Element index
      !-local-variables---------------------------------------------------------
      integer :: i,j,k
      real :: Rmat (PP_nVar,PP_nVar,0:PP_N,0:PP_N,0:PP_N-1) ! Right eigenvalue matrix for the nodes of each interface
      real :: RmatT(PP_nVar,PP_nVar)                        ! Right eigenvalue matrix transposed
      real :: Dmat         (PP_nVar,0:PP_N,0:PP_N,0:PP_N-1) ! Dissipation matrix for the nodes of each interface
      real :: wTildeL      (PP_nVar,0:PP_N,0:PP_N,0:PP_N)   ! Reconstructed scaled entropy vars variables on the left of *subcell*
      real :: wTildeR      (PP_nVar,0:PP_N,0:PP_N,0:PP_N)   ! Reconstructed scaled entropy vars variables on the right of *subcell*
      real :: w_L          (PP_nVar,0:PP_N,0:PP_N,0:PP_N-1) ! Scaled entropy vars variables on the left of interface
      real :: w_R          (PP_nVar,0:PP_N,0:PP_N,0:PP_N-1) ! Scaled entropy vars variables on the right of interface
      real :: U_L          (PP_nVar)                        ! Rotated state on the left of interface
      real :: U_R          (PP_nVar)                        ! Rotated state on the right of interface
      real :: v_L          (PP_nVar)                        ! Rotated entropy vars on the left of interface
      real :: v_R          (PP_nVar)                        ! Rotated entropy vars on the right of interface
      real :: sigma        (PP_nVar,0:PP_N,0:PP_N)          ! Slope (limited)
      real :: DummyF       (PP_nVar)                        ! Dummy variable to store the central flux when only the dissipation matrices are needed
      real :: Rmat0(PP_nVar,PP_nVar)
      real :: Dmat0(PP_nVar)
      real :: w_L0         (PP_nVar,0:PP_N,0:PP_N)
      real :: w_R0         (PP_nVar,0:PP_N,0:PP_N)
      !-------------------------------------------------------------------------
      
      wTildeL = 0.0
      wTildeR = 0.0
    
!     Get the dissipation matrices and compute the scaled entropy vars
!     ****************************************************************
      ! Loop over the inner faces and compute
      do i=0, PP_N-1
        ! Loop over the DOFs of the face
        do k=0, PP_N ; do j=0, PP_N
          ! Rotate state
          U_L = RotateState(U_(:,j,k,i  ),metrics%nv(:,j,k,i),metrics%t1(:,j,k,i), metrics%t2(:,j,k,i))
          U_R = RotateState(U_(:,j,k,i+1),metrics%nv(:,j,k,i),metrics%t1(:,j,k,i), metrics%t2(:,j,k,i))
          ! Compute EC flux and dissipation matrices
          call RiemannVolFluxAndDissipMatrices(U_L,U_R,F_(:,j,k,i),Dmat(:,j,k,i),Rmat(:,:,j,k,i))
          RmatT = transpose(Rmat(:,:,j,k,i))
          ! get entropy vars from rotated state
          v_L = ConsToEntropy(U_L)
          v_R = ConsToEntropy(U_R)
          ! Compute scaled entropy vars
          w_L(:,j,k,i) = matmul(RmatT,v_L)
          w_R(:,j,k,i) = matmul(RmatT,v_R)
        end do       ; end do !j,k
      end do !i
    
!     Do the reconstruction of the scaled entropy variables
!     *****************************************************
    
      ! First DOF
      !----------
      
      select case (ReconsBoundaries)
        case (RECONS_NONE)
          wTildeR(:,:,:,0) = w_L(:,:,:,0)
        case (RECONS_CENTRAL)
          wTildeR(:,:,:,0) = w_L(:,:,:,0) + (w_R(:,:,:,0) - w_L(:,:,:,0))*sdxL(1)*wGP(0)
        case (RECONS_NEIGHBOR)
          do k=0, PP_N ; do j=0, PP_N
            ! Rotate state
            U_L = RotateState(U_ext(:,j,k,sideL,iElem),metrics%nv(:,j,k,-1),metrics%t1(:,j,k,-1), metrics%t2(:,j,k,-1))
            U_R = RotateState(U_   (:,j,k,1)          ,metrics%nv(:,j,k,-1),metrics%t1(:,j,k,-1), metrics%t2(:,j,k,-1))
            ! get entropy vars from rotated state
            v_L = ConsToEntropy(U_L)
            v_R = ConsToEntropy(U_R)
            ! Get dissipation matrices
            call RiemannVolFluxAndDissipMatrices(U_L,U_R,DummyF,Dmat0,Rmat0)
            Rmat0 = transpose(Rmat0)
            ! Compute scaled entropy vars
            w_L0(:,j,k) = matmul(Rmat0,v_L)
            w_R0(:,j,k) = matmul(Rmat0,v_R)
          end do       ; end do !j,k
          
          sigma = minmod(sdxL(1)*(w_R(:,:,:,0) - w_L(:,:,:,0)),sdxL(1)*(w_R0 - w_L0))
          wTildeR(:,:,:,0) = w_L(:,:,:,0) + sigma*wGP(0)
      end select
      
      ! Middle DOFs
      ! -----------
      do i=1, PP_N-1
        sigma = minmod(sdxR(i)*(w_R(:,:,:,i)-w_L(:,:,:,i)),sdxL(i)*(w_R(:,:,:,i-1)-w_L(:,:,:,i-1)))
        
        wTildeR(:,:,:,i) = w_L(:,:,:,i  ) + sigma * rR(i)
        wTildeL(:,:,:,i) = w_R(:,:,:,i-1) + sigma * rL(i)
      end do
      
      ! Last DOF
      ! --------
      
      select case (ReconsBoundaries)
        case (RECONS_NONE)
          wTildeL(:,:,:,PP_N) = w_R(:,:,:,PP_N-1)
        case (RECONS_CENTRAL)
          wTildeL(:,:,:,PP_N) = w_R(:,:,:,PP_N-1) - (w_R(:,:,:,PP_N-1) - w_L(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
        case (RECONS_NEIGHBOR)
          do k=0, PP_N ; do j=0, PP_N
            ! Rotate state
            U_L = RotateState(U_   (:,j,k,PP_N-1 )    ,metrics%nv(:,j,k,PP_N),metrics%t1(:,j,k,PP_N), metrics%t2(:,j,k,PP_N))
            U_R = RotateState(U_ext(:,j,k,sideR,iElem),metrics%nv(:,j,k,PP_N),metrics%t1(:,j,k,PP_N), metrics%t2(:,j,k,PP_N))
            ! get entropy vars from rotated state
            v_L = ConsToEntropy(U_L)
            v_R = ConsToEntropy(U_R)
            ! Get dissipation matrices
            call RiemannVolFluxAndDissipMatrices(U_L,U_R,DummyF,Dmat0,Rmat0)
            Rmat0 = transpose(Rmat0)
            ! Compute scaled entropy vars
            w_L0(:,j,k) = matmul(Rmat0,v_L)
            w_R0(:,j,k) = matmul(Rmat0,v_R)
          end do       ; end do !j,k
          
          sigma = minmod(sdxL(1)*(w_R(:,:,:,PP_N-1) - w_L(:,:,:,PP_N-1)),sdxL(1)*(w_R0 - w_L0))
          wTildeL(:,:,:,PP_N) = w_R(:,:,:,PP_N-1) - sigma*wGP(PP_N)
          
      end select
      
!     Compute fluxes 
!     **************

!     Loop over inner interfaces
!     --------------------------
      
      do i=0, PP_N-1
        ! Loop over DOFs on inner face
        do k=0, PP_N ; do j=0, PP_N
          ! Add the right amount of dissipation
          F_(:,j,k,i) = F_(:,j,k,i) - 0.5*matmul(Rmat(:,:,j,k,i),Dmat(:,j,k,i)*(wTildeL(:,j,k,i+1)-wTildeR(:,j,k,i)))
          
          ! Rotate flux back to the physical frame
          call RotateFluxBack(F_(:,j,k,i),metrics%nv(:,j,k,i),metrics%t1(:,j,k,i), metrics%t2(:,j,k,i))
        end do       ; end do !j,k
        
#if NONCONS
        ! Copy conservative part
        FR_(:,:,:,i) = F_(:,:,:,i)
        
        ! Add nonconservative fluxes
        CALL AddNonConsFlux(F_(:,:,:,i), &
                            U_(:,:,:,i+1), U_(:,:,:,i),&
                            metrics%nv(:,:,:,i),metrics%t1(:,:,:,i), metrics%t2(:,:,:,i))
        CALL AddNonConsFlux(FR_(:,:,:,i), &
                            U_(:,:,:,i  ), U_(:,:,:,i+1),&
                            metrics%nv(:,:,:,i),metrics%t1(:,:,:,i), metrics%t2(:,:,:,i))
        
        ! Scale right flux
        do k=0, PP_N ; do j=0, PP_N
          FR_(:,j,k,i) = FR_(:,j,k,i) * metrics%norm(j,k,i)
        end do       ; end do
#endif /*NONCONS*/
        
        ! Scale flux
        do k=0, PP_N ; do j=0, PP_N
          F_ (:,j,k,i) = F_ (:,j,k,i) * metrics%norm(j,k,i)
        end do       ; end do
      end do
      
    end subroutine Compute_FVFluxes_1D_TVD_Fjordholm
  end subroutine Compute_FVFluxes_TVD_Fjordholm
!===================================================================================================================================
!> Computes the numerical fluxes and the numerical non-conservative terms for all the DOFs of an inner interface
!> ATTENTION: 1) The U_, UL and UR arrays have to be reshaped before calling this routine
!>            2) Works for Compute_FVFluxes_1st_Order and Compute_FVFluxes_TVD
!===================================================================================================================================
  pure subroutine Compute_FVFluxes_1D(UL,UR,F_ , &
#if NONCONS
                                         U_,FR_, &
#endif /*NONCONS*/
                                            metrics)
    use MOD_Riemann      , only: AdvRiemann
    use MOD_Equation_Vars, only: SolveRiemannProblem
    use MOD_NFVSE_Vars   , only: InnerFaceMetrics_t
#if NONCONS
    use MOD_Riemann   , only: AddNonConsFlux
#endif /*NONCONS*/
    implicit none
    !-arguments---------------------------------------------------------------
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: UL       !< The solution on the left of the subcells
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: UR       !< The solution on the right of the subcells... It can be the same as UL for 1st order: there's no aliasing problem because of the intent(in)
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N ), intent(inout) :: F_       !< Left flux
#if NONCONS
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: U_       !< The nodal solution (to compute the non-conservative terms)... It can be the same as UL and/or UR: there's no aliasing problem because of the intent(in)
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N ), intent(inout) :: FR_      !< Right flux
#endif /*NONCONS*/
    type(InnerFaceMetrics_t)                        , intent(in)    :: metrics  !< Sub-cell metric terms
    !-local-variables---------------------------------------------------------
    integer  :: i,j,k
    !-------------------------------------------------------------------------
    
!   Loop over inner interfaces
!   --------------------------
    do i=0, PP_N-1
      
      call AdvRiemann(F_(:,:,:,i),UR(:,:,:,i),UL(:,:,:,i+1), &
                      metrics%nv(:,:,:,i),metrics%t1(:,:,:,i),metrics%t2(:,:,:,i),SolveRiemannProblem)
      
#if NONCONS
      ! Copy conservative part
      FR_(:,:,:,i) = F_(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddNonConsFlux(F_(:,:,:,i), &
                          U_(:,:,:,i+1), U_(:,:,:,i),&
                          metrics%nv(:,:,:,i),metrics%t1(:,:,:,i),metrics%t2(:,:,:,i))
      CALL AddNonConsFlux(FR_(:,:,:,i), &
                          U_(:,:,:,i  ), U_(:,:,:,i+1),&
                          metrics%nv(:,:,:,i),metrics%t1(:,:,:,i),metrics%t2(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        FR_(:,j,k,i) = FR_(:,j,k,i) * metrics%norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * metrics%norm(j,k,i)
      end do       ; end do
    end do
  end subroutine Compute_FVFluxes_1D
!===================================================================================================================================
!> Perdorms a second order reconstruction of the solution (primitive variables)
!> Attention: 1) prim_ must be reshaped before calling this routine
!===================================================================================================================================
  pure subroutine ReconstructSolution_2nd_Order(prim_,primL,primR,sideL,sideR,iElem)
    use MOD_Preproc
    use MOD_NFVSE_Vars        , only: sdxR, sdxL, rL, rR
    use MOD_NFVSE_Vars        , only: U_ext, ReconsBoundaries, RECONS_CENTRAL, RECONS_NONE, RECONS_NEIGHBOR
    use MOD_Equation_Vars     , only: ConsToPrimVec
    use MOD_Interpolation_Vars, only: wGP
    implicit none
    !-arguments---------------------------------------------------------------
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(in)    :: prim_      !< The solution in primitive variables
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(out)   :: primL      !< The reconstructed solution (prim vars) on the left of the subcells
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N ), intent(out)   :: primR      !< The reconstructed solution (prim vars) on the right of the subcells
    integer                                         , intent(in)    :: sideL
    integer                                         , intent(in)    :: sideR
    integer                                         , intent(in)    :: iElem
    !-local-variables---------------------------------------------------------
    integer  :: i
    real :: prim_extL(PP_nVar, 0:PP_N, 0:PP_N)
    real :: prim_extR(PP_nVar, 0:PP_N, 0:PP_N)
    real :: sigma    (PP_nVar, 0:PP_N, 0:PP_N)
    !-------------------------------------------------------------------------
    
    ! Initialize values
    primL = 0.0
    primR = 0.0
    
    if (ReconsBoundaries >= RECONS_NEIGHBOR) call ConsToPrimVec((PP_N+1)**2,prim_extL,U_ext(:,:,:,sideL,iElem) )
    if (ReconsBoundaries >= RECONS_NEIGHBOR) call ConsToPrimVec((PP_N+1)**2,prim_extR,U_ext(:,:,:,sideR,iElem) )
    
    ! First DOF
    !----------
    
    select case (ReconsBoundaries)
      case (RECONS_NONE)
        primR(:,:,:,0) = prim_(:,:,:,0)
      case (RECONS_CENTRAL)
        primR(:,:,:,0) = prim_(:,:,:,0) + (prim_(:,:,:,1) - prim_(:,:,:,0))*sdxL(1)*wGP(0)
      case (RECONS_NEIGHBOR)
        sigma = minmod(sdxL(1)*(prim_(:,:,:,1) - prim_(:,:,:,0)),sdxL(1)*(prim_(:,:,:,1) - prim_extL))
        primR(:,:,:,0) = prim_(:,:,:,0) + sigma*wGP(0)
    end select
    
    ! Middle DOFs
    ! -----------
    do i=1, PP_N-1
      sigma = minmod(sdxR(i)*(prim_(:,:,:,i+1)-prim_(:,:,:,i)),sdxL(i)*(prim_(:,:,:,i)-prim_(:,:,:,i-1)))
      
      primR(:,:,:,i) = prim_(:,:,:,i) + sigma * rR(i)
      primL(:,:,:,i) = prim_(:,:,:,i) + sigma * rL(i)
    end do
    
    ! Last DOF
    ! --------
    
    select case (ReconsBoundaries)
      case (RECONS_NONE)
        primL(:,:,:,PP_N) = prim_(:,:,:,PP_N)
      case (RECONS_CENTRAL)
        primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - (prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
      case (RECONS_NEIGHBOR)
        sigma = minmod(sdxL(1)*(prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1)),sdxL(1)*(prim_extR-prim_(:,:,:,PP_N-1)))
        primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - sigma*wGP(0)
    end select
  end subroutine ReconstructSolution_2nd_Order
!===================================================================================================================================
!> Traditional minmod limiter
!===================================================================================================================================
  elemental real function minmod(a,b)
    implicit none
    !---------------------------
    real, intent(in) :: a
    real, intent(in) :: b
    !---------------------------
    
    if (abs(a) < abs(b) .and. a*b > 0.0) then
      minmod = a
    else if (abs(b) < abs(a) .and. a*b > 0.0) then
      minmod = b
    else
      minmod = 0.0
    endif
    
  end function minmod
!===================================================================================================================================
!> Routines to compute the blending coefficient for NFVSE
!> -> See Hennemann and Gassner (2020). "Entropy stable shock capturing for the discontinuous galerkin spectral element
!>                                          method with native finite volume sub elements"
!> -> This routine computes the sensor, makes the correction (with alpha_min and alpha_max), and sends the information with MPI
!> -> No propagation is done yet (MPI informationmust be received).
!===================================================================================================================================
  subroutine CalcBlendingCoefficient(U)
    use MOD_PreProc
    use MOD_Mesh_Vars          , only: nElems
    use MOD_NFVSE_MPI          , only: ProlongBlendingCoeffToFaces, PropagateBlendingCoeff
    use MOD_NFVSE_Vars         , only: SpacePropSweeps, TimeRelFactor, alpha, alpha_max, alpha_min, ComputeAlpha, ShockBlendCoef, sharpness, threshold
    use MOD_ShockCapturing_Vars, only: Shock_Indicator, ShockIndicatorNum
    ! For reconstruction on boundaries
#if MPI
    use MOD_Mesh_Vars          , only: nSides
    use MOD_NFVSE_Vars         , only: ReconsBoundaries, MPIRequest_Umaster, RECONS_NEIGHBOR
    use MOD_MPI                , only: StartReceiveMPIData,StartSendMPIData
    USE MOD_MPI_Vars
    use MOD_DG_Vars            , only: U_master
#endif /*MPI*/
    use MOD_IDP_Vars , only: IDPneedsUprev_ext
#if (PP_NodeType==1)
    use MOD_NFVSE_Vars         , only: Ut_DGGauss, U_master_FV
#if FV_BLENDSURFACE
    use MOD_NFVSE_Vars         , only: Ut_FVGauss
#endif /*FV_BLENDSURFACE*/
#endif /*(PP_NodeType==1)*/
    implicit none
    ! Arguments
    !---------------------------------------------------------------------------------------------------------------------------------
    real,dimension(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems), intent(in)  :: U
    !---------------------------------------------------------------------------------------------------------------------------------
    ! Local variables
    real ::  eta(nElems), eta0(nElems)
    integer :: eID, ind
    !---------------------------------------------------------------------------------------------------------------------------------
    
! If we do reconstruction on boundaries, we need to send the U_master
! -------------------------------------------------------------------
#if MPI
    if (ReconsBoundaries >= RECONS_NEIGHBOR .or. IDPneedsUprev_ext) then
      ! We send and receive the nodal (~mean) solution of the boundary subcell of the neighbor elements!
      ! * For Gauss nodes, that means U_master_FV and U_slave_FV
      ! * For LGL nodes, that means U_master and U_slave
#if (PP_NodeType==1)
      ! receive the master
      call StartReceiveMPIData(U_master_FV, DataSizeSide, 1, nSides, &
                               MPIRequest_Umaster(:,1), SendID=1) ! Receive YOUR  (sendID=1) 
      
      ! Send the master
      call StartSendMPIData   (U_master_FV, DataSizeSide, 1, nSides, &
                               MPIRequest_Umaster(:,2),SendID=1) 
#else /*(PP_NodeType==2)*/
      ! receive the master
      call StartReceiveMPIData(U_master, DataSizeSide, 1, nSides, &
                               MPIRequest_Umaster(:,1), SendID=1) ! Receive YOUR  (sendID=1) 
      
      ! Send the master
      call StartSendMPIData   (U_master, DataSizeSide, 1, nSides, &
                               MPIRequest_Umaster(:,2),SendID=1) 
#endif /*(PP_NodeType==1)*/
    end if
#endif /*MPI*/
    
!   Compute the blending coefficients
!   ---------------------------------
    select case (ComputeAlpha)
      case(1)   ! Shock indicator with heuristic mapping of Hennemann
        ! Compute first shock indicator
        eta = Shock_Indicator(1) % compute(U)
        
        ! Compute other shock indicators
        do ind=2, ShockIndicatorNum
          eta0 = Shock_Indicator(ind) % compute(U)
          do eID=1, nElems
            eta(eID) = max(eta(eID), eta0(eID) )
          end do
        end do
        
        ! Compute alpha
        do eID=1, nElems
          alpha(eID) = max(TimeRelFactor*alpha(eID), 1.0 / (1.0 + exp(-sharpness * (eta(eID) - threshold)/threshold )) )
        end do
        
      case(11)   ! Shock indicator capped in interval [0,1]
        ! Compute first shock indicator
        eta = Shock_Indicator(1) % compute(U)
        
        ! Compute other shock indicators
        do ind=2, ShockIndicatorNum
          eta0 = Shock_Indicator(ind) % compute(U)
          do eID=1, nElems
            eta(eID) = max(eta(eID), eta0(eID) )
          end do
        end do
        
        ! Compute alpha
        do eID=1, nElems
          alpha(eID) = max(TimeRelFactor*alpha(eID), 0.0, min(eta(eID), 1.0) )
        end do
        
      case(2)   ! Random indicator
        do eID=1, nElems
          call RANDOM_NUMBER(alpha(eID))
        end do
      case(20)   ! Random indicator, fixed over time (using shockBlendCoef as a switch)
        if(shockBlendCoef.GT.-1.)THEN
          call RANDOM_NUMBER(alpha(:))
          alpha=alpha*alpha_max
          shockBlendCoef=-2.
        end if
        
      case(3)   ! Fixed blending coefficients
        alpha = ShockBlendCoef
    end select
    
!   Impose alpha_max and alpha_min
!   ------------------------------
    where (alpha < alpha_min)
      alpha = 0.0
    elsewhere (alpha >= alpha_max)
      alpha = alpha_max
    end where
    
!   For Gauss shock capturing, we reset the Ut_DGGauss and Ut_FVGauss array
!   -----------------------------------------------------------------------
#if (PP_NodeType==1)
    Ut_DGGauss = 0.0
#if FV_BLENDSURFACE
    Ut_FVGauss = 0.0
#endif /*FV_BLENDSURFACE*/
#endif /*(PP_NodeType==1)*/
    
!   Start first space propagation sweep (MPI-optimized)
!   ---------------------------------------------------
    if (SpacePropSweeps > 0) call ProlongBlendingCoeffToFaces()
    
  end subroutine CalcBlendingCoefficient
!===================================================================================================================================
!> Finalizes the NFVSE module
!===================================================================================================================================
  subroutine FinalizeNFVSE()
    use MOD_NFVSE_Vars, only: SubCellMetrics, sWGP, Compute_FVFluxes, alpha, alpha_Master, alpha_Slave
#if NFVSE_CORR
    use MOD_IDP       , only: Finalize_IDP
#endif /*NFVSE_CORR*/
    use MOD_NFVSE_Vars, only: U_ext, sdxR,sdxL,rR,rL
#if MPI
    use MOD_NFVSE_Vars, only: MPIRequest_alpha
#endif /*MPI*/
    implicit none
    
    SDEALLOCATE (alpha)
    SDEALLOCATE (alpha_Master)
    SDEALLOCATE (alpha_Slave)
    SDEALLOCATE (SubCellMetrics)
    SDEALLOCATE (sWGP)
#if MPI
    SDEALLOCATE (MPIRequest_alpha)
    SDEALLOCATE (MPIRequest_Umaster)
#endif /*MPI*/
#if NFVSE_CORR
    call Finalize_IDP()
#endif /*NFVSE_CORR*/
    
    Compute_FVFluxes => null()
    
    ! Reconstruction
    SDEALLOCATE ( U_ext )
    SDEALLOCATE ( sdxR  )
    SDEALLOCATE ( sdxL  )
    SDEALLOCATE ( rR    )
    SDEALLOCATE ( rL    )
    
    ! Gauss methods
#if (PP_NodeType==1)
    SDEALLOCATE (Ut_DGGauss)
#if FV_BLENDSURFACE
    SDEALLOCATE (Ut_FVGauss)
    SDEALLOCATE (U_master_FV)
    SDEALLOCATE (U_slave_FV)
    SDEALLOCATE (Flux_master_FV)
    SDEALLOCATE (Flux_slave_FV)
#if MPI
    SDEALLOCATE (MPIRequest_U_FV)
    SDEALLOCATE (MPIRequest_Flux_FV)
#endif /*MPI*/
#endif /*FV_BLENDSURFACE*/
#endif /*(PP_NodeType==1)*/
    
  end subroutine FinalizeNFVSE
#endif /*SHOCK_NFVSE*/
end module MOD_NFVSE
