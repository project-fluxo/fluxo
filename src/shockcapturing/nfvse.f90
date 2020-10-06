!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
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
module MOD_NFVSE
#if SHOCK_NFVSE
  use MOD_PreProc
  implicit none
!==================================================================================================================================
!> Computes the spatial operator using Native Finite Volume Sub-Elements (NFVSE)
!> Based on: Hennemann et al. (2020). "A provably entropy stable subcell shock capturing approach for high order split form 
!>                                     DG for the compressible Euler Equations"
!> Attention 1: The current strategy is to perform the operations on "inner faces". This can be improved, but the Riemann solver 
!>              routines have to be adapted
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
!(2) FV BC idx:     -1 0     1       2     3 4
!                    └─┴─────┴───────┴─────┴─┘
!(3) FV subcell idx:  0   1      2      3   4
!
!   * (1) matches (3) for consistency
!===================================================================================================================================
  private
  public :: VolInt_NFVSE, InitNFVSE, FinalizeNFVSE
  public :: DefineParametersNFVSE
  public :: InitNFVSEAfterAdaptation
#if NFVSE_CORR
  public :: Apply_NFVSE_Correction
#endif /*NFVSE_CORR*/
contains
!===================================================================================================================================
!> Defines parameters for NFVSE (to be called in shockcapturing.f90)
!===================================================================================================================================
  subroutine DefineParametersNFVSE()
    use MOD_ReadInTools,  only: prms
    implicit none
    
    CALL prms%CreateIntOption     (   "SubFVMethod",  " Specifies subcell Finite-Volume method to be used "//&
                                              "  1: 1st order FV"//&
                                              "  2: 2nd order TVD method (not ES)"//&
                                              "  3: '2nd' order TVD-ES method"&
                                             ,"1")
    call prms%CreateRealOption    (      "alpha_max", "Maximum value for the blending coefficient", "0.5")
    call prms%CreateRealOption    (      "alpha_min", "Minimum value for the blending coefficient (below this, alpha=0)", "0.01")
    call prms%CreateRealOption    ("SpacePropFactor", "Space propagation factor", "0.5")
    call prms%CreateIntOption     ("SpacePropSweeps", "Number of space propagation sweeps (MPI-optimized only for 0 or 1)", "1")
    call prms%CreateRealOption    (  "TimeRelFactor", "Time relaxation factor", "0.0")
    call prms%CreateLogicalOption("ReconsBoundaries", "Use a reconstruction procedure on boundary subcells.. Only for SubFVMethod=1/2","F")
   
  end subroutine DefineParametersNFVSE
!===================================================================================================================================
!> Initializes the NFVSE module
!===================================================================================================================================
  subroutine InitNFVSE()
    USE MOD_Globals
    use MOD_NFVSE_Vars         , only: SubCellMetrics, sWGP, MPIRequest_alpha, Fsafe, Fblen, Compute_FVFluxes, SubFVMethod, ReconsBoundaries, MPIRequest_Umaster
    use MOD_NFVSE_Vars         , only: SpacePropFactor, SpacePropSweeps, TimeRelFactor
    use MOD_ReadInTools        , only: GETINT, GETREAL, GETLOGICAL
    use MOD_NFVSE_Vars         , only: U_ext, sdxR, sdxL, rL, rR
    use MOD_Mesh_Vars          , only: nElems, Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
    use MOD_Interpolation_Vars , only: wGP, xGP
    use MOD_ShockCapturing_Vars, only: alpha_old, alpha_max, alpha_min
    use MOD_DG_Vars            , only: D
#if MPI
    use MOD_MPI_Vars           , only: nNbProcs
#endif /*MPI*/
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    integer :: iElem
    integer :: i,j,k,l,m         !DOF counters
    real :: Metrics_fCont(3,0:PP_N,0:PP_N,0:PP_N) ! Container for the (reshaped) xi metrics
    real :: Metrics_gCont(3,0:PP_N,0:PP_N,0:PP_N) ! Container for the (reshaped) eta metrics
    real :: Metrics_hCont(3,0:PP_N,0:PP_N,0:PP_N) ! Container for the (reshaped) zeta metrics
    real, parameter :: half = 0.5d0
    real :: sumWm1
    !--------------------------------------------------------------------------------------------------------------------------------
     SDEALLOCATE(SubCellMetrics)
    ! Get parameters
    ! --------------
    
    alpha_max        = GETREAL   ('alpha_max','0.5')
    alpha_min        = GETREAL   ('alpha_min','0.01')
    SubFVMethod      = GETINT    ('SubFVMethod','1')
    ReconsBoundaries = GETLOGICAL('ReconsBoundaries','F')
    SpacePropFactor  = GETREAL   ('SpacePropFactor','0.5')
    SpacePropSweeps  = GETINT    ('SpacePropSweeps','1')
    TimeRelFactor    = GETREAL   ('TimeRelFactor'  ,'0.0')
    
    ! Initialize everything
    ! ---------------------
    
    select case(SubFVMethod)
      case(1) ; Compute_FVFluxes => Compute_FVFluxes_1st_Order
      case(2) ; Compute_FVFluxes => Compute_FVFluxes_2ndOrder
      case(3) ; Compute_FVFluxes => Compute_FVFluxes_TVD2ES
    end select
    
    ! Allocate storage
    allocate ( sWGP(0:PP_N) )
    
    allocate ( SubCellMetrics(nElems) )
    call SubCellMetrics % construct(PP_N)
    
    call ComputeSubcellMetrics()
    
    ! Compute the inverse of the quadrature weights (sub-cell dimensions)
    sWGP = 1.d0 / wGP
    
    ! Compute variables for 2nd order FV
    ! **********************************
    select case(SubFVMethod)
    case(2,3)
      allocate ( U_ext(1:PP_nVar,0:PP_N,0:PP_N,6,nElems) )
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
    end select
#if NFVSE_CORR
    allocate ( Fsafe(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) )
    allocate ( Fblen(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) )
    allocate ( alpha_old(nElems) )
    Fsafe = 0.
    Fblen = 0.
    alpha_old = 0.
#endif /*NFVSE_CORR*/
    
#if MPI
    allocate(MPIRequest_alpha(nNbProcs,4)    ) ! 1: send slave, 2: send master, 3: receive slave, 4, receive master
    
    if (ReconsBoundaries) then
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
      
      do i=0, PP_N-1
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
      
      do i=0, PP_N-1
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
      do i=0, PP_N-1
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
!> ATTENTION: The subcell metrics are always recomputed, as the metrics of the high-order DG elements
!===================================================================================================================================
  subroutine InitNFVSEAfterAdaptation(ChangeElem,nElemsOld)
    USE MOD_Globals
    use MOD_NFVSE_Vars         , only: SubCellMetrics, Fsafe, Fblen
    use MOD_ShockCapturing_Vars, only: alpha_old
    use MOD_Mesh_Vars          , only: nElems
#if MPI
    use MOD_MPI_Vars           , only: nNbProcs
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha, MPIRequest_Umaster, ReconsBoundaries
#endif /*MPI*/
    implicit none
    !-arguments-----------------------------------
    integer, intent(in) :: ChangeElem(8,nElems)
    integer, intent(in) :: nElemsOld
    !---------------------------------------------
    
    ! Reallocate storage if needed
    if (nElems /= nElemsOld) then
      SDEALLOCATE(SubCellMetrics)
      allocate ( SubCellMetrics(nElems) )
      call SubCellMetrics % construct(PP_N)
#if NFVSE_CORR
      SDEALLOCATE(Fsafe)
      SDEALLOCATE(Fblen)
      SDEALLOCATE(alpha_old)
      allocate ( Fsafe(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) )
      allocate ( Fblen(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) )
      allocate ( alpha_old(nElems) )
      Fsafe = 0.
      Fblen = 0.
      alpha_old = 0.
#endif /*NFVSE_CORR*/
    end if
    
    ! Compute subcell metrics
    call ComputeSubcellMetrics()
    
    ! Reallocate MPI variables
#if MPI
    SDEALLOCATE(MPIRequest_alpha)
    allocate(MPIRequest_alpha(nNbProcs,4)    ) ! 1: send slave, 2: send master, 3: receive slave, 4, receive master
    
    if (ReconsBoundaries) then
      SDEALLOCATE(MPIRequest_Umaster)
      allocate(MPIRequest_Umaster(nNbProcs,2)) ! 1: send master, 2: receive master
    end if
#endif
  end subroutine InitNFVSEAfterAdaptation
!===================================================================================================================================
!> Computes the "volume integral": Spatial contribution to Ut by the subcell finite volumes
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: This routine has to be called after VolInt_adv_SplitForm, since here Ut is updated with the finite volume contribution
!===================================================================================================================================
  subroutine VolInt_NFVSE(Ut)
  !----------------------------------------------------------------------------------------------------------------------------------
    ! Modules
    use MOD_PreProc
    use MOD_DG_Vars            , only: U
    use MOD_Mesh_Vars          , only: nElems
    use MOD_NFVSE_Vars         , only: SubCellMetrics, sWGP, Fsafe, Fblen, Compute_FVFluxes, ReconsBoundaries, MPIRequest_Umaster, SpacePropSweeps
    use MOD_ShockCapturing_Vars, only: alpha, alpha_max
    use MOD_Basis              , only: ALMOSTEQUAL
    use MOD_NFVSE_MPI          , only: PropagateBlendingCoeff
#if MPI
    use MOD_MPI                , only: FinishExchangeMPIData
    use MOD_MPI_Vars           , only: nNbProcs
#endif /*MPI*/
    ! IMPLICIT VARIABLE HANDLING
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    real,intent(inout)                              :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
    !< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
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
    integer                                         :: i,j,k,iElem
    !===============================================================================================================================
    
    if (SpacePropSweeps > 0) call PropagateBlendingCoeff()
    
    !if reconstruction:
    if (ReconsBoundaries) then
#if MPI
      call FinishExchangeMPIData(2*nNbProcs,MPIRequest_Umaster) 
#endif /*MPI*/
      call Get_externalU()
    end if
    
    do iElem=1,nElems
#if !defined(NFVSE_CORR)
      if ( ALMOSTEQUAL(alpha(iElem),0.d0) ) cycle
#endif /*NFVSE_CORR*/
      
!     Compute the finite volume fluxes
!     --------------------------------
      call Compute_FVFluxes( U(:,:,:,:,iElem), ftilde , gtilde , htilde , &
#if NONCONS
                                               ftildeR, gtildeR, htildeR, &
#endif /*NONCONS*/
                                               SubCellMetrics(iElem), iElem )
      
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
               
#if NFVSE_CORR
        Fsafe(:,i,j,k,iElem) = (1. - alpha_max   ) * Ut(:,i,j,k,iElem) +  alpha_max    * F_FV
        Fblen(:,i,j,k,iElem) = (1. - alpha(iElem)) * Ut(:,i,j,k,iElem) +  alpha(iElem) * F_FV
        Ut(:,i,j,k,iElem)    = Fblen(:,i,j,k,iElem)
#else
        Ut   (:,i,j,k,iElem) = (1. - alpha(iElem)) * Ut(:,i,j,k,iElem) +  alpha(iElem) * F_FV
#endif /*NFVSE_CORR*/
      end do       ; end do       ; end do ! i,j,k
      
    end do ! iElem
  end subroutine VolInt_NFVSE
  
!============================================================================================================================
!> Get Pressure
!============================================================================================================================
  pure subroutine GetPressure(U,p)
#ifdef mhd
    use MOD_Equation_Vars      , only: s2mu_0
#endif /*mhd*/
    use MOD_Equation_Vars      , only: KappaM1
    implicit none
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: p
    
    p = KappaM1*(U(5)-0.5*(SUM(U(2:4)*U(2:4))/U(1)))
#ifdef mhd
    p = p - KappaM1*s2mu_0*SUM(U(6:8)*U(6:8))
#ifdef PP_GLM
    p = p - 0.5*KappaM1*U(9)*U(9)
#endif /*PP_GLM*/
#endif /*mhd*/
  end subroutine GetPressure
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a FV consistent flux
!===================================================================================================================================
  subroutine Compute_FVFluxes_1st_Order (U, F , G , H , &
#if NONCONS
                                            FR, GR, HR, &
#endif /*NONCONS*/
                                            sCM, iElem )
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: SubCellMetrics_t
    use MOD_Riemann   , only: AdvRiemann
#if NONCONS
    USE MOD_Riemann   , only: AddWholeNonConsFlux, AddInnerNonConsFlux
    use MOD_Mesh_Vars , only: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
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
    integer  :: i,j,k
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
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(F_ (:,:,:,-1), U_(:,:,:,0), Metrics_fTilde(:,0,:,:,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
                   sCM % xi   % nv(:,:,:,i),sCM % xi   % t1(:,:,:,i), sCM % xi   % t2(:,:,:,i))
      
#if NONCONS
      ! Copy conservative part
      FR_(:,:,:,i) = F_(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(F_(:,:,:,i), &
                          U_(:,:,:,i+1), U_(:,:,:,i),&
                          sCM % xi   % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(FR_(:,:,:,i), &
                          U_(:,:,:,i  ), U_(:,:,:,i+1),&
                          sCM % xi   % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        FR_(:,j,k,i) = FR_(:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
    end do ! i (xi planes)
    
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(FR_(:,:,:,PP_N), U_(:,:,:,PP_N), Metrics_fTilde(:,PP_N,:,:,iElem))
#endif /*NONCONS*/
    
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
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(F_ (:,:,:,-1), U_(:,:,:,0), Metrics_gTilde(:,:,0,:,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
                   sCM % eta  % nv(:,:,:,i),sCM % eta  % t1(:,:,:,i), sCM % eta  % t2(:,:,:,i))
      
#if NONCONS
      ! Copy conservative part
      FR_(:,:,:,i) = F_(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(F_(:,:,:,i), &
                          U_(:,:,:,i+1), U_(:,:,:,i  ),&
                          sCM % eta  % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(FR_(:,:,:,i), &
                          U_(:,:,:,i  ), U_(:,:,:,i+1),&
                          sCM % eta  % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        FR_(:,j,k,i) = FR_(:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
    end do ! i (eta planes)
   
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(FR_(:,:,:,PP_N), U_(:,:,:,PP_N), Metrics_gTilde(:,:,PP_N,:,iElem))
#endif /*NONCONS*/
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    G  = reshape(F_ , shape(G ), order = [1,2,4,3])
#if NONCONS
    GR = reshape(FR_, shape(GR), order = [1,2,4,3])
#endif /*NONCONS*/

!   ***********    
!   Zeta-planes
!   ***********

!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(H(:,:,:,-1), U(:,:,:,0), Metrics_hTilde(:,:,:,0,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(H(:,:,:,i),U(:,:,:,i),U(:,:,:,i+1), &
                   sCM % zeta % nv(:,:,:,i),sCM % zeta % t1(:,:,:,i), sCM % zeta % t2(:,:,:,i))
      
#if NONCONS
      ! Copy conservative part
      HR(:,:,:,i) = H(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(H(:,:,:,i), &
                          U(:,:,:,i+1), U(:,:,:,i  ),&
                          sCM % zeta % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(HR(:,:,:,i), &
                          U(:,:,:,i  ), U(:,:,:,i+1),&
                          sCM % zeta % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        HR(:,j,k,i) = HR(:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        H (:,j,k,i) = H (:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
    end do ! i (zeta planes)
    
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(HR(:,:,:,PP_N), U(:,:,:,PP_N), Metrics_hTilde(:,:,:,PP_N,iElem))
#endif /*NONCONS*/
    
  end subroutine Compute_FVFluxes_1st_Order
  
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a FV consistent flux
!===================================================================================================================================
  subroutine Compute_FVFluxes_2ndOrder(U, F , G , H , &
#if NONCONS
                                 FR, GR, HR, &
#endif /*NONCONS*/
                                              sCM, iElem )
    use MOD_PreProc
    use MOD_Riemann   , only: AdvRiemann
    use MOD_NFVSE_Vars, only: SubCellMetrics_t, sdxR, sdxL, rL, rR
    use MOD_NFVSE_Vars, only: U_ext, ReconsBoundaries
    use MOD_DG_Vars   , only: nTotal_vol
    use MOD_Equation_Vars, only: ConsToPrimVec,PrimToConsVec
    use MOD_Interpolation_Vars , only: wGP
#if NONCONS
    USE MOD_Riemann   , only: AddWholeNonConsFlux, AddInnerNonConsFlux
    use MOD_Mesh_Vars , only: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
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
    integer  :: i,j,k
    real :: U_ (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)
    real :: UL   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the left
    real :: UR   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the right
    real :: prim (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables
    real :: prim_ext(PP_nVar,0:PP_N,0:PP_N, 6)  ! Primitive variables
    real :: prim_(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables after reshape
    real :: primL(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables left
    real :: primR(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables right
    real :: sigma(PP_nVar,0:PP_N,0:PP_N)
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
    
    call ConsToPrimVec(nTotal_vol,prim,U)
    
    if (ReconsBoundaries) call ConsToPrimVec(6*(PP_N+1)**2,prim_ext,U_ext(:,:,:,:,iElem) )
    
!   *********
!   Xi-planes
!   *********
    F_  = 0.0
    prim_ = reshape(prim , shape(prim_), order = [1,4,2,3])
    U_    = reshape(U    , shape(U_)   , order = [1,4,2,3])
    primL = 0.0
    primR = 0.0
    
!   Do the solution reconstruction
!   ******************************
    
    ! First DOF
    !----------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,1) - prim_(:,:,:,0)),sdxL(1)*(prim_(:,:,:,1) - prim_ext(:,:,:,5)))
      primR(:,:,:,0) = prim_(:,:,:,0) + sigma*wGP(0)
    else
      ! Central scheme
      primR(:,:,:,0) = prim_(:,:,:,0) + (prim_(:,:,:,1) - prim_(:,:,:,0))*sdxL(1)*wGP(0)
    end if
    
    ! Middle DOFs
    ! -----------
    do i=1, PP_N-1
      sigma = minmod(sdxR(i)*(prim_(:,:,:,i+1)-prim_(:,:,:,i)),sdxL(i)*(prim_(:,:,:,i)-prim_(:,:,:,i-1)))
      
      primR(:,:,:,i) = prim_(:,:,:,i) + sigma * rR(i)
      primL(:,:,:,i) = prim_(:,:,:,i) + sigma * rL(i)
    end do
    
    ! Last DOF
    ! --------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1)),sdxL(1)*(prim_ext(:,:,:,3)-prim_(:,:,:,PP_N-1)))
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - sigma*wGP(0)
    else
      ! Central scheme
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - (prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
    end if
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(F_ (:,:,:,-1), U_(:,:,:,0), Metrics_fTilde(:,0,:,:,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(F_(:,:,:,i),UR(:,:,:,i),UL(:,:,:,i+1), &
                   sCM % xi   % nv(:,:,:,i),sCM % xi   % t1(:,:,:,i), sCM % xi   % t2(:,:,:,i))
      
#if NONCONS
      ! Copy conservative part
      FR_(:,:,:,i) = F_(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(F_(:,:,:,i), &
                          U_(:,:,:,i+1), U_(:,:,:,i),&
                          sCM % xi   % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(FR_(:,:,:,i), &
                          U_(:,:,:,i  ), U_(:,:,:,i+1),&
                          sCM % xi   % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        FR_(:,j,k,i) = FR_(:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
    end do ! i (xi planes)
    
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(FR_(:,:,:,PP_N), U_(:,:,:,PP_N), Metrics_fTilde(:,PP_N,:,:,iElem))
#endif /*NONCONS*/
    
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
    primL = 0.0
    primR = 0.0
    
!   Do the solution reconstruction
!   ******************************
    
    ! First DOF
    !----------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,1) - prim_(:,:,:,0)),sdxL(1)*(prim_(:,:,:,1) - prim_ext(:,:,:,2)))
      primR(:,:,:,0) = prim_(:,:,:,0) + sigma*wGP(0)
    else
      ! Central scheme
      primR(:,:,:,0) = prim_(:,:,:,0) + (prim_(:,:,:,1) - prim_(:,:,:,0))*sdxL(1)*wGP(0)
    end if
    
    ! Middle DOFs
    ! -----------
    do i=1, PP_N-1
      sigma = minmod(sdxR(i)*(prim_(:,:,:,i+1)-prim_(:,:,:,i)),sdxL(i)*(prim_(:,:,:,i)-prim_(:,:,:,i-1)))
      
      primR(:,:,:,i) = prim_(:,:,:,i) + sigma * rR(i)
      primL(:,:,:,i) = prim_(:,:,:,i) + sigma * rL(i)
    end do
    
    ! Last DOF
    ! --------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1)),sdxL(1)*(prim_ext(:,:,:,4)-prim_(:,:,:,PP_N-1)))
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - sigma*wGP(0)
    else
      ! Central scheme
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - (prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
    end if
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(F_ (:,:,:,-1), U_(:,:,:,0), Metrics_gTilde(:,:,0,:,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(F_(:,:,:,i),UR(:,:,:,i),UL(:,:,:,i+1), &
                   sCM % eta  % nv(:,:,:,i),sCM % eta  % t1(:,:,:,i), sCM % eta  % t2(:,:,:,i))
      
#if NONCONS
      ! Copy conservative part
      FR_(:,:,:,i) = F_(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(F_(:,:,:,i), &
                          U_(:,:,:,i+1), U_(:,:,:,i  ),&
                          sCM % eta  % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(FR_(:,:,:,i), &
                          U_(:,:,:,i  ), U_(:,:,:,i+1),&
                          sCM % eta  % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        FR_(:,j,k,i) = FR_(:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
    end do ! i (eta planes)
   
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(FR_(:,:,:,PP_N), U_(:,:,:,PP_N), Metrics_gTilde(:,:,PP_N,:,iElem))
#endif /*NONCONS*/
    
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
!   ******************************
    
    ! First DOF
    !----------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim(:,:,:,1) - prim(:,:,:,0)),sdxL(1)*(prim(:,:,:,1) - prim_ext(:,:,:,1)))
      primR(:,:,:,0) = prim(:,:,:,0) + sigma*wGP(0)
    else
      ! Central scheme
      primR(:,:,:,0) = prim(:,:,:,0) + (prim(:,:,:,1) - prim(:,:,:,0))*sdxL(1)*wGP(0)
    end if
    
    ! Middle DOFs
    ! -----------
    do i=1, PP_N-1
      sigma = minmod(sdxR(i)*(prim(:,:,:,i+1)-prim(:,:,:,i)),sdxL(i)*(prim(:,:,:,i)-prim(:,:,:,i-1)))
      
      primR(:,:,:,i) = prim(:,:,:,i) + sigma * rR(i)
      primL(:,:,:,i) = prim(:,:,:,i) + sigma * rL(i)
    end do
    
    ! Last DOF
    ! --------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim(:,:,:,PP_N) - prim(:,:,:,PP_N-1)),sdxL(1)*(prim_ext(:,:,:,6)-prim(:,:,:,PP_N-1)))
      primL(:,:,:,PP_N) = prim(:,:,:,PP_N) - sigma*wGP(0)
    else
      ! Central scheme
      primL(:,:,:,PP_N) = prim(:,:,:,PP_N) - (prim(:,:,:,PP_N) - prim(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
    end if
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(H(:,:,:,-1), U(:,:,:,0), Metrics_hTilde(:,:,:,0,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(H(:,:,:,i),UR(:,:,:,i),UL(:,:,:,i+1), &
                   sCM % zeta % nv(:,:,:,i),sCM % zeta % t1(:,:,:,i), sCM % zeta % t2(:,:,:,i))
      
#if NONCONS
      ! Copy conservative part
      HR(:,:,:,i) = H(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(H(:,:,:,i), &
                          U(:,:,:,i+1), U(:,:,:,i  ),&
                          sCM % zeta % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(HR(:,:,:,i), &
                          U(:,:,:,i  ), U(:,:,:,i+1),&
                          sCM % zeta % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        HR(:,j,k,i) = HR(:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        H (:,j,k,i) = H (:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
    end do ! i (zeta planes)
    
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(HR(:,:,:,PP_N), U(:,:,:,PP_N), Metrics_hTilde(:,:,:,PP_N,iElem))
#endif /*NONCONS*/
    
  end subroutine Compute_FVFluxes_2ndOrder
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a FV consistent flux
!===================================================================================================================================
  subroutine Compute_FVFluxes_TVD2ES(U, F , G , H , &
#if NONCONS
                                        FR, GR, HR, &
#endif /*NONCONS*/
                                              sCM, iElem )
    use MOD_PreProc
    use MOD_Riemann   , only: AdvRiemannRecons
    use MOD_NFVSE_Vars, only: SubCellMetrics_t, sdxR, sdxL, rL, rR
    use MOD_NFVSE_Vars, only: U_ext, ReconsBoundaries
    use MOD_Interpolation_Vars , only: wGP
    use MOD_Equation_Vars , only: ConsToPrimVec, PrimToConsVec
    use MOD_DG_Vars           , only: nTotal_vol
#if NONCONS
    USE MOD_Riemann   , only: AddWholeNonConsFlux, AddInnerNonConsFlux
    use MOD_Mesh_Vars , only: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
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
    integer  :: i,j,k
    real :: U_   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)
    real :: UL   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the left
    real :: UR   (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed solution on the right
    real :: prim (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables
    real :: prim_ext(PP_nVar,0:PP_N,0:PP_N, 6)  ! Primitive variables
    real :: prim_(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Primitive variables after reshape
    real :: primL(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables left
    real :: primR(PP_nVar,0:PP_N,0:PP_N, 0:PP_N)  ! Reconstructed Primitive variables right
    real :: sigma(PP_nVar,0:PP_N,0:PP_N)
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
    
    call ConsToPrimVec(nTotal_vol,prim,U)
    if (ReconsBoundaries) call ConsToPrimVec(6*(PP_N+1)**2,prim_ext,U_ext(:,:,:,:,iElem) )
    
!   *********
!   Xi-planes
!   *********
    F_  = 0.0
    prim_ = reshape(prim , shape(prim_), order = [1,4,2,3])
    U_    = reshape(U    , shape(U_)   , order = [1,4,2,3])
    primL = 0.0
    primR = 0.0
    
!   Do the solution reconstruction
!   ******************************
    
    ! First DOF
    !----------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,1) - prim_(:,:,:,0)),sdxL(1)*(prim_(:,:,:,1) - prim_ext(:,:,:,5)))
      primR(:,:,:,0) = prim_(:,:,:,0) + sigma*wGP(0)
    else
      ! Central scheme
      primR(:,:,:,0) = prim_(:,:,:,0) + (prim_(:,:,:,1) - prim_(:,:,:,0))*sdxL(1)*wGP(0)
    end if
    
    
    ! Middle DOFs
    ! -----------
    do i=1, PP_N-1
      sigma = minmod(sdxR(i)*(prim_(:,:,:,i+1)-prim_(:,:,:,i)),sdxL(i)*(prim_(:,:,:,i)-prim_(:,:,:,i-1)))
      
      primR(:,:,:,i) = prim_(:,:,:,i) + sigma * rR(i)
      primL(:,:,:,i) = prim_(:,:,:,i) + sigma * rL(i)
    end do
    
    ! Last DOF
    ! --------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1)),sdxL(1)*(prim_ext(:,:,:,3)-prim_(:,:,:,PP_N-1)))
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - sigma*wGP(0)
    else
      ! Central scheme
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - (prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
    end if
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)

!   Entropy fix: Fall to first order if entropy condition is not fulfilled  
!   **********************************************************************
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(F_ (:,:,:,-1), U_(:,:,:,0), Metrics_fTilde(:,0,:,:,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemannRecons(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1),UR(:,:,:,i),UL(:,:,:,i+1), &
                   sCM % xi   % nv(:,:,:,i),sCM % xi   % t1(:,:,:,i), sCM % xi   % t2(:,:,:,i) )
#if NONCONS
      ! Copy conservative part
      FR_(:,:,:,i) = F_(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(F_(:,:,:,i), &
                          U_(:,:,:,i+1), U_(:,:,:,i),&
                          sCM % xi   % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(FR_(:,:,:,i), &
                          U_(:,:,:,i  ), U_(:,:,:,i+1),&
                          sCM % xi   % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        FR_(:,j,k,i) = FR_(:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
    end do ! i (xi planes)
    
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(FR_(:,:,:,PP_N), U_(:,:,:,PP_N), Metrics_fTilde(:,PP_N,:,:,iElem))
#endif /*NONCONS*/
    
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
    primL = 0.0
    primR = 0.0
    
!   Do the solution reconstruction
!   ******************************
    
    ! First DOF
    !----------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,1) - prim_(:,:,:,0)),sdxL(1)*(prim_(:,:,:,1) - prim_ext(:,:,:,2)))
      primR(:,:,:,0) = prim_(:,:,:,0) + sigma*wGP(0)
    else
      ! Central scheme
      primR(:,:,:,0) = prim_(:,:,:,0) + (prim_(:,:,:,1) - prim_(:,:,:,0))*sdxL(1)*wGP(0)
    end if
    
    ! Middle DOFs
    ! -----------
    do i=1, PP_N-1
      sigma = minmod(sdxR(i)*(prim_(:,:,:,i+1)-prim_(:,:,:,i)),sdxL(i)*(prim_(:,:,:,i)-prim_(:,:,:,i-1)))
      
      primR(:,:,:,i) = prim_(:,:,:,i) + sigma * rR(i)
      primL(:,:,:,i) = prim_(:,:,:,i) + sigma * rL(i)
    end do
    
    ! Last DOF
    ! --------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1)),sdxL(1)*(prim_ext(:,:,:,4)-prim_(:,:,:,PP_N-1)))
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - sigma*wGP(0)
    else
      ! Central scheme
      primL(:,:,:,PP_N) = prim_(:,:,:,PP_N) - (prim_(:,:,:,PP_N) - prim_(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
    end if
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(F_ (:,:,:,-1), U_(:,:,:,0), Metrics_gTilde(:,:,0,:,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemannRecons(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1),UR(:,:,:,i),UL(:,:,:,i+1), &
                   sCM % eta  % nv(:,:,:,i),sCM % eta  % t1(:,:,:,i), sCM % eta  % t2(:,:,:,i))
#if NONCONS
      ! Copy conservative part
      FR_(:,:,:,i) = F_(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(F_(:,:,:,i), &
                          U_(:,:,:,i+1), U_(:,:,:,i  ),&
                          sCM % eta  % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(FR_(:,:,:,i), &
                          U_(:,:,:,i  ), U_(:,:,:,i+1),&
                          sCM % eta  % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        FR_(:,j,k,i) = FR_(:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
    end do ! i (eta planes)
   
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(FR_(:,:,:,PP_N), U_(:,:,:,PP_N), Metrics_gTilde(:,:,PP_N,:,iElem))
#endif /*NONCONS*/
    
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
!   ******************************
    
    ! First DOF
    !----------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim(:,:,:,1) - prim(:,:,:,0)),sdxL(1)*(prim(:,:,:,1) - prim_ext(:,:,:,1)))
      primR(:,:,:,0) = prim(:,:,:,0) + sigma*wGP(0)
    else
      ! Central scheme
      primR(:,:,:,0) = prim(:,:,:,0) + (prim(:,:,:,1) - prim(:,:,:,0))*sdxL(1)*wGP(0)
    end if
    
    ! Middle DOFs
    ! -----------
    do i=1, PP_N-1
      sigma = minmod(sdxR(i)*(prim(:,:,:,i+1)-prim(:,:,:,i)),sdxL(i)*(prim(:,:,:,i)-prim(:,:,:,i-1)))
      
      primR(:,:,:,i) = prim(:,:,:,i) + sigma * rR(i)
      primL(:,:,:,i) = prim(:,:,:,i) + sigma * rL(i)
    end do
    
    ! Last DOF
    ! --------
    
    if (ReconsBoundaries) then
      ! v3
      sigma = minmod(sdxL(1)*(prim(:,:,:,PP_N) - prim(:,:,:,PP_N-1)),sdxL(1)*(prim_ext(:,:,:,6)-prim(:,:,:,PP_N-1)))
      primL(:,:,:,PP_N) = prim(:,:,:,PP_N) - sigma*wGP(0)
    else
      ! Central scheme
      primL(:,:,:,PP_N) = prim(:,:,:,PP_N) - (prim(:,:,:,PP_N) - prim(:,:,:,PP_N-1))*sdxL(1)*wGP(PP_N)
    end if
    
    call PrimToConsVec(nTotal_vol,primL,UL)
    call PrimToConsVec(nTotal_vol,primR,UR)
    
!   Fill left boundary if non-conservative terms are present
!   --------------------------------------------------------
#if NONCONS
    FR_ = 0.0
    CALL AddInnerNonConsFlux(H(:,:,:,-1), U(:,:,:,0), Metrics_hTilde(:,:,:,0,iElem))
#endif /*NONCONS*/
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemannRecons(H(:,:,:,i),U(:,:,:,i),U(:,:,:,i+1),UR(:,:,:,i),UL(:,:,:,i+1), &
                   sCM % zeta % nv(:,:,:,i),sCM % zeta % t1(:,:,:,i), sCM % zeta % t2(:,:,:,i))
#if NONCONS
      ! Copy conservative part
      HR(:,:,:,i) = H(:,:,:,i)
      
      ! Add nonconservative fluxes
      CALL AddWholeNonConsFlux(H(:,:,:,i), &
                          U(:,:,:,i+1), U(:,:,:,i  ),&
                          sCM % zeta % nv(:,:,:,i))
      CALL AddWholeNonConsFlux(HR(:,:,:,i), &
                          U(:,:,:,i  ), U(:,:,:,i+1),&
                          sCM % zeta % nv(:,:,:,i))
      
      ! Scale right flux
      do k=0, PP_N ; do j=0, PP_N
        HR(:,j,k,i) = HR(:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
#endif /*NONCONS*/
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        H (:,j,k,i) = H (:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
    end do ! i (zeta planes)
    
!   Fill right boundary if non-conservative terms are present
!   ---------------------------------------------------------
#if NONCONS
    CALL AddInnerNonConsFlux(HR(:,:,:,PP_N), U(:,:,:,PP_N), Metrics_hTilde(:,:,:,PP_N,iElem))
#endif /*NONCONS*/
    
  end subroutine Compute_FVFluxes_TVD2ES
!  
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
  
!
!> Gets outer solution for the TVD reconstruction
!> ATTENTION 1: Must be called after FinishExchangeMPIData
!> ATTENTION 2: Mortar faces are not considered (big TODO)
! TODO: We need to send the U_master tooooooooo
!  
  subroutine Get_externalU()
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: U_ext
    use MOD_DG_Vars   , only: U, U_master, U_slave
    use MOD_Mesh_Vars , only: nBCSides, SideToElem, S2V, firstInnerSide, lastMPISide_YOUR
    implicit none
    !-local-variables-----------------------------------------
    integer :: SideID, ElemID, locSide, p, q, ijk(3), flip, nbElemID, nblocSide, nbFlip
    !---------------------------------------------------------
    
    ! First assign the inner value on boundaries (we don't use an "external state")
    do SideID=1, nBCSides
      ElemID  = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side
      locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
      flip    = 0 ! flip is 0 in master faces!
      
      DO q=0,PP_N; DO p=0,PP_N
        ijk(:)=S2V(:,0,p,q,flip,locSide)
        U_ext(:,p,q,locSide,ElemID)=U(:,ijk(1),ijk(2),ijk(3),ElemID)
      END DO; END DO !p,q=0,PP_N
    end do
    
    ! TODO: Include routines for mortar faces here....
    
    ! Now do the rest
    do SideID=firstInnerSide, lastMPISide_YOUR
      
!     Master side
!     -----------
      ElemID    = SideToElem(S2E_ELEM_ID,SideID) !element belonging to master side
      !master sides(ElemID,locSide and flip =-1 if not existing)
      IF(ElemID.NE.-1)THEN ! element belonging to master side is on this processor
        locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
        U_ext(:,:,:,locSide,ElemID)=U_slave(:,:,:,SideID) ! Slave side is force-aligned with master
      end if
      
!     Slave side
!     ----------
      nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID) !element belonging to slave side
      IF(nbElemID.NE.-1)THEN! element belonging to slave side is on this processor
        nblocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
        nbFlip    = SideToElem(S2E_FLIP,SideID)
        
        select case (nblocSide)
          case(1,6) ! Side is normal to zeta
            DO q=0,PP_N; DO p=0,PP_N
              ijk(:)=S2V(:,0,p,q,nbflip,nblocSide)
              U_ext(:,ijk(1),ijk(2),nblocSide,nbElemID) = U_master(:,p,q,SideID)
            END DO; END DO !p,q=0,PP_N
          case(2,4) ! Side is normal to eta
            DO q=0,PP_N; DO p=0,PP_N
              ijk(:)=S2V(:,0,p,q,nbflip,nblocSide)
              U_ext(:,ijk(3),ijk(1),nblocSide,nbElemID) = U_master(:,p,q,SideID)
            END DO; END DO !p,q=0,PP_N
          case(3,5) ! Side is normal to xi
            DO q=0,PP_N; DO p=0,PP_N
              ijk(:)=S2V(:,0,p,q,nbflip,nblocSide)
              U_ext(:,ijk(2),ijk(3),nblocSide,nbElemID) = U_master(:,p,q,SideID)
            END DO; END DO !p,q=0,PP_N
        end select
      end if
    end do
    
    ! TODO: Do mortar MPI sides
    
  end subroutine Get_externalU
!===================================================================================================================================
!> Corrects U and Ut after the Runge-Kutta stage
!===================================================================================================================================
#if NFVSE_CORR
  subroutine Apply_NFVSE_Correction(U,Ut,t,dt)
    use MOD_ShockCapturing_Vars, only: alpha, alpha_max, PositCorrFactor, alpha_old, PositMaxIter
    use MOD_NFVSE_Vars         , only: Fsafe, Fblen
    use MOD_Mesh_Vars          , only: nElems, offsetElem
    use MOD_Basis              , only: ALMOSTEQUAL
    use MOD_Equation_Vars      , only: sKappaM1
    use MOD_Mesh_Vars          , only: sJ
    use MOD_NFVSE_MPI
    USE MOD_Globals
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: t                                         !< Current time (in time step!)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real, parameter :: eps = 1.e-8
    real    :: Usafe        (PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    real    :: Fsafe_m_Fblen(PP_nVar,0:PP_N,0:PP_N,0:PP_N)  ! Fsafe - Fblen
    real    :: FFV_m_FDG    (PP_nVar,0:PP_N,0:PP_N,0:PP_N)  ! Finite Volume Ut minus DG Ut
    real    :: a   ! a  = PositCorrFactor * rho_safe - rho
    real    :: ap  ! ap = (PositCorrFactor * p_safe   - p) / (kappa-1)
    real    :: pres
    real    :: corr, corr1
    real    :: p_safe(0:PP_N,0:PP_N,0:PP_N) ! pressure obtained with alpha_max for an element
    real    :: alphadiff
    real    :: alphacont  !container for alpha
    real    :: sdt
    integer :: eID
    integer :: i,j,k
    integer :: iter
    !--------------------------------------------------------
    
!   Some definitions
!   ****************
    alpha_old = alpha 
    sdt = 1./dt
    
!   Iterate over the elements
!   *************************
    do eID=1, nElems
      
!     ----------------------------------
!     Check if it makes sense correcting
!     ----------------------------------
      alphadiff = alpha(eID) - alpha_max
      if ( abs(alphadiff) < eps ) cycle ! Not much to do for this element...
      
!     ----------------------------------------------------------
!     Get Usafe for each point of the element and check validity
!     ----------------------------------------------------------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Compute Fsafe-Fblend
        Fsafe_m_Fblen(:,i,j,k) = ( Fsafe(:,i,j,k,eID) - Fblen(:,i,j,k,eID) ) * (-sJ(i,j,k,eID)) ! Account for the sign change and the Jacobian division
        
        ! Compute Usafe
        Usafe(:,i,j,k) = U(:,i,j,k,eID) + dt * Fsafe_m_Fblen(:,i,j,k)
        
        ! Check if this is a valid state
        call GetPressure(Usafe(:,i,j,k),p_safe(i,j,k))
        if (p_safe(i,j,k) < 0.) then
          print*, 'ERROR: safe pressure not safe el=', eID+offsetElem, p_safe(i,j,k)
          stop
        end if
        if (Usafe(1,i,j,k) < 0.) then
          print*, 'ERROR: safe dens not safe el=', eID+offsetElem, Usafe(1,i,j,k)
          stop
        end if
      end do       ; end do       ; end do ! i,j,k
      
      ! Compute F_FV-F_DG
      FFV_m_FDG = -Fsafe_m_Fblen/alphadiff
      
!     ----------------------------
!     Iterate to get a valid state
!     ----------------------------
      do iter=1, PositMaxIter
        corr = -eps ! Safe initialization
        
!       Compute correction factors
!       --------------------------
        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          
          ! Density correction
          a = (PositCorrFactor * Usafe(1,i,j,k) - U(1,i,j,k,eID))
          if (a > 0.) then ! This DOF needs a correction
            corr1 = a / FFV_m_FDG(1,i,j,k)
            corr = max(corr,corr1)
          end if
          
          ! Pressure correction
          call GetPressure(U(:,i,j,k,eID),pres)
          
          ap = (PositCorrFactor * p_safe(i,j,k) - pres) * sKappaM1
          if (ap > 0.) then
            corr1 = ap / FFV_m_FDG(5,i,j,k)
            corr = max(corr,corr1)
          end if
          
        end do       ; end do       ; end do ! i,j,k
        
      
!       Do the correction if needed
!       ---------------------------
        if ( corr > 0. ) then
          
          ! Change the alpha for output
          alphacont  = alpha(eID)
          alpha(eID) = alpha(eID) + corr * sdt
          
          ! Change inconsistent alphas
          if (alpha(eID) > alpha_max) then
            alpha(eID) = alpha_max
            corr = (alpha_max - alphacont ) * dt
          end if
          
          ! Correct!
          do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
            ! Correct U
            U (:,i,j,k,eID) = U (:,i,j,k,eID) + corr * FFV_m_FDG(:,i,j,k)
            ! Correct Ut
            Ut(:,i,j,k,eID) = Ut(:,i,j,k,eID) + (alpha(eID)-alphacont) * FFV_m_FDG(:,i,j,k)
          end do       ; end do       ; enddo
          
        else
          exit
        end if
      
      end do !iter
      
      if (iter > PositMaxIter) then
        write(*,'(A,I0,A,I0)') 'WARNING: Not able to perform NFVSE correction within ', PositMaxIter, ' iterations. Elem: ', eID + offsetElem
      end if
      
    end do !eID
    
  end subroutine Apply_NFVSE_Correction
#endif /*NFVSE_CORR*/
!===================================================================================================================================
!> Finalizes the NFVSE module
!===================================================================================================================================
  subroutine FinalizeNFVSE()
    use MOD_NFVSE_Vars, only: SubCellMetrics, sWGP, MPIRequest_alpha, Fsafe, Fblen, Compute_FVFluxes
    use MOD_NFVSE_Vars, only: U_ext, sdxR,sdxL,rR,rL
    implicit none
    
    SDEALLOCATE (SubCellMetrics)
    SDEALLOCATE (sWGP)
    SDEALLOCATE (MPIRequest_alpha)
    SDEALLOCATE (Fsafe)
    SDEALLOCATE (Fblen)
    
    Compute_FVFluxes => null()
    
    ! Reconstruction
    SDEALLOCATE ( U_ext )
    SDEALLOCATE ( sdxR  )
    SDEALLOCATE ( sdxL  )
    SDEALLOCATE ( rR    )
    SDEALLOCATE ( rL    )
    
  end subroutine FinalizeNFVSE
#endif /*SHOCK_NFVSE*/
end module MOD_NFVSE

