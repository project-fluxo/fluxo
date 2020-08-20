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
!> Based on: Hennemann et al. (2020). "Entropy stable shock capturing for the discontinuous galerkin spectral element
!>                                     method with native finite volume sub elements"
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

contains

!===================================================================================================================================
!> Initializes the NFVSE module
!===================================================================================================================================
  subroutine InitNFVSE()
    use MOD_NFVSE_Vars         , only: SubCellMetrics, sWGP
    use MOD_ReadInTools        , only: GETINT
    use MOD_Mesh_Vars          , only: nElems, Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
    use MOD_Interpolation_Vars , only: wGP
    use MOD_DG_Vars            , only: D
    use MOD_Globals            , only: CROSS
#if MPI
    use MOD_NFVSE_Vars         , only: MPIRequest_alpha
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
    !--------------------------------------------------------------------------------------------------------------------------------
    
    ! Allocate storage
    allocate ( SubCellMetrics(nElems) )
    call SubCellMetrics % construct(PP_N)
    
    allocate ( sWGP(0:PP_N) )
    
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
    
    ! Compute the inverse of the quadrature weights (sub-cell dimensions)
    sWGP = 1.d0 / wGP
    
#if MPI
    allocate(MPIRequest_alpha(nNbProcs,4)    ) ! 1: send slave, 2: send master, 3: receive slave, 4, receive master
#endif
  end subroutine InitNFVSE
  
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
    use MOD_NFVSE_Vars         , only: SubCellMetrics, sWGP
    use MOD_ShockCapturing_Vars, only: alpha
    use MOD_Basis              , only: ALMOSTEQUAL
    use MOD_NFVSE_MPI          , only: PropagateBlendingCoeff
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
    real,dimension(PP_nVar)                         :: F_FV
    integer                                         :: i,j,k,iElem
    !===============================================================================================================================
    
    call PropagateBlendingCoeff()
    
    do iElem=1,nElems
      
      if ( ALMOSTEQUAL(alpha(iElem),0.d0) ) cycle
      
!     Compute the finite volume fluxes
!     --------------------------------
      call Compute_FVFluxes( U(:,:,:,:,iElem), ftilde , gtilde , htilde , SubCellMetrics(iElem) )
      
!     Update Ut
!     ---------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Get Finite Volume Ut
        F_FV =   sWGP(i) * ( ftilde(:,i,j,k) - ftilde(:,i-1,j  ,k  ) ) &
               + sWGP(j) * ( gtilde(:,i,j,k) - gtilde(:,i  ,j-1,k  ) ) &
               + sWGP(k) * ( htilde(:,i,j,k) - htilde(:,i  ,j  ,k-1) )

        ! Update Ut
        Ut   (:,i,j,k,iElem) = (1. - alpha(iElem)) * Ut(:,i,j,k,iElem) +  alpha(iElem) * F_FV
        
      end do       ; end do       ; end do ! i,j,k
      
    end do ! iElem
  end subroutine VolInt_NFVSE

!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a FV consistent flux
!===================================================================================================================================
  subroutine Compute_FVFluxes (U, F , G , H , sCM )
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: SubCellMetrics_t
    use MOD_Riemann   , only: AdvRiemann
    implicit none
    !-arguments---------------------------------------------------------------
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), intent(in)    :: U   !< The element solution
    real,dimension(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N), intent(inout) :: F   !< Left flux in xi
    real,dimension(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N), intent(inout) :: G   !< Left flux in eta
    real,dimension(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N), intent(inout) :: H   !< Left flux in zeta
    type(SubCellMetrics_t)                         , intent(in)    :: sCM !< Sub-cell metric terms
    !-local-variables---------------------------------------------------------
    integer  :: i,j,k
    real :: U_ (PP_nVar,0:PP_N,0:PP_N, 0:PP_N)
    real :: F_ (PP_nVar,0:PP_N,0:PP_N,-1:PP_N)
    !-------------------------------------------------------------------------
    
!   *****************
!   Initialize values
!   *****************
    
    F   = 0.0
    G   = 0.0
    H   = 0.0
    
!   *********
!   Xi-planes
!   *********
    F_  = 0.0
    U_ = reshape(U , shape(U_), order = [1,4,2,3])
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
                   sCM % xi   % nv(:,:,:,i),sCM % xi   % t1(:,:,:,i), sCM % xi   % t2(:,:,:,i))
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
    end do ! i (xi planes)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    F  = reshape(F_ , shape(F ), order = [1,3,4,2])

!   **********    
!   Eta-planes
!   **********
    F_ = 0.0
    U_ = reshape(U , shape(U_), order = [1,2,4,3])
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
                   sCM % eta  % nv(:,:,:,i),sCM % eta  % t1(:,:,:,i), sCM % eta  % t2(:,:,:,i))
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_ (:,j,k,i) = F_ (:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
    end do ! i (eta planes)
    
!   Reshape arrays back to original storage structure
!   -------------------------------------------------
    G  = reshape(F_ , shape(G ), order = [1,2,4,3])

!   ***********    
!   Zeta-planes
!   ***********
    
!   Fill inner interfaces
!   ---------------------
    do i=0, PP_N-1
      
      call AdvRiemann(H(:,:,:,i),U(:,:,:,i),U(:,:,:,i+1), &
                   sCM % zeta % nv(:,:,:,i),sCM % zeta % t1(:,:,:,i), sCM % zeta % t2(:,:,:,i))
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        H (:,j,k,i) = H (:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
    end do ! i (zeta planes)
    
  end subroutine Compute_FVFluxes
  

!===================================================================================================================================
!> Finalizes the NFVSE module
!===================================================================================================================================
  subroutine FinalizeNFVSE()
    use MOD_NFVSE_Vars, only: SubCellMetrics, sWGP, MPIRequest_alpha
    implicit none
    
    SDEALLOCATE (SubCellMetrics)
    SDEALLOCATE (sWGP)
    SDEALLOCATE (MPIRequest_alpha)
    
    
  end subroutine FinalizeNFVSE
#endif /*SHOCK_NFVSE*/
end module MOD_NFVSE

