!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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

!==================================================================================================================================
!> Computes the spatial operator using Native Finite Volume Sub-Elements (NFVSE)
!> Based on: Hennemann and Gassner (2020). "Entropy stable shock capturing for the discontinuous galerkin spectral element
!>                                          method with native finite volume sub elements"
!> Attention 1: The current strategy is to perform the operations on "inner faces". This can be improved, but the Riemann solver 
!>              routines have to be adapted
!==================================================================================================================================
module MOD_NFVSE
#if SHOCK_NFVSE
  use MOD_PreProc
  implicit none
  
  private
  public :: VolInt_NFVSE, InitNFVSE, FinalizeNFVSE
  
contains
!===================================================================================================================================
!> Initializes the NFVSE module
!===================================================================================================================================
  subroutine InitNFVSE()
    use MOD_NFVSE_Vars        , only: SubCellMetrics, sWGP, MPIRequest_alpha
    use MOD_MPI_Vars          , only: nNbProcs
    use MOD_Mesh_Vars         , only: nElems, Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
    use MOD_Interpolation_Vars, only: wGP
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    integer :: iElem
    integer :: i,j,k      !DOF counters
    real, parameter :: half = 0.5d0
    !--------------------------------------------------------------------------------------------------------------------------------
    
    ! Allocate storage
    allocate ( SubCellMetrics(nElems) )
    call SubCellMetrics % construct(PP_N)
    
    allocate ( sWGP(0:PP_N) )
    
    ! Compute inner normal and tangent vectors (the storage access to Metrics_*Tilde is not optimized... but this is done only once!)
    do iElem=1, nElems
      
      do i=0, PP_N-1
        ! Compute vectors
        SubCellMetrics(iElem) % xi   % nv(:,:,:,i) = half * ( Metrics_fTilde(:,i,:,:,iElem) + Metrics_fTilde(:,i+1,:,:,iElem) )
        SubCellMetrics(iElem) % xi   % t1(:,:,:,i) = half * ( Metrics_gTilde(:,i,:,:,iElem) + Metrics_gTilde(:,i+1,:,:,iElem) )
        SubCellMetrics(iElem) % xi   % t2(:,:,:,i) = half * ( Metrics_hTilde(:,i,:,:,iElem) + Metrics_hTilde(:,i+1,:,:,iElem) )
        
        SubCellMetrics(iElem) % eta  % nv(:,:,:,i) = half * ( Metrics_gTilde(:,:,i,:,iElem) + Metrics_gTilde(:,:,i+1,:,iElem) )
        SubCellMetrics(iElem) % eta  % t1(:,:,:,i) = half * ( Metrics_hTilde(:,:,i,:,iElem) + Metrics_hTilde(:,:,i+1,:,iElem) )
        SubCellMetrics(iElem) % eta  % t2(:,:,:,i) = half * ( Metrics_fTilde(:,:,i,:,iElem) + Metrics_fTilde(:,:,i+1,:,iElem) )
        
        SubCellMetrics(iElem) % zeta % nv(:,:,:,i) = half * ( Metrics_hTilde(:,:,:,i,iElem) + Metrics_hTilde(:,:,:,i+1,iElem) )
        SubCellMetrics(iElem) % zeta % t1(:,:,:,i) = half * ( Metrics_fTilde(:,:,:,i,iElem) + Metrics_fTilde(:,:,:,i+1,iElem) )
        SubCellMetrics(iElem) % zeta % t2(:,:,:,i) = half * ( Metrics_gTilde(:,:,:,i,iElem) + Metrics_gTilde(:,:,:,i+1,iElem) )
        
        ! Normalize each
        do k=0, PP_N ; do j=0, PP_N
          SubCellMetrics(iElem) % xi   % norm(j,k,i) = norm2 (SubCellMetrics(iElem) % xi   % nv(:,j,k,i))
          SubCellMetrics(iElem) % xi   % nv(:,j,k,i) = SubCellMetrics(iElem) % xi   % nv(:,j,k,i) / SubCellMetrics(iElem) % xi   % norm(j,k,i)
          SubCellMetrics(iElem) % xi   % t1(:,j,k,i) = SubCellMetrics(iElem) % xi   % t1(:,j,k,i) / norm2 (SubCellMetrics(iElem) % xi   % t1(:,j,k,i))
          SubCellMetrics(iElem) % xi   % t2(:,j,k,i) = SubCellMetrics(iElem) % xi   % t2(:,j,k,i) / norm2 (SubCellMetrics(iElem) % xi   % t2(:,j,k,i))
          
          SubCellMetrics(iElem) % eta  % norm(j,k,i) = norm2 (SubCellMetrics(iElem) % eta  % nv(:,j,k,i))
          SubCellMetrics(iElem) % eta  % nv(:,j,k,i) = SubCellMetrics(iElem) % eta  % nv(:,j,k,i) / SubCellMetrics(iElem) % eta  % norm(j,k,i)
          SubCellMetrics(iElem) % eta  % t1(:,j,k,i) = SubCellMetrics(iElem) % eta  % t1(:,j,k,i) / norm2 (SubCellMetrics(iElem) % eta  % t1(:,j,k,i))
          SubCellMetrics(iElem) % eta  % t2(:,j,k,i) = SubCellMetrics(iElem) % eta  % t2(:,j,k,i) / norm2 (SubCellMetrics(iElem) % eta  % t2(:,j,k,i))
          
          SubCellMetrics(iElem) % zeta % norm(j,k,i) = norm2 (SubCellMetrics(iElem) % zeta % nv(:,j,k,i))
          SubCellMetrics(iElem) % zeta % nv(:,j,k,i) = SubCellMetrics(iElem) % zeta % nv(:,j,k,i) / SubCellMetrics(iElem) % zeta % norm(j,k,i)
          SubCellMetrics(iElem) % zeta % t1(:,j,k,i) = SubCellMetrics(iElem) % zeta % t1(:,j,k,i) / norm2 (SubCellMetrics(iElem) % zeta % t1(:,j,k,i))
          SubCellMetrics(iElem) % zeta % t2(:,j,k,i) = SubCellMetrics(iElem) % zeta % t2(:,j,k,i) / norm2 (SubCellMetrics(iElem) % zeta % t2(:,j,k,i))
        end do       ; end do
      end do
    end do
    
    ! Compute the inverse of the quadrature weights (sub-cell dimensions)
    sWGP = 1.d0 / wGP
    
#if MPI
    allocate(MPIRequest_alpha(nNbProcs,4)    ) ! 1: send slave, 2: send master, 3: receive slave, 4, receive master
#endif
  end subroutine InitNFVSE
  
!===================================================================================================================================
!> Computes the "volume integral": Spatial contribution to Ut by the subcell finite volumes
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!> Attention 3: The element boundaries are neglected in this routine
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
    real,intent(inout)                                 :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
    !< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    real,dimension(PP_nVar,-1:PP_N,-1:PP_N,-1:PP_N) :: ftilde, gtilde, htilde ! transformed flux inter-subcell fluxes (with ghost cells)
    integer                                            :: i,j,k,iElem
    !===============================================================================================================================
    
    call PropagateBlendingCoeff()
    
    ftilde = 0.d0
    gtilde = 0.d0
    htilde = 0.d0
    
    do iElem=1,nElems
      if ( ALMOSTEQUAL(alpha(iElem),0.d0) ) cycle
      !compute inner Riemann solutions (TODO: Not optimal.. Improve)
      call Compute_VolFluxes( U(:,:,:,:,iElem), SubCellMetrics(iElem), ftilde(:,0:PP_N-1, 0:PP_N  , 0:PP_N  ), &
                                                                       gtilde(:,0:PP_N  , 0:PP_N-1, 0:PP_N  ), &
                                                                       htilde(:,0:PP_N  , 0:PP_N  , 0:PP_N-1)  )
      
      ! Update the inner nodes
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        Ut(:,i,j,k,iElem) = (1.d0 - alpha(iElem)) * Ut(:,i,j,k,iElem) &
                          +  alpha(iElem) * (  sWGP(i) * ( ftilde(:,i,j,k) - ftilde(:,i-1,j  ,k  ) ) &
                                             + sWGP(j) * ( gtilde(:,i,j,k) - gtilde(:,i  ,j-1,k  ) ) &
                                             + sWGP(k) * ( htilde(:,i,j,k) - htilde(:,i  ,j  ,k-1) ) )
      end do       ; end do       ; end do ! i,j,k
      
    end do ! iElem
  end subroutine VolInt_NFVSE
  
!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a FV consistent flux
!> Attention 1: It's not optimal to reshape the vectors (TODO: Improve)
! TODO: add parabolic part...
!===================================================================================================================================
  subroutine Compute_VolFluxes(U,sCM,F,G,H)
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: SubCellMetrics_t
    use MOD_Riemann   , only: Riemann
    implicit none
    !
    real,dimension(PP_nVar,0:PP_N  ,0:PP_N  ,0:PP_N  ), intent(in)    :: U         !< The element solution
    type(SubCellMetrics_t)                            , intent(in)    :: sCM       !< Sub-cell metric terms
    real,dimension(PP_nVar,0:PP_N-1,0:PP_N  ,0:PP_N  ), intent(inout) :: F
    real,dimension(PP_nVar,0:PP_N  ,0:PP_N-1,0:PP_N  ), intent(inout) :: G
    real,dimension(PP_nVar,0:PP_N  ,0:PP_N  ,0:PP_N-1), intent(inout) :: H
    !
    integer  :: i,j,k
    real :: U_(PP_nVar,0:PP_N,0:PP_N,0:PP_N  )
    real :: F_(PP_nVar,0:PP_N,0:PP_N,0:PP_N-1)
    
    !--
    
!   Xi-planes
!   ---------
    U_ = reshape(U , shape(U_), order = [1,4,2,3])
    do i=0, PP_N-1
      
      call Riemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif
                   sCM % xi   % nv(:,:,:,i),sCM % xi   % t1(:,:,:,i), sCM % xi   % t2(:,:,:,i))
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_(:,j,k,i) = F_(:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
    end do
    F  = reshape(F_, shape(F ), order = [1,3,4,2])
    
!   Eta-planes
!   ----------
    U_ = reshape(U , shape(U_), order = [1,2,4,3])
    do i=0, PP_N-1
      
      call Riemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif
                   sCM % eta  % nv(:,:,:,i),sCM % eta  % t1(:,:,:,i), sCM % eta  % t2(:,:,:,i))
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_(:,j,k,i) = F_(:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
    end do
    G  = reshape(F_, shape(G ), order = [1,2,4,3])
    
!   Zeta-planes
!   -----------
    do i=0, PP_N-1
      
      call Riemann(H(:,:,:,i),U(:,:,:,i),U(:,:,:,i+1), &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif
                   sCM % zeta % nv(:,:,:,i),sCM % zeta % t1(:,:,:,i), sCM % zeta % t2(:,:,:,i))
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        H (:,j,k,i) = H (:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
    end do
    
  end subroutine Compute_VolFluxes
  
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

