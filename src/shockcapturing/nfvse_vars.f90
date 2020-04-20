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

!===================================================================================================================================
!> summary: Contains global variables for the NFVSE routines
!===================================================================================================================================
module MOD_NFVSE_Vars
  implicit none
  
  private
  public :: SubCellMetrics, SubCellMetrics_t, sWGP, MPIRequest_alpha, Fsafe, Fblen
  
!-----------------------------------------------------------------------------------------------------------------------------------
! New types
!-----------------------------------------------------------------------------------------------------------------------------------
  type :: InnerFaceMetrics_t
    real, allocatable :: nv(:,:,:,:) !< Sub-cell normal vectors            : (1:3,0:N,0:N,0:N-1)
    real, allocatable :: t1(:,:,:,:) !< Sub-cell tangent vectors (1)       : (1:3,0:N,0:N,0:N-1)
    real, allocatable :: t2(:,:,:,:) !< Sub-cell tangent vectors (2)       : (1:3,0:N,0:N,0:N-1)
    real, allocatable :: norm(:,:,:) !< Norm of the sub-cell normal vectors:     (0:N,0:N,0:N-1)
!                             | | |_ Inner face index (FV BC idx)
!                             | |___ (FV subcell idx)
!                             |_____ (FV subcell idx)
  end type InnerFaceMetrics_t
  
  type :: SubCellMetrics_t
    type(InnerFaceMetrics_t) :: xi
    type(InnerFaceMetrics_t) :: eta
    type(InnerFaceMetrics_t) :: zeta
    contains
      procedure :: construct => SubCellMetrics_construct
      procedure :: TestMetricIdentities => SubCellMetrics_TestMetricIdentities
  end type SubCellMetrics_t
  
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables
!-----------------------------------------------------------------------------------------------------------------------------------
  type(SubCellMetrics_t)    , allocatable :: SubCellMetrics(:)      !< Metric terms for the native sub-cell finite volumes
  real                      , allocatable :: sWGP(:)                !< Inverse of the Gauss quadrature weights
  integer                   , allocatable :: MPIRequest_alpha(:,:)  !< MPI request for the transfer of the blending coefficient
                                                                    !< (nNbProcs,4)... 1: send slave, 2: send master, 3: receive slave, 4, receive master
  real                      , allocatable :: Fsafe(:,:,:,:,:)
  real                      , allocatable :: Fblen(:,:,:,:,:)
!===================================================================================================================================
  contains
    elemental subroutine SubCellMetrics_construct(this,N)
      implicit none
      class(SubCellMetrics_t), intent(inout) :: this
      integer                , intent(in)    :: N
      
      allocate ( this % xi   % nv (1:3,0:N,0:N,0:N-1) )
      allocate ( this % xi   % t1 (1:3,0:N,0:N,0:N-1) )
      allocate ( this % xi   % t2 (1:3,0:N,0:N,0:N-1) )
      allocate ( this % xi   % norm   (0:N,0:N,0:N-1) )
      
      allocate ( this % eta  % nv (1:3,0:N,0:N,0:N-1) )
      allocate ( this % eta  % t1 (1:3,0:N,0:N,0:N-1) )
      allocate ( this % eta  % t2 (1:3,0:N,0:N,0:N-1) )
      allocate ( this % eta   % norm   (0:N,0:N,0:N-1) )
      
      allocate ( this % zeta % nv (1:3,0:N,0:N,0:N-1) )
      allocate ( this % zeta % t1 (1:3,0:N,0:N,0:N-1) )
      allocate ( this % zeta % t2 (1:3,0:N,0:N,0:N-1) )
      allocate ( this % zeta % norm   (0:N,0:N,0:N-1) )
      
    end subroutine SubCellMetrics_construct
    
    subroutine SubCellMetrics_TestMetricIdentities(this,iSC,jSC,kSC,div)
      use MOD_Interpolation_Vars , only: wGP
      implicit none
      class(SubCellMetrics_t), intent(inout) :: this
      integer                , intent(in)    :: iSC, jSC, kSC !< Subcell index: 1...N-1
      real, dimension(3)     , intent(out)   :: div
      !----------------------------------------------
      real, dimension(3) :: a, b, c, d, e, f
      real               :: anorm, bnorm, cnorm, dnorm, enorm, fnorm
      !----------------------------------------------
        
      a     = this % xi   % nv(:,jSC,kSC,iSC)
      anorm = this % xi   % norm(jSC,kSC,iSC)
      b     = this % xi   % nv(:,jSC,kSC,iSC-1)
      bnorm = this % xi   % norm(jSC,kSC,iSC-1)
      c     = this % eta  % nv(:,iSC,kSC,jSC)
      cnorm = this % eta  % norm(iSC,kSC,jSC)
      d     = this % eta  % nv(:,iSC,kSC,jSC-1)
      dnorm = this % eta  % norm(iSC,kSC,jSC-1)
      e     = this % zeta % nv(:,iSC,jSC,kSC)
      enorm = this % zeta % norm(iSC,jSC,kSC)
      f     = this % zeta % nv(:,iSC,jSC,kSC-1)
      fnorm = this % zeta % norm(iSC,jSC,kSC-1)
      
      div = (a * anorm - b * bnorm) / wGP(iSC) + (c * cnorm - d * dnorm) / wGP(jSC) + (e * enorm - f * fnorm) / wGP(kSC)
    end subroutine SubCellMetrics_TestMetricIdentities
end module MOD_NFVSE_Vars
