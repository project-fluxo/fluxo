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

!==================================================================================================================================
!> Contains global variables used by the shock capturing modules.
!==================================================================================================================================
MODULE MOD_ShockCapturing_Vars
!===================================================================================================================================
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: ShockCapturingInitIsDone=.FALSE.     !< marks whether the shock capturing init routine is complete
REAL,ALLOCATABLE :: sVdm_Leg(:,:)                        !< 1D inverse Vandermondematrix to Legendre polynomials

REAL             :: nu_max                               !< Globally maximum artificial viscosity
#if SHOCK_ARTVISC
REAL,ALLOCATABLE :: nu(:),nu_Master(:),nu_Slave(:)       !< Cellwise (and Facewise) artificial viscosity
#endif /*SHOCK_ARTVISC*/

#if SHOCK_LOC_ARTVISC
real, allocatable :: FilterMat(:,:)
real, allocatable :: artVisc     (:,:,:,:,:)   ! Artificial viscosities in each Gauss point
real, allocatable :: physSens    (:,:,:,:,:)   ! Physics based sensors
real, allocatable :: artVisc_master(:,:,:,:)   ! Artificial viscosities on the master side of a face
real, allocatable :: artVisc_slave (:,:,:,:)   ! Artificial viscosities on the slave side of a face
real, allocatable :: Mh_inv    (:,:,:,:,:,:)   ! Inverse of metric tensor in each DOF. Accordung to Möller and Hansbo: "On advancing front mesh generation in three dimensions" --> ( (d\vec{x}/d\vec{xi})^T (d\vec{x}/d\vec{xi}) )^-1
real, allocatable :: covMetrics(:,:,:,:,:,:)   ! Covariant metric tensor: (d\vec{x}/d\vec{xi})^T = ( (d\vec{xi}/d\vec{x})^T )^{-1}
real, allocatable :: sigmamin_Mh   (:,:,:,:)   ! Minimum eigenvalue of Möller's metric tensor
#endif /*SHOCK_LOC_ARTVISC*/

! For NFVSE
real, allocatable :: alpha(:)                             !< Element-wise blending function (modified every time Ut is computed)
real, allocatable :: alpha_old(:)                         !< Element-wise blending function (before correction)
real, allocatable :: alpha_Master(:)                      !< Blending function on master sides
real, allocatable :: alpha_Slave(:)                       !< Blending function on slave sides
real              :: threshold
real              :: alpha_max        ! Maximum blending factor
real              :: alpha_min        ! Minimum blending factor
integer           :: ModalThreshold
#if NFVSE_CORR
real              :: PositCorrFactor  ! Limiting factor for NFVSE correction
integer           :: PositMaxIter
#endif /*NFVSE_CORR*/
real, parameter   :: sharpness = log((1.d0-1.d-4)/1.d-4)

!===================================================================================================================================

END MODULE MOD_ShockCapturing_Vars
