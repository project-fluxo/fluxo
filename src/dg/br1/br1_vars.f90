!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
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

#if PARABOLIC
!==================================================================================================================================
!> Contains global variables used by the BR1 lifting.
!==================================================================================================================================
MODULE MOD_Lifting_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! interior face gradients for all elements
REAL,ALLOCATABLE :: gradPx_slave(:,:,:,:)         !< slave side gradients in x-dir
REAL,ALLOCATABLE :: gradPy_slave(:,:,:,:)         !< slave side gradients in x-dir
REAL,ALLOCATABLE :: gradPz_slave(:,:,:,:)         !< slave side gradients in x-dir
REAL,ALLOCATABLE :: gradPx_master(:,:,:,:)        !< master side gradients in x-dir
REAL,ALLOCATABLE :: gradPy_master(:,:,:,:)        !< master side gradients in x-dir
REAL,ALLOCATABLE :: gradPz_master(:,:,:,:)        !< master side gradients in x-dir
REAL,ALLOCATABLE :: gradPx(:,:,:,:,:)             !< gradients in x-dir at degree N
REAL,ALLOCATABLE :: gradPy(:,:,:,:,:)             !< gradients in y-dir at degree N
REAL,ALLOCATABLE :: gradPz(:,:,:,:,:)             !< gradients in z-dir at degree N
REAL,ALLOCATABLE :: FluxX(:,:,:,:)                !< = 1/2(U_outer-U_inner)*shat*nvec_x
REAL,ALLOCATABLE :: FluxY(:,:,:,:)                !< = 1/2(U_outer-U_inner)*shat*nvec_y
REAL,ALLOCATABLE :: FluxZ(:,:,:,:)                !< = 1/2(U_outer-U_inner)*shat*nvec_z
LOGICAL          :: LiftingInitIsDone=.FALSE.     !< marks whether the lifting init routines are complete
!==================================================================================================================================
END MODULE MOD_Lifting_Vars
#endif /*PARABOLIC*/
