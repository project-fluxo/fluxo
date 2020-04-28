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
#if PARABOLIC
!===================================================================================================================================
!>Contains global variables used by the DG modules.
!===================================================================================================================================
MODULE MOD_Lifting_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! interior face gradients for all elements
REAL,ALLOCATABLE                      :: gradPx_Slave(:,:,:,:),gradPx_Master(:,:,:,:)
REAL,ALLOCATABLE                      :: gradPy_Slave(:,:,:,:),gradPy_Master(:,:,:,:)
REAL,ALLOCATABLE                      :: gradPz_Slave(:,:,:,:),gradPz_Master(:,:,:,:)
REAL,ALLOCATABLE                      :: FluxX(:,:,:,:)
REAL,ALLOCATABLE                      :: FluxY(:,:,:,:)
REAL,ALLOCATABLE                      :: FluxZ(:,:,:,:)
! the lifted gradients needed for NSE
REAL,ALLOCATABLE                      :: gradPx(:,:,:,:,:)
REAL,ALLOCATABLE                      :: gradPy(:,:,:,:,:)
REAL,ALLOCATABLE                      :: gradPz(:,:,:,:,:)
REAL                                  :: etaBR2              !<stabilization constant, the higher the stiffer (timestep!)
LOGICAL                               :: LiftingInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Lifting_Vars
#endif /*PARABOLIC*/
