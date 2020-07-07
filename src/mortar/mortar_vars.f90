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
!==================================================================================================================================
!> Variables used for mortars: mortar interpolation and projection matrices
!==================================================================================================================================
MODULE MOD_Mortar_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

REAL,ALLOCATABLE :: Mint(:,:,:)          !< 1D-Mortar Operator: interpolation full interval 0: [-1,1] to left interval 1: [-1,0] 
                                         !< and right intervall 2: [0,1], size is (0:N,0:N,1:2), 
                                        
REAL,ALLOCATABLE :: Mint_h(:,:,:)        !< Mint*0.5
REAL,ALLOCATABLE :: Mproj(:,:,:)         !< 1D-Mortar Operator: projection left interval 1: [-1,0] 
                                         !< and right intervall 2: [0,1] to full intervall 0: [-1,1],size is (0:N,0:N,1:2)
REAL,ALLOCATABLE :: Mproj_h(:,:,:)       !< Mproj*0.5
REAL,ALLOCATABLE :: U_small(:,:,:,:,:)   !< interpolated solution for the big mortar side to small.
                                         !< two sides for 2-1 mortars, and 4+2 sides for  4-1 mortars
#ifdef JESSE_MORTAR
REAL,ALLOCATABLE :: delta_flux_jesse(:,:,:,:)  !< contribution of the two-point flux correction to the big mortar side 
REAL,ALLOCATABLE :: Ns_small(:,:,:,:,:)  !< Ja normal of big face interpolated to small faces (scaled normal vector nvec*surfelem
                                         !< with factor 4 (4-1) or factor 2 (2-1).), also with intermediate solution for 4-1 mortar
#endif /*JESSE_MORTAR*/
LOGICAL          :: MortarInitIsDone=.FALSE. !< marks whether mortar init routines are complete
!==================================================================================================================================
END MODULE MOD_Mortar_Vars
