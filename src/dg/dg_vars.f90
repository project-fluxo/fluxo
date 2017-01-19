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

!==================================================================================================================================
!> Contains global variables provided by the DG modules.
!>
!> This module contains the building blocks for the DGSEM operator, i.e. the interpolation and differentiation 
!> operators and the variables for solution storage as well as auxilliary variables. 
!==================================================================================================================================
MODULE MOD_DG_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! DG basis, contains the differentiation and interpolation operators
REAL,ALLOCATABLE                      :: D(:,:)                 !< Differentiation matrix of size [0..N,0..N], contains the 
                                                                !< first derivative of each Lagrange polynomial at each node. 

REAL,ALLOCATABLE                      :: D_T(:,:)               !< Transpose of differentiation matrix, size [0..N,0..N]. 

REAL,ALLOCATABLE                      :: D_Hat(:,:)             !< Differentiation matrix premultiplied by
                                                                !< mass matrix, \f$ \hat{D} = M^{-1} D^T M \f$, size [0..N,0..N].  

REAL,ALLOCATABLE                      :: D_Hat_T(:,:)           !< Transpose of differentiation matrix premultiplied by 
                                                                !< mass matrix, size [0..N,0..N].
 
REAL,ALLOCATABLE                      :: L_HatMinus(:)          !< Lagrange polynomials evaluated at \f$\xi=-1\f$
                                                                !< premultiplied by mass matrix

REAL,ALLOCATABLE                      :: L_HatPlus(:)           !< Lagrange polynomials evaluated at \f$\xi=+1\f$
                                                                !< premultiplied by mass matrix
!----------------------------------------------------------------------------------------------------------------------------------
! DG solution (JU or U) vectors)
REAL,ALLOCATABLE,TARGET               :: U(:,:,:,:,:)           !< Solution variable for each equation, node and element, 
                                                                !< size [1..NVar,0..N,0..N,0..N,nElems]. 

!----------------------------------------------------------------------------------------------------------------------------------
! DG time derivative or Residual U_t
REAL,ALLOCATABLE                      :: Ut(:,:,:,:,:)          !< Residual/time derivative, size [1..NVar,0..N,0..N,0..N,nElems]. 

!----------------------------------------------------------------------------------------------------------------------------------
! auxilliary counters: number of entries in U, Ut, gradUx, gradUy, gradUz, used of optimization 
INTEGER                               :: nTotalU                !< Total number of entries in U / size of U. 

INTEGER                               :: nDOFElem               !< Degrees of freedom in single element(per equation) 
                                                                !< $ nDOFElem=(N+1)^3 $.  

!----------------------------------------------------------------------------------------------------------------------------------
! interior face values for all elements
REAL,ALLOCATABLE                      :: U_master(:,:,:,:)      !< 2D Solution on face nodes for the master sides, 
                                                                !< size [1..NVar,0..N,0..N,all_master_sides] 

REAL,ALLOCATABLE                      :: U_slave(:,:,:,:)       !< 2D Solution on face nodes for the slave sides, 
REAL,ALLOCATABLE                      :: Flux_master(:,:,:,:)    !< Fluxes on face, size [1..NVar,0..N,0..N,allsides]. 
REAL,ALLOCATABLE                      :: Flux_slave (:,:,:,:)    !< Fluxes on face, size [1..NVar,0..N,0..N,allsides]. 
!----------------------------------------------------------------------------------------------------------------------------------
! Auxilliary variables
LOGICAL                               :: DGInitIsDone=.FALSE.   !< Switch to check DGInit status
!==================================================================================================================================
END MODULE MOD_DG_Vars
