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
!> Contains global variables provided by the CT modules.
!>
!> This module contains the variables for the constraint transport implementation, where div-free magnetic fields are constructed 
!>
!==================================================================================================================================
MODULE MOD_CT_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! CT solution 
REAL,ALLOCATABLE                      :: curlA(:,:,:,:,:)       !< div-free B computeed from curl(A)
                                                                !< size \([1..3,0..N,0..N,0..N,nElems]\). 
!----------------------------------------------------------------------------------------------------------------------------------
! curl of time derivative of A
REAL,ALLOCATABLE                      :: curlAt(:,:,:,:,:)      !< Residual/time derivative, size \([1..3,0..N,0..N,0..N,nElems]\). 

!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                               :: nTotalcurlA            !< Total number of entries in curlA / size of curlA. 
                                                                !< ntotalcurlA\(=3(N+1)^3nElems\).  

REAL,ALLOCATABLE                      :: JacMat(:,:,:,:,:,:)    !< inverse of inverse jacobian matrix from (f,g,h)_metrics 
!----------------------------------------------------------------------------------------------------------------------------------
! for projection
INTEGER                               :: Nover
INTEGER                               :: Mint
REAL,ALLOCATABLE                      :: x_int(:),w_int(:)    !< integration points and weights
REAL,ALLOCATABLE                      :: Vdm_GLN_int(:,:)     !< vandermonde from Gauss-Lob. degree N to integration points
REAL,ALLOCATABLE                      :: Vdm_e_GLN(:,:)       !< vandermonde form edge basis functions (degree N-1) to Gauss-Lob. degree N
REAL,ALLOCATABLE                      :: covVec_1(:,:,:,:,:)  !<
REAL,ALLOCATABLE                      :: covVec_2(:,:,:,:,:)  !<
REAL,ALLOCATABLE                      :: covVec_3(:,:,:,:,:)  !<
REAL,ALLOCATABLE                      :: elem_xint_1(:,:,:,:,:)  !<
REAL,ALLOCATABLE                      :: elem_xint_2(:,:,:,:,:)  !<
REAL,ALLOCATABLE                      :: elem_xint_3(:,:,:,:,:)  !<
!----------------------------------------------------------------------------------------------------------------------------------
! Auxilliary variables
INTEGER                               :: use_CT
LOGICAL                               :: doEvalSource_A=.TRUE.  !< switch is set to false if no source is used 
LOGICAL                               :: CTInitIsDone=.FALSE.   !< Switch to check CTInit status
!==================================================================================================================================
END MODULE MOD_CT_Vars
