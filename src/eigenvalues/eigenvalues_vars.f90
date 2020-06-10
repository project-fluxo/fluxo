!==================================================================================================================================
! Copyright (c) 2012 - 2019 Florian Hindenlang
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

!===================================================================================================================================
!> Contains global variables used by the eigenvalues module.
!>
!===================================================================================================================================
MODULE MOD_EigenValues_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

LOGICAL           :: calcEigVals    !< switch on/off calculation of eigenvalues 
INTEGER           :: stableDtfromEV !< =-1: off,
                                    !< =0: find stable timestep from Eigenvalues & charact. polynomial of time integator, 
                                    !< =1: apply new timestep
INTEGER           :: FD_degree      !< polynomial degree for approximation of Finite Diff. of Matvec (d/dU R(U))*v,
                                    !< exact for R(U) ~ U^FD_degree,  can be 1,2,4
INTEGER           :: Nx,Ny,Nz       !< max. index (either 0/PP_N) in each direction 
INTEGER           :: localDOF       !< total number of DOF on proc
INTEGER           :: globalDOF      !< total number of DOF over all procs
INTEGER           :: EV_maxIter     !< Numer of iterations of Arnoldi, can be smaller than globalDOF
REAL,ALLOCATABLE  :: U0(:,:,:,:,:)  !< Linearization state for matrix-Vector computation
REAL,ALLOCATABLE  :: Hmat(:,:)      !< upper Hesseberg matrix,size totalDOF^2, only saved on MPIroot
LOGICAL           :: InitEigenValuesDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_EigenValues_Vars
