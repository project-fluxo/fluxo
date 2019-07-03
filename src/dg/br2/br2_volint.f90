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
! Containes the DG volume integral for BR2 lifting
!===================================================================================================================================
MODULE MOD_Lifting_VolInt
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE Lifting_VolInt
  MODULE PROCEDURE Lifting_VolInt
END INTERFACE

PUBLIC::Lifting_VolInt
!===================================================================================================================================
CONTAINS


!===================================================================================================================================
!> Computes the Volume integral of the BR2 scheme, in strong form so that first gradients in xi,eta,zeta are computed 
!> and then transformed to the x,y,z gradients all together, including the Jacobian!
!> the transpose of the D matrix  (D_T) is only used because of memory access
!===================================================================================================================================
SUBROUTINE Lifting_VolInt(gradPx,gradPy,gradPz)
!-----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_DG_Vars            ,ONLY:D_T
USE MOD_Mesh_Vars          ,ONLY:nElems,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,sJ   ! metrics
USE MOD_Flux               ,ONLY: EvalLiftingVolumeFlux
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                           :: gradPx(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradient of U in x
REAL,INTENT(INOUT)                           :: gradPy(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradient of U in y
REAL,INTENT(INOUT)                           :: gradPz(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< gradient of U in z
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)                      :: gradPxi,gradPeta,gradPzeta
REAL                                         :: Flux(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
INTEGER                                      :: iElem,i,j,k
INTEGER                                      :: l
!===================================================================================================================================
! volume integral, gradients should already be initialized =0
!gradPx=0.
!gradPy=0.
!gradPz=0.

DO iElem=1,nElems
  CALL EvalLiftingVolumeFlux(iElem,Flux)
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        gradPxi(:)   = D_T(0,i)*Flux(:,0,j,k)
        gradPeta(:)  = D_T(0,j)*Flux(:,i,0,k)
        gradPzeta(:) = D_T(0,k)*Flux(:,i,j,0)
        DO l=1,PP_N
          gradPxi(:)   = gradPxi(:)   +D_T(l,i)*Flux(:,l,j,k)
          gradPeta(:)  = gradPeta(:)  +D_T(l,j)*Flux(:,i,l,k)
          gradPzeta(:) = gradPzeta(:) +D_T(l,k)*Flux(:,i,j,l)
        END DO !i
        gradPx(:,i,j,k,iElem) = gradPx(:,i,j,k,iElem) + ( Metrics_fTilde(1,i,j,k,iElem)*gradPxi(:)    &
                                                         +Metrics_gTilde(1,i,j,k,iElem)*gradPeta(:)   &
                                                         +Metrics_hTilde(1,i,j,k,iElem)*gradPzeta(:) )*sJ(i,j,k,iElem)
        gradPy(:,i,j,k,iElem) = gradPy(:,i,j,k,iElem) + ( Metrics_fTilde(2,i,j,k,iElem)*gradPxi(:)    &
                                                         +Metrics_gTilde(2,i,j,k,iElem)*gradPeta(:)   &
                                                         +Metrics_hTilde(2,i,j,k,iElem)*gradPzeta(:) )*sJ(i,j,k,iElem)
        gradPz(:,i,j,k,iElem) = gradPz(:,i,j,k,iElem) + ( Metrics_fTilde(3,i,j,k,iElem)*gradPxi(:)    &
                                                         +Metrics_gTilde(3,i,j,k,iElem)*gradPeta(:)   &
                                                         +Metrics_hTilde(3,i,j,k,iElem)*gradPzeta(:) )*sJ(i,j,k,iElem)
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem=1,nElems
END SUBROUTINE Lifting_VolInt

END MODULE MOD_Lifting_VolInt
#endif /*PARABOLIC*/
