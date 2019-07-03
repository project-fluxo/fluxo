!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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
!> Computes two-point average fluxes for the volint when using the split-form (DiscType=2)
!==================================================================================================================================
MODULE MOD_Flux_Average
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE EvalEulerFluxTilde3D
  MODULE PROCEDURE EvalEulerFluxTilde3D
END INTERFACE

INTERFACE EvalUaux
  MODULE PROCEDURE EvalUaux
END INTERFACE

INTERFACE StandardDGFlux
  MODULE PROCEDURE StandardDGFlux
END INTERFACE

INTERFACE StandardDGFluxVec
  MODULE PROCEDURE StandardDGFluxVec
END INTERFACE

INTERFACE StandardDGFluxDealiasedMetricVec
  MODULE PROCEDURE StandardDGFluxDealiasedMetricVec
END INTERFACE


PUBLIC::EvalEulerFluxTilde3D
PUBLIC::EvalUaux
PUBLIC::StandardDGFlux
PUBLIC::StandardDGFluxVec
PUBLIC::StandardDGFluxDealiasedMetricVec
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute linadvdiff transformed fluxes using conservative variables and derivatives for every volume Gauss point.
!> directly apply metrics and output the tranformed flux 
!==================================================================================================================================
SUBROUTINE EvalEulerFluxTilde3D(iElem,ftilde,gtilde,htilde,Uaux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars ,ONLY:nAuxVar
USE MOD_Equation_Vars ,ONLY:AdvVel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem !< current element index in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed flux f(iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed flux g(iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed flux h(iVar,i,j,k)
REAL,DIMENSION(nAuxVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: Uaux                      !< auxiliary variables, not needed here
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Uaux=0.
ftilde(1,:,:,:) = ( AdvVel(1)*Metrics_fTilde(1,:,:,:,iElem) &
                   +AdvVel(2)*Metrics_fTilde(2,:,:,:,iElem) &
                   +AdvVel(3)*Metrics_fTilde(3,:,:,:,iElem)) *U(1,:,:,:,iElem)
gtilde(1,:,:,:) = ( AdvVel(1)*Metrics_gTilde(1,:,:,:,iElem) &
                   +AdvVel(2)*Metrics_gTilde(2,:,:,:,iElem) &
                   +AdvVel(3)*Metrics_gTilde(3,:,:,:,iElem)) *U(1,:,:,:,iElem)
htilde(1,:,:,:) = ( AdvVel(1)*Metrics_hTilde(1,:,:,:,iElem) &
                   +AdvVel(2)*Metrics_hTilde(2,:,:,:,iElem) &
                   +AdvVel(3)*Metrics_hTilde(3,:,:,:,iElem)) *U(1,:,:,:,iElem)
END SUBROUTINE EvalEulerFluxTilde3D


!==================================================================================================================================
!> computes auxiliary nodal variables (1/rho,v_1,v_2,v_3,p,|v|^2) 
!==================================================================================================================================
SUBROUTINE EvalUaux(iElem,Uaux)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:nAuxVar
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem !< current element index in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(nAuxVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: Uaux  !< auxiliary variables:(srho,v1,v2,v3,p,|v|^2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Uaux=0.
END SUBROUTINE EvalUaux


!==================================================================================================================================
!> Computes the standard flux in normal-direction
!==================================================================================================================================
SUBROUTINE StandardDGFlux(Fstar,UL,UR,normal)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL,UR   !< left and right state
REAL,DIMENSION(3),INTENT(IN)        :: normal  !< normal vector
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Fstar=SUM(AdvVel(1:3)*normal(1:3))*0.5*(UL+UR)
END SUBROUTINE StandardDGFlux


!==================================================================================================================================
!> Computes the standard DG flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the scalar advection eq.
!> for curved metrics, no dealiasing is done (exactly = standard DG )!
!==================================================================================================================================
SUBROUTINE StandardDGFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:AdvVel
USE MOD_Equation_Vars ,ONLY:nAuxVar
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
REAL,DIMENSION(nAuxVar),INTENT(IN)  :: UauxL          !< left auxiliary variables
REAL,DIMENSION(nAuxVar),INTENT(IN)  :: UauxR          !< right auxiliary variables
REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

!without metric dealising (standard DG flux on curved meshes)
Fstar=0.5*(  AdvVel(1)*(metric_L(1)*UL+metric_R(1)*UR) &
           + AdvVel(2)*(metric_L(2)*UL+metric_R(2)*UR) &
           + AdvVel(3)*(metric_L(3)*UL+metric_R(3)*UR) )

END SUBROUTINE StandardDGFluxVec

!==================================================================================================================================
!> Computes the standard DG flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the scalar adv. equation
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE StandardDGFluxDealiasedMetricVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
USE MOD_Equation_Vars ,ONLY:nAuxVar
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
REAL,DIMENSION(nAuxVar),INTENT(IN)  :: UauxL          !< left auxiliary variables
REAL,DIMENSION(nAuxVar),INTENT(IN)  :: UauxR          !< right auxiliary variables
REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
REAL,INTENT(IN)                     :: metric_R(3)    !< right mertric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
!Fstar=SUM(AdvVel(1:3)*metric(1:3))*0.5*(UL+UR)

Fstar=0.5*(  AdvVel(1)*metric(1)*(UL+UR) &
           + AdvVel(2)*metric(2)*(UL+UR) &
           + AdvVel(3)*metric(3)*(UL+UR) )

END SUBROUTINE StandardDGFluxDealiasedMetricVec

END MODULE MOD_Flux_Average
