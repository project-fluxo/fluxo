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

#if (PP_DiscType==2)
INTERFACE EvalAdvFluxAverage3D
  MODULE PROCEDURE EvalAdvFluxAverage3D
END INTERFACE
#endif /*PP_DiscType==2*/

INTERFACE StandardDGFlux
  MODULE PROCEDURE StandardDGFlux
END INTERFACE

INTERFACE StandardDGFluxVec
  MODULE PROCEDURE StandardDGFluxVec
END INTERFACE

INTERFACE StandardDGFluxDealiasedMetricVec
  MODULE PROCEDURE StandardDGFluxDealiasedMetricVec
END INTERFACE

#if (PP_DiscType==2)
PUBLIC:: EvalAdvFluxAverage3D
PUBLIC:: EvalAdvFluxAverage
#endif /*PP_DiscType==2*/
PUBLIC::StandardDGFlux
PUBLIC::StandardDGFluxVec
PUBLIC::StandardDGFluxDealiasedMetricVec
!==================================================================================================================================
! local definitions for inlining / optimizing routines, DEFAULT=-1: USE POINTER defined at runtime!
#if PP_VolFlux==-1
#  define PP_VolumeFluxAverageVec VolumeFluxAverageVec
#elif PP_VolFlux==0
#  define PP_VolumeFluxAverageVec StandardDGFluxVec 
#elif PP_VolFlux==1
#  define PP_VolumeFluxAverageVec StandardDGFluxDealiasedMetricVec
#endif


CONTAINS

#if (PP_DiscType==2)
!==================================================================================================================================
!> Compute flux differences in 3D, making use of the symmetry and appling also directly the metrics  
!==================================================================================================================================
SUBROUTINE EvalAdvFluxAverage3D(U_in,M_f,M_g,M_h,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
#if PP_VolFlux==-1
USE MOD_Equation_Vars  ,ONLY:VolumeFluxAverageVec !pointer to flux averaging routine
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1   ,0:PP_N,0:PP_N,0:PP_N),INTENT(IN ) :: U_in        !< solution
REAL,DIMENSION(1:3 ,0:PP_N,0:PP_N,0:PP_N),INTENT(IN ) :: M_f,M_g,M_h !< metrics
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde,gtilde,htilde !< 4D transformed fluxes (iVar,i,,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER        :: i,j,k,l
!==================================================================================================================================


!opt_v1
!Uaux not needed for linadv 
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part not needed since diagonal of DvolSurfMat is zero!
  !ftilde(:,i,i,j,k)=ftilde_c(:,i,j,k) 
  ftilde(:,i,i,j,k)=0.
  DO l=i+1,PP_N
    CALL PP_VolumeFluxAverageVec(U_in(:,i,j,k),U_in(:,l,j,k), &
                                  M_f(:,i,j,k), M_f(:,l,j,k), &
                             ftilde(:,l,i,j,k)                )
    ftilde(:,i,l,j,k)=ftilde(:,l,i,j,k) !symmetric
  END DO!l=i+1,N
END DO; END DO; END DO ! i,j,k
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part not needed since diagonal of DvolSurfMat is zero!
  !gtilde(:,j,i,j,k)=gtilde_c(:,i,j,k) 
  gtilde(:,j,i,j,k)=0.
  DO l=j+1,PP_N
    CALL PP_VolumeFluxAverageVec(U_in(:,i,j,k),U_in(:,i,l,k), &
                                  M_g(:,i,j,k), M_g(:,i,l,k), &
                             gtilde(:,l,i,j,k)                )
    gtilde(:,j,i,l,k)=gtilde(:,l,i,j,k) !symmetric
  END DO!l=j+1,N
END DO; END DO; END DO ! i,j,k
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part not needed since diagonal of DvolSurfMat is zero!
  !htilde(:,k,i,j,k)=htilde_c(:,i,j,k) 
  htilde(:,k,i,j,k)=0.
  DO l=k+1,PP_N
    CALL PP_VolumeFluxAverageVec(U_in(:,i,j,k),U_in(:,i,j,l), &
                                  M_h(:,i,j,k), M_h(:,i,j,l), &
                             htilde(:,l,i,j,k)                )
    htilde(:,k,i,j,l)=htilde(:,l,i,j,k) !symmetric
  END DO!l=k+1,N
END DO; END DO; END DO ! i,j,k

END SUBROUTINE EvalAdvFluxAverage3D


!==================================================================================================================================
!> Compute flux differences between two points appling also directly the metrics
!==================================================================================================================================
SUBROUTINE EvalAdvFluxAverage(UL,UR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
#if PP_VolFlux==-1
USE MOD_Equation_Vars  ,ONLY:VolumeFluxAverageVec !pointer to flux averaging routine
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar          !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

CALL PP_VolumeFluxAverageVec(UL,UR,metric_L,metric_R,Fstar)

END SUBROUTINE EvalAdvFluxAverage
#endif /*PP_DiscType==2*/

!==================================================================================================================================
!> Computes the standard flux in normal-direction
!==================================================================================================================================
PURE SUBROUTINE StandardDGFlux(Fstar,UL,UR,normal)
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
PURE SUBROUTINE StandardDGFluxVec(UL,UR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:AdvVel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
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
PURE SUBROUTINE StandardDGFluxDealiasedMetricVec(UL,UR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
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

#undef PP_VolumeFluxAverageVec

END MODULE MOD_Flux_Average
