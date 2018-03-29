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
#include "defines.h"

!==================================================================================================================================
!> Contains the routine EvalFlux3D which computes the complete flux f,g,h for all DOFs in one Element: used in volume integral
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE EvalFluxTilde3D
  MODULE PROCEDURE EvalFluxTilde3D
END INTERFACE


#if PARABOLIC
INTERFACE EvalDiffFluxTilde3D
  MODULE PROCEDURE EvalDiffFluxTilde3D
END INTERFACE

INTERFACE EvalLiftingVolumeFlux
  MODULE PROCEDURE EvalLiftingVolumeFlux
END INTERFACE

INTERFACE EvalLiftingSurfFlux
  MODULE PROCEDURE EvalLiftingSurfFlux
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalFluxTilde3D
#if PARABOLIC
PUBLIC::EvalDiffFluxTilde3D
PUBLIC::EvalLiftingVolumeFlux
PUBLIC::EvalLiftingSurfFlux
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute linear scalar advection & diffusion fluxes with using the solution and its gradient for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Equation_Vars ,ONLY:AdvVel
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
#if PARABOLIC
USE MOD_Equation_Vars ,ONLY:DiffC
USE MOD_Lifting_Vars  ,ONLY:gradPx,gradPy,gradPz
#endif
#ifdef OPTIMIZED
USE MOD_DG_Vars       ,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem !< Determines the actual element
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed flux f(iVar,i,j,k)
REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed flux g(iVar,i,j,k)
REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i 
#ifndef OPTIMIZED
INTEGER             :: j,k
#endif
REAL                :: f,g,h
!==================================================================================================================================

#ifdef OPTIMIZED
DO i=0,nTotal_vol-1
#else /*OPTIMIZED*/
DO k=0,PP_N;  DO j=0,PP_N; DO i=0,PP_N
#endif /*OPTIMIZED*/

  f=AdvVel(1)*U(1,PP_IJK,iElem)
  g=AdvVel(2)*U(1,PP_IJK,iElem)
  h=AdvVel(3)*U(1,PP_IJK,iElem)
#if PARABOLIC
  f=f-DiffC*gradPx(1,PP_IJK,iElem)
  g=g-DiffC*gradPy(1,PP_IJK,iElem)
  h=h-DiffC*gradPz(1,PP_IJK,iElem)
#endif /*PARABOLIC*/

  ftilde(1,PP_IJK) =   f*Metrics_fTilde(1,PP_IJK,iElem) &
                     + g*Metrics_fTilde(2,PP_IJK,iElem) &
                     + h*Metrics_fTilde(3,PP_IJK,iElem)
  gtilde(1,PP_IJK) =   f*Metrics_gTilde(1,PP_IJK,iElem) &
                     + g*Metrics_gTilde(2,PP_IJK,iElem) &
                     + h*Metrics_gTilde(3,PP_IJK,iElem)
  htilde(1,PP_IJK) =   f*Metrics_hTilde(1,PP_IJK,iElem) &
                     + g*Metrics_hTilde(2,PP_IJK,iElem) &
                     + h*Metrics_hTilde(3,PP_IJK,iElem)

#ifdef OPTIMIZED
END DO ! i
#else /*OPTIMIZED*/
END DO; END DO; END DO ! i,j,k
#endif /*OPTIMIZED*/

END SUBROUTINE EvalFluxTilde3D


#if PARABOLIC
!==================================================================================================================================
!> Compute linear scalar  diffusion fluxes with using the gradient of the solution  for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFluxTilde3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Lifting_Vars  ,ONLY:gradPx,gradPy,gradPz
USE MOD_Equation_Vars ,ONLY:DiffC
#ifdef OPTIMIZED
USE MOD_DG_Vars,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem !< Determines the actual element
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed flux f(iVar,i,j,k)
REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed flux g(iVar,i,j,k)
REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! general curved metrics
ftilde(1,:,:,:) = -DiffC*(  Metrics_fTilde(1,:,:,:,iElem)*gradPx(1,:,:,:,iElem) &
                          + Metrics_fTilde(2,:,:,:,iElem)*gradPy(1,:,:,:,iElem) &
                          + Metrics_fTilde(3,:,:,:,iElem)*gradPz(1,:,:,:,iElem) )
gtilde(1,:,:,:) = -DiffC*(  Metrics_gTilde(1,:,:,:,iElem)*gradPx(1,:,:,:,iElem) &
                          + Metrics_gTilde(2,:,:,:,iElem)*gradPy(1,:,:,:,iElem) &
                          + Metrics_gTilde(3,:,:,:,iElem)*gradPz(1,:,:,:,iElem) )
htilde(1,:,:,:) = -DiffC*(  Metrics_hTilde(1,:,:,:,iElem)*gradPx(1,:,:,:,iElem) &
                          + Metrics_hTilde(2,:,:,:,iElem)*gradPy(1,:,:,:,iElem) &
                          + Metrics_hTilde(3,:,:,:,iElem)*gradPz(1,:,:,:,iElem) )
END SUBROUTINE EvalDiffFluxTilde3D


!==================================================================================================================================
!> Compute the lifting flux depending on the variable to be used for the gradient
!> for linadv, only U as variable
!==================================================================================================================================
SUBROUTINE EvalLiftingVolumeFlux(iElem,Flux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:nTotal_vol,U
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)     :: iElem       !< current element
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)       :: flux(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< lifting flux, depending on lifting_var
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
  Flux=U(:,:,:,:,iElem)
END SUBROUTINE EvalLiftingVolumeFlux


!==================================================================================================================================
!> Compute the average lifting surface flux in strong form U*=1/2(U_s+U_m)  (BR1/BR2),
!> careful, we use the strong form, so that the surface flux becomes: F=U^*-U_m = 1/2(U_s-U_m)
!> for linadv, only U is variable
!==================================================================================================================================
SUBROUTINE EvalLiftingSurfFlux(SideID,Flux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY: U_master,U_slave
USE MOD_Mesh_Vars ,ONLY: SurfElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)     :: SideID       !< current side index
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)       :: flux(PP_nVar,0:PP_N,0:PP_N) !< lifting surface flux, depending on lifting_var
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q 
!==================================================================================================================================
  !strong lifting flux: 
  DO q=0,PP_N; DO p=0,PP_N
      Flux(:,p,q)=0.5*(U_slave(:,p,q,SideID)-U_master(:,p,q,SideID))*SurfElem(p,q,SideID)
  END DO; END DO
END SUBROUTINE EvalLiftingSurfFlux

#endif /*PARABOLIC*/

END MODULE MOD_Flux
