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

!INTERFACE EvalFluxTilde3D
!  MODULE PROCEDURE EvalFluxTilde3D
!END INTERFACE


#if PARABOLIC
!INTERFACE EvalDiffFluxTilde3D
!  MODULE PROCEDURE EvalDiffFluxTilde3D
!END INTERFACE

!INTERFACE EvalLiftingVolumeFlux
!  MODULE PROCEDURE EvalLiftingVolumeFlux
!END INTERFACE

INTERFACE EvalLiftingSurfFlux
  MODULE PROCEDURE EvalLiftingSurfFlux
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalAdvFluxTilde3D
PUBLIC::EvalFluxTilde3D
#if PARABOLIC
PUBLIC::EvalDiffFluxTilde3D
PUBLIC::EvalLiftingVolumeFlux
PUBLIC::EvalLiftingSurfFlux
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Compute linadvdiff transformed fluxes using conservative variables and derivatives for every volume Gauss point.
!> directly apply metrics and output the tranformed flux 
!==================================================================================================================================
SUBROUTINE EvalAdvFluxTilde3D(iElem,U_in,M_f,M_g,M_h,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:AdvVel
USE MOD_DG_Vars       ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN ):: iElem                !< element number
REAL,INTENT(IN )   :: U_in(1,1:nTotal_vol) !< solution state (conservative vars)
REAL,INTENT(IN )   :: M_f( 3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g( 3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h( 3,1:nTotal_vol) !< metrics for htilde                 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(1,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(1,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(1,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ftilde(1,i) =   (AdvVel(1)*M_f(1,i) + AdvVel(2)*M_f(2,i) + AdvVel(3)*M_f(3,i))*U_in(1,i)
  gtilde(1,i) =   (AdvVel(1)*M_g(1,i) + AdvVel(2)*M_g(2,i) + AdvVel(3)*M_g(3,i))*U_in(1,i)
  htilde(1,i) =   (AdvVel(1)*M_h(1,i) + AdvVel(2)*M_h(2,i) + AdvVel(3)*M_h(3,i))*U_in(1,i)
END DO ! i
END SUBROUTINE EvalAdvFluxTilde3D


!==================================================================================================================================
!> Compute linear scalar advection & diffusion fluxes with using the solution and its gradient for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(iElem,U_in,M_f,M_g,M_h, &
#if PARABOLIC
                           gradPx_in,gradPy_in,gradPz_in,&
#endif /*PARABOLIC*/
                           ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:AdvVel
#if PARABOLIC
USE MOD_Equation_Vars ,ONLY:DiffC
#endif
USE MOD_DG_Vars       ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN ):: iElem                !< element number
REAL,INTENT(IN )   :: U_in(1,1:nTotal_vol) !< solution state (conservative vars)
REAL,INTENT(IN )   :: M_f( 3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g( 3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h( 3,1:nTotal_vol) !< metrics for htilde                 
#if PARABOLIC
REAL,INTENT(IN )   :: gradPx_in(1,1:nTotal_vol) !< gradient x in primitive variables 
REAL,INTENT(IN )   :: gradPy_in(1,1:nTotal_vol) !< gradient y in primitive variables 
REAL,INTENT(IN )   :: gradPz_in(1,1:nTotal_vol) !< gradient z in primitive variables 
#endif /*PARABOLIC*/
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(1,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(1,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(1,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i 
REAL                :: f,g,h
!==================================================================================================================================

DO i=1,nTotal_vol

#if PARABOLIC
  f=AdvVel(1)*U_in(1,i)-DiffC*gradPx_in(1,i)
  g=AdvVel(2)*U_in(1,i)-DiffC*gradPy_in(1,i)
  h=AdvVel(3)*U_in(1,i)-DiffC*gradPz_in(1,i)
#else
  f=AdvVel(1)*U_in(1,i)
  g=AdvVel(2)*U_in(1,i)
  h=AdvVel(3)*U_in(1,i)
#endif /*PARABOLIC*/
  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(1,i) =   f*M_f(1,i) + g*M_f(2,i) + h*M_f(3,i)
  gtilde(1,i) =   f*M_g(1,i) + g*M_g(2,i) + h*M_g(3,i)
  htilde(1,i) =   f*M_h(1,i) + g*M_h(2,i) + h*M_h(3,i)
END DO ! i
END SUBROUTINE EvalFluxTilde3D


#if PARABOLIC
!==================================================================================================================================
!> Compute linear scalar  diffusion fluxes with using the gradient of the solution  for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFluxTilde3D(iElem,U_in,M_f,M_g,M_h,gradPx_in,gradPy_in,gradPz_in,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:DiffC
USE MOD_DG_Vars,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN ):: iElem                     !< element number
REAL,INTENT(IN )   :: U_in(     1,1:nTotal_vol) !< state in conservative variables
REAL,INTENT(IN )   :: M_f(      3,1:nTotal_vol) !< metrics for ftilde  
REAL,INTENT(IN )   :: M_g(      3,1:nTotal_vol) !< metrics for gtilde  
REAL,INTENT(IN )   :: M_h(      3,1:nTotal_vol) !< metrics for htilde  
REAL,INTENT(IN )   :: gradPx_in(1,1:nTotal_vol) !< gradient x in primitive variables 
REAL,INTENT(IN )   :: gradPy_in(1,1:nTotal_vol) !< gradient y in primitive variables 
REAL,INTENT(IN )   :: gradPz_in(1,1:nTotal_vol) !< gradient z in primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(1,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(1,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(1,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ftilde(1,i) = -DiffC*(gradPx_in(1,i)*M_f(1,i) + gradPy_in(1,i)*M_f(2,i) + gradPz_in(1,i)*M_f(3,i))
  gtilde(1,i) = -DiffC*(gradPx_in(1,i)*M_g(1,i) + gradPy_in(1,i)*M_g(2,i) + gradPz_in(1,i)*M_g(3,i))
  htilde(1,i) = -DiffC*(gradPx_in(1,i)*M_h(1,i) + gradPy_in(1,i)*M_h(2,i) + gradPz_in(1,i)*M_h(3,i))
END DO ! i
END SUBROUTINE EvalDiffFluxTilde3D


!==================================================================================================================================
!> Compute the lifting flux depending on the variable to be used for the gradient
!> for linadv, only U as variable
!==================================================================================================================================
SUBROUTINE EvalLiftingVolumeFlux(U_in,Flux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN )   :: U_in(PP_nVar,1:nTotal_vol) !< state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: flux(PP_nVar,1:nTotal_vol) !< lifting flux, depending on lifting_var
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
  Flux=U_in
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
