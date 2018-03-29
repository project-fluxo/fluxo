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

PUBLIC::EvalFluxTilde3D
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute advection flux of Maxwell's equations for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc       ! PP_N
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars ,ONLY:c2,c_corr,c_corr_c2
USE MOD_DG_Vars       ,ONLY:U
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem  !< current element treated in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(8,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed flux f(iVar,i,j,k)
REAL,DIMENSION(8,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed flux g(iVar,i,j,k)
REAL,DIMENSION(8,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f(8),g(8),h(8)    ! Cartesian fluxes 
INTEGER             :: i,j,k 
!==================================================================================================================================
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ASSOCIATE(Uin=>U(:,i,j,k,iElem))
  !A
  f(1) = Uin(8)*c_corr_c2    ! phi*chi*c^2
  f(2) = Uin(6)*c2           ! B3*c^2
  f(3) =-Uin(5)*c2           ! -B2*c^2
  f(4) = Uin(7)*c_corr       ! psi*c_corr
  f(5) =-Uin(3)              ! -E3
  f(6) = Uin(2)              ! E2
  f(7) = Uin(4)*c_corr_c2    ! B1*c_corr*c^2
  f(8) = Uin(1)*c_corr       ! E1*c_corr
  !B
  g(1) =-f(2)                ! -B3*c^2
  g(2) = f(1)                ! phi*c_corr*c^2
  g(3) = Uin(4)*c2           ! B1*c^2
  g(4) = Uin(3)              ! E3
  g(5) = f(4)                ! psi*c_corr
  g(6) =-Uin(1)              ! -E1
  g(7) = Uin(5)*c_corr_c2    ! B2*c_corr*c^2
  g(8) = Uin(2)*c_corr       ! E2*c_corr
  !C                              
  h(1) =-f(3)                ! B2*c^2
  h(2) =-g(3)                ! -B1*c^2
  h(3) = f(1)                ! phi*c_corr*c^2
  h(4) =-Uin(2)              ! -E2
  h(5) = Uin(1)              ! E1
  h(6) = f(4)                ! psi*c_corr
  h(7) = Uin(6)*c_corr_c2    ! B3*c_corr*c^2
  h(8) = Uin(3)*c_corr       ! E3*c_corr
  END ASSOCIATE !Uin
  ! general curved metrics
  ftilde(:,i,j,k) =   f(:)*Metrics_fTilde(1,i,j,k,iElem)  &
                    + g(:)*Metrics_fTilde(2,i,j,k,iElem)  &
                    + h(:)*Metrics_fTilde(3,i,j,k,iElem)
  gtilde(:,i,j,k) =   f(:)*Metrics_gTilde(1,i,j,k,iElem)  &
                    + g(:)*Metrics_gTilde(2,i,j,k,iElem)  &
                    + h(:)*Metrics_gTilde(3,i,j,k,iElem)
  htilde(:,i,j,k) =   f(:)*Metrics_hTilde(1,i,j,k,iElem)  &
                    + g(:)*Metrics_hTilde(2,i,j,k,iElem)  &
                    + h(:)*Metrics_hTilde(3,i,j,k,iElem)
      
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalFluxTilde3D

END MODULE MOD_Flux
