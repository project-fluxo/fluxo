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
ABSTRACT INTERFACE
  subroutine i_sub_EvalFluxTilde3D(iElem,U_in,M_f,M_g,M_h,ftilde,gtilde,htilde)
    USE MOD_DG_Vars       ,ONLY:nTotal_vol
    INTEGER,INTENT(IN ):: iElem                !< element number
    REAL,INTENT(IN )   :: U_in(8,1:nTotal_vol) !< solution state (conservative vars)
    REAL,INTENT(IN )   :: M_f( 3,1:nTotal_vol) !< metrics for ftilde                 
    REAL,INTENT(IN )   :: M_g( 3,1:nTotal_vol) !< metrics for gtilde                 
    REAL,INTENT(IN )   :: M_h( 3,1:nTotal_vol) !< metrics for htilde                 
    REAL,INTENT(OUT)   :: ftilde(8,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
    REAL,INTENT(OUT)   :: gtilde(8,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
    REAL,INTENT(OUT)   :: htilde(8,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
  end subroutine i_sub_EvalFluxTilde3D
END INTERFACE

procedure(i_sub_EvalFluxTilde3D), pointer :: EvalAdvFluxTilde3D => EvalFluxTilde3D

PUBLIC::EvalFluxTilde3D
PUBLIC::EvalAdvFluxTilde3D
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute advection flux of Maxwell's equations for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(iElem,U_in,M_f,M_g,M_h, &
                           ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc       ! PP_N
USE MOD_Equation_Vars ,ONLY:c2,c_corr,c_corr_c2
USE MOD_DG_Vars       ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN ):: iElem                !< element number
REAL,INTENT(IN )   :: U_in(8,1:nTotal_vol) !< solution state (conservative vars)
REAL,INTENT(IN )   :: M_f( 3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g( 3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h( 3,1:nTotal_vol) !< metrics for htilde                 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(8,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(8,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(8,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f(8),g(8),h(8)    ! Cartesian fluxes 
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  !A
  f(1) = U_in(8,i)*c_corr_c2    ! phi*chi*c^2
  f(2) = U_in(6,i)*c2           ! B3*c^2
  f(3) =-U_in(5,i)*c2           ! -B2*c^2
  f(4) = U_in(7,i)*c_corr       ! psi*c_corr
  f(5) =-U_in(3,i)              ! -E3
  f(6) = U_in(2,i)              ! E2
  f(7) = U_in(4,i)*c_corr_c2    ! B1*c_corr*c^2
  f(8) = U_in(1,i)*c_corr       ! E1*c_corr
  !B
  g(1) =-f(2)                ! -B3*c^2
  g(2) = f(1)                ! phi*c_corr*c^2
  g(3) = U_in(4,i)*c2           ! B1*c^2
  g(4) = U_in(3,i)              ! E3
  g(5) = f(4)                ! psi*c_corr
  g(6) =-U_in(1,i)              ! -E1
  g(7) = U_in(5,i)*c_corr_c2    ! B2*c_corr*c^2
  g(8) = U_in(2,i)*c_corr       ! E2*c_corr
  !C                              
  h(1) =-f(3)                ! B2*c^2
  h(2) =-g(3)                ! -B1*c^2
  h(3) = f(1)                ! phi*c_corr*c^2
  h(4) =-U_in(2,i)              ! -E2
  h(5) = U_in(1,i)              ! E1
  h(6) = f(4)                ! psi*c_corr
  h(7) = U_in(6,i)*c_corr_c2    ! B3*c_corr*c^2
  h(8) = U_in(3,i)*c_corr       ! E3*c_corr
  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,i) =   f(:)*M_f(1,i) + g(:)*M_f(2,i) + h(:)*M_f(3,i)
  gtilde(:,i) =   f(:)*M_g(1,i) + g(:)*M_g(2,i) + h(:)*M_g(3,i)
  htilde(:,i) =   f(:)*M_h(1,i) + g(:)*M_h(2,i) + h(:)*M_h(3,i)
      
END DO ! i
END SUBROUTINE EvalFluxTilde3D

END MODULE MOD_Flux
