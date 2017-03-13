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
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

PUBLIC::Riemann
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Computes the numerical flux of the linear scalar advection equation on a given face, and also adds the 
!> diffusion part 1/2*(DiffC*grad U_L + DiffC*gradU_R) as needed for the BR1 / BR2 type of lifting.
!==================================================================================================================================
SUBROUTINE Riemann(F,U_L,U_R,&
#if PARABOLIC
                  gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif /*PARABOLIC*/
                  nv,t1,t2)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:AdvVel,DiffC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N),INTENT(IN) :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
#endif /*PARABOLIC*/
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N),t1(3,0:PP_N,0:PP_N),t2(3,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                                             :: LambdaMax(0:PP_N,0:PP_N)
!==================================================================================================================================
LambdaMax = AdvVel(1)*nv(1,:,:) +  AdvVel(2)*nv(2,:,:) + AdvVel(3)*nv(3,:,:)
! Compute the classic upwind flux into normal direction for each face GP
F(1,:,:) = 0.5*( (LambdaMax + ABS(LambdaMax))*U_L(1,:,:) + (LambdaMax-ABS(LambdaMax))*U_R(1,:,: ))
#if PARABOLIC
! Diffusion flux
F(1,:,:) = F(1,:,:)-DiffC*0.5*(  (gradUx_L(1,:,:)+gradUx_R(1,:,:))*nv(1,:,:) &
                               + (gradUy_L(1,:,:)+gradUy_R(1,:,:))*nv(2,:,:) & 
                               + (gradUz_L(1,:,:)+gradUz_R(1,:,:))*nv(3,:,:))
#endif /*PARABOLIC*/
END SUBROUTINE Riemann

END MODULE MOD_Riemann
