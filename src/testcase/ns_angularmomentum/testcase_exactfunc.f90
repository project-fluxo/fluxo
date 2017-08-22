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
!> Subroutines defining one specific testcase with all necessary variables
!==================================================================================================================================
MODULE MOD_Testcase_ExactFunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE TestcaseExactFunc
  MODULE PROCEDURE TestcaseExactFunc
END INTERFACE

PUBLIC:: TestcaseExactFunc

CONTAINS

!==================================================================================================================================
!> Specifies all the initial conditions.
!==================================================================================================================================
SUBROUTINE TestcaseExactFunc(ExactFunction,t,x,resu,resu_t,resu_tt) 
! MODULES
USE MOD_Globals,      ONLY: Abort,CROSS
USE MOD_Equation_Vars,ONLY: Kappa,KappaM1,sKappaM1
USE MOD_Equation_Vars,ONLY: PrimToCons
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction    !< determines the exact function
REAL,INTENT(IN)                 :: t                !< current simulation time
REAL,INTENT(IN)                 :: x(3)             !< position in physical coordinates
REAL,INTENT(OUT)                :: resu(PP_nVar)    !< exact fuction evaluated at tIn, returning state in conservative variables
REAL,INTENT(OUT)                :: resu_t(PP_nVar)  !< first time deriv of exact fuction
REAL,INTENT(OUT)                :: resu_tt(PP_nVar) !< second time deriv of exact fuction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: du,xi(3),r(3),axis(3),xa,dT,prim(5)
!==================================================================================================================================
SELECT CASE(ExactFunction)
CASE(10001) !vortex with deplaced density bump,2d
  ASSOCIATE( rho_0 => 1.0, p_0=>1.0, A=>0.3, H=>0.1,xc=>[0.0,0.0,0.0] ,dx=>0.7 , &
             rho =>prim(1),vel=>prim(2:4),pres=>prim(5) )
  xi = (x(:)-xc(:))/H

  du     = A*EXP(0.5*(1.0-SUM(xi(1:2)*xi(1:2))))
  vel(1) = - du *xi(2)
  vel(2) =   du *xi(1)
  vel(3) = 0.0
  du     = 3.0*A*EXP(0.5*(1.0-SUM((xi(1:2)-dx)*(xi(1:2)-dx))))
  dT     = (KappaM1)*rho_0/(2.0*Kappa*p_0)*du*du
  rho    = (1.0-dT)**(sKappaM1)
  pres   = p_0*rho**(Kappa)
  rho    = rho_0*rho
  END ASSOCIATE
  CALL PrimToCons(Prim(:),Resu(:))
  Resu_t=0.
  Resu_tt=0.
CASE(10002) !vortex with deplaced density bump
  ASSOCIATE( rho_0 => 1.0, p_0=>1.0, A=>0.3, H=>0.1,xc=>[0.0,0.0,0.0] ,dx=>0.6*[-0.5,-0.5,0.9] ,  &
             rho =>prim(1),vel=>prim(2:4),pres=>prim(5) )
  axis=[1.,1.,1.]
  axis = axis/SQRT(SUM(axis*axis)) !normalize
  xi = (x(:)-xc(:))/H
  xa = SUM(xi*axis)
  r  = xi-xa*axis

  du     = A*EXP(0.5*(1.0-SUM(r*r)-xa*xa*4.))
!  vel(1) = - du *xi(2)
!  vel(2) =   du *xi(1)
!  vel(3) = 0.0
  vel(:) = du*CROSS(r,axis)
  du     = 3.0*A*EXP(0.5*(1.0-SUM((xi-dx)*(xi-dx))))
  dT     = (KappaM1)*rho_0/(2.0*Kappa*p_0)*du*du
  rho    = (1.0-dT)**(sKappaM1)
  pres   = p_0*rho**(Kappa)
  rho    = rho_0*rho
  END ASSOCIATE
  CALL PrimToCons(Prim(:),Resu(:))
  Resu_t=0.
  Resu_tt=0.
CASE DEFAULT
  CALL abort(__STAMP__,'Exactfunction not specified!')
END SELECT !ExactFunction
END SUBROUTINE TestcaseExactFunc

END MODULE MOD_Testcase_ExactFunc
