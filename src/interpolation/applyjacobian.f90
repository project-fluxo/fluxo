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
!> TODO
!==================================================================================================================================

MODULE MOD_ApplyJacobianCons
IMPLICIT NONE
PRIVATE

INTERFACE ApplyJacobianCons
   MODULE PROCEDURE ApplyJacobianCons
END INTERFACE

PUBLIC::ApplyJacobianCons

CONTAINS
!==================================================================================================================================
!> Convert solution between physical <-> reference space (DG elements only), input will be overwritten with transformed solution 
!==================================================================================================================================
SUBROUTINE ApplyJacobianCons(U,toPhysical)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars ,ONLY: sJ,nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
LOGICAL,INTENT(IN) :: toPhysical
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
!==================================================================================================================================
IF(toPhysical)THEN
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      U(:,i,j,k,iElem)=U(:,i,j,k,iElem)*sJ(i,j,k,iElem)
    END DO; END DO; END DO
  END DO
ELSE
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      U(:,i,j,k,iElem)=U(:,i,j,k,iElem)/sJ(i,j,k,iElem)
    END DO; END DO; END DO
  END DO
END IF
END SUBROUTINE ApplyJacobianCons

END MODULE MOD_ApplyJacobianCons

