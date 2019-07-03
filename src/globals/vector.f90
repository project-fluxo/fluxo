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
!> Routines to map simple operations (copy,scalar multiply) from multi-dimensional arrays onto 1D arrays for better performance.
!==================================================================================================================================
MODULE MOD_Vector
! MODULES
IMPLICIT NONE
PUBLIC

!==================================================================================================================================
CONTAINS


!==================================================================================================================================
!> Y=0
!==================================================================================================================================
SUBROUTINE VNullify(nTotal,Vec)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal                                               !< vector length
REAL,INTENT(OUT)      :: Vec(nTotal)                                          !< input vector
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Vec=0.
END SUBROUTINE VNullify


!==================================================================================================================================
!> Y=a
!==================================================================================================================================
SUBROUTINE VSetConst(nTotal,Vec,Const)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal                                               !< vector length
REAL,INTENT(OUT)      :: Vec(nTotal)                                          !< input vector
REAL,INTENT(IN)       :: Const                                                !< constant to set
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
Vec=Const
END SUBROUTINE VSetConst


!==================================================================================================================================
!> Y=X
!==================================================================================================================================
SUBROUTINE VCopy(nTotal,VecOut,VecIn)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal                                               !< vector length
REAL,INTENT(IN)       :: VecIn(nTotal)                                        !< input vector
REAL,INTENT(OUT)      :: VecOut(nTotal)                                       !< output vector
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
VecOut=VecIn
END SUBROUTINE VCopy


!==================================================================================================================================
!> Y=AY+BX
!==================================================================================================================================
SUBROUTINE VAXPBY(nTotal,VecOut,VecIn,ConstOut,ConstIn)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)    :: nTotal                                               !< vector length
REAL,INTENT(IN)       :: VecIn(nTotal)                                        !< input vector
REAL,INTENT(INOUT)    :: VecOut(nTotal)                                       !< output vector
REAL,INTENT(IN),OPTIONAL :: ConstIn                                           !< constant to multiply with input vec
REAL,INTENT(IN),OPTIONAL :: ConstOut                                          !< constant to multiply with output vec
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i
!==================================================================================================================================
IF(    PRESENT(ConstIn).AND.(.NOT.PRESENT(ConstOut)))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)+VecIn(i)*ConstIn
  END DO
ELSEIF(PRESENT(ConstOut).AND.(.NOT.PRESENT(ConstIn)))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)*ConstOut+VecIn(i)
  END DO
ELSEIF(PRESENT(ConstIn).AND.PRESENT(ConstOut))THEN
  DO i=1,nTotal
    VecOut(i)=VecOut(i)*ConstOut+VecIn(i)*ConstIn
  END DO
ELSE
  DO i=1,nTotal
    VecOut(i)=VecOut(i)+VecIn(i)
  END DO
END IF
END SUBROUTINE VAXPBY


SUBROUTINE V2D_M_V1D(dim1,nTotal,Y,X)
!===================================================================================================================================
! Y(1:dim1,i)=Y(1:dim1,i)*X(i)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: dim1
INTEGER,INTENT(IN)    :: nTotal
REAL,INTENT(IN)       :: X(nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: Y(dim1,nTotal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER               :: i
!===================================================================================================================================
  DO i=1,nTotal
    Y(:,i)=Y(:,i)*X(i)
  END DO
END SUBROUTINE V2D_M_V1D


END MODULE MOD_Vector
