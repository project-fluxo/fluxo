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
!> Module for the positivity preservation routines
!==================================================================================================================================


MODULE MOD_PositivityPreservation
! MODULES
IMPLICIT NONE
PRIVATE
! ----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersPositivityPreservation
   MODULE PROCEDURE DefineParametersPositivityPreservation
END INTERFACE

INTERFACE InitPositivityPreservation
   MODULE PROCEDURE InitPositivityPreservation
END INTERFACE

#if POSITVITYPRES
INTERFACE MakeSolutionPositive
   MODULE PROCEDURE MakeSolutionPositive
END INTERFACE
#endif /*POSITVITYPRES*/

INTERFACE FinalizePositivityPreservation
   MODULE PROCEDURE FinalizePositivityPreservation
END INTERFACE

PUBLIC::DefineParametersPositivityPreservation
PUBLIC::InitPositivityPreservation
#if POSITIVITYPRES
PUBLIC::MakeSolutionPositive
#endif /*POSITIVITYPRES*/
PUBLIC::FinalizePositivityPreservation


real, allocatable :: Jac(:,:,:,:)  ! det of jacobian.
real, allocatable :: vol      (:)  ! Volume of cell.

!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersPositivityPreservation()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("PositivityPreservation")
END SUBROUTINE DefineParametersPositivityPreservation


SUBROUTINE InitPositivityPreservation()
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PositivityPreservation_Vars
USE MOD_ReadInTools
USE MOD_mesh_Vars                  , ONLY: sJ,nElems
USE MOD_Interpolation_Vars         , ONLY: wGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
integer :: i,j,k, l
!============================================================================================================================
IF (PositivityPreservationInitIsDone) THEN
  SWRITE(*,*) "InitPositivityPreservation already called."
  RETURN
END IF

allocate ( Jac(0:PP_N,0:PP_N,0:PP_N,nElems) )
allocate ( vol(nElems) )
vol = 0.
do l=1, nElems
  do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
    Jac(i,j,k,l) = 1./sJ(i,j,k,l)
    vol(l) = vol(l) + Jac(i,j,k,l) * wGP(i)*wGP(j)*wGP(k)
  end do       ; end do       ; end do
end do

PositivityPreservationInitIsDone = .TRUE.

if (PositivityPreservationInitFirst) then
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' INIT POSITIVITYPRESERVATION...'
  SWRITE(UNIT_stdOut,'(A)')' INIT POSITIVITYPRESERVATION DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')
  PositivityPreservationInitFirst = .FALSE.
end if

END SUBROUTINE InitPositivityPreservation

#if POSITIVITYPRES
SUBROUTINE MakeSolutionPositive(U)
!===================================================================================================================================
! Use the framework of Shu et. al. and enforce positivity of density and pressure
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Interpolation_Vars         , ONLY: wGP
use MOD_Equation_Vars              , only: Get_Pressure
USE MOD_mesh_Vars                  , ONLY: nElems
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems),INTENT(INOUT) :: U
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
INTEGER                                      :: i,j,k,l
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N) :: Uloc
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_N)         :: p
REAL,DIMENSION(PP_nVar)                      :: Umean
REAL                                         :: delta,rho_min,p_min,p_mean

DO l = 1,nElems
  Uloc    = U(:,:,:,:,l)
  Umean   = 0.
! Get minimum value of the density
  rho_min = HUGE(1.)
  DO k = 0,PP_N
    DO j = 0,PP_N
      DO i = 0,PP_N
        Umean   = Umean + Uloc(:,i,j,k)*wGP(i)*wGP(j)*wGP(k)*Jac(i,j,k,l)
        rho_min = MIN(rho_min,Uloc(1,i,j,k))
      END DO ! i
    END DO ! j
  END DO ! k
! Compute the cell average
  Umean = Umean/vol(l)
! Limit the density
  IF(rho_min.LT.0.) THEN
    delta = Umean(1)/(Umean(1) - rho_min)
    delta = 0.999*delta ! make delta > 0
    DO k = 0,PP_N
      DO j = 0,PP_N
        DO i = 0,PP_N
          Uloc(1,i,j,k) = Umean(1) + delta*(Uloc(1,i,j,k) - Umean(1))
        END DO ! i
      END DO ! j
    END DO ! k
  END IF
! Get the minimum value of the pressure and the new Umean
  Umean = 0.0
  p_min = HUGE(1.)
  DO k = 0,PP_N
    DO j = 0,PP_N
      DO i = 0,PP_N
        Umean   = Umean + Uloc(:,i,j,k)*wGP(i)*wGP(j)*wGP(k)*Jac(i,j,k,l)
        call Get_Pressure(Uloc(:,i,j,k),p(i,j,k))
        p_min = MIN(p_min,p(i,j,k))
      END DO ! i
    END DO ! j
  END DO ! k
  ! Compute the cell average
  Umean = Umean/vol(l)
  
! Limit the pressure
  IF (p_min.LT.0.) THEN
    call Get_Pressure(Umean,p_mean) ! We compute the pressure with the mean value, as we assume that Jensen's inequality holds
    delta  = p_mean/(p_mean - p_min)
    delta  = 0.999*delta ! make delta > 0
    DO k = 0,PP_N
      DO j = 0,PP_N
        DO i = 0,PP_N
          U(:,i,j,k,l) = Umean + delta*(Uloc(:,i,j,k) - Umean)
        END DO ! i
      END DO ! j
    END DO ! k
  END IF
END DO ! l

END SUBROUTINE MakeSolutionPositive
#endif /*POSITIVITYPRES*/


SUBROUTINE FinalizePositivityPreservation()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PositivityPreservation_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!============================================================================================================================

SDEALLOCATE (vol)
SDEALLOCATE (Jac)

IF (.NOT.PositivityPreservationInitIsDone) THEN
  WRITE(UNIT_stdOut,*) "InitPositivityPreservation was not called before."
  RETURN
END IF
PositivityPreservationInitIsDone = .FALSE.
END SUBROUTINE FinalizePositivityPreservation

END MODULE MOD_PositivityPreservation
