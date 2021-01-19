!==================================================================================================================================
! Copyright (c) 2018 - 2020 Alexander Astanin
! Copyright (c) 2020 - 2020 Andrés Rueda
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
#include "amr_f.h"

#if USE_AMR
!==================================================================================================================================
!> Contains routines for AMR as trecking function
!==================================================================================================================================
MODULE MOD_AMR_tracking

    INTERFACE PerformAMR
        MODULE PROCEDURE PerformAMR
    END INTERFACE


    ! IMPLICIT NONE
    INTEGER :: doLBalance = 0
    INTEGER :: Count = 0
CONTAINS

! FUNCTION GetShockCapturing(Uin) result(eta_dof)
!     USE MOD_PreProc
!     USE MOD_ChangeBasis,            ONLY : ChangeBasis3D
!     IMPLICIT NONE
!     REAL, DIMENSION(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), INTENT(IN) :: Uin
!     REAL, DIMENSION(1:1, 0:PP_N, 0:PP_N, 0:PP_N) :: Umod
!     REAL, DIMENSION(0:PP_N, 0:PP_N) :: Vdm_Leg, sVdm_Leg
!     REAL :: LU, LUM1, LUM2, LU_N, LU_NM1!, eta_dof, eta_min, eta_max, eps0, RhoInf, Pinf, RhoMax, RhoMin, Xmin(3), Xmax(3), Abst
!     REAL eta_dof
!     CALL ChangeBasis3D(1, PP_N, PP_N, sVdm_Leg, Uin(1:1,:,:,:), Umod)
!     LU = SUM(Umod(1, :, :, :)**2)
!     LUM1 = SUM(Umod(1, 0:PP_N - 1, 0:PP_N - 1, 0:PP_N - 1)**2)
!     LUM2 = SUM(Umod(1, 0:PP_N - 2, 0:PP_N - 2, 0:PP_N - 2)**2)
!     LU_N = LU - LUM1
!     LU_NM1 = LUM1 - LUM2
!     ! DOF energy indicator
!     eta_dof = LOG10(MAX(LU_N / LU, LU_NM1 / LUM1, TINY(1.0)))
    
! END FUNCTION 

!==================================================================================================================================
!> Specifies the initial AMR refinement
!==================================================================================================================================
subroutine InitialAMRRefinement()
  USE MOD_PreProc     , only: PP_N
  use MOD_AMR_Vars    , only: InitialRefinement, UseAMR, MaxLevel, MinLevel, IniHalfwidthAMR
  use MOD_AMR         , only: RunAMR
  use MOD_Mesh_Vars   , only: nElems, Elem_xGP
  use MOD_DG_Vars , only: U
  use MOD_Globals, only: MPIRoot
  implicit none
  !-local-variables-----------------------------------------
  real    :: r
  logical :: RefineElem
  integer :: iter
  integer :: iElem, i,j,k
  integer, allocatable :: ElemToRefineAndCoarse(:)
  !---------------------------------------------------------
  
  if (.not. UseAMR) return
  
  select case (InitialRefinement)  
    case default ! Use the default indicator up to max-level
      do iter = 1,MaxLevel
        call PerformAMR()
        call InitData()
      end do
    
    case(1) ! Refine any element containing a node in the spherewith radius r=IniHalfwidthAMR to the MaxLevel
      do iter = 1,MaxLevel
        allocate (ElemToRefineAndCoarse(1:nElems))!
        
        ! Fill ElemToRefineAndCoarse
        ! --------------------------
        do iElem=1, nElems  
          ! Check if the element has a node on the desired region
          RefineElem = .FALSE.
          outer: do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
            r = sqrt(sum(Elem_xGP(:,i,j,k,iElem)**2))
            if (r <= IniHalfwidthAMR) then
              RefineElem = .TRUE.
              exit outer
            end if
          end do       ; end do       ; end do outer
          
          ! Set that element for refinement
          if (RefineElem) then
            ElemToRefineAndCoarse(iElem) = MaxLevel
          else
            ElemToRefineAndCoarse(iElem) = MinLevel
          end if
          
        end do
        
        ! Refine
        ! ------
        CALL RunAMR(ElemToRefineAndCoarse)
        
        deallocate (ElemToRefineAndCoarse)
      end do
      call InitData()
  end select
  

end subroutine


    SUBROUTINE PerformAMR()
        !   USE MOD_AMR_vars,            ONLY: P4EST_PTR, CONNECTIVITY_PTR
        USE MOD_PreProc
        USE MOD_Globals,                ONLY : MPIroot
        USE MOD_DG_Vars,                ONLY : U
        USE MOD_AMR,                    ONLY : RunAMR, LoadBalancingAMR, SaveMesh;
        USE MOD_Mesh_Vars,              ONLY : nElems
        USE MOD_Basis,                  ONLY : BuildLegendreVdm
        USE MOD_AMR_Vars,               ONLY : MinLevel, MaxLevel, RefineVal, CoarseVal, AMR_Indicator
        USE MOD_P4EST,                  ONLY: SaveP4est
#if SHOCK_NFVSE
        use MOD_NFVSE_Vars,             only: alpha, alpha_max, SpacePropFactor
#endif /*SHOCK_NFVSE*/
        ! USE MOD_Equation_Vars,      ONLY: kappaM1, RefStatePrim, IniRefState
        IMPLICIT NONE
        ! SAVE
        !Local variables
        INTEGER, ALLOCATABLE, TARGET :: ElemToRefineAndCoarse(:) ! positive Number - refine, negative - coarse, 0 - do nothing
        INTEGER :: iElem
        
        REAL    :: eta_dof(nElems)
            
        ALLOCATE(ElemToRefineAndCoarse(1:nElems))!
        ElemToRefineAndCoarse = MinLevel;
        
        eta_dof = AMR_Indicator % compute(U)
        DO iElem = 1, nElems
          IF (eta_dof(iElem) .GE. RefineVal) THEN
            ElemToRefineAndCoarse(iElem) = MaxLevel
          ELSE IF (eta_dof(iElem) .LE. CoarseVal) THEN
            ElemToRefineAndCoarse(iElem) = -MinLevel - 1
          ELSE
            ElemToRefineAndCoarse(iElem) = MinLevel
          END IF
#if SHOCK_NFVSE
          ! Always refine if the shock capturing is firing
          if ( alpha(iElem) > min(0.1,alpha_max*max(0.5,SpacePropFactor))) ElemToRefineAndCoarse(iElem) = MaxLevel
#endif /*SHOCK_NFVSE*/
        ENDDO

        CALL RunAMR(ElemToRefineAndCoarse);
        
        Deallocate(ElemToRefineAndCoarse)
    
    END SUBROUTINE PerformAMR

    SUBROUTINE InitData()
        USE MOD_PreProc
        USE MOD_DG_Vars, ONLY : U
        USE MOD_Equation_Vars, ONLY :IniExactFunc
        USE MOD_Equation, ONLY : FillIni
        USE MOD_Restart_Vars, only: DoRestart
        USE MOD_AMR, ONLY : RunAMR
        IMPLICIT NONE

        if (.not. DoRestart) call FillIni(IniExactFunc,U)

    END SUBROUTINE InitData

END MODULE MOD_AMR_tracking
#endif
