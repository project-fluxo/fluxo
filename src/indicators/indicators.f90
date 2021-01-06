!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
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
!> This module only initializes the equation specific parameters and computes analytical functions and the evaluation of the source
!==================================================================================================================================
module MOD_Indicators
  use MOD_Equation_Vars, only: i_indicatorFunction
  implicit none
  
  private
  public :: Indicator_PerssonPeraire
  
!============================================================================================================================
!> Class for the modal shock sensor of Persson and Peraire
!> Persson, P. O.; Peraire, J. (2006). "Sub-cell shock capturing for discontinuous Galerkin methods". In 44th AIAA Aerospace Sciences Meeting and Exhibit (p. 112).
!============================================================================================================================
  type Indicator_PerssonPeraire
    logical            :: InitIsDone =.FALSE. ! True if indicator has been initialized
    logical            :: ReturnLog =.FALSE.  ! True if the
    character(len=255) :: IndicatorQuantity   ! Name of the indicator quantity
    procedure(i_indicatorFunction), nopass, pointer :: GetIndicatorQuantity => null()
    contains
      procedure :: construct  => PerssonPeraire_Construct
      procedure :: compute    => PerssonPeraire_Compute
      procedure :: destruct   => PerssonPeraire_Destruct
      procedure, private :: GetIndicator_3D  => PerssonPeraire_GetIndicator_3D
  end type Indicator_PerssonPeraire

contains
!===================================================================================================================================
!> Constructor for Indicator_PerssonPeraire 
!===================================================================================================================================
  subroutine PerssonPeraire_Construct(this,IndicatorQuantity,ReturnLog)
    USE MOD_Globals
    USE MOD_PreProc
    USE MOD_Indicators_vars   , only: NumberOfPerssonInd
    use MOD_Equation_Vars     , only: SetIndicatorFunction
    USE MOD_Interpolation_Vars, only: xGP, InterpolationInitIsDone
   
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    !-arguments------------------------------------------------------------------------------------------------------------------
    class(Indicator_PerssonPeraire), intent(inout) :: this
    character(len=255)             , intent(in)    :: IndicatorQuantity
    logical, optional              , intent(in)    :: ReturnLog
    !============================================================================================================================
    
    ! Initial checks
    IF (this % InitIsDone .OR. (.NOT.InterpolationInitIsDone)) THEN
      SWRITE(*,*) "InitIndicator not ready to be called or already called."
      RETURN
    END IF
    IF (PP_N.LT.2) THEN
      CALL abort(__STAMP__,'Polynomial Degree too small for Indicator!',999,999.)
      RETURN
    END IF
    
    ! Compute Vandermonde matrix to change basis
    CALL InitBasisTrans(PP_N,xGP)
    
    ! Set ReturnLog if necessary
    if ( present(ReturnLog) ) this % ReturnLog = ReturnLog
    
    ! Assign indicator quantity
    this % IndicatorQuantity = IndicatorQuantity
    call SetIndicatorFunction(IndicatorQuantity,this % GetIndicatorQuantity)
    
    ! All is done
    this % InitIsDone = .TRUE.
    NumberOfPerssonInd = NumberOfPerssonInd+1
    
  end subroutine PerssonPeraire_Construct
!===================================================================================================================================
!> Destructor for Indicator_PerssonPeraire 
!===================================================================================================================================
  subroutine PerssonPeraire_Destruct(this)
    use MOD_Indicators_vars, only: NumberOfPerssonInd, sVdm_Leg
    implicit none
    class(Indicator_PerssonPeraire), intent(inout) :: this
    
    if (.not. this % InitIsDone) return
    
    NumberOfPerssonInd = NumberOfPerssonInd-1
    if (NumberOfPerssonInd<1) deallocate(sVdm_Leg)
    
    this % GetIndicatorQuantity => null()
    this % InitIsDone =.FALSE.
    this % ReturnLog =.FALSE.
    this % IndicatorQuantity = ''
  end subroutine PerssonPeraire_Destruct
    
!============================================================================================================================
!> Compute the indicator of Persson and Peraire
!============================================================================================================================
  function PerssonPeraire_Compute(this,U) RESULT (eta)
    use MOD_PreProc
    use MOD_ChangeBasis                     ,only: ChangeBasis3D
    use MOD_Mesh_Vars                       ,only: nElems
    use MOD_Indicators_vars                 ,only: sVdm_Leg
    implicit none
    ! Arguments
    !---------------------------------------------------------------------------------------------------------------------------------
    class(Indicator_PerssonPeraire), intent(in) :: this
    real                           , intent(in) :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
    real                                        :: eta(nElems)
    ! Local variables
    !---------------------------------------------------------------------------------------------------------------------------------
    real, dimension(1:1,0:PP_N,0:PP_N,0:PP_N) :: Uind,Umod
    real                                      :: LU,LUM1,LUM2,LU_N,LU_NM1
    integer                                   :: eID
    !---------------------------------------------------------------------------------------------------------------------------------
    
    do eID=1, nElems
      ! Get the indicator variable (Uind)
      call this % GetIndicator_3D(Uind,U(:,:,:,:,eID))
      
      ! Transform Uind into modal Legendre interpolant Umod
      CALL ChangeBasis3D(1,PP_N,PP_N,sVdm_Leg,Uind,Umod)
      
      ! Compute (truncated) error norms
      LU     = SUM(Umod(1,0:PP_N  ,0:PP_N  ,0:PP_N  )**2)
      LUM1   = SUM(Umod(1,0:PP_N-1,0:PP_N-1,0:PP_N-1)**2)
      LUM2   = SUM(Umod(1,0:PP_N-2,0:PP_N-2,0:PP_N-2)**2)
      LU_N   = LU-LUM1
      LU_NM1 = LUM1-LUM2
  
      ! DOF energy indicator
      eta(eID) = MAX(LU_N/LU,LU_NM1/LUM1)
      
    end do
    
    where (eta < epsilon(1.0)) eta = epsilon(eta(eID))
    
    if (this % ReturnLog) eta = log10(eta)
    
  end function PerssonPeraire_Compute

  
!============================================================================================================================
!> Get the shock indicator quantity for all degrees of freedom of an element
!============================================================================================================================
  pure subroutine PerssonPeraire_GetIndicator_3D(this,Uind,U)
    use MOD_PreProc
    implicit none
    !-arguments-------------------------------------------
    class(Indicator_PerssonPeraire), intent(in) :: this
    real, intent(in)  :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    real, intent(out) :: Uind (1:1,0:PP_N,0:PP_N,0:PP_N)
    !-local-variables-------------------------------------
    integer :: i,j,k
    !-----------------------------------------------------
    
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      call this % GetIndicatorQuantity(U(:,i,j,k),Uind(1,i,j,k))
    end do       ; end do       ; end do
  end subroutine PerssonPeraire_GetIndicator_3D
  
!===================================================================================================================================
!> Initialize Vandermodematrix for basis transformation
!===================================================================================================================================
  subroutine InitBasisTrans(N_in,xGP)
    ! MODULES
    use MOD_Indicators_vars         ,ONLY: sVdm_Leg, NumberOfPerssonInd
    use MOD_Basis                   ,ONLY: BuildLegendreVdm
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    INTEGER,INTENT(IN)                         :: N_in
    REAL,INTENT(IN),DIMENSION(0:N_in)          :: xGP
    REAL,DIMENSION(0:N_in,0:N_in)              :: Vdm_Leg
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !===================================================================================================================================
    
    if (NumberOfPerssonInd>0) return
    
    ! Compute the 1D Vandermondematrix, needed to tranform the nodal basis into a modal (Legendre) basis
    allocate(sVdm_Leg(0:N_in,0:N_in))
    call BuildLegendreVdm(N_in,xGP,Vdm_Leg,sVdm_Leg)
  end subroutine InitBasisTrans

end module MOD_Indicators
