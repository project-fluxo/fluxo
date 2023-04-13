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

!==================================================================================================================================
!> This module only initializes the equation specific parameters and computes analytical functions and the evaluation of the source
!==================================================================================================================================
MODULE MOD_Indicators_vars
    use MOD_Equation_Vars, only: i_indicatorFunction
    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: Indicator_Generic
    PUBLIC :: sVdm_Leg
    PUBLIC :: NumberOfPerssonInd
    SAVE
!============================================================================================================================
!> Generic Indicator class
!============================================================================================================================
    type Indicator_Generic
      logical            :: InitIsDone =.FALSE. ! True if indicator has been initialized
      logical            :: ReturnLog =.FALSE.  ! True if the
      character(len=255) :: IndicatorQuantity   ! Name of the indicator quantity
      procedure(i_indicatorFunction), nopass, pointer :: GetIndicatorQuantity => null()
      contains
        procedure :: construct  => Generic_Construct
        procedure :: compute    => Generic_Compute
        procedure :: destruct   => Generic_Destruct
        procedure :: GetIndicator_3D  => Generic_GetIndicator_3D
    end type Indicator_Generic
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
    REAL,ALLOCATABLE :: sVdm_Leg(:,:)                        !< 1D inverse Vandermondematrix to Legendre polynomials
    integer          :: NumberOfPerssonInd =0

 CONTAINS

!===================================================================================================================================
!> Constructor for Indicator_Generic 
!===================================================================================================================================
  subroutine Generic_Construct(this,IndicatorQuantity,ReturnLog)
    ! IMPLICIT VARIABLE HANDLING
    USE MOD_Globals
    use MOD_Equation_Vars     , only: SetIndicatorFunction
    IMPLICIT NONE
    !-arguments------------------------------------------------------------------------------------------------------------------
    class(Indicator_Generic), intent(inout) :: this
    character(len=255)             , intent(in)    :: IndicatorQuantity
    logical, optional              , intent(in)    :: ReturnLog
    !============================================================================================================================
    
    ! Initial checks
    IF (this % InitIsDone) THEN
      SWRITE(*,*) "GenericIndicator_Construct already called."
      RETURN
    END IF
    
    ! Set ReturnLog if necessary
    if ( present(ReturnLog) ) this % ReturnLog = ReturnLog
    
    ! Assign indicator quantity
    this % IndicatorQuantity = IndicatorQuantity
    call SetIndicatorFunction(IndicatorQuantity,this % GetIndicatorQuantity)
    
  end subroutine Generic_Construct
!===================================================================================================================================
!> Destructor for Indicator_Generic 
!===================================================================================================================================
  subroutine Generic_Destruct(this)
    implicit none
    class(Indicator_Generic), intent(inout) :: this
    
    this % GetIndicatorQuantity => null()
    this % ReturnLog =.FALSE.
    this % IndicatorQuantity = ''
    
  end subroutine Generic_Destruct
    
!============================================================================================================================
!> Compute the Indicator_Generic
!============================================================================================================================
  function Generic_Compute(this,U) RESULT (eta)
    use MOD_PreProc
    use MOD_Mesh_Vars                       ,only: nElems
    implicit none
    ! Arguments
    !---------------------------------------------------------------------------------------------------------------------------------
    class(Indicator_Generic), intent(in) :: this
    real                           , intent(in) :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
    real                                        :: eta(nElems)
    ! Local variables
    !---------------------------------------------------------------------------------------------------------------------------------
    
    stop 'Compute not defined'
    
  end function Generic_Compute
!============================================================================================================================
!> Get the shock indicator quantity for all degrees of freedom of an element
!============================================================================================================================
  pure subroutine Generic_GetIndicator_3D(this,Uind,U)
    use MOD_PreProc
    implicit none
    !-arguments-------------------------------------------
    class(Indicator_Generic), intent(in) :: this
    real, intent(in)  :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    real, intent(out) :: Uind (1:1,0:PP_N,0:PP_N,0:PP_N)
    !-local-variables-------------------------------------
    integer :: i,j,k
    !-----------------------------------------------------
    
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      call this % GetIndicatorQuantity(U(:,i,j,k),Uind(1,i,j,k))
    end do       ; end do       ; end do
  end subroutine Generic_GetIndicator_3D
  
END MODULE MOD_Indicators_vars
